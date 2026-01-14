import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.special import sph_harm_y, spherical_jn
from CuO_do import build_DO, DO_coeffs_b1_L, DO_coeffs_b2_L, DO_coeffs_a1_L, DO_coeffs_b1_R, DO_coeffs_b2_R, DO_coeffs_a1_R, b1_L_norm, b1_R_norm, b2_L_norm, b2_R_norm, a1_L_norm, a1_R_norm
from build_mode import build_continuum_mode
# ---------------------------------------------
# Define Grid
L = 19
n_pts = 101
x = np.linspace(-L, L, n_pts)
y = np.linspace(-L, L, n_pts)
z = np.linspace(-L, L, n_pts) 
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]
dV = dx * dy * dz
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

def integrate_3d(f, dV):
    """3D integration using Simpson's rule along each axis."""
    from scipy.integrate import simpson
    # Integrate along each axis sequentially
    # Note: dV = dx * dy * dz, so dx = dV^(1/3) for uniform grids
    dx = dV ** (1/3)
    result = simpson(simpson(simpson(f, dx=dx, axis=0), dx=dx, axis=0), dx=dx, axis=0)
    return result

# Constants
hartree = 27.2114  # eV
c = 137.036        # a.u. speed of light

# ---------------------------------------------
# Define Point Dipole Continuum
# ---------------------------------------------
def r_theta_phi(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(np.clip(z / (r + 1e-12), -1.0, 1.0))
    phi = np.arctan2(y, x)
    return r, theta, phi

def planewave_expansion(k, l, m, x, y, z):
    r, theta, phi = r_theta_phi(x, y, z)
    R = spherical_jn(l, k * r)
    Y = sph_harm_y(l, m, theta, phi)
    return Y * R

# ---------------------------------------------
# Cross section calculation
# ---------------------------------------------
def calculate_cross_sections(E_photon_grid, Trans_E, FC_factors, x, y, z, DO_coeffs_L, DO_coeffs_R, L_max):
    DO_L = build_DO(DO_coeffs_L, x, y, z, dV)
    DO_R = build_DO(DO_coeffs_R, x, y, z, dV)
    polarizations = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]
    n_vib = len(Trans_E)
    rel_cross_sections = np.zeros((len(E_photon_grid), n_vib))

    for j, E_photon in tqdm(enumerate(E_photon_grid)):
        weights = []
        for i in range(n_vib):
            E_bind = Trans_E[i]
            eKE = (E_photon - E_bind) / hartree
            if eKE <= 0:
                weights.append(0.0)
                continue
            k_mag = np.sqrt(2 * eKE)
            fc2 = FC_factors[i]**2
            prefactor = fc2 * (8 * np.pi * k_mag * (E_photon / hartree) / c)

            total_A = 0
            for pol in polarizations:
                for l in range(L_max + 1):
                    for m in range(-l, l + 1):
                        dipole = pol[0]*x + pol[1]*y + pol[2]*z
                        psi_el = planewave_expansion(k_mag, l, m, x, y, z)
                        integrand_L = np.conj(psi_el) * dipole * DO_L
                        integrand_R = np.conj(DO_R) * dipole * psi_el
                        amp_L = integrate_3d(integrand_L, dV)
                        amp_R = integrate_3d(integrand_R, dV)
                        total_A += amp_L * amp_R

            sigma = prefactor * total_A
            weights.append(sigma)

        weights = np.array(weights)
        if np.sum(weights) > 0:
            rel_cross_sections[j, :] = weights / np.sum(weights)

    return rel_cross_sections

def calculate_total_cross_sections(E_photon_grid, Trans_E, x, y, z, DO_coeffs_L, DO_coeffs_R, DO_L_norm, DO_R_norm, L_max, D):
    DO_L = build_DO(DO_coeffs_L, x, y, z, dV)
    DO_R = build_DO(DO_coeffs_R, x, y, z, dV)
    total_cross_sections = []

    for E_photon in tqdm(E_photon_grid):
        eKE = (E_photon - Trans_E) / hartree
        k_mag = np.sqrt(2 * eKE)
        prefactor = 8 * np.pi * k_mag * (E_photon / hartree) / c
        total_A = 0.0
        
        for m in range(-L_max, L_max + 1):
            for n in range(L_max + 1 - abs(m)):
                Cklm_squared = 0.0
                for pol in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
                    psi_el = build_continuum_mode(D, m, n, eKE, x, y, z, l_max=L_max)
                    if psi_el is None:
                        continue
                    dipole = pol[0]*x + pol[1]*y + pol[2]*z
                    integrand_L = np.conj(psi_el) * dipole * DO_L
                    integrand_R = np.conj(DO_R) * dipole * psi_el
                    amp_L = integrate_3d(integrand_L, dV)
                    amp_R = integrate_3d(integrand_R, dV)
                    
                    Cklm_squared += amp_L * amp_R

                total_A += Cklm_squared / 3 * (2 * DO_L_norm * DO_R_norm)

        sigma = prefactor * total_A
        total_cross_sections.append(sigma)

    return np.real(total_cross_sections)


# ---------------------------------------------
# Calculate and Plot
# ---------------------------------------------
vib_transitions = np.array([
    [1.7780, 8.492010e-01],
    [1.8574, 5.006826e-01],
    [1.9367, 1.656282e-01],
    [2.0161, 2.703144e-02],
])

Trans_E = vib_transitions[:, 0]
FC_factors = vib_transitions[:, 1]
E_grid = np.linspace(0.1, 10.1, 21)
E_photon_grid = Trans_E[0] + E_grid
L_max = 4

# D values from 0.0 to 0.6 (dipole strength = 2*D)
D_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]

# Calculate cross sections for each D value
results = {'E_photon': E_photon_grid}
for D in D_values:
    print(f"\nCalculating for D = {2*D}...")
    cross_sections = calculate_total_cross_sections(
        E_photon_grid, Trans_E[0], X, Y, Z,
        DO_coeffs_a1_L, DO_coeffs_a1_R,
        a1_L_norm, a1_R_norm, L_max, D
    )
    results[str(2*D)] = cross_sections

# Save to CSV
df = pd.DataFrame(results)
df.to_csv('total_cross_sections_a1.csv', index=False)
print("\nSaved results to 'total_cross_sections_a1.csv'")
