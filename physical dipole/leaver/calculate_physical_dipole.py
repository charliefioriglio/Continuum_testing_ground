import numpy as np
from tqdm import tqdm
from build_mode_annie import psi_el
import sys
import os

# Add the parent directory to the path to import CuO_do
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'point dipole'))
from CuO_do import build_DO, DO_coeffs_b1_L, DO_coeffs_b2_L, DO_coeffs_a1_L, DO_coeffs_b1_R, DO_coeffs_b2_R, DO_coeffs_a1_R, b2_L_norm, b2_R_norm, b1_L_norm, b1_R_norm, a1_L_norm, a1_R_norm, a

import matplotlib.pyplot as plt

# ---------------------------------------------
# Define Grid
# ---------------------------------------------
L = 19
n_pts = 50
x = np.linspace(-L, L, n_pts)
y = np.linspace(-L, L, n_pts)
z = np.linspace(-L, L, n_pts) 
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]
dV = dx * dy * dz
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

def integrate_3d(f, dV):
    return np.sum(f) * dV

# Constants
hartree = 27.2114  # eV
c = 137.036        # a.u. speed of light

# ---------------------------------------------
# Build continuum mode for physical dipole case
# ---------------------------------------------
def build_continuum_mode_physical_dipole(D, m, n, eKE, x, y, z, a):
    try:
        psi_continuum, _ = psi_el(n, m, eKE, a, D, x, y, z)
        
        if psi_continuum is None:
            return None
            
        return psi_continuum
        
    except Exception as e:
        print(f"Warning: Failed to build continuum mode (D={D}, m={m}, n={n}, E={eKE}): {e}")
        return None

# ---------------------------------------------
# Cross section calculation
# ---------------------------------------------
def calculate_total_cross_sections_physical_dipole(E_photon_grid, Trans_E, x, y, z, DO_coeffs_L, DO_coeffs_R, DO_L_norm, DO_R_norm, L_max, D, a):
    DO_L = build_DO(DO_coeffs_L, x, y, z, dV)
    DO_R = build_DO(DO_coeffs_R, x, y, z, dV)
    total_cross_sections = []

    for E_photon in tqdm(E_photon_grid):
        eKE = (E_photon - Trans_E) / hartree
        k_mag = np.sqrt(2 * eKE)
        prefactor = 8 * np.pi * k_mag * (E_photon / hartree) / c
        total_A = 0.0
        
        for m in range(-L_max, L_max + 1):
            for N in range(L_max + 1 - abs(m)):
                Cklm_squared = 0.0
                for pol in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
                    
                    # Build the continuum wavefunction
                    psi_el = build_continuum_mode_physical_dipole(D, m, N, eKE, x, y, z, a)
                    
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
def main():
    # Vibrational transitions data (same as in original)
    vib_transitions = np.array([
        [1.7780, 8.492010e-01],
        [1.8574, 5.006826e-01],
        [1.9367, 1.656282e-01],
        [2.0161, 2.703144e-02],
    ])

    Trans_E = vib_transitions[:, 0]
    FC_factors = vib_transitions[:, 1]
    E_grid = np.linspace(0.01, 5.01, 10)
    E_photon_grid = Trans_E[0] + E_grid
    L_max = 3
    
    # Physical dipole parameters
    D_values = np.array([0.0])
    a_focal = a

    # Calculate for different dipole moments
    for D in D_values:
        total_cross_sections = calculate_total_cross_sections_physical_dipole(E_photon_grid, Trans_E[0], X, Y, Z, DO_coeffs_b1_L, DO_coeffs_b1_R, b1_L_norm, b1_R_norm, L_max, D, a_focal)
            
        print(f"{'E_photon':>12} {'s_total':>15}")
        print("-" * 28)
        for E, s in zip(E_photon_grid, total_cross_sections):
            print(f"{E:12.6f} {s:15.6f}")

if __name__ == "__main__":
    main()
