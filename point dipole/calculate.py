from __future__ import annotations

import numpy as np
from tqdm import tqdm
from build_mode import build_directional_wavefunction
from CuO_do import (
    build_DO,
    DO_coeffs_b1_L,
    DO_coeffs_b2_L,
    DO_coeffs_a1_L,
    DO_coeffs_b1_R,
    DO_coeffs_b2_R,
    DO_coeffs_a1_R,
    b2_L_norm,
    b2_R_norm,
    b1_L_norm,
    b1_R_norm,
    a1_L_norm,
    a1_R_norm,
)
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------
# Define Grid
L = 10
n_pts = 100
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
POLARIZATION_DIRECTIONS = (
    np.array([1.0, 0.0, 0.0]),
    np.array([0.0, 1.0, 0.0]),
    np.array([0.0, 0.0, 1.0]),
)

# ---------------------------------------------
# Cross section calculation
# ---------------------------------------------
def calculate_cross_sections(
    E_photon_grid,
    Trans_E,
    FC_factors,
    x,
    y,
    z,
    DO_coeffs_L,
    DO_coeffs_R,
    DO_L_norm,
    DO_R_norm,
    L_max,
    D,
):
    DO_L = build_DO(DO_coeffs_L, x, y, z, dV)
    DO_R = build_DO(DO_coeffs_R, x, y, z, dV)
    n_vib = len(Trans_E)
    rel_cross_sections = np.zeros((len(E_photon_grid), n_vib))
    Xg, Yg, Zg = np.meshgrid(x, y, z, indexing="ij")

    continuum_cache: dict[tuple[int, float], np.ndarray] = {}

    def get_continuum(direction_idx: int, energy_au: float) -> np.ndarray:
        key = (direction_idx, float(np.round(energy_au, decimals=12)))
        if key not in continuum_cache:
            continuum_cache[key] = build_directional_wavefunction(
                D,
                POLARIZATION_DIRECTIONS[direction_idx],
                energy_au,
                Xg,
                Yg,
                Zg,
                l_max=L_max,
            )
        return continuum_cache[key]

    for j, E_photon in enumerate(E_photon_grid):
        weights = []
        for i in range(n_vib):
            E_bind = Trans_E[i]
            eKE_ev = E_photon - E_bind
            eKE = eKE_ev / hartree
            if eKE <= 0.0:
                weights.append(0.0)
                continue

            k_mag = np.sqrt(2.0 * eKE)
            fc2 = FC_factors[i] ** 2
            prefactor = fc2 * (8.0 * np.pi * k_mag * (E_photon / hartree) / c)

            channel_A = 0.0
            for idx, pol in enumerate(POLARIZATION_DIRECTIONS):
                continuum_field = get_continuum(idx, eKE)
                dipole = pol[0] * Xg + pol[1] * Yg + pol[2] * Zg
                integrand_L = np.conj(continuum_field) * dipole * DO_L
                integrand_R = np.conj(DO_R) * dipole * continuum_field
                amp_L = integrate_3d(integrand_L, dV)
                amp_R = integrate_3d(integrand_R, dV)
                channel_A += amp_L * amp_R

            channel_A /= len(POLARIZATION_DIRECTIONS)
            channel_A *= 2.0 * DO_L_norm * DO_R_norm

            sigma = prefactor * channel_A
            weights.append(float(np.real_if_close(sigma)))

        weights = np.array(weights)
        total = weights.sum()
        if total > 0.0:
            rel_cross_sections[j, :] = weights / total

    return rel_cross_sections

def calculate_total_cross_sections(
    E_photon_grid,
    Trans_E,
    x,
    y,
    z,
    DO_coeffs_L,
    DO_coeffs_R,
    DO_L_norm,
    DO_R_norm,
    L_max,
    D,
):
    DO_L = build_DO(DO_coeffs_L, x, y, z, dV)
    DO_R = build_DO(DO_coeffs_R, x, y, z, dV)
    total_cross_sections = []
    Xg, Yg, Zg = np.meshgrid(x, y, z, indexing="ij")

    continuum_cache: dict[tuple[int, float], np.ndarray] = {}

    def get_continuum(direction_idx: int, energy_au: float) -> np.ndarray:
        key = (direction_idx, float(np.round(energy_au, decimals=12)))
        if key not in continuum_cache:
            continuum_cache[key] = build_directional_wavefunction(
                D,
                POLARIZATION_DIRECTIONS[direction_idx],
                energy_au,
                Xg,
                Yg,
                Zg,
                l_max=L_max,
            )
        return continuum_cache[key]

    for E_photon in tqdm(E_photon_grid):
        eKE_ev = E_photon - Trans_E
        eKE = eKE_ev / hartree
        if eKE <= 0.0:
            total_cross_sections.append(0.0)
            continue

        k_mag = np.sqrt(2.0 * eKE)
        prefactor = 8.0 * np.pi * k_mag * (E_photon / hartree) / c

        channel_A = 0.0
        for idx, pol in enumerate(POLARIZATION_DIRECTIONS):
            continuum_field = get_continuum(idx, eKE)
            dipole = pol[0] * Xg + pol[1] * Yg + pol[2] * Zg
            integrand_L = np.conj(continuum_field) * dipole * DO_L
            integrand_R = np.conj(DO_R) * dipole * continuum_field
            amp_L = integrate_3d(integrand_L, dV)
            amp_R = integrate_3d(integrand_R, dV)
            channel_A += amp_L * amp_R

        channel_A /= len(POLARIZATION_DIRECTIONS)
        channel_A *= 2.0 * DO_L_norm * DO_R_norm

        sigma = prefactor * channel_A
        total_cross_sections.append(float(np.real_if_close(sigma)))

    return np.array(total_cross_sections)

# ---------------------------------------------
# Calculate and Plot
# ---------------------------------------------
vib_transitions = np.array([
    [1.7780, 8.492010e-01],
    [1.8574, 5.006826e-01],
    [1.9367, 2.056282e-01],
])

Trans_E = vib_transitions[:, 0]
FC_factors = vib_transitions[:, 1]
E_photon_grid = np.linspace(1.8, 2.3, 20)
L_max = 5
D = 0.3


rel_cross_sections = calculate_cross_sections(
    E_photon_grid,
    Trans_E,
    FC_factors,
    x,
    y,
    z,
    DO_coeffs_b1_L,
    DO_coeffs_b1_R,
    a1_L_norm,
    a1_R_norm,
    L_max,
    D,
)

# Save results to CSV
data_dict = {'Photon_Energy_eV': E_photon_grid}
for i in range(len(Trans_E)):
    data_dict[f'v={i}_RelCrossSection'] = rel_cross_sections[:, i]

df = pd.DataFrame(data_dict)
csv_filename = f'relative_cross_sections_D_{2*D:.1f}.csv'
df.to_csv(csv_filename, index=False)
print(f"\nResults saved to {csv_filename}")

# Plot results
for i in range(len(Trans_E)):
    plt.plot(E_photon_grid, rel_cross_sections[:, i], label=f'v={i}, D={2*D}')
plt.xlabel('Photon Energy (eV)')
plt.ylabel('Relative Cross Section (a.u.)')
plt.title('Relative Vibrational Cross Sections vs Photon Energy')
plt.legend()
plt.show()

