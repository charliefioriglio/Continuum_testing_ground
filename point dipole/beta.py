import numpy as np
from tqdm import tqdm
from build_mode import build_directional_wavefunction
from CuO_do import (
    build_DO,
    DO_coeffs_b1_L,
    DO_coeffs_b1_R,
    DO_coeffs_b2_L,
    DO_coeffs_b2_R,
    DO_coeffs_a1_R,
    DO_coeffs_a1_L,
)

import matplotlib.pyplot as plt
from ang_grid import repulsion_orientations

# ---------------------------------------------
# Define Grid
N_ang = 20
L = 10
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

# ---------------------------------------------
# Build continuum
# ---------------------------------------------
def continuum(k_hat, D, L_max, eKE, X, Y, Z):
    """Build the summed continuum for a specific emission direction."""

    return build_directional_wavefunction(
        D,
        k_hat,
        eKE,
        X,
        Y,
        Z,
        l_max=L_max,
    )

# ---------------------------------------------
# Build rotation matrix from Euler angles (ZYZ convention)
# ---------------------------------------------
def rotation_matrix(alpha, beta, gamma):
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta),  np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)
    Rz1 = np.array([[ ca, -sa, 0], [ sa,  ca, 0], [  0,   0, 1]])
    Ry  = np.array([[ cb,   0, sb], [  0,   1,  0], [-sb,  0, cb]])
    Rz2 = np.array([[ cg, -sg, 0], [ sg,  cg, 0], [  0,   0, 1]])
    return Rz1 @ Ry @ Rz2

def rotate_grid(X, Y, Z, alpha, beta, gamma):
    R = rotation_matrix(alpha, beta, gamma)

    # Flatten grid arrays
    shape = X.shape
    coords = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=0)  # shape (3, N)

    # Apply rotation: rotated_coords = R @ coords
    rotated_coords = R @ coords

    # Reshape to original grid shape
    X_rot = rotated_coords[0].reshape(shape)
    Y_rot = rotated_coords[1].reshape(shape)
    Z_rot = rotated_coords[2].reshape(shape)

    return X_rot, Y_rot, Z_rot

# ---------------------------------------------
# Orientational averaging routine
# ---------------------------------------------
def compute_average_amplitudes(
    X,
    Y,
    Z,
    DO_coeffs_L,
    DO_coeffs_R,
    ang_grid,
    mu_field,
    continuum_field,
):
    """Rotate the molecule/DOs while keeping continuum + polarization fixed in lab frame."""

    N = ang_grid.shape[0]
    amps = np.zeros(N, dtype=complex)

    for i in range(N):
        alpha, beta, gamma = ang_grid[i]
        X_rot, Y_rot, Z_rot = rotate_grid(X, Y, Z, alpha, beta, gamma)

        DO_L = build_DO(DO_coeffs_L, X_rot, Y_rot, Z_rot, dV)
        DO_R = build_DO(DO_coeffs_R, X_rot, Y_rot, Z_rot, dV)

        integrand_L = np.conj(DO_L) * mu_field * continuum_field
        integrand_R = np.conj(continuum_field) * mu_field * DO_R
        A_L = integrate_3d(integrand_L, dV)
        A_R = integrate_3d(integrand_R, dV)

        amps[i] = A_L * A_R

    return np.real(amps)

# ---------------------------------------------
# Main script
# ---------------------------------------------
if __name__ == '__main__':

    directions = {
        "par": np.array([0.0, 0.0, 1.0]),
        "perp1": np.array([1.0, 0.0, 0.0]),
        "perp2": np.array([0.0, 1.0, 0.0]),
    }

    ang_grid = repulsion_orientations(n_orientations=N_ang)

    E_eV = np.linspace(0.01, 1.51, 20)
    L_max = 5

    D_values = np.array([0.0, 0.3])
    degeneracies = [
        (DO_coeffs_b1_L, DO_coeffs_b1_R),
        (DO_coeffs_b2_L, DO_coeffs_b2_R)
    ]

    pol_lab = np.array([0.0, 0.0, 1.0])
    mu_field = pol_lab[0] * X + pol_lab[1] * Y + pol_lab[2] * Z

    results = {}

    for D in D_values:
        print(f"\nCalculating for D = {D}")
        results[D] = []

        for E in tqdm(E_eV, desc=f"D={D}"):
            E_au = E / hartree

            continuum_cache = {}
            for key, direction in directions.items():
                continuum_cache[key] = continuum(direction, D, L_max, E_au, X, Y, Z)

            sigma_par = 0.0
            sigma_perp_accum = 0.0

            for DO_L, DO_R in degeneracies:
                A_par = compute_average_amplitudes(
                    X,
                    Y,
                    Z,
                    DO_L,
                    DO_R,
                    ang_grid,
                    mu_field,
                    continuum_cache["par"],
                )
                A_perp1 = compute_average_amplitudes(
                    X,
                    Y,
                    Z,
                    DO_L,
                    DO_R,
                    ang_grid,
                    mu_field,
                    continuum_cache["perp1"],
                )
                A_perp2 = compute_average_amplitudes(
                    X,
                    Y,
                    Z,
                    DO_L,
                    DO_R,
                    ang_grid,
                    mu_field,
                    continuum_cache["perp2"],
                )

                sigma_par += np.mean(A_par)
                sigma_perp_accum += 0.5 * (np.mean(A_perp1) + np.mean(A_perp2))

            sigma_perp = sigma_perp_accum
            beta = 2 * (sigma_par - sigma_perp) / (sigma_par + 2 * sigma_perp)

            results[D].append((E, sigma_par, sigma_perp, beta))

# ---------------------------------------------
# Save Results as CSV
# ---------------------------------------------
    
    for D in D_values:
        E_vals = np.array([row[0] for row in results[D]])
        beta_vals = np.array([row[3] for row in results[D]])
        
        # Create CSV with eKE and beta columns
        csv_data = np.column_stack((E_vals, beta_vals))
        csv_filename = f"CuO_point_dipole_beta_D_{2*D:.2f}.csv"
        np.savetxt(csv_filename, csv_data, delimiter=",", 
                   header="eKE_eV,beta", comments="")
        print(f"Saved {csv_filename}")

# ---------------------------------------------
# Print Results and Create Comparative Plot
# ---------------------------------------------
    
    for D in D_values:
        print(f"\nD = {2*D}")
        print(f"{'E_KE (eV)':<10} {'sigma_par':>15} {'sigma_perp':>15} {'beta':>10}")
        for E, s_par, s_perp, beta in results[D]:
            print(f"{E:10.2f} {s_par:15.6e} {s_perp:15.6e} {beta:10.4f}")

    # Create comparative plot
    plt.figure(figsize=(10, 6))
    
    colors = ['blue', 'red', 'green', 'purple']
    markers = ['o', 's', '^', 'd']
    
    for i, D in enumerate(D_values):
        E_vals = [row[0] for row in results[D]]
        beta_vals = [row[3] for row in results[D]]
        
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]

        plt.plot(E_vals, beta_vals, 
            color=color, 
            marker=marker, 
                linewidth=2, 
                markersize=6,
                label=f'D = {2*D}')
    
    plt.xlabel('eKE (eV)', fontsize=15)
    plt.ylabel('β', fontsize=15)
    plt.ylim(-1, 2)
    plt.legend(fontsize=11)
    
    plt.show()

