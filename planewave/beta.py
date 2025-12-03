import time
import numpy as np
import matplotlib.pyplot as plt

try:
    from ang_grid import repulsion_orientations
except ImportError:  # optional helper, so keep going if absent
    repulsion_orientations = None

from CuO_do import (
    DO_coeffs_a1_L,
    DO_coeffs_a1_R,
    build_DO
)

# ---------------------------------------------
# Define spatial grid (kept identical to CN.py defaults)
# ---------------------------------------------
L = 10
n_pts = 100
x = np.linspace(-L, L, n_pts)
y = np.linspace(-L, L, n_pts)
z = np.linspace(-L, L, n_pts)
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]
dV = dx * dy * dz
X, Y, Z = np.meshgrid(x, y, z, indexing="ij")


# ---------------------------------------------
# Basic utilities
# ---------------------------------------------
def integrate_3d(f, dV):
    return np.sum(f) * dV


def plane_wave(k_vec, X, Y, Z):
    k_vec = np.asarray(k_vec, dtype=float)
    k_dot_r = k_vec[0] * X + k_vec[1] * Y + k_vec[2] * Z
    return np.exp(1j * k_dot_r)


def rotation_matrix(alpha, beta, gamma):
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)
    Rz1 = np.array([[ca, -sa, 0], [sa, ca, 0], [0, 0, 1]])
    Ry = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])
    Rz2 = np.array([[cg, -sg, 0], [sg, cg, 0], [0, 0, 1]])
    return Rz1 @ Ry @ Rz2


def rotate_grid(X, Y, Z, alpha, beta, gamma):
    R = rotation_matrix(alpha, beta, gamma)
    shape = X.shape
    coords = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=0)
    rotated_coords = R @ coords
    X_rot = rotated_coords[0].reshape(shape)
    Y_rot = rotated_coords[1].reshape(shape)
    Z_rot = rotated_coords[2].reshape(shape)
    return X_rot, Y_rot, Z_rot


# ---------------------------------------------
# Orientational averaging with rotating Dyson orbitals
# ---------------------------------------------
def compute_average_amplitudes_rotating_DO(
    X,
    Y,
    Z,
    DO_coeffs_L,
    DO_coeffs_R,
    ang_grid,
    mu_field,
    continuum_field,
):
    """Rotate the molecule/DOs while keeping the continuum fixed in the lab frame."""

    if ang_grid.shape[0] == 3 and ang_grid.shape[1] != 3:
        ang_grid_use = ang_grid.T
    else:
        ang_grid_use = ang_grid

    N = ang_grid_use.shape[0]
    amps = np.zeros(N, dtype=complex)

    for i in range(N):
        alpha, beta, gamma = ang_grid_use[i]
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
# Main calculation wrapper
# ---------------------------------------------
def calculate_beta(
    E_eV,
    ang_grid,
    degeneracies=None,
    pol_lab=None,
):
    """Compute beta(E) using a plane-wave continuum and rotating DOs."""

    hartree = 27.2114

    if degeneracies is None:
        degeneracies = [
            (DO_coeffs_a1_L, DO_coeffs_a1_R)
        ]

    if pol_lab is None:
        pol_lab = np.array([0.0, 0.0, 1.0])

    pol_lab = np.asarray(pol_lab, dtype=float)
    mu_field = pol_lab[0] * X + pol_lab[1] * Y + pol_lab[2] * Z

    directions = {
        "par": np.array([0.0, 0.0, 1.0]),
        "perp1": np.array([1.0, 0.0, 0.0]),
        "perp2": np.array([0.0, 1.0, 0.0]),
    }

    results = []

    for E in E_eV:
        E_au = E / hartree
        k_mag = np.sqrt(2.0 * E_au)

        continuum_cache = {}
        for key, direction in directions.items():
            k_vec = k_mag * direction
            continuum_cache[key] = plane_wave(k_vec, X, Y, Z)

        sigma_par = 0.0
        sigma_perp_accum = 0.0

        for DO_L, DO_R in degeneracies:
            A_par = compute_average_amplitudes_rotating_DO(
                X,
                Y,
                Z,
                DO_L,
                DO_R,
                ang_grid,
                mu_field,
                continuum_cache["par"],
            )
            A_perp1 = compute_average_amplitudes_rotating_DO(
                X,
                Y,
                Z,
                DO_L,
                DO_R,
                ang_grid,
                mu_field,
                continuum_cache["perp1"],
            )
            A_perp2 = compute_average_amplitudes_rotating_DO(
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
        results.append((E, sigma_par, sigma_perp, beta))

    return results


# ---------------------------------------------
# Example usage / plotting
# ---------------------------------------------
if __name__ == "__main__":
    start_time = time.time()

    # Configure angular grid. Use a repulsion grid when available, otherwise
    # fall back to the CN.py benchmark grid to keep parity with older runs.
    USE_REPULSION_GRID = True
    N_ORIENTATIONS = 20

    if USE_REPULSION_GRID and repulsion_orientations is not None:
        ang_grid = repulsion_orientations(n_orientations=N_ORIENTATIONS)
    else:
        ang_grid = np.array([
            [1.07393, 3.31052, 1.74962, 1.64604, 6.10750, 1.75961, 5.53027, 5.52842, 2.28479, 3.58752],
            [2.37487, 2.01802, 0.19298, 1.42732, 0.94148, 2.87768, 0.43346, 2.23837, 2.38175, 1.86742],
            [0.00649606, 0.00677567, 0.00640311, 0.00666473, 0.00669922, 0.00663394, 0.00674343, 0.00667127, 0.00668776, 0.00676721]
        ]).T

    # Photoelectron kinetic energies
    E_eV = np.array([0.1, 0.3])

    print(
        f"Calculating beta (plane-wave continuum, rotating DOs) for {len(E_eV)} energies "
        f"with {ang_grid.shape[0]} orientations..."
    )

    results = calculate_beta(E_eV, ang_grid)

    print("\n" + "=" * 50)
    print(f"{'E_KE (eV)':<12} {'sigma_par':<12} {'sigma_perp':<12} {'beta':<12}")
    print("=" * 50)
    for E, s_par, s_perp, beta in results:
        print(f"{E:12.4f} {s_par:12.6e} {s_perp:12.6e} {beta:12.6f}")
    print("=" * 50)

    end_time = time.time()
    print(f"\nCalculation time: {end_time - start_time:.2f} seconds")

    results_array = np.array(results)

    plt.figure(figsize=(8, 6))
    plt.plot(results_array[:, 0], results_array[:, 3], marker="o", linewidth=2, markersize=6)
    plt.xlabel("Photoelectron KE (eV)", fontsize=12)
    plt.ylabel(r"$\beta$", fontsize=12)
    plt.title("Anisotropy Parameter (plane-wave continuum)", fontsize=13)
    plt.ylim(-1, 2)
    plt.axhline(y=0, color="k", linestyle="--", linewidth=0.8, alpha=0.5)
    plt.tight_layout()
    plt.show()
