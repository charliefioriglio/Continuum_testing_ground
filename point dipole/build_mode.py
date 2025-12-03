from __future__ import annotations

import numpy as np
from scipy.special import sph_harm

from angular import build_dipole_angular_matrix
from radial import radial_function


_SPHERICAL_CACHE: dict[tuple, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
_OMEGA_CACHE: dict[tuple, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
_EIGEN_CACHE: dict[tuple, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}


def _ensure_cartesian_grid(X, Y, Z):
    """Return 3D mesh arrays even if 1D axes are provided."""

    X_arr = np.asarray(X)
    Y_arr = np.asarray(Y)
    Z_arr = np.asarray(Z)
    if X_arr.ndim == Y_arr.ndim == Z_arr.ndim == 1:
        return np.meshgrid(X_arr, Y_arr, Z_arr, indexing="ij")
    return X_arr, Y_arr, Z_arr


def _grid_key(X, Y, Z):
    return (
        id(X),
        id(Y),
        id(Z),
        X.shape,
        X.strides,
        Y.shape,
        Z.shape,
    )


def cartesian_to_spherical(X, Y, Z):
    Xg, Yg, Zg = _ensure_cartesian_grid(X, Y, Z)
    key = _grid_key(Xg, Yg, Zg)
    if key in _SPHERICAL_CACHE:
        return _SPHERICAL_CACHE[key]

    r = np.sqrt(Xg**2 + Yg**2 + Zg**2)
    theta = np.zeros_like(r)
    phi = np.zeros_like(r)

    mask = r > 1.0e-12
    theta[mask] = np.arccos(np.clip(Zg[mask] / r[mask], -1.0, 1.0))
    phi[mask] = np.mod(np.arctan2(Yg[mask], Xg[mask]), 2 * np.pi)

    _SPHERICAL_CACHE[key] = (r, theta, phi)
    return r, theta, phi


def _get_eigensystem(D, lam, l_max):
    key = (float(D), int(lam), int(l_max))
    if key in _EIGEN_CACHE:
        return _EIGEN_CACHE[key]

    H, l_vals = build_dipole_angular_matrix(D, lam, l_max)
    eigvals, eigvecs = np.linalg.eigh(H)
    data = (eigvals, eigvecs, l_vals)
    _EIGEN_CACHE[key] = data
    return data


def _omega_fields(D, lam, l_max, theta_grid, phi_grid):
    key = (float(D), int(lam), int(l_max), id(theta_grid), id(phi_grid))
    if key in _OMEGA_CACHE:
        return _OMEGA_CACHE[key]

    eigvals, eigvecs, l_vals = _get_eigensystem(D, lam, l_max)
    harmonics = np.array([sph_harm(lam, l, phi_grid, theta_grid) for l in l_vals])
    omega_funcs = np.tensordot(eigvecs.T, harmonics, axes=([1], [0]))
    data = (omega_funcs, eigvals, eigvecs, l_vals)
    _OMEGA_CACHE[key] = data
    return data


def build_continuum_mode(D, lam, N, E, X, Y, Z, l_max=10):
    r_grid, theta_grid, phi_grid = cartesian_to_spherical(X, Y, Z)
    if E <= 0.0:
        return np.zeros_like(r_grid, dtype=np.complex128)

    omega_funcs, eigvals, _, _ = _omega_fields(D, lam, l_max, theta_grid, phi_grid)

    if N >= len(eigvals):
        raise ValueError(f"Mode index N={N} exceeds available eigenvalues ({len(eigvals)})")

    if eigvals[N] < -0.25 or np.isnan(eigvals[N]):
        return None

    L_N = 0.5 * (-1.0 + np.sqrt(1.0 + 4.0 * eigvals[N]))
    k = np.sqrt(2.0 * E)
    R_N = radial_function(L_N, k, r_grid)

    psi = R_N * omega_funcs[N]
    return psi


def build_directional_wavefunction(
    D,
    k_vec,
    energy_au,
    X,
    Y,
    Z,
    *,
    l_max=10,
    max_modes_per_m=None,
    coefficient_threshold=1.0e-12,
):
    """Assemble the summed continuum wavefunction for a specific emission direction."""

    r_grid, theta_grid, phi_grid = cartesian_to_spherical(X, Y, Z)
    if energy_au <= 0.0:
        return np.zeros_like(r_grid, dtype=np.complex128)

    k_vec = np.asarray(k_vec, dtype=float)
    if k_vec.shape != (3,):
        raise ValueError("k_vec must be a 3-element vector")
    norm = np.linalg.norm(k_vec)
    if norm == 0.0:
        raise ValueError("k_vec must define a non-zero direction")
    k_hat = k_vec / norm

    theta_k = np.arccos(np.clip(k_hat[2], -1.0, 1.0))
    phi_k = np.mod(np.arctan2(k_hat[1], k_hat[0]), 2.0 * np.pi)
    k_mag = np.sqrt(2.0 * energy_au)

    total = np.zeros_like(r_grid, dtype=np.complex128)
    for lam in range(-l_max, l_max + 1):
        omega_funcs, eigvals, eigvecs, l_vals = _omega_fields(D, lam, l_max, theta_grid, phi_grid)
        if eigvals.size == 0:
            continue

        harmonics_dir = np.array([sph_harm(lam, l, phi_k, theta_k) for l in l_vals])
        omega_dir = eigvecs.T @ harmonics_dir

        mode_limit = eigvals.size if max_modes_per_m is None else min(int(max_modes_per_m), eigvals.size)
        for mode_idx in range(mode_limit):
            eigval = eigvals[mode_idx]
            if eigval < -0.25 or not np.isfinite(eigval):
                continue

            coeff = omega_dir[mode_idx]
            if abs(coeff) < coefficient_threshold:
                continue

            L_eff = 0.5 * (-1.0 + np.sqrt(1.0 + 4.0 * eigval))
            radial = radial_function(L_eff, k_mag, r_grid)
            phase = np.exp(0.5j * np.pi * L_eff)
            weight = 4.0 * np.pi * phase * np.conjugate(coeff)
            total += weight * radial * omega_funcs[mode_idx]

    return total


def coefficient(D, lam, N, k_hat, l_max=10):
    k_hat = np.asarray(k_hat, dtype=float)
    norm = np.linalg.norm(k_hat)
    if norm == 0.0:
        raise ValueError("k_hat must be non-zero")
    k_hat = k_hat / norm

    kx, ky, kz = k_hat
    theta_k = np.arccos(np.clip(kz, -1.0, 1.0))
    phi_k = np.mod(np.arctan2(ky, kx), 2 * np.pi)

    eigvals, eigvecs, l_vals = _get_eigensystem(D, lam, l_max)
    if N >= len(eigvals):
        raise ValueError(f"Mode index N={N} exceeds available eigenvalues ({len(eigvals)})")

    harmonics = np.array([sph_harm(lam, l, phi_k, theta_k) for l in l_vals])
    omega_vals = eigvecs.T @ harmonics
    L_N = 0.5 * (-1.0 + np.sqrt(1.0 + 4.0 * eigvals[N]))
    phase = np.exp(0.5j * np.pi * L_N)
    coeff = 4.0 * np.pi * phase * np.conjugate(omega_vals[N])
    return coeff

