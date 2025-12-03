import numpy as np
from scipy.special import sph_harm

from angular import build_dipole_angular_matrix
from radial import radial_function


_SPHERICAL_CACHE = {}
_OMEGA_CACHE = {}
_EIGEN_CACHE = {}


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
    key = _grid_key(X, Y, Z)
    if key in _SPHERICAL_CACHE:
        return _SPHERICAL_CACHE[key]

    r = np.sqrt(X**2 + Y**2 + Z**2)
    theta = np.zeros_like(r)
    phi = np.zeros_like(r)

    mask = r > 1.0e-12
    theta[mask] = np.arccos(np.clip(Z[mask] / r[mask], -1.0, 1.0))
    phi[mask] = np.mod(np.arctan2(Y[mask], X[mask]), 2 * np.pi)

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
    _OMEGA_CACHE[key] = (omega_funcs, eigvals)
    return omega_funcs, eigvals


def build_continuum_mode(D, lam, N, E, X, Y, Z, l_max=10):
    r_grid, theta_grid, phi_grid = cartesian_to_spherical(X, Y, Z)
    omega_funcs, eigvals = _omega_fields(D, lam, l_max, theta_grid, phi_grid)

    if N >= len(eigvals):
        raise ValueError(f"Mode index N={N} exceeds available eigenvalues ({len(eigvals)})")

    if eigvals[N] < -0.25 or np.isnan(eigvals[N]):
        return None

    L_N = 0.5 * (-1.0 + np.sqrt(1.0 + 4.0 * eigvals[N]))
    k = np.sqrt(2.0 * E)
    R_N = radial_function(L_N, k, r_grid)

    psi = R_N * omega_funcs[N]
    return psi


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

