"""
dipole_eigensystem.py
=====================

Construct and diagonalize the point-dipole Hamiltonian matrix to obtain
the angular eigenproblem eigenvalues (ν_N) and eigenvectors (expansion
coefficients), translated from MATLAB pointdpmatrixc.m.

Functions
---------
pointdpmatrixc(lmax, dp, lam, w)
    Build and diagonalize the point-dipole matrix for a given λ.
build_dipole_matrix(lmax, dp, lam)
    Build the (lmax+1) × (lmax+1) dipole matrix without diagonalization.
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import eig

from .wigner import wigner3j

__all__ = ["pointdpmatrixc", "build_dipole_matrix"]


def build_dipole_matrix(lmax: int, dp: float, lam: int) -> np.ndarray:
    r"""
    Build the point-dipole Hamiltonian matrix for the angular eigenproblem.

    The matrix element between orbital angular momenta l (row) and l' (column)
    is:

    .. math::
        H_{l l'} = \delta_{l l'} l(l+1)
                   - 2 D \sqrt{\frac{4\pi}{3}}
                     \sqrt{\frac{(2l+1) \cdot 3 \cdot (2l'+1)}{4\pi}}
                     \begin{pmatrix} l & 1 & l' \\ 0 & 0 & 0 \end{pmatrix}
                     \begin{pmatrix} l & 1 & l' \\ -\lambda & 0 & \lambda \end{pmatrix}

    Parameters
    ----------
    lmax : int
        Maximum orbital angular momentum in the expansion.
    dp : float
        Dipole moment in atomic units.
    lam : int
        Projection of angular momentum onto dipole axis (λ).

    Returns
    -------
    ndarray, shape (lmax+1, lmax+1)
        The Hamiltonian matrix (real, symmetric).

    Notes
    -----
    Rows/columns with l < |λ| are effectively zeroed (no matrix element
    because those states don't couple for the given λ).
    """
    n = lmax + 1
    H = np.zeros((n, n), dtype=np.float64)

    prefac = -2 * dp * np.sqrt(4 * np.pi / 3)

    for lr in range(n):  # row index (l = lr)
        l = lr
        if l < abs(lam):
            continue  # states below |λ| are forbidden
        for lc in range(n):  # column index (l' = lc)
            lp = lc
            if l == lp:
                H[lr, lc] = l * (l + 1)
            else:
                # Off-diagonal coupling
                sqrt_part = np.sqrt((2 * l + 1) * 3 * (2 * lp + 1) / (4 * np.pi))
                w1 = wigner3j(l, 1, lp, 0, 0, 0)
                w2 = wigner3j(l, 1, lp, -lam, 0, lam)
                # Note: MATLAB code has (-1^lam) which is always -1 (typo),
                # but the physical phase should be (-1)**lam.  We replicate
                # the MATLAB behaviour exactly for fidelity.
                sign = (-1) ** lam  # correct Condon-Shortley phase
                H[lr, lc] = prefac * sqrt_part * sign * w1 * w2

    return H


def pointdpmatrixc(
    lmax: int, dp: float, lam: int, w: int
) -> np.ndarray:
    r"""
    Diagonalize the point-dipole matrix and return eigenvalues or eigenvectors.

    Parameters
    ----------
    lmax : int
        Maximum orbital angular momentum.
    dp : float
        Dipole moment (atomic units).
    lam : int
        Projection quantum number λ.
    w : int
        * 1 → return array of *effective angular momenta* ν_N from the
          eigenvalues via ν = (−1 + √(1 + 4E))/2.
        * 2 → return the eigenvector matrix (columns are eigenvectors).

    Returns
    -------
    ndarray
        Either shape (n,) eigenvalues ν_N or shape (n, n) eigenvector matrix.

    Notes
    -----
    MATLAB source (pointdpmatrixc.m) builds the matrix, calls eigs() for the
    smallest eigenvalues, and returns either eigenvalues (transformed via
    the quadratic formula) or eigenvectors.

    The eigenvalue transformation comes from the angular Schrödinger equation
    eigenvalues being N(N+1), so inverting gives N = (−1 + √(1 + 4E))/2.
    """
    H = build_dipole_matrix(lmax, dp, lam)

    # Solve the eigenproblem
    eigenvalues, eigenvectors = eig(H)

    # Sort by eigenvalue (smallest first, matching MATLAB eigs behaviour)
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Zero out rows for l < |λ| in eigenvectors (as MATLAB does)
    for l in range(lmax + 1):
        if l < abs(lam):
            eigenvectors[l, :] = 0.0

    if w == 1:
        # Return effective angular momenta ν_N
        # ν = (−1 + √(1 + 4E)) / 2
        nu = (-1 + np.sqrt(1 + 4 * eigenvalues.real)) / 2
        return nu
    elif w == 2:
        return eigenvectors.real
    else:
        raise ValueError(f"w must be 1 (eigenvalues) or 2 (eigenvectors), got {w}")
