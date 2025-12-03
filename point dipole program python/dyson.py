"""
dyson.py
========

Build Dyson orbitals from a Gaussian basis set expansion, translated from
MATLAB Dyson.m.

Functions
---------
dyson(molecule_data, ps, step, X, Y, Z, tol, Dystol)
    Construct left and right Dyson orbitals on a Cartesian grid.
"""

from __future__ import annotations

from typing import Callable, Tuple

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .harmonics import Ylm
from .integration import simp3D
from .utilities import dfac

__all__ = ["dyson"]

# Bohr radius in Å
A0_ANGSTROM = 0.5291772109


def _double_factorial_odd(n: int) -> int:
    """Return (2n-1)!! = 1 * 3 * 5 * ... * (2n-1)."""
    return dfac(2 * n - 1)


def _double_factorial_odd_plus(n: int) -> int:
    """Return (2n+1)!! = 1 * 3 * 5 * ... * (2n+1)."""
    return dfac(2 * n + 1)


def dyson(
    molecule_data: Callable[[], Tuple],
    ps: int,
    step: float,
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    tol: int = 5,
    Dystol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, int]:
    r"""
    Build left and right Dyson orbitals on a Cartesian grid.

    Parameters
    ----------
    molecule_data : callable
        A function that returns (A, N, CZ, wL, wR, llim) where:
        - A: ndarray (natoms × 3) atom coordinates in Ångströms
        - N: int (unused, present for MATLAB compatibility)
        - CZ: list of basis set specifications per atom
        - wL: 1-D array of left Dyson coefficients for each AO
        - wR: 1-D array of right Dyson coefficients for each AO
        - llim: int, maximum angular momentum in basis
    ps : int
        Grid parameter; npts = ps + 1.
    step : float
        Grid spacing in atomic units.
    X, Y, Z : ndarray, shape (npts, npts, npts)
        Cartesian coordinate meshes.
    tol : int, optional
        Decimal places for rounding centroid (default 5).
    Dystol : float, optional
        Threshold for including AO contributions (default 1e-8).

    Returns
    -------
    orbL : ndarray (npts, npts, npts)
        Normalized, centered left Dyson orbital.
    orbR : ndarray (npts, npts, npts)
        Normalized, centered right Dyson orbital.
    llim : int
        Maximum angular momentum.

    Notes
    -----
    Translates MATLAB Dyson.m which:

    1. Loads atomic coordinates and converts from Ångströms → atomic units.
    2. Reads Gaussian primitive exponents and contraction coefficients.
    3. Normalizes primitives and contracted functions.
    4. Sums over AOs with left/right Dyson coefficients.
    5. Centers each orbital at the centroid of its density.
    """
    # ---------- 1. Load molecule data ----------
    A, _N, CZ, wL, wR, llim = molecule_data()

    # Convert coordinates to atomic units
    A = np.asarray(A, dtype=np.float64) / A0_ANGSTROM
    natoms = A.shape[0]
    npts = ps + 1

    # Precompute double factorials
    f2 = np.array([_double_factorial_odd(ell) for ell in range(llim + 1)])  # (2l-1)!!
    f2p = np.array([_double_factorial_odd_plus(ell) for ell in range(llim + 1)])  # (2l+1)!!

    # ---------- 2. Compute r² from each atom centre ----------
    rmat = np.zeros((npts, npts, npts, natoms), dtype=np.float64)
    for a in range(natoms):
        dx = X - A[a, 0]
        dy = Y - A[a, 1]
        dz = Z - A[a, 2]
        rmat[..., a] = dx**2 + dy**2 + dz**2

    # ---------- 3. Build real harmonics centred on each atom ----------
    n_m = 2 * llim + 1
    S = np.zeros((npts, npts, npts, llim + 1, n_m, natoms), dtype=np.float64)
    for a in range(natoms):
        dX = X - A[a, 0]
        dY = Y - A[a, 1]
        dZ = Z - A[a, 2]
        for ell in range(llim + 1):
            for m in range(-ell, ell + 1):
                S[..., ell, m + llim, a] = Ylm(ell, m, dX, dY, dZ)

    # ---------- 4. Parse basis set: l values, exponents, coeffs ----------
    # Determine number of AOs and max primitives per atom
    o = []  # number of AOs per atom
    p = []  # max primitives per AO per atom
    for a in range(natoms):
        basis = CZ[a]
        o.append(int(basis[0, 0]))
        p.append(int(basis[0, 1]))

    nmax = max(o)
    mmax = max(p)

    # l values for each AO on each atom
    l_arr = np.zeros((nmax, natoms), dtype=int)
    # count of AOs for each l on each atom
    v = np.zeros((llim + 1, natoms), dtype=int)

    for a in range(natoms):
        basis = CZ[a]
        row = 1
        for i in range(int(basis[0, 0])):
            l_val = int(basis[row, 0])
            l_arr[i, a] = l_val
            v[l_val, a] += 1
            n_prim = int(basis[row, 1])
            row += 1 + n_prim

    # ---------- 5. Exponents / contraction coefficients + primitive norm ----------
    sl = int(v.max())
    c = np.zeros((mmax, sl, llim + 1, natoms), dtype=np.float64)
    ze = np.zeros((mmax, sl, llim + 1, natoms), dtype=np.float64)
    m_arr = np.zeros((llim + 1, sl, natoms), dtype=int)
    v2 = np.zeros((llim + 1, natoms), dtype=int)

    for a in range(natoms):
        basis = CZ[a]
        row = 1
        for i in range(int(basis[0, 0])):
            l_val = int(basis[row, 0])
            n_prim = int(basis[row, 1])
            idx = v2[l_val, a]
            m_arr[l_val, idx, a] = n_prim
            row += 1
            for k in range(n_prim):
                zeta = basis[row, 0]
                coeff = basis[row, 1]
                ze[k, idx, l_val, a] = zeta
                # Primitive normalization (Eq. from Helgaker / Szabo-Ostlund)
                prim_norm = np.sqrt(
                    (2 * zeta / np.pi) ** 0.5
                    * (2 * zeta) ** (l_val + 1)
                    * 2 ** (l_val + 2)
                    / f2p[l_val]
                )
                c[k, idx, l_val, a] = coeff * prim_norm
                row += 1
            v2[l_val, a] += 1

    # ---------- 6. Normalize contracted basis functions ----------
    nm = np.zeros((llim + 1, sl, natoms), dtype=np.float64)
    for a in range(natoms):
        for ell in range(llim + 1):
            for j in range(v2[ell, a]):
                s_overlap = 0.0
                for k in range(m_arr[ell, j, a]):
                    for kd in range(m_arr[ell, j, a]):
                        zsum = ze[k, j, ell, a] + ze[kd, j, ell, a]
                        s_overlap += (
                            c[k, j, ell, a]
                            * c[kd, j, ell, a]
                            * f2p[ell]
                            * np.sqrt(np.pi)
                            / 2 ** (ell + 2)
                            / zsum ** (ell + 1.5)
                        )
                nm[ell, j, a] = 1.0 / np.sqrt(s_overlap) if s_overlap > 0 else 0.0

    # ---------- 7. Build Dyson orbitals ----------
    AOSL = np.zeros((npts, npts, npts), dtype=np.float64)
    AOSR = np.zeros((npts, npts, npts), dtype=np.float64)
    AOC = 0  # AO counter

    for a in range(natoms):
        for ell in range(llim + 1):
            for j in range(v2[ell, a]):
                # Radial part (contracted Gaussian)
                addmat = np.zeros((npts, npts, npts), dtype=np.float64)
                for k in range(m_arr[ell, j, a]):
                    addmat += nm[ell, j, a] * c[k, j, ell, a] * np.exp(
                        -ze[k, j, ell, a] * rmat[..., a]
                    )
                # Loop over magnetic quantum numbers
                for ml in range(-ell, ell + 1):
                    wL_val = wL[AOC]
                    wR_val = wR[AOC]
                    AOC += 1
                    if abs(wL_val) > Dystol or abs(wR_val) > Dystol:
                        ang = S[..., ell, ml + llim, a]
                        AOSL += wL_val * addmat * ang
                        AOSR += wR_val * addmat * ang

    # ---------- 8. Normalize orbitals ----------
    print("\t normalizing and centering")
    HsqL = AOSL**2
    HsqR = AOSR**2
    normL = 1.0 / np.sqrt(simp3D(HsqL, step, step, step))
    normR = 1.0 / np.sqrt(simp3D(HsqR, step, step, step))
    print(f"Dyson Left norm {normL:.4f}, \t Dyson Right norm {normR:.4f}")

    # Centroid of each orbital
    def _centroid(psi_sq: np.ndarray, norm: float) -> Tuple[float, float, float]:
        psi2_normed = psi_sq * norm**2
        xc = round(simp3D(X * psi2_normed, step, step, step), tol)
        yc = round(simp3D(Y * psi2_normed, step, step, step), tol)
        zc = round(simp3D(Z * psi2_normed, step, step, step), tol)
        return xc, yc, zc

    xcl, ycl, zcl = _centroid(HsqL, normL)
    xcr, ycr, zcr = _centroid(HsqR, normR)
    print(f"Initial Left Dyson centroid  x = {xcl:.5f}  y = {ycl:.5f}  z = {zcl:.5f}")
    print(f"Initial Right Dyson centroid x = {xcr:.5f}  y = {ycr:.5f}  z = {zcr:.5f}")

    # ---------- 9. Recenter orbitals via interpolation ----------
    print("Recentering The Left and Right Dyson Orbitals")

    # Build coordinate vectors (assumed uniform)
    x_vec = X[:, 0, 0]
    y_vec = Y[0, :, 0]
    z_vec = Z[0, 0, :]

    interpL = RegularGridInterpolator(
        (x_vec, y_vec, z_vec), AOSL, method="linear", bounds_error=False, fill_value=0.0
    )
    interpR = RegularGridInterpolator(
        (x_vec, y_vec, z_vec), AOSR, method="linear", bounds_error=False, fill_value=0.0
    )

    # Shifted query points
    pts_L = np.stack([X + xcl, Y + ycl, Z + zcl], axis=-1)
    pts_R = np.stack([X + xcr, Y + ycr, Z + zcr], axis=-1)

    orbL = normL * interpL(pts_L)
    orbR = normR * interpR(pts_R)

    return orbL, orbR, llim
