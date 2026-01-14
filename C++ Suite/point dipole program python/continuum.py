"""
continuum.py
============

Continuum (photoelectron) wavefunction construction translated from
MATLAB continuum.m.

Functions
---------
continuum(N, lam, lmax, eign, eigv, q, r, S, npts)
    Build the point-dipole continuum wavefunction on a 3-D grid.
build_S_grid(lmax, X, Y, Z)
    Pre-compute the real harmonic grid S[l, λ] for all l, λ.
"""

from __future__ import annotations

import numpy as np
from scipy.special import jv  # Bessel J_nu

from .harmonics import Ylm

__all__ = ["continuum", "build_S_grid"]


def build_S_grid(
    lmax: int,
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
) -> np.ndarray:
    r"""
    Pre-compute the real spherical harmonic grids S_{l,m}(X, Y, Z).

    Parameters
    ----------
    lmax : int
        Maximum angular momentum.
    X, Y, Z : ndarray, shape (npts, npts, npts)
        Cartesian coordinate meshes.

    Returns
    -------
    S : ndarray, shape (npts, npts, npts, lmax+1, 2*lmax+1)
        Indexed as S[..., l, m + lmax].

    Notes
    -----
    MATLAB stores S as S(:,:,:, l+1, lam+lmax+1) with λ (lam) spanning
    -lmax to +lmax.  We follow the same convention in Python (0-indexed).
    """
    npts = X.shape[0]
    n_m = 2 * lmax + 1
    S = np.zeros((npts, npts, npts, lmax + 1, n_m), dtype=np.float64)

    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            S[..., l, m + lmax] = Ylm(l, m, X, Y, Z)

    return S


def continuum(
    N: int,
    lam: int,
    lmax: int,
    eign: np.ndarray,
    eigv: np.ndarray,
    q: float,
    r: np.ndarray,
    S: np.ndarray,
    npts: int,
) -> np.ndarray:
    r"""
    Compute the point-dipole continuum wavefunction on a 3-D grid.

    Parameters
    ----------
    N : int
        Continuum mode index (column of eigenvector matrix).
    lam : int
        Projection quantum number λ.
    lmax : int
        Maximum orbital angular momentum in the expansion.
    eign : ndarray, shape (lmax+1, 2*lmax+1)
        Eigenvalue array indexed as eign[N, lam + lmax].
    eigv : ndarray, shape (lmax+1, lmax+1, 2*lmax+1)
        Eigenvector array indexed as eigv[l, N, lam + lmax].
    q : float
        Photoelectron momentum (wavenumber) in atomic units.
    r : ndarray, shape (npts, npts, npts)
        Radial grid sqrt(X² + Y² + Z²).
    S : ndarray, shape (npts, npts, npts, lmax+1, 2*lmax+1)
        Pre-computed real harmonic grids (from `build_S_grid`).
    npts : int
        Number of grid points per dimension.

    Returns
    -------
    CFN : ndarray, shape (npts, npts, npts)
        The continuum wavefunction on the grid.

    Notes
    -----
    MATLAB source (continuum.m)::

        Omega = 0;
        for l = abs(lam):lmax
            if abs(eigv(l+1, N+1, lmax+lam+1)) > 1E-5
                Omega = Omega + eigv(l+1, N+1, lmax+lam+1) * S(:,:,:, l+1, lmax+lam+1);
            end
        end
        CFN = sqrt(r) * sqrt(pi/(2*q)) * besselj(eign(N+1,lam+lmax+1)+0.5, q*r) * Omega;
        CFN(cen,cen,cen) = 0;

    The wavefunction normalizes as described in Gallup (PRA 23, 632, 1981):

    .. math::
        \psi_{N,\lambda}(\mathbf{r}) = \sqrt{r}\,\sqrt{\frac{\pi}{2 q}}\,
            J_{\nu_N + 1/2}(q r)\,\Omega_{N,\lambda}(\theta, \phi)

    where Ω is the angular part expanded in real harmonics.
    """
    # Build angular part Ω
    Omega = np.zeros((npts, npts, npts), dtype=np.float64)
    lam_idx = lam + lmax  # index into eign / eigv / S

    for l in range(abs(lam), lmax + 1):
        coeff = eigv[l, N, lam_idx]
        if abs(coeff) > 1e-5:
            Omega += coeff * S[..., l, lam_idx]

    # Radial part: sqrt(r) * sqrt(π / 2q) * J_{ν + 0.5}(qr)
    nu = eign[N, lam_idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        radial = np.sqrt(r) * np.sqrt(np.pi / (2 * q)) * jv(nu + 0.5, q * r)

    CFN = radial * Omega

    # Zero out the origin to avoid singularity artefacts
    cen = (npts - 1) // 2
    CFN[cen, cen, cen] = 0.0

    return CFN
