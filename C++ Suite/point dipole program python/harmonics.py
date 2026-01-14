"""
harmonics.py
============

Spherical harmonic functions and angular coordinate utilities translated
from the MATLAB point dipole program.

Functions
---------
prefactor(l, m)
    Normalization prefactor for real (cubic) harmonics  (MATLAB: prefactor.m)
thph(X, Y, Z)
    Convert Cartesian grid to (r, theta, phi)  (MATLAB: thph.m)
Ylmc(l, m, theta, phi)
    Complex spherical harmonics via Legendre generating function (MATLAB: Ylmc.m)
Ylm(l, m, X, Y, Z)
    Real (cubic) spherical harmonics on a Cartesian grid  (MATLAB: Ylm.m)
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

__all__ = ["prefactor", "thph", "Ylmc", "Ylm"]


# ---------------------------------------------------------------------------
# prefactor.m  – normalization prefactor for real harmonics
# ---------------------------------------------------------------------------
@lru_cache(maxsize=256)
def prefactor(l: int, m: int) -> float:
    r"""
    Compute the normalization prefactor for real (cubic) spherical harmonics.

    Parameters
    ----------
    l : int
        Angular momentum quantum number.
    m : int
        Magnetic quantum number (non-negative, since m<0 is handled via
        the sine/cosine angular factor).

    Returns
    -------
    float
        :math:`\sqrt{\frac{(2l+1)(l-m)!}{2(l+m)!}} / \sqrt{\pi (1 + \delta_{m,0})}`

    Notes
    -----
    MATLAB source (prefactor.m)::

        function p = prefactor(l, m)
        p = sqrt( (2*l+1) * factorial(l-m) / 2 / factorial(l+m) ) ...
            / sqrt( pi * (1 + (m==0)) );
    """
    delta = 1 if m == 0 else 0
    return math.sqrt(
        (2 * l + 1) * math.factorial(l - m) / 2 / math.factorial(l + m)
    ) / math.sqrt(math.pi * (1 + delta))


# ---------------------------------------------------------------------------
# thph.m  – Cartesian to spherical coordinate conversion
# ---------------------------------------------------------------------------
def thph(
    X: np.ndarray, Y: np.ndarray, Z: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert Cartesian coordinates to spherical (r, theta, phi).

    Parameters
    ----------
    X, Y, Z : ndarray
        Cartesian coordinate arrays (same shape).

    Returns
    -------
    r : ndarray
        Radial coordinate sqrt(X² + Y² + Z²).
    theta : ndarray
        Polar angle arccos(Z/r), in [0, π].
    phi : ndarray
        Azimuthal angle atan2(Y, X), in [-π, π].

    Notes
    -----
    MATLAB source (thph.m)::

        function [th, ph] = thph(X, Y, Z)
        r = sqrt(X.^2 + Y.^2 + Z.^2);
        r(r == 0) = eps;   % avoid division by zero
        th = acos(Z ./ r);
        ph = atan2(Y, X);
    """
    r = np.sqrt(X**2 + Y**2 + Z**2)
    r_safe = np.where(r == 0, np.finfo(float).eps, r)
    theta = np.arccos(Z / r_safe)
    phi = np.arctan2(Y, X)
    return r, theta, phi


# ---------------------------------------------------------------------------
# Ylmc.m  – complex spherical harmonics
# ---------------------------------------------------------------------------
def Ylmc(l: int, m: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    r"""
    Compute complex spherical harmonics via Legendre generating function.

    Parameters
    ----------
    l : int
        Degree of the harmonic.
    m : int
        Order of the harmonic (can be negative).
    theta : ndarray
        Polar angle array.
    phi : ndarray
        Azimuthal angle array.

    Returns
    -------
    ndarray (complex)
        :math:`Y_l^m(\theta, \phi)` with Condon-Shortley phase convention.

    Notes
    -----
    Uses the recursion relation for associated Legendre polynomials:

    .. math::
        P_l^l(\cos\theta) &= -(2l-1) \sin\theta \; P_{l-1}^{l-1}(\cos\theta) \\
        P_l^m(\cos\theta) &= \frac{\cos\theta (2l-1) P_{l-1}^m
                                   - (l+m-1) P_{l-2}^m}{l - m}

    MATLAB source (Ylmc.m)::

        function Y = Ylmc(l, m, theta, phi)
        ...
        % Uses generating function for P_l^m
    """
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    abs_m = abs(m)

    # Build P_{|m|}^{|m|} (diagonal recursion)
    if abs_m == 0:
        P_prev = np.ones_like(theta, dtype=np.float64)
    else:
        P_prev = np.ones_like(theta, dtype=np.float64)
        for k in range(1, abs_m + 1):
            P_prev = -(2 * k - 1) * sin_theta * P_prev

    # Now raise l from |m| to l using standard recursion
    if l == abs_m:
        P_lm = P_prev
    else:
        P_prev2 = np.zeros_like(theta, dtype=np.float64)  # P_{|m|-1}^{|m|} = 0
        P_prev1 = P_prev.copy()  # P_{|m|}^{|m|}
        for ell in range(abs_m + 1, l + 1):
            P_new = (
                cos_theta * (2 * ell - 1) * P_prev1 - (ell + abs_m - 1) * P_prev2
            ) / (ell - abs_m)
            P_prev2 = P_prev1
            P_prev1 = P_new
        P_lm = P_prev1

    # Normalization factor with Condon-Shortley phase
    norm = math.sqrt(
        (2 * l + 1) / (4 * math.pi)
        * math.factorial(l - abs_m) / math.factorial(l + abs_m)
    )
    # Condon-Shortley phase is built into the diagonal recursion above
    if m < 0:
        # Y_l^{-|m|} = (-1)^m (Y_l^{|m|})^*
        phase = (-1) ** abs_m
        Y = phase * norm * P_lm * np.exp(-1j * abs_m * phi)
    else:
        Y = norm * P_lm * np.exp(1j * m * phi)

    return Y


# ---------------------------------------------------------------------------
# Ylm.m  – real (cubic) spherical harmonics
# ---------------------------------------------------------------------------
def Ylm(l: int, m: int, X: np.ndarray, Y: np.ndarray, Z: np.ndarray) -> np.ndarray:
    r"""
    Compute real (cubic) spherical harmonics on a Cartesian grid.

    The real harmonics are defined as combinations of complex harmonics:

    .. math::
        S_l^m = \begin{cases}
            \frac{1}{\sqrt{2}} (Y_l^{-|m|} + (-1)^m Y_l^{|m|}) & m > 0 \text{ (cosine)} \\
            Y_l^0 & m = 0 \\
            \frac{i}{\sqrt{2}} (Y_l^{-|m|} - (-1)^m Y_l^{|m|}) & m < 0 \text{ (sine)}
        \end{cases}

    The MATLAB code uses an explicit generating function approach.

    Parameters
    ----------
    l : int
        Degree of the harmonic.
    m : int
        Order of the harmonic (any sign).
    X, Y, Z : ndarray
        Cartesian coordinate arrays.

    Returns
    -------
    ndarray (real)
        Real spherical harmonic evaluated on the grid.

    Notes
    -----
    MATLAB source (Ylm.m) is ~175 lines implementing an explicit generating
    function for the associated Legendre polynomial.  We translate the core
    algorithm here.
    """
    _, theta, phi_ang = thph(X, Y, Z)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    abs_m = abs(m)

    # ---------- Build P_l^{|m|}(cos θ) via generating function ----------
    if abs_m == 0:
        P_mm = np.ones_like(theta)
    else:
        P_mm = np.ones_like(theta)
        for k in range(1, abs_m + 1):
            P_mm = -(2 * k - 1) * sin_theta * P_mm

    if l == abs_m:
        P_lm = P_mm
    else:
        P_prev2 = np.zeros_like(theta)
        P_prev1 = P_mm.copy()
        for ell in range(abs_m + 1, l + 1):
            P_new = (
                cos_theta * (2 * ell - 1) * P_prev1 - (ell + abs_m - 1) * P_prev2
            ) / (ell - abs_m)
            P_prev2 = P_prev1
            P_prev1 = P_new
        P_lm = P_prev1

    # ---------- Normalization ----------
    pf = prefactor(l, abs_m)

    # ---------- Angular factor (sine/cosine) ----------
    if m > 0:
        angular = np.cos(m * phi_ang)
    elif m < 0:
        angular = np.sin(abs_m * phi_ang)
    else:
        angular = np.ones_like(phi_ang)

    return pf * P_lm * angular
