"""
parper.py
=========

Precompute P_l^m(±1) values needed for parallel/perpendicular polarization
in the β calculation, translated from MATLAB parper.m.

Functions
---------
parper(lmax)
    Compute P_l^m(cos θ) at θ = 0 and θ = π.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = ["parper"]


def parper(lmax: int) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Compute P_l^m(cos θ) at θ = 0 (parallel) and θ = π (perpendicular).

    At θ = 0, cos θ = +1:
        P_l^0(1) = 1, P_l^m(1) = 0 for m ≠ 0.

    At θ = π, cos θ = -1:
        P_l^0(-1) = (-1)^l, P_l^m(-1) = 0 for m ≠ 0.

    This function returns arrays containing these values plus the
    normalization prefactors used in the real harmonics.

    Parameters
    ----------
    lmax : int
        Maximum angular momentum to compute.

    Returns
    -------
    ppar : ndarray, shape (lmax+1, 2*lmax+3)
        Values at cos θ = +1 (parallel polarization, θ = 0).
        Index as ppar[l, m + lmax + 1].
    pper : ndarray, shape (lmax+1, 2*lmax+3)
        Values at cos θ = -1 (perpendicular polarization, θ = π).
        Index as pper[l, m + lmax + 1].

    Notes
    -----
    MATLAB source (parper.m)::

        function [ppar, pper] = parper(lmax)
        ppar = zeros(lmax+1, 2*lmax+3);
        pper = zeros(lmax+1, 2*lmax+3);
        for l = 0:lmax
            for m = -l:l
                pf = sqrt( (2*l+1) * factorial(l-abs(m)) / 2 / factorial(l+abs(m)) ) ...
                     / sqrt( pi * (1 + (m==0)) );
                if m == 0
                    ppar(l+1, m+lmax+2) = pf;
                    pper(l+1, m+lmax+2) = pf * (-1)^l;
                end
            end
        end

    The physical interpretation:

    * Parallel (θ = 0): Light polarization along the molecular axis.
    * Perpendicular (θ = π): Light polarization perpendicular to the axis.

    Since P_l^m(±1) = 0 for m ≠ 0, only the m = 0 entries are non-zero.
    """
    n_m = 2 * lmax + 3  # m index range: -lmax-1 to lmax+1
    ppar = np.zeros((lmax + 1, n_m), dtype=np.float64)
    pper = np.zeros((lmax + 1, n_m), dtype=np.float64)

    for l in range(lmax + 1):
        # Only m = 0 gives non-zero P_l^m(±1)
        m = 0
        delta = 1  # since m == 0
        pf = math.sqrt(
            (2 * l + 1) * math.factorial(l - abs(m)) / 2 / math.factorial(l + abs(m))
        ) / math.sqrt(math.pi * (1 + delta))

        idx_m = m + lmax + 1  # Python 0-indexed offset
        ppar[l, idx_m] = pf * 1.0           # P_l^0(+1) = 1
        pper[l, idx_m] = pf * ((-1) ** l)   # P_l^0(-1) = (-1)^l

    return ppar, pper
