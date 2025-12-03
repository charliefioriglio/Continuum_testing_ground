"""
clebsch_array.py
================

Pre-computation of Clebsch-Gordan coefficient arrays for efficient β
evaluation, translated from MATLAB clebarray1.m.

Functions
---------
clebarray1(l, el, lmax)
    Build 3-D array of CG coefficients indexed by (m2, m1, m3).
"""

from __future__ import annotations

import numpy as np

from .wigner import clebschgordan

__all__ = ["clebarray1"]


def clebarray1(l: int, el: int, lmax: int) -> np.ndarray:
    r"""
    Pre-compute Clebsch-Gordan array for a given (l, L) pair.

    Stores the coefficients ⟨l m₂; 1 m₁ | L m₃⟩ in an array for rapid
    lookup during the β calculation.

    Parameters
    ----------
    l : int
        First angular momentum quantum number (orbital angular momentum).
    el : int
        Second angular momentum quantum number (total angular momentum L).
    lmax : int
        Maximum angular momentum to allocate array indexing.

    Returns
    -------
    ndarray, shape (3, 2*lmax+1, 2*lmax+3)
        The CG coefficient array CG[m2+2, m1+lmax+1, m3+lmax+2] where
        * m2 ∈ {-1, 0, +1}      →  index 1, 2, 3  (MATLAB) / 0, 1, 2 (Python)
        * m1 ∈ {-lmax, …, lmax} →  index 1..2*lmax+1 (MATLAB) / 0..2*lmax (Python)
        * m3 ∈ {-lmax-1, …, lmax+1} → index 1..2*lmax+3 (MATLAB) / 0..2*lmax+2 (Python)

    Notes
    -----
    MATLAB source (clebarray1.m)::

        function CG = clebarray1(l, el, lmax)
        CG = zeros(l+1, el+1, 3, 2*lmax+1, 2*lmax+3);
        for m2 = -1:1
            for m1 = -lmax:lmax
                for m3 = -lmax-1:lmax+1
                    CG(l+1, el+1, m2+2, m1+lmax+1, m3+lmax+2) = ...
                        clebschgordan(l, 1, el, m1, m2, m3);
                end
            end
        end

    The MATLAB version allocates a 5-D array.  Here we return a 3-D slice
    for a specific (l, L) pair since Python indexing is more natural and
    we can call this function inside a loop.
    """
    # Array dimensions
    n_m2 = 3                 # m2 ∈ {-1, 0, 1}
    n_m1 = 2 * lmax + 1      # m1 ∈ {-lmax, ..., lmax}
    n_m3 = 2 * lmax + 3      # m3 ∈ {-lmax-1, ..., lmax+1}

    CG = np.zeros((n_m2, n_m1, n_m3), dtype=np.float64)

    for m2 in range(-1, 2):        # -1, 0, 1
        idx_m2 = m2 + 1            # 0, 1, 2
        for m1 in range(-lmax, lmax + 1):
            idx_m1 = m1 + lmax     # 0 .. 2*lmax
            for m3 in range(-lmax - 1, lmax + 2):
                idx_m3 = m3 + lmax + 1  # 0 .. 2*lmax+2
                CG[idx_m2, idx_m1, idx_m3] = clebschgordan(l, 1, el, m1, m2, m3)

    return CG
