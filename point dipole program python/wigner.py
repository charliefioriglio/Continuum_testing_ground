"""
wigner.py
=========

Wigner 3-j symbol and Clebsch-Gordan coefficient routines translated
from the MATLAB point dipole program.

Functions
---------
wigner3j(j1, j2, j3, m1, m2, m3)
    Wigner 3-j symbol  (MATLAB: wigner3j.m)
clebschgordan(j1, j2, j3, m1, m2, m3)
    Clebsch-Gordan coefficient  (MATLAB: clebschgordan.m)
"""

from __future__ import annotations

import math
from functools import lru_cache

from .utilities import tricoeff, x

__all__ = ["wigner3j", "clebschgordan"]


# ---------------------------------------------------------------------------
# wigner3j.m  – Wigner 3-j symbol
# ---------------------------------------------------------------------------
@lru_cache(maxsize=4096)
def wigner3j(
    j1: int, j2: int, j3: int, m1: int, m2: int, m3: int
) -> float:
    r"""
    Compute the Wigner 3-j symbol.

    Parameters
    ----------
    j1, j2, j3 : int
        Angular momentum quantum numbers (non-negative integers).
    m1, m2, m3 : int
        Magnetic quantum numbers satisfying |mi| ≤ ji.

    Returns
    -------
    float
        The Wigner 3-j symbol

        .. math::
            \begin{pmatrix}
              j_1 & j_2 & j_3 \\
              m_1 & m_2 & m_3
            \end{pmatrix}

        or 0 if the selection rules are violated.

    Notes
    -----
    MATLAB source (wigner3j.m)::

        function w = wigner3j(j1, j2, j3, m1, m2, m3)
        if abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3
            w = 0;
        elseif m1 + m2 + m3 ~= 0
            w = 0;
        elseif j3 > j1+j2 || j3 < abs(j1-j2)
            w = 0;
        else
            w = (-1)^(j1-j2-m3) * tricoeff(j1,j2,j3) ...
                * sqrt( factorial(j1+m1)*factorial(j1-m1)...
                       *factorial(j2+m2)*factorial(j2-m2)...
                       *factorial(j3+m3)*factorial(j3-m3) ) ...
                * x(j1, j2, j3, m1, m2);
        end
    """
    # Selection rules
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if m1 + m2 + m3 != 0:
        return 0.0
    if j3 > j1 + j2 or j3 < abs(j1 - j2):
        return 0.0

    sign = (-1) ** (j1 - j2 - m3)
    tri = tricoeff(j1, j2, j3)
    fac_part = math.sqrt(
        math.factorial(j1 + m1)
        * math.factorial(j1 - m1)
        * math.factorial(j2 + m2)
        * math.factorial(j2 - m2)
        * math.factorial(j3 + m3)
        * math.factorial(j3 - m3)
    )
    sum_part = x(j1, j2, j3, m1, m2)

    return sign * tri * fac_part * sum_part


# ---------------------------------------------------------------------------
# clebschgordan.m  – Clebsch-Gordan coefficient
# ---------------------------------------------------------------------------
@lru_cache(maxsize=4096)
def clebschgordan(
    j1: int, j2: int, j3: int, m1: int, m2: int, m3: int
) -> float:
    r"""
    Compute the Clebsch-Gordan coefficient ⟨j1 m1; j2 m2 | j3 m3⟩.

    Parameters
    ----------
    j1, j2, j3 : int
        Angular momentum quantum numbers.
    m1, m2, m3 : int
        Magnetic quantum numbers.

    Returns
    -------
    float
        The Clebsch-Gordan coefficient (real-valued by convention).

    Notes
    -----
    MATLAB source (clebschgordan.m)::

        function c = clebschgordan(j1, j2, j3, m1, m2, m3)
        c = (-1)^(j1-j2+m3) * sqrt(2*j3+1) * wigner3j(j1, j2, j3, m1, m2, -m3);

    Relation to Wigner 3-j::

        ⟨j1 m1; j2 m2 | j3 m3⟩ = (-1)^{j1-j2+m3} * sqrt(2*j3+1)
                                * ( j1  j2   j3  )
                                  ( m1  m2  -m3  )
    """
    sign = (-1) ** (j1 - j2 + m3)
    return sign * math.sqrt(2 * j3 + 1) * wigner3j(j1, j2, j3, m1, m2, -m3)
