"""
utilities.py
============

Low-level helper functions translated from the MATLAB point dipole program.

Functions
---------
dfac(n)
    Double factorial  (MATLAB: dfac.m)
tricoeff(a, b, c)
    Triangle coefficient for Wigner symbols  (MATLAB: tricoeff.m)
x(j1, j2, j3, m1, m2)
    Summation helper for Wigner 3-j symbols  (MATLAB: x.m)
"""

from __future__ import annotations

import math
from functools import lru_cache

__all__ = ["dfac", "tricoeff", "x"]


# ---------------------------------------------------------------------------
# dfac.m  – double factorial
# ---------------------------------------------------------------------------
@lru_cache(maxsize=128)
def dfac(n: int) -> int:
    """
    Compute the double factorial n!!.

    Parameters
    ----------
    n : int
        Non-negative integer.

    Returns
    -------
    int
        n!! = n * (n-2) * (n-4) * ... * (1 or 2).

    Notes
    -----
    MATLAB source (dfac.m)::

        function d = dfac(n)
        d = 1;
        for i = n:-2:1
            d = d*i;
        end
    """
    if n < 0:
        raise ValueError("Double factorial not defined for negative integers")
    d = 1
    for i in range(n, 0, -2):
        d *= i
    return d


# ---------------------------------------------------------------------------
# tricoeff.m  – triangle coefficient
# ---------------------------------------------------------------------------
@lru_cache(maxsize=256)
def tricoeff(a: int, b: int, c: int) -> float:
    """
    Compute the triangle coefficient for Wigner 3-j / CG symbols.

    Parameters
    ----------
    a, b, c : int
        Angular momentum quantum numbers satisfying the triangle inequality.

    Returns
    -------
    float
        sqrt( (a+b-c)! (a-b+c)! (-a+b+c)! / (a+b+c+1)! )

    Notes
    -----
    MATLAB source (tricoeff.m)::

        function t = tricoeff(a, b, c)
        t = sqrt( factorial(a+b-c) * factorial(a-b+c) * factorial(-a+b+c) ...
                  / factorial(a+b+c+1) );
    """
    f1 = math.factorial(a + b - c)
    f2 = math.factorial(a - b + c)
    f3 = math.factorial(-a + b + c)
    f4 = math.factorial(a + b + c + 1)
    return math.sqrt(f1 * f2 * f3 / f4)


# ---------------------------------------------------------------------------
# x.m  – summation helper for Wigner 3-j
# ---------------------------------------------------------------------------
@lru_cache(maxsize=512)
def x(j1: int, j2: int, j3: int, m1: int, m2: int) -> float:
    r"""
    Summation helper used in the Wigner 3-j symbol calculation.

    This computes the sum over index *t* that appears in the standard
    Racah formula.

    Parameters
    ----------
    j1, j2, j3 : int
        Angular momentum quantum numbers.
    m1, m2 : int
        Magnetic quantum numbers for j1 and j2.

    Returns
    -------
    float
        The value of

        .. math::
            \sum_t \frac{(-1)^t}{t! (j_3 - j_2 + t + m_1)! (j_3 - j_1 + t - m_2)!
                         (j_1 + j_2 - j_3 - t)! (j_1 - t - m_1)! (j_2 - t + m_2)!}

        summed over all *t* for which every factorial argument is ≥ 0.

    Notes
    -----
    MATLAB source (x.m)::

        function s = x(j1, j2, j3, m1, m2)
        s = 0;
        for t = 0:100
            f1 = j3-j2+t+m1; f2 = j3-j1+t-m2;
            f3 = j1+j2-j3-t; f4 = j1-t-m1; f5 = j2-t+m2;
            if min([f1 f2 f3 f4 f5]) >= 0
                s = s + (-1)^t / ( factorial(t)*factorial(f1)*factorial(f2)...
                                  *factorial(f3)*factorial(f4)*factorial(f5) );
            end
        end
    """
    s = 0.0
    for t in range(101):
        f1 = j3 - j2 + t + m1
        f2 = j3 - j1 + t - m2
        f3 = j1 + j2 - j3 - t
        f4 = j1 - t - m1
        f5 = j2 - t + m2
        if min(f1, f2, f3, f4, f5) >= 0:
            denom = (
                math.factorial(t)
                * math.factorial(f1)
                * math.factorial(f2)
                * math.factorial(f3)
                * math.factorial(f4)
                * math.factorial(f5)
            )
            s += ((-1) ** t) / denom
    return s
