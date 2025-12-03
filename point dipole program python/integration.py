"""
integration.py
==============

3-D Simpson's rule integration translated from MATLAB simp3D.m.

Functions
---------
simp3D(f, hx, hy, hz)
    Integrate a 3-D array using Simpson's rule with uniform spacing.
"""

from __future__ import annotations

import numpy as np

__all__ = ["simp3D"]


def simp3D(f: np.ndarray, hx: float, hy: float, hz: float) -> float:
    r"""
    Integrate a 3-D function using Simpson's rule on a uniform grid.

    Parameters
    ----------
    f : ndarray, shape (nx, ny, nz)
        Function values on a 3-D Cartesian grid.
        Each dimension must have an *odd* number of points (n = 2k+1)
        to apply Simpson's rule correctly.
    hx, hy, hz : float
        Grid spacing in x, y, z directions.

    Returns
    -------
    float
        The approximate integral ∫∫∫ f dx dy dz.

    Notes
    -----
    Simpson's rule in 1-D with n = 2k+1 points:

    .. math::
        \int f \,dx \approx \frac{h}{3}
          \left[ f_0 + 4\sum_{i \text{ odd}} f_i + 2\sum_{i \text{ even}>0}^{n-2} f_i + f_{n-1} \right]

    We apply this factorized over each dimension.

    MATLAB source (simp3D.m)::

        function I = simp3D(f, hx, hy, hz)
        % Build 1D Simpson weight vectors
        nx = size(f,1); ny = size(f,2); nz = size(f,3);
        wx = simpweights(nx) * hx/3;
        wy = simpweights(ny) * hy/3;
        wz = simpweights(nz) * hz/3;
        % Contract over each dimension
        I = 0;
        for i=1:nx
            for j=1:ny
                for k=1:nz
                    I = I + wx(i)*wy(j)*wz(k)*f(i,j,k);
                end
            end
        end

    With Python broadcasting we can achieve this more efficiently.
    """
    nx, ny, nz = f.shape

    # Build Simpson weight vectors (1, 4, 2, 4, 2, ..., 4, 1)
    wx = _simpson_weights(nx) * (hx / 3.0)
    wy = _simpson_weights(ny) * (hy / 3.0)
    wz = _simpson_weights(nz) * (hz / 3.0)

    # Contract using outer product of weights
    # W[i,j,k] = wx[i] * wy[j] * wz[k]
    W = np.einsum("i,j,k->ijk", wx, wy, wz)
    return float(np.sum(W * f))


def _simpson_weights(n: int) -> np.ndarray:
    """Return Simpson's rule weight vector of length n (assumed odd)."""
    if n < 3 or n % 2 == 0:
        raise ValueError(f"Simpson's rule requires odd n ≥ 3, got {n}")
    w = np.ones(n, dtype=np.float64)
    w[1:-1:2] = 4.0  # odd indices
    w[2:-1:2] = 2.0  # even indices > 0
    return w
