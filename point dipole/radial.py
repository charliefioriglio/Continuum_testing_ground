"""Radial helpers for the point-dipole continuum."""

from __future__ import annotations

import numpy as np
from scipy.special import jv

_RADIAL_EPS = 1.0e-10


def _safe_argument(values: np.ndarray) -> np.ndarray:
    """Return a version of *values* with a small floor to avoid division by zero."""

    return np.where(values > _RADIAL_EPS, values, _RADIAL_EPS)


def radial_function(order: float, k_mag: float, r_grid: np.ndarray) -> np.ndarray:
    """Return the Gallup-style radial factor for a single continuum mode.

    The legacy ``continuum.m`` implementation constructs the radial profile as
    ``sqrt(r) * sqrt(pi/(2*k)) * J_{order+1/2}(k*r)`` before contracting with
    the angular fields. Reproducing that expression verbatim keeps the Python
    port consistent with the benchmark MATLAB workflow while remaining finite
    at the nucleus.
    """

    if not np.isfinite(k_mag) or k_mag <= 0.0:
        return np.zeros_like(r_grid, dtype=np.complex128)

    r_safe = _safe_argument(np.asarray(r_grid, dtype=float))
    kr = k_mag * r_safe
    prefactor = np.sqrt(np.pi / (2.0 * k_mag))
    bessel = jv(order + 0.5, kr)
    radial = prefactor * np.sqrt(r_safe) * bessel
    result = np.array(radial, dtype=np.complex128)
    result = np.where(r_grid > _RADIAL_EPS, result, 0.0)
    return result
