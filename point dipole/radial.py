import numpy as np
from scipy.special import jv

def spherical_bessel_general_order(L, kr):
    kr_safe = np.maximum(kr, 1e-6)  # avoid kr → 0
    return np.sqrt(np.pi / (2 * kr_safe)) * jv(L + 0.5, kr_safe)

def radial_function(L_N, k, r_grid):
    kr = k * r_grid
    j_L = spherical_bessel_general_order(L_N, kr)
    R = j_L #* r_grid * np.sqrt(k)
    return R
