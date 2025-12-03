import numpy as np
from scipy.special import jv

def spherical_bessel_general_order(L, x):
    return np.sqrt(np.pi / (2 * x)) * jv(L + 0.5, x)

def radial_function(L_N, c, xi_grid):
    c_xi = c * xi_grid
    j_L = spherical_bessel_general_order(L_N, c_xi)
    
    # Raw radial function (unnormalized for now)
    R = xi_grid * np.sqrt(c) * j_L

    return R
