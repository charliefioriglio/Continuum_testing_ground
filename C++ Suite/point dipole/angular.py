import numpy as np
from sympy.physics.wigner import wigner_3j
from scipy.special import sph_harm

def build_dipole_angular_matrix(D, lam, l_max):
    l_min = max(abs(lam), 0)
    l_vals = np.arange(l_min, l_max + 1)
    n_l = len(l_vals)
    H = np.zeros((n_l, n_l), dtype=np.float64)

    wigner_cache = {}
    for i, l in enumerate(l_vals):
        for dl in [-1, 1]:
            lp = l + dl
            if lp < l_min or lp > l_max:
                continue
            key = (lp, l)
            if key not in wigner_cache:
                tj1 = float(wigner_3j(lp, 1, l, 0, 0, 0).evalf())
                tj2 = float(wigner_3j(lp, 1, l, -lam, 0, lam).evalf())
                wigner_cache[key] = tj1 * tj2

    for i, l in enumerate(l_vals):
        H[i, i] = l * (l + 1)
        for dl in [-1, 1]:
            lp = l + dl
            if lp < l_min or lp > l_max:
                continue
            j = np.where(l_vals == lp)[0][0]
            coeff = -2 * D * np.sqrt((2 * l + 1) * (2 * lp + 1))
            H[i, j] += coeff * wigner_cache[(lp, l)]
            H[j, i] = H[i, j]

    return H, l_vals

def build_angular_functions(D, lam, l_max, theta_grid, phi_grid):
    H, l_vals = build_dipole_angular_matrix(D, lam, l_max)
    eigvals, eigvecs = np.linalg.eigh(H)
    Y_stack = np.array([sph_harm(lam, l, phi_grid, theta_grid) for l in l_vals])
    omega_funcs = np.tensordot(eigvecs.T, Y_stack, axes=([1], [0]))

    return omega_funcs, eigvals

