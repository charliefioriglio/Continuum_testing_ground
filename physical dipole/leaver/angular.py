import numpy as np
from scipy.special import factorial
from scipy.linalg import eigh

def analytic(m, L_max, E, a, D):
    m = np.abs(m)
    c = np.sqrt(2 * E * a**2)
    ell_vals = np.arange(m, L_max + 1)
    N = len(ell_vals)
    H = np.zeros((N, N))
    S_diag = np.array([2 * factorial(l + m) / ((2 * l + 1) * factorial(l - m)) for l in ell_vals])
    S = np.diag(S_diag)

    for i, l in enumerate(ell_vals):
        if l - m >= 0:
            f1 = -(l * (l + 1)) * S_diag[i]
            f2 = -c**2 * (
                (l + m) * (l - m) / ((2 * l + 1) * (2 * l - 1)) +
                (l - m + 1) * (l + m + 1) / ((2 * l + 1) * (2 * l + 3))
            ) * S_diag[i]
            H[i, i] += f1 + f2

        if i + 1 < N:
            f = (-2 * D / (2 * l + 1)) * (l + 1 - m)
            f *= 2 * factorial(l + m + 1) / ((2 * l + 3) * factorial(l - m + 1))
            H[i + 1, i] += f
            H[i, i + 1] += f
        if i - 1 >= 0:
            f = (-2 * D / (2 * l + 1)) * (l + m)
            f *= 2 * factorial(l + m - 1) / ((2 * l - 1) * factorial(l - m - 1))
            H[i - 1, i] += f
            H[i, i - 1] += f

        if i + 2 < N:
            f = -c**2 * (l - m + 1) * (l - m + 2)
            f *= 2 * factorial(l + m + 2) / (
                (2 * l + 1) * (2 * l + 3) * (2 * l + 5) * factorial(l - m + 2)
            )
            H[i + 2, i] += f
            H[i, i + 2] += f
        if i - 2 >= 0:
            f = -c**2 * (l + m) * (l + m - 1)
            f *= 2 * factorial(l + m - 2) / (
                (2 * l + 1) * (2 * l - 1) * (2 * l - 3) * factorial(l - m - 2)
            )
            H[i - 2, i] += f
            H[i, i - 2] += f

    eigvals, eigvecs = eigh(H, S)
    eigvals = eigvals[::-1]
    eigvecs = eigvecs[:, ::-1]
    return eigvals, eigvecs, ell_vals