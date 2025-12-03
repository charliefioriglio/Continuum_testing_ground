import numpy as np
from scipy.special import spherical_jn, factorial
from scipy.linalg import eig, null_space

def build_radial_coefficient_matrix(m, lam_mn, c, R_max, parity):
    m_abs = abs(m)
    if parity == 'even':
        r_vals = np.arange(0, 2 * R_max, 2)
    elif parity == 'odd':
        r_vals = np.arange(1, 2 * R_max, 2)
    else:
        raise ValueError("parity must be 'even' or 'odd'")

    N = len(r_vals)
    A = np.zeros((N, N))

    for i, r in enumerate(r_vals):
        l = m_abs + r

        # Diagonal term
        A[i, i] = l * (l + 1) - lam_mn + ((2 * l * (l + 1) - 2 * m_abs**2 - 1) * c**2) / ((2 * l - 1) * (2 * l + 3))

        # Off-diagonal terms (r ↔ r ± 2)
        if i + 1 < N:
            A[i, i + 1] = ((2 * m_abs + r + 2) * (2 * m_abs + r + 1) * c**2) / (
                (2 * m_abs + 2 * r + 3) * (2 * m_abs + 2 * r + 5)
            )

        if i - 1 >= 0:
            A[i, i - 1] = (r * (r - 1) * c**2) / (
                (2 * m_abs + 2 * r - 3) * (2 * m_abs + 2 * r - 1)
            )

    return A, r_vals

def compute_radial_function(m, n, lam_mn, c, R_max, xi_vals):
    """
    Compute the radial prolate spheroidal function S_{mn}(ξ)
    by solving the eigenvalue problem from the recurrence relation.
    """
    m_abs = abs(m)
    parity = 'even' if abs(n - m_abs) % 2 == 0 else 'odd'

    A, r_vals = build_radial_coefficient_matrix(m, lam_mn, c, R_max, parity)

    # Solve eigenproblem
    eigvals, eigvecs = eig(A)
    d_coeffs = eigvecs[:, n]

    # Normalize (Flammer-style)
    norm_sum = np.sum([
        d * factorial(2 * m_abs + r) / factorial(r)
        for d, r in zip(d_coeffs, r_vals)
    ])
    norm = 1.0 / norm_sum

    # Build radial function S_{mn}(ξ)
    prefactor = ((xi_vals**2 - 1) / xi_vals**2)**(m_abs / 2)
    sum_series = np.zeros_like(xi_vals, dtype=complex)

    for d, r in zip(d_coeffs, r_vals):
        l = m_abs + r
        coeff = factorial(2 * m_abs + r) / factorial(r)
        j_l = spherical_jn(l, c * xi_vals)
        sum_series += d * coeff * j_l

    S_mn_xi = norm * prefactor * sum_series
    return S_mn_xi