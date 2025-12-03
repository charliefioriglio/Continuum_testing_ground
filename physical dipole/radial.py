import numpy as np
from scipy.integrate import solve_ivp

def compute_bj_coefficients(m, J_max):
    """Compute the coefficients in the series expansion of -m^2/(xi^2 - 1) about xi=1."""
    b = np.zeros(J_max)
    for j in range(J_max):
        b[j] = -m**2 * (-1)**(j+1) * (j + 1)
    return b

def compute_general_series_initial_conditions(c, lam, m, N, xi0, J_max):
    """Compute initial S(xi0), dS(xi0) using generalized recursion relation for all m."""
    a = np.zeros(N)
    a[0] = 1.0
    b = compute_bj_coefficients(m, J_max)

    for k in range(1, N):
        sum_b = 0.0
        for j in range(J_max):
            idx = k - (j + 1)
            if idx >= 0:
                sum_b += a[idx] * b[j]

        term_1 = 2 * k**2 * a[k-1] if k >= 1 else 0
        term_2 = (k * (k - 1) + c**2 - lam) * a[k-1] if k >= 1 else 0
        term_3 = 2 * c**2 * a[k-2] if k >= 2 else 0
        term_4 = c**2 * a[k-3] if k >= 3 else 0

        a[k] = - (term_2 + term_3 + term_4 + sum_b) / (2 * k**2)

    dx = xi0 - 1
    S = sum(a[k] * dx**k for k in range(N))
    dS = sum(k * a[k] * dx**(k - 1) for k in range(1, N))
    return S, dS

def radial_rhs_full(xi, y, c, lam, m):
    S, dS = y
    denom = xi**2 - 1
    d2S = (lam * S - c**2 * xi**2 * S + (m**2 / denom) * S -
           (2 * xi * dS)) / denom
    return [dS, d2S]

def solve_radial_full(c, lam, xi_start, xi_end, N_terms, m, J_max):
    S0, dS0 = compute_general_series_initial_conditions(
        c, lam, m, N_terms, xi_start, J_max=J_max
    )
    sol = solve_ivp(
        lambda xi, y: radial_rhs_full(xi, y, c, lam, m),
        [xi_start, xi_end],
        [S0, dS0],
        method='LSODA',
        rtol=1e-10,
        atol=1e-12,
        dense_output=True
    )
    return sol

def compute_radial_function_full(m, n, lam_mn, c, xi_vals, N_terms=50,
                                 xi_start=1+1e-10, J_max=20):
    sol = solve_radial_full(c, lam_mn, xi_start=xi_start,
                            xi_end=xi_vals[-1], N_terms=N_terms,
                            m=m, J_max=J_max)
    S_vals = sol.sol(xi_vals)[0]
    return S_vals


