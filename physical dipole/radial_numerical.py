import numpy as np
from scipy.integrate import solve_ivp

# For D = 1.7 au
alpha = 0.95234
b = 0.03638
cat = 1.03991 
d = -0.46058
def compute_series_initial_conditions(c, lam, N, xi0):
    a = np.zeros(N)
    a[0] = 1 / (alpha + b * np.cos(2 * cat * np.log(c) + d))
    a[1] = -0.5 * (c**2 - lam) * a[0]
    for k in range(2, N):
        a[k] = -1 / (2 * k**2) * (
            (c**2 - lam + k * (k - 1)) * a[k-1] +
            2 * c**2 * a[k-2] +
            (c**2 * a[k-3] if k >= 3 else 0)
        )
    dx = xi0 - 1
    S = sum(a[k] * dx**k for k in range(N))
    dS = sum((k) * a[k] * dx**(k - 1) for k in range(1, N))
    return S, dS

def radial_rhs(xi, y, c, lam, m):
    S, dS = y
    denom = xi**2 - 1
    d2S = (lam * S - c**2 * xi**2 * S + (m**2 / denom) * S -
           (2 * xi * dS)) / denom
    return [dS, d2S]

def solve_radial(c, lam, xi_start, xi_end, N_terms, m):
    S0, dS0 = compute_series_initial_conditions(c, lam, N=N_terms, xi0=xi_start)
    sol = solve_ivp(
        lambda xi, y: radial_rhs(xi, y, c, lam, m),
        [xi_start, xi_end],
        [S0, dS0],
        method='LSODA',
        rtol=1e-10,
        atol=1e-12,
        dense_output=True
    )
    return sol

def compute_radial_function(m, n, lam_mn, c, xi_vals, N_terms=50,
                            xi_start=1+1e-10):
    sol = solve_radial(c, lam_mn, xi_start=xi_start,
                       xi_end=xi_vals[-1], N_terms=N_terms, m=m)
    S_vals = sol.sol(xi_vals)[0]

    return S_vals
