import numpy as np
from angular import analytic
from radial_point import radial_function
from scipy.special import lpmv
def build_single_mode_wavefunction_on_xyz_grid(m, n, E, a, D, X, Y, Z, L_max=40.0, xi_vals=None, eta_vals=None, phi_vals=None):

    # Convert Cartesian (X,Y,Z) → prolate spheroidal (ξ, η, φ)
    Z_A = a
    Z_B = -a
    rA = np.sqrt((X)**2 + (Y)**2 + (Z - Z_A)**2)
    rB = np.sqrt((X)**2 + (Y)**2 + (Z - Z_B)**2)

    xi = (rA + rB) / (2 * a)
    eta = (rA - rB) / (2 * a)
    phi = np.arctan2(Y, X)

    # Precompute 1D coordinate grids if not provided
    xi_max = np.max(xi) * 1.05
    if xi_vals is None:
        xi_vals = np.linspace(1.0, xi_max, 200)
    if eta_vals is None:
        eta_vals = np.linspace(-1.0, 1.0, 200)
    if phi_vals is None:
        phi_vals = np.linspace(-np.pi, np.pi, 100)

    # Solve angular and radial equations
    m_abs = abs(m)
    c = np.sqrt(2 * E * a**2)

    eigvals_ang, eigvecs_ang, ell_vals = analytic(m_abs, L_max, E, a, D)
    v_ang = eigvecs_ang[:, n]

    # Angular function T_{m,n}(η)
    P_basis = np.array([lpmv(m_abs, l, eta_vals) for l in ell_vals])
    T_eta_vals = P_basis.T @ v_ang
    if m < 0:
        T_eta_vals = T_eta_vals * (-1)**m_abs

    # Interpolate T_eta onto grid η
    from scipy.interpolate import interp1d
    T_eta_func = interp1d(eta_vals, T_eta_vals, kind='cubic')
    T_eta = T_eta_func(eta)

    # === Compute λ_n and check if valid (λ < 0)
    lam_n = -eigvals_ang[n]
    if lam_n < -0.25 or np.isnan(lam_n):
        return None

    # === Compute effective angular momentum N
    N_eff = (-1 + np.sqrt(1 + 4 * lam_n)) / 2

    # === Build radial function R(ξ) = j_N(c ξ) on 1D grid
    R_xi_vals = radial_function(N_eff, c, xi_vals)

    # Interpolate onto full 3D grid
    R_xi_func = interp1d(xi_vals, R_xi_vals, kind='cubic', bounds_error=False, fill_value=0.0)
    R_xi = R_xi_func(xi)

    # === Azimuthal function Φ(φ)
    Phi_phi = np.exp(1j * m * phi) / np.sqrt(2 * np.pi)

    # === Total wavefunction Ψ(ξ, η, φ) on 3D grid
    Psi = R_xi * T_eta * Phi_phi

    return Psi
