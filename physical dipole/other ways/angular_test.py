import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lpmv, factorial
from scipy.linalg import eig

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

    eigvals, eigvecs = eig(H, S)
    return eigvals, eigvecs, ell_vals


# ---------------------
# Parameters
L_max = 50
m = 0
n = 0  # lowest mode
E = 0.1 / 27.2
a = 1.67 / (52.9e-2 * 2)
eta = np.linspace(-1, 1, 300)
Dipoles = np.array([-1, -0.5, 0, 0.5])

plt.figure(figsize=(10, 6))
for D in Dipoles:
    eigvals, eigvecs, ell_vals = analytic(m, L_max, E, a, D)
    coeffs = eigvecs[:, n]  # n-th eigenvector
    print(D, -eigvals[n])
    # Construct T_{m,n}(η) = ∑_l c_l * P_l^m(η)
    P_basis = np.array([lpmv(m, l, eta) for l in ell_vals])
    T_eta = coeffs @ P_basis

    plt.plot(eta, T_eta, label=f"D = {D:.2f}")

plt.xlabel(r'$\eta$')
plt.ylabel(r'$T_{m,n}(\eta)$')
plt.ylim(-2, 2)
plt.title(rf'Angular functions $T_{{m,n}}(\eta)$ for $m={m}$, $n={n}$ at varying dipole moments')
plt.legend(loc='best', fontsize='small')
plt.grid(True)
plt.tight_layout()
plt.show()
