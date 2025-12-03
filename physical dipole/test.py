import numpy as np
import matplotlib.pyplot as plt
from radial_numerical import compute_radial_function
from radial import compute_radial_function_full
from angular import analytic

# Physical/system parameters
m = 0
n = 1   
L_max = 50    
R_max = 50
E = 0.05 / 27.2
a = 1.67 / (2 * 52.9e-2) 
Dipoles = np.array([0.85, 0])
xi_vals = np.linspace(1+1e-10, 12, 200)
c = np.sqrt(2 * E * a**2)

# Plotting
plt.figure(figsize=(8, 5))
for D in Dipoles:
    eigvals, eigvecs, ell_vals = analytic(m, L_max, E, a, D)
    lam_mn = -eigvals[n]
    R_vals = compute_radial_function_full(m, n, lam_mn, c, xi_vals)
    plt.plot(xi_vals, R_vals, label=fr"$D={D}$")
plt.xlabel(r"$\xi$")
plt.ylabel(r"$R_{m,n}(\xi)$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

