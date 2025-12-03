import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm, spherical_jn, eval_legendre

from angular import build_angular_functions
from radial import radial_function
from build_mode import coefficient

# --- Parameters ---
D = 0.1
m = 0
N = 0
l = abs(m) + N
l_max = 40
E = 1
k = np.sqrt(2 * E)

# --- Angular plot: Ω_Nλ(θ) ---
theta = np.linspace(0, np.pi, 500)
phi = np.zeros_like(theta)  # fix φ = 0

# --- Angular plot: Ω_Nλ(θ) vs Y_lm(θ) ---
plt.figure(figsize=(6, 4))

omega_funcs, eigvals = build_angular_functions(D, m, l_max, theta, phi)
omega_theta = omega_funcs[N]
plt.plot(theta, omega_theta**2, label=r'Point dipole $|\Omega_{N\lambda}(\theta)|^2$')

# Plane wave equivalent: |Y_lm(θ, φ=0)|^2 with l=0, m=0
Y00 = sph_harm(m, l, phi, theta)
plt.plot(theta, np.abs(Y00)**2, '--', label=r'Plane wave $|Y_{00}(\theta)|^2$')

plt.xlabel(r'$\theta$ (rad)')
plt.ylabel(r'Angular Density')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# --- Radial plot: R_Nλ(kr) ---
r = np.linspace(1e-4, 10, 1000)
kr = k * r

plt.figure(figsize=(6, 4))
_, eigvals = build_angular_functions(D, m, l_max, theta, phi)
L_N = 0.5 * (-1 + np.sqrt(1 + 4 * eigvals[N]))
print(L_N)

# Point dipole radial
R_dipole = radial_function(L_N, k, r)
plt.plot(r, R_dipole, label=r'Point dipole $R_{N\lambda}(kr)$')

# Plane wave radial: j_0(kr) for l=0
j0 = spherical_jn(l, kr)
plt.plot(r, j0, '--', label=r'Plane wave $j_0(kr)$')

plt.xlabel(r'$r$')
plt.ylabel('Radial Function')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# --- 2D yZ slice vs plane-wave expansion (k || z) ---
slice_extent = 6.0
slice_points = 301
y_vals = np.linspace(-slice_extent, slice_extent, slice_points)
z_vals = np.linspace(-slice_extent, slice_extent, slice_points)
Y, Z = np.meshgrid(y_vals, z_vals, indexing='ij')
X = np.zeros_like(Y)

r_slice = np.sqrt(X**2 + Y**2 + Z**2)
theta_slice = np.zeros_like(r_slice)
phi_slice = np.zeros_like(r_slice)
nonzero_mask = r_slice > 1e-8
theta_slice[nonzero_mask] = np.arccos(np.clip(Z[nonzero_mask] / r_slice[nonzero_mask], -1.0, 1.0))
phi_slice[nonzero_mask] = np.mod(np.arctan2(Y[nonzero_mask], X[nonzero_mask]), 2 * np.pi)

omega_grid, eigvals_slice = build_angular_functions(D, m, l_max, theta_slice, phi_slice)
omega_slice = omega_grid[N]
R_slice = radial_function(L_N, k, r_slice)
psi_dipole = R_slice * omega_slice

cos_theta = np.ones_like(r_slice)
cos_theta[nonzero_mask] = np.clip(Z[nonzero_mask] / r_slice[nonzero_mask], -1.0, 1.0)
kr_slice = k * r_slice
plane_wave = np.zeros_like(r_slice, dtype=np.complex128)
l_expansion = l_max
for ell in range(l_expansion + 1):
	jl = spherical_jn(ell, kr_slice)
	Pl = eval_legendre(ell, cos_theta)
	plane_wave += (2 * ell + 1) * (1j ** ell) * jl * Pl
plane_wave[~nonzero_mask] = 1.0

# Point-dipole expansion using continuum coefficients
k_hat = np.array([0.0, 0.0, 1.0])
psi_expansion = np.zeros_like(r_slice, dtype=np.complex128)
num_modes = omega_grid.shape[0]
for mode_index in range(num_modes):
	L_mode = 0.5 * (-1.0 + np.sqrt(1.0 + 4.0 * eigvals_slice[mode_index]))
	if eigvals_slice[mode_index] < -0.25 or np.isnan(eigvals_slice[mode_index]):
		continue
	R_mode = radial_function(L_mode, k, r_slice)
	coeff = coefficient(D, m, mode_index, k_hat, l_max=l_max)
	psi_expansion += coeff * R_mode * omega_grid[mode_index]

fields = [np.real(psi_dipole), np.real(plane_wave), np.real(psi_expansion)]
titles = [
	'Point dipole Re[$\\psi_N$] (yZ slice)',
	'Plane-wave expansion Re[$\\psi$]',
	'Point-dipole expansion Re[$\\Psi$]',
]
vmax = max(np.max(np.abs(field)) for field in fields if np.any(np.isfinite(field)))

fig, axes = plt.subplots(1, 3, figsize=(16, 4), sharex=True, sharey=True)
extent = [z_vals.min(), z_vals.max(), y_vals.min(), y_vals.max()]

images = []
for ax, field, title in zip(axes, fields, titles):
	im = ax.imshow(field, extent=extent, origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
	ax.set_title(title)
	ax.set_xlabel('z')
	images.append(im)

axes[0].set_ylabel('y')

for ax, im in zip(axes, images):
	plt.colorbar(im, ax=ax, shrink=0.85)

plt.tight_layout()
plt.show()

