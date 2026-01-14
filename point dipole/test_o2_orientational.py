"""
O2- Beta Test Script for Orientational Averaging Code
Runs beta calculations for D=0 and D=0.3 at 10 energies from 0.01 to 5.01 eV
"""

import numpy as np
import sys
sys.path.insert(0, '/Users/charliefioriglio/Desktop/research/code/continuum_testing_ground/point dipole')

from build_mode import build_directional_wavefunction
from o2_sto_do import build_DO, DO_coeffs_b3g, DO_coeffs_b2g
from ang_grid import repulsion_orientations

# Energy and dipole parameters
E_values = np.linspace(0.01, 5.01, 10)  # eV
D_values = [0.0, 0.3]  # Physical D values (not 2D)
hartree = 27.2114

# Grid parameters
L_max = 4  # Angular momentum for continuum
N_ang = 20  # Number of orientations for averaging
L = 10  # Grid extent
n_pts = 50  # Grid points
x = np.linspace(-L, L, n_pts)
y = np.linspace(-L, L, n_pts)
z = np.linspace(-L, L, n_pts)
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]
dV = dx * dy * dz
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

def integrate_3d(f, dV):
    return np.sum(f) * dV

def rotation_matrix(alpha, beta, gamma):
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta),  np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)
    Rz1 = np.array([[ ca, -sa, 0], [ sa,  ca, 0], [  0,   0, 1]])
    Ry  = np.array([[ cb,   0, sb], [  0,   1,  0], [-sb,  0, cb]])
    Rz2 = np.array([[ cg, -sg, 0], [ sg,  cg, 0], [  0,   0, 1]])
    return Rz1 @ Ry @ Rz2

def rotate_grid(X, Y, Z, alpha, beta, gamma):
    R = rotation_matrix(alpha, beta, gamma)
    shape = X.shape
    coords = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=0)
    rotated_coords = R @ coords
    X_rot = rotated_coords[0].reshape(shape)
    Y_rot = rotated_coords[1].reshape(shape)
    Z_rot = rotated_coords[2].reshape(shape)
    return X_rot, Y_rot, Z_rot

def continuum(k_hat, D, L_max, eKE, X, Y, Z):
    return build_directional_wavefunction(D, k_hat, eKE, X, Y, Z, l_max=L_max)

def compute_average_amplitudes(X, Y, Z, DO_coeffs, ang_grid, mu_field, continuum_field, dV):
    N = ang_grid.shape[0]
    amps = np.zeros(N, dtype=complex)
    
    for i in range(N):
        alpha, beta, gamma = ang_grid[i]
        X_rot, Y_rot, Z_rot = rotate_grid(X, Y, Z, alpha, beta, gamma)
        
        DO = build_DO(DO_coeffs, X_rot, Y_rot, Z_rot, dV)
        
        integrand_L = np.conj(DO) * mu_field * continuum_field
        integrand_R = np.conj(continuum_field) * mu_field * DO
        A_L = integrate_3d(integrand_L, dV)
        A_R = integrate_3d(integrand_R, dV)
        
        amps[i] = A_L * A_R
    
    return np.real(amps)

print("O2- Beta Calculation (Orientational Averaging)")
print("=" * 55)
print(f"Grid: {n_pts}^3, extent = ±{L} a.u.")
print(f"L_max = {L_max}, N_orientations = {N_ang}")
print(f"Energies: {len(E_values)} points from {E_values[0]} to {E_values[-1]} eV")
print()

# Get orientation grid
ang_grid = repulsion_orientations(n_orientations=N_ang)

directions = {
    "par": np.array([0.0, 0.0, 1.0]),
    "perp1": np.array([1.0, 0.0, 0.0]),
    "perp2": np.array([0.0, 1.0, 0.0]),
}

degeneracies = [DO_coeffs_b3g, DO_coeffs_b2g]

pol_lab = np.array([0.0, 0.0, 1.0])
mu_field = pol_lab[0] * X + pol_lab[1] * Y + pol_lab[2] * Z

results = {}

for D in D_values:
    print(f"\n=== D = {D:.2f} ===")
    results[D] = []
    
    for E in E_values:
        E_au = E / hartree
        
        # Build continuum for each direction
        continuum_cache = {}
        for key, direction in directions.items():
            continuum_cache[key] = continuum(direction, D, L_max, E_au, X, Y, Z)
        
        sigma_par = 0.0
        sigma_perp_accum = 0.0
        
        for DO_coeffs in degeneracies:
            A_par = compute_average_amplitudes(X, Y, Z, DO_coeffs, ang_grid, 
                                               mu_field, continuum_cache["par"], dV)
            A_perp1 = compute_average_amplitudes(X, Y, Z, DO_coeffs, ang_grid,
                                                 mu_field, continuum_cache["perp1"], dV)
            A_perp2 = compute_average_amplitudes(X, Y, Z, DO_coeffs, ang_grid,
                                                 mu_field, continuum_cache["perp2"], dV)
            
            sigma_par += np.mean(A_par)
            sigma_perp_accum += 0.5 * (np.mean(A_perp1) + np.mean(A_perp2))
        
        sigma_perp = sigma_perp_accum
        beta = 2 * (sigma_par - sigma_perp) / (sigma_par + 2 * sigma_perp)
        
        results[D].append((E, beta))
        print(f"  E = {E:.2f} eV, beta = {beta:.4f}")

# Print summary
print("\n" + "=" * 60)
print("SUMMARY - O2- Orientational Averaging")
print("=" * 60)
print(f"{'E (eV)':<10}", end="")
for D in D_values:
    print(f"{'D='+str(D):<15}", end="")
print()
print("-" * 40)

for i, E in enumerate(E_values):
    print(f"{E:<10.2f}", end="")
    for D in D_values:
        beta = results[D][i][1]
        print(f"{beta:<15.4f}", end="")
    print()

np.savez("o2_orientational_results.npz", 
         E_values=E_values, 
         D_values=D_values,
         results=results)
print("\nResults saved to o2_orientational_results.npz")
