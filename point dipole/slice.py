import numpy as np
from build_mode import build_continuum_mode, coefficient
import matplotlib.pyplot as plt

def plane_wave(k_hat, D, L_max, eKE, x, y, z):
    psi = np.zeros_like(x, dtype=complex)
    for m in range(-L_max, L_max + 1):
        for N in range(L_max + 1 - abs(m)):
            mode = build_continuum_mode(D, m, N, eKE, x, y, z)
            C = coefficient(D, m, N, k_hat)
            psi += C * mode
    return psi

def run_example(D=0.0, L_max=5, eKE=1.0, grid_size=101, L=20.0, k_hat=(0, 1, 0)):
    """
    Compute and plot continuum wavefunction slice for given parameters.

    Parameters
    ----------
    D : float
        Dipole strength (internuclear half-separation).
    L_max : int
        Angular momentum cutoff.
    eKE : float
        Electron kinetic energy.
    grid_size : int
        Number of grid points per axis.
    L : float
        Extent of grid in each dimension.
    k_hat : tuple
        Unit vector for momentum direction.
    """
    # Cartesian grid
    x = np.linspace(-L, L, grid_size)
    y = np.linspace(-L, L, grid_size)
    z = np.linspace(-L, L, grid_size)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Build plane wave on grid
    psi = plane_wave(k_hat, D, L_max, eKE, X, Y, Z)

    # Take a 2D slice at z = 0 plane
    mid = grid_size // 2
    psi_slice = np.real(psi[:, :, mid])  # real part

    # Plot
    plt.figure(figsize=(6,5))
    plt.pcolormesh(x, y, psi_slice.T, shading='auto', cmap='viridis')
    plt.colorbar(label='Re[ψ(x,y,z=0)]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Continuum Wavefunction Slice\nD={D}, L_max={L_max}, E={eKE}, k̂={k_hat}')
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_example()