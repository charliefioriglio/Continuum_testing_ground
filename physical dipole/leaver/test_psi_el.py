"""
Minimal test script for psi_el - creates a single yz-plane slice plot.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import spherical_jn
from scipy.special import sph_harm_y as sph_harm
from build_mode_annie import psi_el

def test_psi_el_yz_slice():
    """Test psi_el and create a yz-plane slice plot with spherical comparison."""
    
    # ====== SPECIFY YOUR PARAMETERS HERE ======
    n = 0        # Principal quantum number
    m = -1        # Magnetic quantum number  
    E = 0.1/ 27.2     # Energy parameter
    a = 1.57     # Dipole separation
    D = 0.0      # Dipole strength
    # ===========================================
    
    print(f"Testing psi_el with parameters: n={n}, m={m}, E={E}, a={a}, D={D}")
    
    # Create 3D grid
    N = 50
    L = 20.0  # Grid extends from -L to L
    coords = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')
    
    try:
        # Compute prolate spheroidal wavefunction
        Psi_prolate, info = psi_el(n, m, E, a, D, X, Y, Z)
       
        # Compute single-center spherical wavefunction j_l(kr) * Y_l^m(θ,φ)
        l = n + np.abs(m)  # For comparison, use l = n + |m|
        k = np.sqrt(2 * E)  # Wave number from energy
        
        # Convert to spherical coordinates
        r = np.sqrt(X**2 + Y**2 + Z**2)
        theta = np.arccos(np.divide(Z, r, out=np.zeros_like(Z), where=(r!=0)))
        phi = np.arctan2(Y, X)
        
        # Compute spherical Bessel function and spherical harmonic
        j_l = spherical_jn(l, k * r)
        Y_lm = sph_harm(l, m, theta, phi)
        Psi_spherical = j_l * Y_lm
        
        # Create yz-plane slice (x=0) for 2D visualization
        mid_x = X.shape[0] // 2
        y_slice = Y[mid_x, :, :]
        z_slice = Z[mid_x, :, :]
        psi_prolate_slice = Psi_prolate[mid_x, :, :]
        psi_spherical_slice = Psi_spherical[mid_x, :, :]
        
        # Create angular and radial function comparisons
        
        # Angular comparison: Y_l^m(θ, φ=0) vs prolate angular function
        eta_vals = np.linspace(-1, 1, 100)  # η = cos(θ) from -1 to 1
        theta_vals = np.arccos(eta_vals)  # Convert back to θ for spherical harmonics
        phi_val = 0.0  # Fix φ = 0 for comparison
        
        # Spherical harmonics at φ=0
        Y_spherical_angular = sph_harm(l, m, theta_vals, phi_val)
        
        # Get the prolate spheroidal angular function: eigenvec[n] * P_l^m(η)
        from angular import analytic
        eigvals, eigvecs, *_ = analytic(abs(m), 20, E, a, D)
        
        if n < eigvecs.shape[1]:
            # Build the angular function as sum over l: eigenvec[n][i] * P_l^m(η)
            from scipy.special import lpmv
            Y_prolate_angular = np.zeros_like(eta_vals, dtype=complex)
            
            # Get the l values used in the expansion (starts from |m|)
            l_vals = np.arange(abs(m), abs(m) + eigvecs.shape[0])
            eigvec_n = eigvecs[:, n]  # Get the n-th eigenvector
            
            for i, l_val in enumerate(l_vals):
                if i < len(eigvec_n):
                    # P_l^m(η) where η = cos(θ)
                    P_lm = lpmv(abs(m), l_val, eta_vals)
                    Y_prolate_angular += eigvec_n[i] * P_lm * 1/np.sqrt(2*np.pi)  # Normalization factor
        else:
            Y_prolate_angular = np.zeros_like(eta_vals)
        
        # Radial comparison: j_l(kr) vs prolate radial function  
        r_vals = np.linspace(0.1, 10.0, 100)  # Avoid r=0
        xi_vals = r_vals / a  # Convert to prolate spheroidal coordinate ξ = r/a
        
        # Spherical Bessel function
        j_l_radial = spherical_jn(l, xi_vals * np.sqrt(2*E*a**2))
        
        # Get the prolate spheroidal radial function
        from radial_annie import radial_function_leaver
        lambda_mn = abs(eigvals[n])
        R_prolate_radial, v = radial_function_leaver(xi_vals, np.sqrt(2*E*a**2), abs(m), lambda_mn, n)
        
        # Create comprehensive comparison plot
        fig = plt.figure(figsize=(14, 10))
        
        # 1. Full wavefunction comparison (yz-plane slice)
        ax1 = plt.subplot(2, 2, 1)
        im1 = ax1.contourf(y_slice, z_slice, np.abs(psi_prolate_slice), levels=20, cmap='viridis')
        ax1.set_title(f'Prolate Spheroidal |Ψ_{{{n},{m}}}|\nD={D}, a={a}, E={E}')
        ax1.set_xlabel('y')
        ax1.set_ylabel('z')
        ax1.set_aspect('equal')
        plt.colorbar(im1, ax=ax1, label='|Ψ|')
        
        ax2 = plt.subplot(2, 2, 2)
        im2 = ax2.contourf(y_slice, z_slice, np.abs(psi_spherical_slice), levels=20, cmap='viridis')
        ax2.set_title(f'Spherical j_{{{l}}}(kr)Y_{{{l}}}^{{{m}}}\nk={k:.3f}, E={E}')
        ax2.set_xlabel('y')
        ax2.set_ylabel('z')
        ax2.set_aspect('equal')
        plt.colorbar(im2, ax=ax2, label='|Ψ|')
        
        # 2. Angular function comparison vs η = cos(θ)
        ax3 = plt.subplot(2, 2, 3)
        ax3.plot(eta_vals, np.real(Y_prolate_angular), 'b-', linewidth=2, label=f'Prolate Y_{{{n},{abs(m)}}}(η)')
        ax3.plot(eta_vals, np.real(Y_spherical_angular), 'r--', linewidth=2, label=f'Spherical Y_{{{l}}}^{{{m}}}(θ)')
        ax3.set_xlabel('η = cos(θ)')
        ax3.set_ylabel('Re[Y]')
        ax3.set_xlim(-1, 1)
        ax3.set_ylim(-1, 1)
        ax3.set_title(f'Angular Functions (φ=0)\nProlate vs Spherical')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 3. Radial function comparison
        ax4 = plt.subplot(2, 2, 4)
        ax4.plot(xi_vals, np.real(R_prolate_radial), 'b-', linewidth=2, label=f'Prolate R_{{{n},{abs(m)}}}(ξ)')
        ax4.plot(r_vals/a, j_l_radial, 'r--', linewidth=2, label=f'Spherical j_{{{l}}}(kr)')
        ax4.set_xlabel('ξ = r/a')
        ax4.set_ylabel('R(ξ) / j(kr)')
        ax4.set_title(f'Radial Functions\nProlate vs Spherical')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = f'psi_yz_slice_n{n}_m{m}_D{D:.3f}_a{a:.1f}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"Plot saved as: {filename}")
        plt.show()
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_psi_el_yz_slice()
