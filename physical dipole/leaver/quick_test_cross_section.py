#!/usr/bin/env python3
"""
Quick test run of the physical dipole cross section calculation.
This uses smaller grids and fewer parameters for faster testing.
"""

import numpy as np
from tqdm import tqdm
from calculate_physical_dipole import (
    build_continuum_mode_physical_dipole,
    calculate_total_cross_sections_physical_dipole,
    integrate_3d
)
import sys
import os

# Add the parent directory to the path to import CuO_do
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'point dipole'))
from CuO_do import (build_DO, DO_coeffs_b1_L, DO_coeffs_b1_R, 
                   b1_L_norm, b1_R_norm)

def quick_test():
    """Quick test with minimal parameters."""
    
    print("Physical Dipole Cross Section - Quick Test")
    print("=" * 45)
    
    # Small grid for testing
    L = 8.0  # Reduced from 19
    n_pts = 15  # Reduced from 50
    x = np.linspace(-L, L, n_pts)
    y = np.linspace(-L, L, n_pts)  
    z = np.linspace(-L, L, n_pts)
    dx = x[1] - x[0]
    dV = dx**3
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    print(f"Grid: {n_pts}³ points over ±{L} a.u.")
    print(f"Volume element: {dV:.6f} a.u.³")
    
    # Test parameters
    Trans_E = 1.778  # Single transition energy
    E_photon_grid = np.array([Trans_E + 2.7, Trans_E + 5.4])  # Higher photon energies for better convergence
    L_max = 2  # Reduced from 3
    D = 0.0   # No dipole for first test
    a_focal = 1.0
    
    print(f"Transition energy: {Trans_E} eV")
    print(f"Photon energies: {E_photon_grid}")
    print(f"L_max: {L_max}")
    print(f"Dipole strength D: {D}")
    print(f"Semi-focal distance a: {a_focal}")
    print()
    
    # Constants
    hartree = 27.2114  # eV
    
    # Test single continuum mode first
    print("Testing single continuum mode...")
    eKE = (E_photon_grid[0] - Trans_E) / hartree
    print(f"Electron kinetic energy: {eKE:.6f} a.u.")
    
    try:
        psi_test = build_continuum_mode_physical_dipole(D, 0, 0, eKE, X, Y, Z, a_focal)
        if psi_test is not None:
            print(f"✓ Continuum mode (m=0, n=0): |ψ|_max = {np.max(np.abs(psi_test)):.3e}")
            
            # Test dipole matrix element
            dipole_z = Z
            DO_L = build_DO(DO_coeffs_b1_L, X, Y, Z, dV)
            
            integrand = np.conj(psi_test) * dipole_z * DO_L
            matrix_element = integrate_3d(integrand, dV)
            
            print(f"✓ Sample matrix element: {matrix_element:.6e}")
        else:
            print("✗ Failed to build continuum mode")
            return False
    except Exception as e:
        print(f"✗ Error building continuum mode: {e}")
        return False
    
    # Test full cross section calculation  
    print("\nTesting full cross section calculation...")
    try:
        total_cross_sections = calculate_total_cross_sections_physical_dipole(
            E_photon_grid, Trans_E, X, Y, Z,
            DO_coeffs_b1_L, DO_coeffs_b1_R, b1_L_norm, b1_R_norm,
            L_max, D, a_focal
        )
        
        print("✓ Cross section calculation completed!")
        print(f"   Results: {total_cross_sections}")
        
        # Check if results are reasonable
        if all(isinstance(x, (int, float, complex)) and not np.isnan(x) for x in total_cross_sections):
            print("✓ Results appear valid (no NaN values)")
            return True
        else:
            print("⚠ Warning: Some results are NaN or invalid")
            return False
            
    except Exception as e:
        print(f"✗ Cross section calculation failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = quick_test()
    
    print("\n" + "="*45)
    if success:
        print("🎉 Quick test PASSED! Cross section calculation works.")
        print("Ready to run full calculation with complete parameters.")
    else:
        print("❌ Quick test FAILED! Please debug before full calculation.")
    print("="*45)
