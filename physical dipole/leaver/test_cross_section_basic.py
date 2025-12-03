#!/usr/bin/env python3
"""
Simple test script for the physical dipole cross section calculation.

This script tests the basic functionality of the physical dipole cross section 
calculation before running the full version.
"""

import numpy as np
from build_mode_annie import psi_el
import sys
import os

def test_physical_dipole_wavefunction():
    """Test if we can build a basic physical dipole wavefunction."""
    
    print("Testing physical dipole wavefunction construction...")
    print("="*55)
    
    # Simple grid for testing
    L = 10.0  # Smaller box for testing
    n_pts = 20  # Fewer points for speed
    coords = np.linspace(-L, L, n_pts)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')
    
    # Physical parameters
    n = 0  # Principal quantum number
    m = 0  # Azimuthal quantum number
    E = 0.1  # Energy (in atomic units)
    a = 1.0  # Semi-focal distance
    D = 0.0  # No dipole for initial test
    
    print(f"Grid: {n_pts}³ points over ±{L} a.u.")
    print(f"Quantum numbers: n={n}, m={m}")
    print(f"Parameters: E={E}, a={a}, D={D}")
    print()
    
    try:
        # Build the wavefunction
        psi, info = psi_el(n, m, E, a, D, X, Y, Z)
        
        print("✓ Wavefunction calculated successfully!")
        print(f"  Shape: {psi.shape}")
        print(f"  Max |ψ|: {np.max(np.abs(psi)):.6e}")
        print(f"  Characteristic exponent ν: {info['eigenvalues']['nu']:.6f}")
        print(f"  Angular eigenvalue λ: {info['eigenvalues']['lambda']:.6f}")
        
        # Check if it's mostly non-zero
        non_zero_fraction = np.sum(np.abs(psi) > 1e-10) / psi.size
        print(f"  Non-zero elements: {non_zero_fraction:.1%}")
        
        return True
        
    except Exception as e:
        print("✗ Wavefunction calculation failed!")
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_multiple_quantum_numbers():
    """Test building wavefunctions for multiple quantum numbers."""
    
    print("\n\nTesting multiple quantum numbers...")
    print("="*35)
    
    # Simple grid
    L = 8.0
    n_pts = 15
    coords = np.linspace(-L, L, n_pts)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')
    
    # Test parameters
    E = 0.1
    a = 1.0
    D = 0.0  # Start with no dipole
    
    # Test various (n, m) combinations
    quantum_numbers = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0)]
    
    successful = 0
    total = len(quantum_numbers)
    
    for n, m in quantum_numbers:
        try:
            psi, info = psi_el(n, m, E, a, D, X, Y, Z)
            max_val = np.max(np.abs(psi))
            nu = info['eigenvalues']['nu']
            
            print(f"  (n={n}, m={m}): ✓ |ψ|_max = {max_val:.3e}, ν = {nu:.3f}")
            successful += 1
            
        except Exception as e:
            print(f"  (n={n}, m={m}): ✗ Failed - {e}")
    
    print(f"\nSuccess rate: {successful}/{total} ({100*successful/total:.0f}%)")
    return successful == total


def test_with_dipole():
    """Test with non-zero dipole moment."""
    
    print("\n\nTesting with non-zero dipole moment...")
    print("="*38)
    
    # Simple grid
    L = 8.0
    n_pts = 15
    coords = np.linspace(-L, L, n_pts)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')
    
    # Test with dipole
    n, m = 0, 0
    E = 0.1
    a = 1.0
    D_values = [0.0, 0.5, 1.0]
    
    success_count = 0
    for D in D_values:
        try:
            psi, info = psi_el(n, m, E, a, D, X, Y, Z)
            max_val = np.max(np.abs(psi))
            nu = info['eigenvalues']['nu']
            
            print(f"  D={D}: ✓ |ψ|_max = {max_val:.3e}, ν = {nu:.3f}")
            success_count += 1
            
        except Exception as e:
            print(f"  D={D}: ✗ Failed - {e}")
    
    return success_count == len(D_values)


def test_cross_section_ingredients():
    """Test the basic ingredients needed for cross section calculation."""
    
    print("\n\nTesting cross section calculation ingredients...")
    print("="*48)
    
    # Try to import the molecular orbital coefficients
    try:
        sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'point dipole'))
        from CuO_do import build_DO, DO_coeffs_b1_L, DO_coeffs_b1_R, b1_L_norm, b1_R_norm
        print("✓ Successfully imported molecular orbital data")
        
        # Test building molecular orbitals
        L = 10.0
        n_pts = 20
        coords = np.linspace(-L, L, n_pts)
        X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')
        dx = coords[1] - coords[0]
        dV = dx**3
        
        DO_L = build_DO(DO_coeffs_b1_L, X, Y, Z, dV)
        DO_R = build_DO(DO_coeffs_b1_R, X, Y, Z, dV)
        
        print(f"✓ Built molecular orbitals: |DO_L|_max = {np.max(np.abs(DO_L)):.3e}")
        print(f"                          |DO_R|_max = {np.max(np.abs(DO_R)):.3e}")
        
    except Exception as e:
        print(f"✗ Failed to import/build molecular orbitals: {e}")
        return False
        
    return True


def main():
    """Run all tests."""
    
    print("Physical Dipole Cross Section - Basic Tests")
    print("=" * 50)
    print()
    
    tests = [
        test_physical_dipole_wavefunction,
        test_multiple_quantum_numbers, 
        test_with_dipole,
        test_cross_section_ingredients
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"\n✗ Test failed with exception: {e}")
            results.append(False)
    
    print("\n" + "="*50)
    print("TEST SUMMARY")
    print("="*50)
    
    test_names = [
        "Basic wavefunction",
        "Multiple quantum numbers", 
        "Non-zero dipole",
        "Cross section ingredients"
    ]
    
    passed = sum(results)
    total = len(results)
    
    for name, result in zip(test_names, results):
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{name:.<30} {status}")
    
    print("-" * 50)
    print(f"Overall: {passed}/{total} tests passed ({100*passed/total:.0f}%)")
    
    if passed == total:
        print("\n🎉 All tests passed! Ready to run full cross section calculation.")
    else:
        print(f"\n⚠️  {total-passed} test(s) failed. Please fix before running full calculation.")
    
    return passed == total


if __name__ == "__main__":
    main()
