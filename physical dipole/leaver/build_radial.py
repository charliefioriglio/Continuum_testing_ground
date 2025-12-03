"""
Minimal script to calculate radial functions for wavefunction construction.

This provides a clean interface for building complete wavefunctions by
computing the radial part R(xi) for given quantum numbers and coordinates.
"""

import numpy as np
from radial_annie import radial_function_leaver
from angular import analytic

def build_radial(n, m, E, D, a, xi_coordinates):
    """
    Build radial function R(xi) for given quantum numbers and coordinates.
    
    This is the minimal interface for wavefunction construction.
    
    Parameters:
    -----------
    n : int
        Principal quantum number (radial node count)
    m : int  
        Azimuthal quantum number
    D : float
        Dipole strength parameter
    a : float
        Semi-focal distance parameter
    xi_coordinates : array_like
        Xi coordinate values where to evaluate R(xi)
        
    Returns:
    --------
    R_xi : array
        Radial function values R(xi)
    info : dict
        Additional information: {'v': characteristic_exponent, 'c': parameter_c, 
                                'E': energy, 'lambda': angular_eigenvalue}
    """
    
    c = np.sqrt(2 * E * a**2)
    
    # Get angular eigenvalue
    L_max = max(20, n + abs(m) + 5)
    eigvals, *_ = analytic(m, L_max, E, a, D)
    
    if n >= len(eigvals):
        raise ValueError(f"Not enough eigenvalues: need index {n}, got {len(eigvals)}")
    
    lambda_mn = abs(eigvals[n])
    
    # Compute radial function using Leaver's method
    xi_coordinates = np.asarray(xi_coordinates)
    R_xi, v = radial_function_leaver(xi_coordinates, c, m, lambda_mn, n)
    
    # Package additional information
    info = {
        'v': v,
        'c': c, 
        'E': E,
        'lambda': lambda_mn,
        'n': n,
        'm': m,
        'D': D,
        'a': a
    }
    
    return R_xi, info

def build_radial_simple(n, m, D, a, xi_min=1.0, xi_max=10.0, num_points=100):
    """
    Simplified interface with automatic coordinate generation.
    
    Parameters:
    -----------
    n, m : int
        Quantum numbers
    D : float
        Dipole strength
    a : float
        Semi-focal distance
    xi_min, xi_max : float
        Range of xi coordinates (default: 1.0 to 10.0)
    num_points : int
        Number of coordinate points (default: 100)
        
    Returns:
    --------
    xi_coords : array
        Xi coordinate values
    R_xi : array
        Radial function values
    info : dict
        Additional information
    """
    
    # Generate coordinate grid
    xi_coords = np.linspace(xi_min, xi_max, num_points)
    
    # Compute radial function
    R_xi, info = build_radial(n, m, D, a, xi_coords)
    
    return xi_coords, R_xi, info