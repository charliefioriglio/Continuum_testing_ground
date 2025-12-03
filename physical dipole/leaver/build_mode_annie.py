import numpy as np
from angular import analytic
from build_radial import build_radial
from phi import compute_Phi_m
from scipy.special import lpmv
from scipy.interpolate import interp1d

def build_single_mode_wavefunction_on_xyz_grid(m, n, E, a, D, X, Y, Z, L_max=20, xi_max=None, n_xi=200, n_eta=200, n_phi=100):
    """
    Build a single (m,n) wavefunction on a Cartesian (X,Y,Z) grid using prolate spheroidal coordinates.
    
    The wavefunction is:
    Ψ_{m,n}(X,Y,Z) = S_{m,n}(ξ) × T_{m,n}(η) × Φ_m(φ)
    
    where (ξ,η,φ) are prolate spheroidal coordinates.
    
    Parameters:
    -----------
    m, n : int
        Quantum numbers (azimuthal and radial)
    E : float
        Energy parameter  
    a : float
        Semi-focal distance
    D : float
        Dipole moment parameter
    X, Y, Z : ndarray
        Cartesian coordinate grids
    L_max : int, optional
        Maximum L for angular calculation (default: 20)
    xi_max : float, optional
        Maximum ξ value for radial grid (default: auto from X,Y,Z)
    n_xi, n_eta, n_phi : int, optional
        Grid sizes for interpolation (default: 200, 200, 100)
        
    Returns:
    --------
    Psi_mn : ndarray
        Wavefunction values on the (X,Y,Z) grid
    info : dict
        Information about the calculation including ν value
    """

    
    # Convert Cartesian (X,Y,Z) → prolate spheroidal (ξ, η, φ)
    Z_A = a
    Z_B = -a
    rA = np.sqrt(X**2 + Y**2 + (Z - Z_A)**2)
    rB = np.sqrt(X**2 + Y**2 + (Z - Z_B)**2)

    xi = (rA + rB) / (2 * a)
    eta = (rA - rB) / (2 * a)
    phi = np.arctan2(Y, X)
    
    # Determine coordinate ranges for interpolation grids
    if xi_max is None:
        xi_max = np.max(xi) * 1.1  # Add 10% buffer
    
    # Create 1D coordinate grids for function evaluation
    xi_vals = np.linspace(1.0, xi_max, n_xi)
    eta_vals = np.linspace(-1.0, 1.0, n_eta) 
    phi_vals = np.linspace(-np.pi, np.pi, n_phi)

    # Physical parameters
    m_abs = abs(m)
    c = np.sqrt(2 * E * a**2)

    # ========================================
    # 1. ANGULAR PART: T_{m,n}(η)
    # ========================================
    
    # Get angular eigenvalues and eigenvectors
    eigvals_ang, eigvecs_ang, ell_vals = analytic(m_abs, L_max, E, a, D)
    if n >= len(eigvals_ang):
        raise ValueError(f"Not enough angular modes. Have {len(eigvals_ang)}, need n={n}")
    
    v_ang = eigvecs_ang[:, n]  # Angular eigenvector for mode n
    
    # Build angular function on η grid using Legendre polynomials
    P_basis = np.array([lpmv(m_abs, l, eta_vals) for l in ell_vals])
    T_eta_vals = P_basis.T @ v_ang
    
    # Apply phase factor for negative m
    if m < 0:
        T_eta_vals = T_eta_vals * (-1)**m_abs
    
    # Create interpolation function for T(η)
    T_eta_func = interp1d(eta_vals, T_eta_vals, kind='cubic', 
                          bounds_error=False, fill_value=0.0)

    # ========================================  
    # 2. RADIAL PART: S_{m,n}(ξ)
    # ========================================
    
    try:
        # Use our new radial function builder
        S_xi_vals, radial_info = build_radial(n, m_abs, E, D, a, xi_vals)
    
    except Exception as e:
        raise RuntimeError(f"Failed to compute radial function: {e}")
    
    # Create interpolation function for S(ξ)
    S_xi_func = interp1d(xi_vals, S_xi_vals, kind='cubic',
                         bounds_error=False, fill_value=0.0)

    # ========================================
    # 3. AZIMUTHAL PART: Φ_m(φ)  
    # ========================================
    
    # Φ_m(φ) = exp(i m φ) / sqrt(2π)
    Phi_m = compute_Phi_m(m, phi)

    # ========================================
    # 4. COMBINE ALL PARTS
    # ========================================
    
    # Evaluate functions on the actual coordinate grids
    T_eta = T_eta_func(eta)
    S_xi = S_xi_func(xi)
    
    # Final wavefunction: Ψ_{m,n}(X,Y,Z) = S_{m,n}(ξ) × T_{m,n}(η) × Φ_m(φ)
    Psi_mn = S_xi * T_eta * Phi_m
    
    # Compile information about the calculation
    info = {
        'quantum_numbers': (m, n),
        'parameters': {'E': E, 'a': a, 'D': D, 'c': c},
        'characteristic_exponent': radial_info['v'],
        'angular_eigenvalue': abs(eigvals_ang[n]),
        'coordinate_ranges': {
            'xi_range': [np.min(xi), np.max(xi)],
            'eta_range': [np.min(eta), np.max(eta)],
            'phi_range': [np.min(phi), np.max(phi)]
        },
        'grid_sizes': {'n_xi': n_xi, 'n_eta': n_eta, 'n_phi': n_phi},
        'radial_info': radial_info
    }

    
    return Psi_mn, info


def build_multiple_mode_wavefunctions(quantum_numbers, E, a, D, X, Y, Z, **kwargs):
    """
    Build multiple (m,n) wavefunctions on the same Cartesian grid.
    
    Parameters:
    -----------
    quantum_numbers : list of tuples
        List of (m, n) quantum number pairs
    E, a, D : float
        Physical parameters
    X, Y, Z : ndarray
        Cartesian coordinate grids
    **kwargs : dict
        Additional arguments passed to build_single_mode_wavefunction_on_xyz_grid
        
    Returns:
    --------
    wavefunctions : dict
        Dictionary with keys (m,n) and values (Psi_mn, info)
    """
    
    
    wavefunctions = {}
    
    for i, (m, n) in enumerate(quantum_numbers):
        
        try:
            Psi_mn, info = build_single_mode_wavefunction_on_xyz_grid(
                m, n, E, a, D, X, Y, Z, **kwargs
            )
            
            wavefunctions[(m, n)] = {'wavefunction': Psi_mn, 'info': info}
            
        except Exception as e:
            print(f"ERROR: Failed to build (m={m}, n={n}): {e}")
            wavefunctions[(m, n)] = None
    
    successful = sum(1 for v in wavefunctions.values() if v is not None)
    
    return wavefunctions


def psi_el(n, m, E, a, D, X, Y, Z):
    """
    Construct a single continuum mode Ψ_{nm} = R_{nm}(ξ) × Y_{nm}(η,φ) on an xyz grid.
    
    This is the main interface for calculating specific continuum modes on a 3D grid.
    
    Parameters:
    -----------
    n : int
        Principal quantum number (radial node count)
    m : int  
        Azimuthal quantum number (can be negative)
    E : float
        Energy parameter
    a : float
        Semi-focal distance parameter
    D : float
        Dipole strength parameter
    X, Y, Z : ndarray
        Cartesian coordinate grids (must have same shape)
        
    Returns:
    --------
    Psi : ndarray (complex)
        Wavefunction values on the (X,Y,Z) grid
    info : dict
        Calculation information including ν, E, c, λ, etc.
    """

    
    # Convert Cartesian to prolate spheroidal coordinates
    rA = np.sqrt(X**2 + Y**2 + (Z - a)**2)
    rB = np.sqrt(X**2 + Y**2 + (Z + a)**2) 
    
    xi = (rA + rB) / (2 * a)
    eta = (rA - rB) / (2 * a)
    phi = np.arctan2(Y, X)
    
    # Calculate parameter c
    c = np.sqrt(2 * E * a**2)
    m_abs = abs(m)
    
    # ========================================
    # 1. RADIAL PART: R_{nm}(ξ)
    # ========================================
    
    # Get unique ξ values for efficient calculation
    xi_unique = np.unique(xi.flatten())
    xi_unique = xi_unique[xi_unique >= 1.0]  # ξ ≥ 1 for prolate spheroidal
    
    # Compute radial function using our optimized interface
    R_unique, radial_info = build_radial(n, m_abs, E, D, a, xi_unique)
    
    # Interpolate back to full grid
    from scipy.interpolate import interp1d
    R_interp = interp1d(xi_unique, R_unique, kind='cubic', 
                       bounds_error=False, fill_value=0.0)
    R_xi = R_interp(xi)
    
    # ========================================
    # 2. ANGULAR PART: Y_{nm}(η,φ) = T_{nm}(η) × Φ_m(φ)
    # ========================================
    
    # Get angular eigenvalues and eigenvectors
    L_max = max(20, n + m_abs + 5)
    eigvals_ang, eigvecs_ang, ell_vals = analytic(m_abs, L_max, E, a, D)
    
    if n >= len(eigvals_ang):
        raise ValueError(f"Not enough angular modes. Have {len(eigvals_ang)}, need n={n}")
    
    v_ang = eigvecs_ang[:, n]  # Angular eigenvector for mode n
    lambda_mn = abs(eigvals_ang[n])
    
    # Build T_{nm}(η) using Legendre polynomials
    from scipy.special import lpmv
    
    # Get unique η values for efficient calculation  
    eta_unique = np.unique(eta.flatten())
    eta_unique = eta_unique[(eta_unique >= -1.0) & (eta_unique <= 1.0)]
    
    # Compute Legendre polynomial basis on unique η values
    P_basis = np.array([lpmv(m_abs, l, eta_unique) for l in ell_vals])
    T_eta_unique = P_basis.T @ v_ang
    
    # Apply phase factor for negative m
    if m < 0:
        T_eta_unique = T_eta_unique * (-1)**m_abs
    
    # Interpolate back to full grid
    T_interp = interp1d(eta_unique, T_eta_unique, kind='cubic',
                       bounds_error=False, fill_value=0.0)
    T_eta = T_interp(eta)
    
    # Compute Φ_m(φ) = exp(im φ) / sqrt(2π)  
    Phi_m = compute_Phi_m(m, phi)
    
    
    # ========================================
    # 3. COMBINE: Ψ_{nm} = R_{nm}(ξ) × Y_{nm}(η,φ)
    # ========================================
    
    # Full wavefunction
    Psi = R_xi * T_eta * Phi_m
    
    # Handle regions where ξ < 1 (outside physical domain)
    mask = xi < 1.0
    Psi[mask] = 0.0
    
    # Compilation information
    info = {
        'quantum_numbers': (n, m),
        'parameters': {
            'E': E,
            'a': a, 
            'D': D,
            'c': c
        },
        'eigenvalues': {
            'nu': radial_info['v'],
            'lambda': lambda_mn
        },
        'coordinate_ranges': {
            'xi_min': np.min(xi),
            'xi_max': np.max(xi),
            'eta_min': np.min(eta), 
            'eta_max': np.max(eta)
        },
        'grid_shape': Psi.shape,
        'radial_info': radial_info
    }
    
    max_val = np.max(np.abs(Psi))
    
    return Psi, info

def psi_el_simple(n, m, E, a, D, grid_size=50, box_size=10.0):
    """
    Simplified interface that creates its own coordinate grid.
    
    Parameters:
    -----------
    n, m : int
        Quantum numbers
    E : float
        Energy parameter
    a : float
        Semi-focal distance parameter  
    D : float
        Dipole strength parameter
    grid_size : int
        Number of grid points per dimension (default: 50)
    box_size : float
        Half-size of cubic grid in units of 'a' (default: 10.0)
        
    Returns:
    --------
    X, Y, Z : ndarray
        Coordinate grids
    Psi : ndarray
        Wavefunction on grid
    info : dict
        Calculation information
    """
    
    # Create coordinate grids
    L = box_size * a  # Box size in physical units
    coords = np.linspace(-L, L, grid_size)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')
    
    # Compute wavefunction
    Psi, info = psi_el(n, m, E, a, D, X, Y, Z)
    
    return X, Y, Z, Psi, info

