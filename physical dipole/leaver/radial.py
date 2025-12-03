import numpy as np
import matplotlib.pyplot as plt
from radial_annie import radial_function_leaver
from angular import analytic

def compute_radial_function(n, m, D, a, c_xi_values, v_guess=None):
    """
    Compute the radial function R(xi) for given quantum numbers.
    
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
    c_xi_values : array_like
        Values of c*xi where to evaluate the function
    v_guess : float, optional
        Initial guess for characteristic exponent v
        
    Returns:
    --------
    R_values : array
        Radial function values at the given c*xi points
    v : float
        Characteristic exponent found
    xi_values : array
        Corresponding xi values (c_xi_values / c)
    c : float
        Parameter c used in calculation
    """
    
    # Calculate energy and c parameter from D and a
    # For the physical dipole case, E = D^2 / (8 * a^4)
    if D == 0:
        # Special case: no dipole field, use small energy
        E = 0.01
    else:
        E = D**2 / (8 * a**4)
    
    c = np.sqrt(2 * E * a**2)
    print(f"Calculated: E = {E:.6e}, c = {c:.6e}")
    
    # Convert c*xi values to xi values
    c_xi_values = np.asarray(c_xi_values)
    xi_values = c_xi_values / c
    
    # Get angular eigenvalue
    L_max = max(20, n + abs(m) + 5)  # Ensure enough eigenvalues
    eigvals, *_ = analytic(m, L_max, E, a, D)
    
    if n >= len(eigvals):
        raise ValueError(f"Not enough eigenvalues: need {n+1}, got {len(eigvals)}")
    
    lambda_mn = abs(eigvals[n])
    print(f"Angular eigenvalue λ_{m},{n} = {lambda_mn:.6f}")
    
    # Compute radial function using Leaver's method
    R_values, v = radial_function_leaver(xi_values, c, m, lambda_mn, n, v_guess)
    
    return R_values, v, xi_values, c

def plot_radial_function(n, m, D, a, c_xi_max=10.0, num_points=200, save_plot=True):
    """
    Plot the radial function for given quantum numbers.
    
    Parameters:
    -----------
    n, m : int
        Quantum numbers
    D : float
        Dipole strength
    a : float
        Semi-focal distance
    c_xi_max : float
        Maximum value of c*xi to plot
    num_points : int
        Number of points for plotting
    save_plot : bool
        Whether to save the plot to file
    """
    
    # Generate c*xi values
    c_xi_values = np.linspace(0.1, c_xi_max, num_points)
    
    try:
        # Compute radial function
        R_values, v, xi_values, c = compute_radial_function(n, m, D, a, c_xi_values)
        
        # Create the plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        
        # Plot 1: R(xi) vs c*xi
        ax1.plot(c_xi_values, np.real(R_values), 'b-', linewidth=2, label='Real part')
        if np.any(np.imag(R_values) != 0):
            ax1.plot(c_xi_values, np.imag(R_values), 'r--', linewidth=2, label='Imaginary part')
        ax1.set_xlabel('c·ξ')
        ax1.set_ylabel('R(ξ)')
        ax1.set_title(f'Radial Function: n={n}, m={m}, D={D}, a={a}\n' + 
                      f'ν = {v:.6f}, c = {c:.6e}')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Plot 2: R(xi) vs xi 
        ax2.plot(xi_values, np.real(R_values), 'b-', linewidth=2, label='Real part')
        if np.any(np.imag(R_values) != 0):
            ax2.plot(xi_values, np.imag(R_values), 'r--', linewidth=2, label='Imaginary part')
        ax2.set_xlabel('ξ')
        ax2.set_ylabel('R(ξ)')
        ax2.set_title(f'Radial Function vs ξ coordinate')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        
        if save_plot:
            filename = f'radial_n{n}_m{m}_D{D:.3f}_a{a:.3f}.png'
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Plot saved as: {filename}")
        
        plt.show()
        
        return R_values, v, xi_values, c
        
    except Exception as e:
        print(f"Error computing radial function: {e}")
        return None, None, None, None

def compare_radial_modes(m, D, a, n_values=[0, 1, 2], c_xi_max=15.0, num_points=200):
    """
    Compare radial functions for different n values (same m).
    
    Parameters:
    -----------
    m : int
        Azimuthal quantum number (fixed)
    D : float
        Dipole strength
    a : float  
        Semi-focal distance
    n_values : list
        List of n values to compare
    c_xi_max : float
        Maximum c*xi for plotting
    num_points : int
        Number of plotting points
    """
    
    c_xi_values = np.linspace(0.1, c_xi_max, num_points)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
    
    for i, n in enumerate(n_values):
        try:
            R_values, v, xi_values, c = compute_radial_function(n, m, D, a, c_xi_values)
            
            color = colors[i % len(colors)]
            ax.plot(c_xi_values, np.real(R_values), color=color, linewidth=2, 
                   label=f'n={n} (ν={v:.4f})')
            
            print(f"Mode n={n}: ν = {v:.6f}")
            
        except Exception as e:
            print(f"Error for n={n}: {e}")
            continue
    
    ax.set_xlabel('c·ξ')
    ax.set_ylabel('R(ξ)')
    ax.set_title(f'Radial Functions Comparison: m={m}, D={D}, a={a}')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    filename = f'radial_comparison_m{m}_D{D:.3f}_a{a:.3f}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Comparison plot saved as: {filename}")
    plt.show()

def compare_dipole_at_fixed_energy(n, m, E_fixed, a, D_values, c_xi_max=15.0, num_points=200):
    """
    Compare radial functions for different dipole strengths D at a fixed energy E.
    
    This is useful for studying how the dipole field affects the radial structure
    when the total energy is held constant.
    
    Parameters:
    -----------
    n, m : int
        Quantum numbers
    E_fixed : float
        Fixed energy value
    a : float  
        Semi-focal distance
    D_values : list
        List of dipole strength values to compare
    c_xi_max : float
        Maximum c*xi for plotting
    num_points : int
        Number of plotting points
    """
    
    c_xi_values = np.linspace(0.1, c_xi_max, num_points)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(D_values)))
    
    print(f"\nComparing dipole strengths at fixed energy E = {E_fixed:.6f}")
    print("="*60)
    
    for i, D in enumerate(D_values):
        try:
            # For fixed energy, calculate the corresponding c parameter
            c = np.sqrt(2 * E_fixed * a**2)
            xi_values = c_xi_values / c
            
            # Get angular eigenvalue for this D
            L_max = max(20, n + abs(m) + 5)
            eigvals, *_ = analytic(m, L_max, E_fixed, a, D)
            
            if n >= len(eigvals):
                print(f"D={D}: Not enough eigenvalues")
                continue
                
            lambda_mn = abs(eigvals[n])
            
            # Compute radial function using Leaver's method
            R_values, v = radial_function_leaver(xi_values, c, m, lambda_mn, n)
            
            # Plot on both axes
            ax1.plot(c_xi_values, np.real(R_values), color=colors[i], linewidth=2, 
                    label=f'D={D:.1f} (ν={v:.4f})')
            ax2.plot(xi_values, np.real(R_values), color=colors[i], linewidth=2, 
                    label=f'D={D:.1f}')
            
            print(f"D={D:.1f}: ν = {v:.6f}, λ = {lambda_mn:.6f}, c = {c:.6f}")
            
        except Exception as e:
            print(f"Error for D={D}: {e}")
            continue
    
    # Format plots
    ax1.set_xlabel('c·ξ')
    ax1.set_ylabel('R(ξ)')
    ax1.set_title(f'Radial Functions vs c·ξ: n={n}, m={m}, E={E_fixed:.4f}, a={a}')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_xlabel('ξ')
    ax2.set_ylabel('R(ξ)')
    ax2.set_title(f'Radial Functions vs ξ coordinate')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    filename = f'radial_dipole_comparison_n{n}_m{m}_E{E_fixed:.4f}_a{a:.3f}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Dipole comparison plot saved as: {filename}")
    plt.show()

def compare_dipole_strengths(n, m, a, D_values, c_xi_max=15.0, num_points=200):
    """
    Compare radial functions for different dipole strengths D.
    
    Each dipole strength will have its own energy E = D²/(8a⁴).
    This shows how both the energy scale and dipole field affect the function.
    
    Parameters:
    -----------
    n, m : int
        Quantum numbers
    a : float  
        Semi-focal distance
    D_values : list
        List of dipole strength values to compare
    c_xi_max : float
        Maximum c*xi for plotting
    num_points : int
        Number of plotting points
    """
    
    c_xi_values = np.linspace(0.1, c_xi_max, num_points)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    colors = plt.cm.plasma(np.linspace(0, 1, len(D_values)))
    
    print(f"\nComparing dipole strengths with corresponding energies")
    print("="*60)
    
    for i, D in enumerate(D_values):
        try:
            # Calculate energy for this dipole strength
            if D == 0:
                E = 0.01  # Small reference energy for D=0 case
            else:
                E = D**2 / (8 * a**4)
            
            R_values, v, xi_values, c = compute_radial_function(n, m, D, a, c_xi_values)
            
            # Plot on both axes
            ax1.plot(c_xi_values, np.real(R_values), color=colors[i], linewidth=2, 
                    label=f'D={D:.1f} (E={E:.4f}, ν={v:.4f})')
            ax2.plot(xi_values, np.real(R_values), color=colors[i], linewidth=2, 
                    label=f'D={D:.1f}')
            
            print(f"D={D:.1f}: E = {E:.6f}, ν = {v:.6f}, c = {c:.6f}")
            
        except Exception as e:
            print(f"Error for D={D}: {e}")
            continue
    
    # Format plots
    ax1.set_xlabel('c·ξ')
    ax1.set_ylabel('R(ξ)')
    ax1.set_title(f'Radial Functions vs c·ξ: n={n}, m={m}, a={a}')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_xlabel('ξ')
    ax2.set_ylabel('R(ξ)')
    ax2.set_title(f'Radial Functions vs ξ coordinate')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    filename = f'radial_dipole_strengths_n{n}_m{m}_a{a:.3f}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Dipole strengths comparison plot saved as: {filename}")
    plt.show()

def test_radial_functions():
    """
    Test the radial function calculation with various parameters.
    """
    print("="*70)
    print("TESTING RADIAL FUNCTIONS")
    print("="*70)
    
    # Test parameters
    test_cases = [
        # (n, m, D, a, description)
        (0, 0, 0.0, 0.1, "Ground state, no dipole"),
        (1, 0, 0.0, 0.1, "First excited state, no dipole"),
        (0, 1, 0.0, 0.1, "m=1 ground state, no dipole"),
        (0, 0, 0.01, 0.1, "Ground state, weak dipole"),
        (1, 0, 0.01, 0.1, "First excited state, weak dipole"),
    ]
    
    for n, m, D, a, description in test_cases:
        print(f"\nTesting: {description}")
        print(f"Parameters: n={n}, m={m}, D={D}, a={a}")
        
        try:
            # Compute for a small range of c*xi values
            c_xi_test = np.array([0.5, 1.0, 2.0, 5.0])
            R_values, v, xi_values, c = compute_radial_function(n, m, D, a, c_xi_test)
            
            print(f"  ν = {v:.6f}")
            print(f"  c = {c:.6e}")
            print("  Sample values:")
            for i, (cxi, R) in enumerate(zip(c_xi_test, R_values)):
                print(f"    c·ξ = {cxi:.1f}: R = {R:.6e}")
            
            print("  ✓ SUCCESS")
            
        except Exception as e:
            print(f"  ✗ ERROR: {e}")

if __name__ == "__main__":
    # Run tests first
    test_radial_functions()
    
    print("\n" + "="*70)
    print("GENERATING SAMPLE PLOTS")
    print("="*70)
    
    # Generate some sample plots
    try:
        # Plot ground state
        print("\n1. Ground state (n=0, m=0) with no dipole:")
        plot_radial_function(0, 0, 0.0, 0.1, c_xi_max=10.0)
        
        # Plot first excited state
        print("\n2. First excited state (n=1, m=0) with no dipole:")
        plot_radial_function(1, 0, 0.0, 0.1, c_xi_max=10.0)
        
        # Compare different n values
        print("\n3. Comparison of radial modes (m=0):")
        compare_radial_modes(0, 0.0, 0.1, n_values=[0, 1, 2])
        
        print("\n4. Ground state with weak dipole field:")
        plot_radial_function(0, 0, 0.01, 0.1, c_xi_max=8.0)
        
        print("\n5. Dipole strength comparison (natural energies):")
        compare_dipole_strengths(0, 0, 0.1, [0, 0.1, 0.2, 0.5])
        
        print("\n6. Dipole strength comparison (fixed energy E=0.1):")
        compare_dipole_at_fixed_energy(0, 0, 0.1, 0.1, [0, 0.1, 0.2, 0.5])
        
    except Exception as e:
        print(f"Error in plotting: {e}")
    
    print("\nRadial function calculations complete!")
    print("Use plot_radial_function(n, m, D, a) to generate custom plots.")
    print("Use compare_radial_modes(m, D, a, n_values) to compare different n modes.")
    print("Use compare_dipole_at_fixed_energy(n, m, E_fixed, a, D_values) to compare")
    print("  radial functions at fixed energy but varying dipole strength.")
    print("Use compare_dipole_strengths(n, m, a, D_values) to compare radial functions")
    print("  for different dipole strengths with corresponding energies.")
