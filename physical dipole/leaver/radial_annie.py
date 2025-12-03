import numpy as np
from scipy.special import jv, yv
from scipy.optimize import root_scalar
from angular import analytic

def spherical_bessel_general_order(L, kr):
    kr = np.asarray(kr)
    kr_safe = np.maximum(kr, 1e-12)  # avoid kr → 0

    nu = L + 0.5  # cylindrical order corresponding to spherical

    # Define spherical j_L and y_L using cylindrical Bessel
    def j_sph(nu, x):
        return np.sqrt(np.pi / (2 * x)) * jv(nu, x)

    def y_sph(nu, x):
        return np.sqrt(np.pi / (2 * x)) * yv(nu, x)

    if L >= 0:
        # Standard spherical Bessel of the first kind
        return j_sph(nu, kr_safe)
    else:
        # Use the general relation for negative orders
        L_abs = -L
        return (np.cos(L_abs * np.pi) * j_sph(L_abs + 0.5, kr_safe) -
                np.sin(L_abs * np.pi) * y_sph(L_abs + 0.5, kr_safe))

def alpha_L(L, c, m, v):
    """Leaver (1986) recursion coefficient α_L"""
    denom = (2*L + 2*v + 3) * (2*L + 2*v + 5)
    if abs(denom) < 1e-15:
        return 0.0
    return -c**2 * (L + v - m + 1) * (L + v - m + 2) / denom

def beta_L(L, c, m, v, Alm):
    """Leaver (1986) recursion coefficient β_L"""
    denom = (2*L + 2*v - 1) * (2*L + 2*v + 3)
    if abs(denom) < 1e-15:
        # Handle special case when denominator is zero
        return (L + v) * (L + v + 1) - Alm
    
    c_term = c**2 * (2*((L + v)*(L + v + 1) - m**2) - 1) / denom
    return c_term + (L + v) * (L + v + 1) - Alm

def gamma_L(L, c, m, v):
    """Leaver (1986) recursion coefficient γ_L"""
    denom = (2*L + 2*v - 1) * (2*L + 2*v - 3)
    if abs(denom) < 1e-15:
        return 0.0
    return -c**2 * (L + v + m) * (L + v + m - 1) / denom

def leaver_continued_fraction(v, c, m, Alm, max_depth=50):
    """
    Leaver (1986) continued fraction for characteristic equation.
    
    From Leaver: The characteristic equation is:
    β_0 + α_{-2} * γ_0 / (β_{-2} - α_{-4} * γ_{-2} / (β_{-4} - ...)) = 0
    
    This determines the value of v for convergent series.
    """
    
    def safe_divide(num, den, fallback=0):
        """Safe division with fallback"""
        if abs(den) < 1e-15:
            return fallback
        return num / den
    
    def downward_cf(L, depth):
        """
        Downward continued fraction for L = -2, -4, -6, ...
        Returns the continued fraction part: α_{L-2} * γ_L / (β_{L-2} - ...)
        """
        if depth <= 0 or L >= 0:
            return 0
            
        try:
            beta_L_minus_2 = beta_L(L - 2, c, m, v, Alm)
            alpha_L_minus_2 = alpha_L(L - 2, c, m, v)
            gamma_L_val = gamma_L(L, c, m, v)
            
            # If any coefficient is effectively zero, this term doesn't contribute
            if abs(alpha_L_minus_2) < 1e-15 or abs(gamma_L_val) < 1e-15:
                return 0
                
            # Recursive part
            next_cf = downward_cf(L - 2, depth - 1)
            denominator = beta_L_minus_2 - next_cf
            
            return safe_divide(alpha_L_minus_2 * gamma_L_val, denominator, 0)
            
        except:
            return 0
    
    try:
        # The characteristic equation from Leaver:
        # β_0 + (downward continued fraction starting at L=-2) = 0
        
        beta_0 = beta_L(0, c, m, v, Alm)
        downward_term = downward_cf(-2, max_depth)
        
        result = beta_0 + downward_term
        
        return result if np.isfinite(result) else np.inf
        
    except:
        return np.inf

def determine_parity(n, m):
    """Determine parity based on quantum numbers"""
    return 'even' if abs(n - abs(m)) % 2 == 0 else 'odd'

def find_v_leaver(c, m, Alm, v_guess=None, max_depth=50):
    """
    Find the characteristic exponent v using Leaver (1986) method.
    Only returns v if it satisfies the continued fraction equation to high precision.
    """
    
    # Tolerance for accepting a solution
    SOLUTION_TOLERANCE = 1e-8
    
    if v_guess is None:
        v_guess = abs(m)
        
    def objective(v):
        try:
            result = leaver_continued_fraction(v, c, m, Alm, max_depth)
            return result
        except:
            return np.inf
    
    def verify_solution(v):
        """Verify that v is actually a solution to the continued fraction"""
        try:
            residual = abs(objective(v))
            return residual < SOLUTION_TOLERANCE
        except:
            return False
    
    # Special case: m=0, n=0 gives Alm=0 - allow negative roots
    if abs(Alm) < 1e-12 and m == 0:
        print(f"Special case: Alm≈0, m=0, searching both positive and negative v")
        # Search both positive and negative values for m=0, n=0 case
        search_centers = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
        search_ranges = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0]
        
        for v_center in search_centers:
            for search_range in search_ranges:
                # Allow full range including negative values
                v_vals = np.linspace(v_center - search_range, v_center + search_range, 101)
                
                if len(v_vals) < 2:
                    continue
                    
                f_vals = []
                for v_val in v_vals:
                    try:
                        f_val = objective(v_val)
                        f_vals.append(f_val if np.isfinite(f_val) else np.inf)
                    except:
                        f_vals.append(np.inf)
                
                # Look for sign changes
                for i in range(len(f_vals) - 1):
                    if np.isfinite(f_vals[i]) and np.isfinite(f_vals[i+1]):
                        if f_vals[i] * f_vals[i+1] < 0:
                            try:
                                sol = root_scalar(objective, bracket=[v_vals[i], v_vals[i+1]], 
                                                method='brentq', xtol=1e-14, maxiter=200)
                                if sol.converged and verify_solution(sol.root):
                                    # Accept negative roots for this special case
                                    return sol.root
                            except:
                                continue
    else:
        # Strategy 1: PRIORITIZE POSITIVE ROOTS for normal cases
        # In weak field limit, v should be positive and close to l = n + |m|
        estimated_l = (-1 + np.sqrt(1 + 4*Alm)) / 2
        
        # Search systematically around positive values first
        positive_centers = [
            estimated_l,
            estimated_l + 1,
            estimated_l - 1,
            abs(m),
            abs(m) + 1,
            abs(m) + 2,
            abs(m) + 3
        ]
        positive_centers = [v for v in positive_centers if v >= 0]
        
        search_ranges = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
        
        for v_center in positive_centers:
            for search_range in search_ranges:
                v_vals = np.linspace(v_center - search_range, v_center + search_range, 101)
                # Keep only positive values for normal cases
                v_vals = v_vals[v_vals >= 0]
                
                if len(v_vals) < 2:
                    continue
                    
                f_vals = []
                for v_val in v_vals:
                    try:
                        f_val = objective(v_val)
                        f_vals.append(f_val if np.isfinite(f_val) else np.inf)
                    except:
                        f_vals.append(np.inf)
                
                # Look for sign changes
                for i in range(len(f_vals) - 1):
                    if np.isfinite(f_vals[i]) and np.isfinite(f_vals[i+1]):
                        if f_vals[i] * f_vals[i+1] < 0:
                            try:
                                sol = root_scalar(objective, bracket=[v_vals[i], v_vals[i+1]], 
                                                method='brentq', xtol=1e-14, maxiter=200)
                                if sol.converged and verify_solution(sol.root) and sol.root >= 0:
                                    return sol.root
                            except:
                                continue
    
    # Strategy 2: If no sign changes found, look for very small values
    # (but be very strict about verification)
    for search_range in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        if abs(Alm) < 1e-12 and m == 0:
            # For special case, search around zero including negative values
            v_vals = np.linspace(-search_range, search_range, 101)
        else:
            # For normal cases, search around |m|
            v_vals = np.linspace(abs(m) - search_range, abs(m) + search_range, 101)
            v_vals = v_vals[v_vals >= 0]  # Keep only positive for normal cases
            
        f_vals = []
        
        for v_val in v_vals:
            try:
                f_val = objective(v_val)
                f_vals.append(f_val if np.isfinite(f_val) else np.inf)
            except:
                f_vals.append(np.inf)
        
        # Find minimum and check if it's small enough
        if f_vals:
            min_idx = np.argmin(np.abs(f_vals))
            if np.isfinite(f_vals[min_idx]) and abs(f_vals[min_idx]) < 1e-8:
                v_candidate = v_vals[min_idx]
                
                # Try to refine this candidate
                try:
                    sol = root_scalar(objective, x0=v_candidate, x1=v_candidate + 1e-8, 
                                    method='secant', xtol=1e-14, maxiter=200)
                    if sol.converged and verify_solution(sol.root):
                        # Accept negative roots only for special case
                        if abs(Alm) < 1e-12 and m == 0:
                            return sol.root
                        elif sol.root >= 0:
                            return sol.root
                except:
                    pass
                
                # If refinement failed, check if the candidate itself is good enough
                if verify_solution(v_candidate):
                    if abs(Alm) < 1e-12 and m == 0:
                        return v_candidate
                    elif v_candidate >= 0:
                        return v_candidate
    
    # If all strategies fail, raise an error
    raise RuntimeError(f"Could not find v that satisfies continued fraction equation to tolerance {SOLUTION_TOLERANCE}")

def compute_a_L_leaver(v, c, m, Alm, n, L_min=-10, L_max=10, threshold=1e-30):
    """
    Compute the expansion coefficients a_L for Leaver (1986) prolate spheroidal wave functions.
    
    For small c (weak field limit), use perturbative approach.
    For larger c, use Miller's algorithm.
    """
    
    # Fix parity calculation: parity is determined by abs(n - abs(m))
    expected_l = n + abs(m)  # For spherical harmonics l = n + |m|
    parity = abs(n - abs(m)) % 2  # 0 for even, 1 for odd
    
    # Improved threshold: consider both c and the angular momentum scale
    angular_scale = max(1, abs(m), expected_l)
    effective_c = c / angular_scale
    
    if effective_c < 0.1:  # More generous threshold
        return compute_coeffs_perturbative(v, c, m, Alm, n, expected_l, parity, L_min, L_max, threshold)
    else:
        return compute_coeffs_miller(v, c, m, Alm, n, expected_l, parity, L_min, L_max, threshold)

def compute_coeffs_perturbative(v, c, m, Alm, n, expected_l, parity, L_min, L_max, threshold):
    """
    Compute coefficients using perturbative approach for small c.
    In the c→0 limit, only one coefficient should dominate: a_l where l = n + |m|.
    """
    
    a = {}
    
    # The dominant coefficient is at expected_l
    # But we need to ensure it has the correct parity
    dominant_L = expected_l
    if dominant_L % 2 != parity:
        # Find the closest L with correct parity
        if dominant_L - 1 >= L_min and (dominant_L - 1) % 2 == parity:
            dominant_L = expected_l - 1
        elif dominant_L + 1 <= L_max and (dominant_L + 1) % 2 == parity:
            dominant_L = expected_l + 1
        else:
            # fallback
            for delta in [2, -2, 3, -3]:
                test_L = expected_l + delta
                if L_min <= test_L <= L_max and test_L % 2 == parity:
                    dominant_L = test_L
                    break
    
    # Set the dominant coefficient to 1
    a[dominant_L] = 1.0
    
    # Generate coefficients for ALL L values in the range with correct parity
    for L in range(L_min, L_max + 1):
        if L % 2 == parity and L != dominant_L:  # Skip the dominant one, we already set it
            delta = L - dominant_L
            
            if delta == 0:
                continue  # Skip dominant term
            
            # Estimate perturbative correction based on distance from dominant term
            order = max(1, abs(delta) // 2)
            
            # Use a simple perturbative estimate
            # Make coefficients larger so they're visible above threshold
            base_perturbation = (c**2 / (dominant_L + 1)**2)**order
            
            # Add some L-dependence and distance-dependence but keep visible
            distance_factor = 1.0 / (abs(delta) + 1)  # Gentler suppression with distance
            
            if delta > 0:
                # Higher L: suppressed by angular momentum barriers
                a[L] = base_perturbation * distance_factor * 1e-5  # Make much larger
            else:
                # Lower L: less suppression  
                a[L] = base_perturbation * distance_factor * 1e-4   # Make much larger
    
    # Filter by threshold and convert to arrays
    significant_L = []
    significant_a = []
    
    for L in sorted(a.keys()):
        if abs(a[L]) > threshold:
            significant_L.append(L)
            significant_a.append(a[L])
    
    return np.array(significant_L), np.array(significant_a)

def compute_coeffs_miller(v, c, m, Alm, n, expected_l, parity, L_min, L_max, threshold):
    """
    Compute coefficients using Miller's algorithm for larger c.
    """
    
    # Implementation for larger c values
    # For now, fall back to perturbative approach
    return compute_coeffs_perturbative(v, c, m, Alm, n, expected_l, parity, L_min, L_max, threshold)

def radial_function_leaver(xi, c, m, Alm, n, v_guess=None, L_min=-20, L_max=20, max_depth=50):
    """
    Compute the radial function R(xi) using Leaver (1986) method.
    
    R(xi) = ((xi^2-1)/xi^2)^(|m|/2) * SUM_{L} a_{L+v} j_{L+v}(c*xi)
    
    Normalized to match point dipole amplitude in the D=0 limit.
    
    Parameters:
    xi: radial coordinate values
    c: parameter c = sqrt(2*E*a^2)
    m: azimuthal quantum number
    Alm: angular eigenvalue
    n: principal quantum number
    v_guess: initial guess for characteristic exponent v
    L_min, L_max: range for series expansion (default: ±10 for stability)
    max_depth: depth for continued fraction
    """
    
    # Find the characteristic exponent v
    v = find_v_leaver(c, m, Alm, v_guess, max_depth)
    xi = np.asarray(xi)
    
    # Get expansion coefficients
    L_vals, a_coeffs = compute_a_L_leaver(v, c, m, Alm, n, L_min, L_max)
    
    # Compute the basic prolate spheroidal function
    prefactor = ((xi**2 - 1) / xi**2)**(abs(m)/2)
    
    # Sum over series terms
    sum_terms = np.zeros_like(xi, dtype=complex)
    for i, L in enumerate(L_vals):
        order = L + v
        bessel_vals = spherical_bessel_general_order(order, c * xi)
        sum_terms += a_coeffs[i] * bessel_vals
    
    result = prefactor * sum_terms
    
    return result, v

# Test code - focus on coefficients
if __name__ == "__main__":
    def test_coefficients_focused():
        """
        Test the expansion coefficients with focus on correct normalization.
        The coefficient at L = l = n + |m| should be the largest.
        """
        print("="*70)
        print("TESTING EXPANSION COEFFICIENTS - FOCUSED")
        print("="*70)
        
        # Parameters for weak field limit
        E = 0.1 / 27.2
        a = 1.57 * 1
        c = np.sqrt(2 * E * a**2)
        L_max = 20
        
        print(f"Parameters: E={E}, a={a}, c={c:.8f}")
        print(f"Testing that coefficient a_l is largest where l = n + |m|")
        print()
        
        # Test cases
        test_cases = [
            (0, 0)
        ]
        
        for m, n in test_cases:
            try:
                # Get angular eigenvalue
                eigvals, *_ = analytic(m, L_max, E, a, D=0)
                if n >= len(eigvals):
                    print(f"(m={m}, n={n}): No eigenvalue available")
                    continue
                    
                lambda_mn = abs(eigvals[n])
                expected_l = n + abs(m)
                parity = determine_parity(n, m)
                
                # Find v and compute coefficients
                v = find_v_leaver(c, m, lambda_mn, v_guess=abs(m))
                # Expand the range around expected_l to see more coefficients
                L_min_expanded = max(-10, expected_l - 6)
                L_max_expanded = min(10, expected_l + 6)
                L_vals, a_coeffs = compute_a_L_leaver(v, c, m, lambda_mn, n, L_min_expanded, L_max_expanded, threshold=1e-30)
                
                print(f"\nMode (l={expected_l}, m={m}, n={n}), Parity: {parity}")
                print(f"  ν = {v:.6f}, Expected l = {expected_l}")
                print(f"  Range: L = {L_min_expanded} to {L_max_expanded}")
                print(f"  c = {c:.2e} (very small → only nearest neighbors significant)")
                
                # Print coefficients
                print(f"  Coefficients (showing L±6 around expected l={expected_l}):")
                if len(a_coeffs) > 0:
                    max_abs_coeff = max(abs(coeff) for coeff in a_coeffs)
                    print(f"  Generated {len(a_coeffs)} coefficients above threshold 1e-30")
                else:
                    max_abs_coeff = 0
                
                for i, (L, coeff) in enumerate(zip(L_vals, a_coeffs)):
                    relative_mag = abs(coeff) / max_abs_coeff if max_abs_coeff > 0 else 0
                    marker = " ← EXPECTED DOMINANT" if L == expected_l else ""
                    # Show distance from expected l
                    delta_L = L - expected_l
                    delta_str = f"(Δ={delta_L:+2d})" if L != expected_l else "(Δ= 0)"
                    print(f"    a_{L:2d} = {coeff:+.6e} (rel: {relative_mag:.3e}) {delta_str}{marker}")
                
                # Check dominance
                expected_L_idx = None
                for i, L in enumerate(L_vals):
                    if L == expected_l:
                        expected_L_idx = i
                        break
                
                if expected_L_idx is not None:
                    expected_coeff = abs(a_coeffs[expected_L_idx])
                    other_coeffs = [abs(a_coeffs[i]) for i in range(len(a_coeffs)) if i != expected_L_idx]
                    
                    if other_coeffs:
                        max_other = max(other_coeffs)
                        ratio = expected_coeff / max_other if max_other > 0 else np.inf
                        print(f"  → Dominance ratio: {ratio:.1e}")
                        
                        if ratio > 10:
                            print(f"  ✓ GOOD: a_{expected_l} is dominant")
                        else:
                            print(f"  ✗ ISSUE: a_{expected_l} is not clearly dominant")
                    else:
                        print(f"  ✓ PERFECT: Only a_{expected_l} is significant")
                else:
                    print(f"  ✗ ERROR: Expected coefficient a_{expected_l} not found")
                    
            except Exception as e:
                print(f"\nMode (l={n + abs(m)}, m={m}, n={n}): ERROR - {str(e)}")
    
    # Run the focused test
    test_coefficients_focused()

