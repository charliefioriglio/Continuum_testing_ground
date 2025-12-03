import numpy as np

def integrate_3d(f, dV):
    return np.sum(f) * dV

# ---------------------------------------------
# Double factorial
# ---------------------------------------------
def double_factorial(n):
    if n <= 0:
        return 1
    else:
        return n * double_factorial(n - 2)

# ---------------------------------------------
# Normalization factor for Cartesian Gaussians
# ---------------------------------------------
def norm_cartesian_gaussian(alpha, a, b, c):
    l = a + b + c
    prefactor = (2 * alpha / np.pi)**(3/4)
    numerator = (4 * alpha)**l
    denom = double_factorial(2*a - 1) * double_factorial(2*b - 1) * double_factorial(2*c - 1)
    return prefactor * np.sqrt(numerator / denom)

# ---------------------------------------------
# Primitive Gaussian with given (a, b, c)
# ---------------------------------------------
def gaussian_primitive(alpha, x, y, z, a, b, c, center):
    x0, y0, z0 = center
    xs, ys, zs = x-x0, y-y0, z-z0
    r2 = xs**2 + ys**2 + zs**2
    norm = norm_cartesian_gaussian(alpha, a, b, c)
    return norm * (xs**a) * (ys**b) * (zs**c) * np.exp(-alpha * r2)

# ---------------------------------------------
# Build AOs
# ---------------------------------------------
def AO_norm(primatives, coeffs, dV):
    norm = 0.0
    n = len(coeffs)
    for i in range(n):
        for j in range(n):
            overlap = np.sum(primatives[i] * primatives[j]) * dV
            norm += coeffs[i] * coeffs[j] * overlap
    return 1.0 / np.sqrt(norm)

def AO(alphas, coeffs, a, b, c, center, x, y, z, dV):
    primitives = [gaussian_primitive(alpha, x, y, z, a, b, c, center) for alpha in alphas]
    AO = sum(c * p for c, p in zip(coeffs, primitives))
    AO_norm_const = AO_norm(primitives, coeffs, dV)
    AO_normalized = AO_norm_const * AO
    return AO_normalized

# ---------------------------------------------
# Build DO
# ---------------------------------------------
bl_A = 1.35  # Bond length in angstroms
bl_au = bl_A / 52.9e-2  # bond length in au
R_A = np.array([0, 0, (bl_au / 2)])
R_B = np.array([0, 0, -(bl_au / 2)])

alpha_1S = np.array([1.30709320e+2, 2.38088610e+1, 6.44360830])
coeffs_1S = np.array([1.54328970e-1, 5.35328140e-1, 4.44634540e-1])
alpha_2SP = np.array([5.03315130, 1.16959610, 3.80389000e-1])
coeffs_2S = np.array([-9.99672300e-2, 3.99512830e-1, 7.00115470e-1])
coeffs_2P  = np.array([1.55916270e-1, 6.07683720e-1, 3.9195739e-1])

DO_coeffs_b3g = np.array([0, 0, 0, -7.47775656e-1, 0, 0, 0, 0, 7.47775656e-1, 0])
DO_coeffs_b2g = np.array([0, 0, -7.47775656e-1, 0, 0, 0, 0, 7.47775656e-1, 0, 0])

def build_DO(DO_coeffs, x, y, z, dV):
    basis_info = [
        # Format: (alphas, coeffs, a, b, c, center, DO_coeff index)
        (alpha_1S,  coeffs_1S,  0, 0, 0, R_A, 0),
        (alpha_2SP, coeffs_2S,  0, 0, 0, R_A, 1),
        (alpha_2SP, coeffs_2P,  1, 0, 0, R_A, 2),
        (alpha_2SP, coeffs_2P,  0, 1, 0, R_A, 3),
        (alpha_2SP, coeffs_2P,  0, 0, 1, R_A, 4),
        (alpha_1S,  coeffs_1S,  0, 0, 0, R_B, 5),
        (alpha_2SP, coeffs_2S,  0, 0, 0, R_B, 6),
        (alpha_2SP, coeffs_2P,  1, 0, 0, R_B, 7),
        (alpha_2SP, coeffs_2P,  0, 1, 0, R_B, 8),
        (alpha_2SP, coeffs_2P,  0, 0, 1, R_B, 9),
    ]

    DO = np.zeros_like(x)
    for alphas, coeffs, a, b, c, center, i in basis_info:
        if DO_coeffs[i] != 0.0:  # Skip zero coefficients
            ao = AO(alphas, coeffs, a, b, c, center, x, y, z, dV)
            DO += DO_coeffs[i] * ao

    DO /= np.sqrt(integrate_3d(abs(DO)**2, dV))

    return DO
