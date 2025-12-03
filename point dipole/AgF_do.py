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
R_Ag = np.array([0, 0, 0.872093 / 52.9e-2])
R_F = np.array([0, 0, -1.077907 / 52.9e-2])
def recenter_DO(DO, x, y, z, dV):
    density = np.abs(DO)**2
    norm = np.sum(density) * dV
    x_c = np.sum(density * x) * dV / norm
    y_c = np.sum(density * y) * dV / norm
    z_c = np.sum(density * z) * dV / norm
    return np.array([x_c, y_c, z_c])


alpha_1S_Ag = np.array([4.74452163e+3, 8.64220538e+2, 2.33891805e+2])
coeffs_1S_Ag = np.array([1.54328970e-1, 5.35328140-1, 4.44634540e-1])
alpha_2SP_Ag = np.array([4.14965207e+2, 9.64289900e+1, 3.13617003e+1])
coeffs_2S_Ag = np.array([-9.99672300e-2, 3.99512830e-1, 7.001154702-1])
coeffs_2P_Ag = np.array([1.55916280e-1, 6.07683720e-1, 3.91957390e-1])
alpha_3SP_Ag = np.array([4.94104861e+1, 1.50717731e+1, 5.81515863])
coeffs_3S_Ag = np.array([-2.27763500e-1, 2.17543600e-1, 9.16676960e-1])
coeffs_3P_Ag = np.array([4.95151000e-3, 5.77766470e-1, 4.84646040e-1])
alpha_4D_Ag = np.array([4.94104861e+1, 1.50717731e+1, 5.81515863e+1])
coeffs_4D_Ag = np.array([2.19767950e-1, 6.55547360e-1, 2.86573260e-1])
alpha_5SP_Ag = np.array([5.29023045, 2.05998832,  9.06811930e-1])
coeffs_5S_Ag = np.array([-3.30610060e-1, 5.76109500e-2, 1.15578745])
coeffs_5P_Ag = np.array([-1.28392760e-1, 5.85204760e-1, 5.43944200e-1])
alpha_6D_Ag = np.array([3.28339567, 1.27853725, 5.62815250e-1])
coeffs_6D_Ag = np.array([1.25066210e-1, 6.68678560e-1, 3.05246820e-1])
alpha_7SP_Ag = np.array([4.37080480e-1, 2.35340820e-1, 1.03954180e-1])
coeffs_7S_Ag = np.array([-3.8426426e-1, -1.97256740e-1, 1.37549551])
coeffs_7P_Ag = np.array([-3.48169150e-1, 6.29032370e-1, 6.66283270e-1])
alpha_1S_F = np.array([1.66679130e+2, 3.03608120e+1, 8.21682070])
coeffs_1S_F = np.array([1.54328970e-1, 5.35328140e-1, 4.44634540e-1])
alpha_2SP_F = np.array([6.46480320, 1.50228120, 4.88588500e-1])
coeffs_2S_F = np.array([-9.99672300e-2, 3.99512830e-1, 7.00115470e-1])
coeffs_2P_F = np.array([1.55916270e-1, 6.07683720e-1, 3.919573902-1])

DO_coeffs_AgF = np.array([9.14529141e-03, -2.69537956e-02, 0, 0, 3.02739714e-02, 5.83796178e-02, 0, 0, -8.04092502e-02, 0, 0, -1.50307315e-01, 0, 0, -1.17098135e-01, 0, 0, 1.90945673e-01, 0, 0, 5.51622500e-01, 0, 0, 5.06467116e-01, 0, 0, -5.96318801e-01, -5.97134365e-03, 7.04717812e-02, 0, 0, 3.39488329e-01])

def build_DO(DO_coeffs, x, y, z, dV, recenter=True, print_geom=True):
    global R_Ag, R_F
    basis_info = [
        # Format: (alphas, coeffs, a, b, c, center, DO_coeff index)
        (alpha_1S_Ag,  coeffs_1S_Ag,  0, 0, 0, R_Ag, 0),
        (alpha_2SP_Ag, coeffs_2S_Ag,  0, 0, 0, R_Ag, 1),
        (alpha_2SP_Ag, coeffs_2P_Ag,  1, 0, 0, R_Ag, 2),
        (alpha_2SP_Ag, coeffs_2P_Ag,  0, 1, 0, R_Ag, 3),
        (alpha_2SP_Ag, coeffs_2P_Ag,  0, 0, 1, R_Ag, 4),
        (alpha_3SP_Ag, coeffs_3S_Ag,  0, 0, 0, R_Ag, 5),
        (alpha_3SP_Ag, coeffs_3P_Ag,  1, 0, 0, R_Ag, 6),
        (alpha_3SP_Ag, coeffs_3P_Ag,  0, 1, 0, R_Ag, 7),
        (alpha_3SP_Ag, coeffs_3P_Ag,  0, 0, 1, R_Ag, 8),
        (alpha_4D_Ag, coeffs_4D_Ag, 1, 1, 0, R_Ag, 9),
        (alpha_4D_Ag, coeffs_4D_Ag, 0, 1, 1, R_Ag, 10),
        #d z^2 done manually
        (alpha_4D_Ag, coeffs_4D_Ag, 1, 0, 1, R_Ag, 12),
        #d x^2-y^2 done manually
        (alpha_5SP_Ag, coeffs_5S_Ag,  0, 0, 0, R_Ag, 14),
        (alpha_5SP_Ag, coeffs_5P_Ag,  1, 0, 0, R_Ag, 15),
        (alpha_5SP_Ag, coeffs_5P_Ag,  0, 1, 0, R_Ag, 16),
        (alpha_5SP_Ag, coeffs_5P_Ag,  0, 0, 1, R_Ag, 17),
        (alpha_6D_Ag, coeffs_6D_Ag, 1, 1, 0, R_Ag, 18),
        (alpha_6D_Ag, coeffs_6D_Ag, 0, 1, 1, R_Ag, 19),
        #d z^2 done manually
        (alpha_6D_Ag, coeffs_6D_Ag, 1, 0, 1, R_Ag, 21),
        #d x^2 - y^2 done manually
        (alpha_7SP_Ag, coeffs_7S_Ag,  0, 0, 0, R_Ag, 23),
        (alpha_7SP_Ag, coeffs_7P_Ag,  1, 0, 0, R_Ag, 24),
        (alpha_7SP_Ag, coeffs_7P_Ag,  0, 1, 0, R_Ag, 25),
        (alpha_7SP_Ag, coeffs_7P_Ag,  0, 0, 1, R_Ag, 26),
        (alpha_1S_F,  coeffs_1S_F,  0, 0, 0, R_F, 27),
        (alpha_2SP_F, coeffs_2S_F,  0, 0, 0, R_F, 28),
        (alpha_2SP_F, coeffs_2P_F,  1, 0, 0, R_F, 29),
        (alpha_2SP_F, coeffs_2P_F,  0, 1, 0, R_F, 30),
        (alpha_2SP_F, coeffs_2P_F,  0, 0, 1, R_F, 31),
    ]

    DO = np.zeros_like(x)
    for alphas, coeffs, a, b, c, center, i in basis_info:
        if DO_coeffs[i] != 0.0:  # Skip zero coefficients
            ao = AO(alphas, coeffs, a, b, c, center, x, y, z, dV)
            DO += DO_coeffs[i] * ao

    # Handle weird d orbitals manually
    # d_z2 = 1/2 * [2zz - xx - yy]
    if DO_coeffs[11] != 0.0:
        ao_zz = AO(alpha_4D_Ag, coeffs_4D_Ag, 0, 0, 2, R_Ag, x, y, z, dV)
        ao_xx = AO(alpha_4D_Ag, coeffs_4D_Ag, 2, 0, 0, R_Ag, x, y, z, dV)
        ao_yy = AO(alpha_4D_Ag, coeffs_4D_Ag, 0, 2, 0, R_Ag, x, y, z, dV)
        DO += DO_coeffs[11] * 0.5 * (2 * ao_zz - ao_xx - ao_yy)

    # d_x2_y2 = sqrt(3)/2 * (xx - yy)
    if DO_coeffs[13] != 0.0:
        ao_xx = AO(alpha_4D_Ag, coeffs_4D_Ag, 2, 0, 0, R_Ag, x, y, z, dV)
        ao_yy = AO(alpha_4D_Ag, coeffs_4D_Ag, 0, 2, 0, R_Ag, x, y, z, dV)
        DO += DO_coeffs[13] * (np.sqrt(3)/2) * (ao_xx - ao_yy)

    if DO_coeffs[20] != 0.0:
        ao_zz = AO(alpha_6D_Ag, coeffs_6D_Ag, 0, 0, 2, R_Ag, x, y, z, dV)
        ao_xx = AO(alpha_6D_Ag, coeffs_6D_Ag, 2, 0, 0, R_Ag, x, y, z, dV)
        ao_yy = AO(alpha_6D_Ag, coeffs_6D_Ag, 0, 2, 0, R_Ag, x, y, z, dV)
        DO += DO_coeffs[20] * 0.5 * (2 * ao_zz - ao_xx - ao_yy)

    if DO_coeffs[22] != 0.0:
        ao_xx = AO(alpha_6D_Ag, coeffs_6D_Ag, 2, 0, 0, R_Ag, x, y, z, dV)
        ao_yy = AO(alpha_6D_Ag, coeffs_6D_Ag, 0, 2, 0, R_Ag, x, y, z, dV)
        DO += DO_coeffs[22] * (np.sqrt(3)/2) * (ao_xx - ao_yy)


        DO /= np.sqrt(integrate_3d(abs(DO)**2, dV))

    if recenter:
        centroid = recenter_DO(DO, x, y, z, dV)
        print("Dyson orbital centroid: {:.6e} {:.6e} {:.6e}".format(*centroid))
        print("Origin of the expansion: {:.6f} {:.6f} {:.6f}".format(*centroid))

        # Shift grid
        x_new = x - centroid[0]
        y_new = y - centroid[1]
        z_new = z - centroid[2]

        # Shift molecular geometry
        R_Ag = R_Ag - centroid
        R_F = R_F - centroid

        if print_geom:
            bohr_to_angstrom = 0.529177
            R_A_ang = R_Ag * bohr_to_angstrom
            R_B_ang = R_F * bohr_to_angstrom

            print("\nShifting molecular geometry to the new center")
            print("New molecular geometry is (in Ångströms):")
            print(" atom         X             Y             Z")
            print("   1      {: .6f}      {: .6f}      {: .6f}".format(*R_A_ang))
            print("   2      {: .6f}      {: .6f}      {: .6f}".format(*R_B_ang))

        # Rebuild on shifted grid
        return build_DO(DO_coeffs, x_new, y_new, z_new, dV, recenter=False)

    return DO


