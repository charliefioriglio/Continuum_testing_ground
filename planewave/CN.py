import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------
# Define Grid
# ---------------------------------------------
L = 10
n_pts = 100
x = np.linspace(-L, L, n_pts)
y = np.linspace(-L, L, n_pts)
z = np.linspace(-L, L, n_pts)
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]
dV = dx * dy * dz
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

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
def AO_norm(primitives, coeffs, dV):
    norm = 0.0
    n = len(coeffs)
    for i in range(n):
        for j in range(n):
            overlap = np.sum(primitives[i] * primitives[j]) * dV
            norm += coeffs[i] * coeffs[j] * overlap
    return 1.0 / np.sqrt(norm)

def AO(alphas, coeffs, a, b, c, center, x, y, z, dV):
    primitives = [gaussian_primitive(alpha, x, y, z, a, b, c, center) for alpha in alphas]
    ao = sum(c * p for c, p in zip(coeffs, primitives))
    ao_norm_const = AO_norm(primitives, coeffs, dV)
    return ao_norm_const * ao

# ---------------------------------------------
# Build DO
# ---------------------------------------------
bl_A = 1.16  # Bond length in angstroms
bl_au = bl_A / 52.9e-2  # bond length in au
a = bl_au / 2
R_C = np.array([0, 0, a])
R_N = np.array([0, 0, -a])

def recenter_DO(DO, x, y, z, dV):
    density = np.abs(DO)**2
    x_c = np.sum(density * x) * dV
    y_c = np.sum(density * y) * dV
    z_c = np.sum(density * z) * dV
    return np.array([x_c, y_c, z_c])

# --- C S shells ---
c_alpha_1S = np.array([8236.0, 1235.0, 280.8, 79.27, 25.59, 8.997, 3.319, 0.3643])
c_coeffs_1S = np.array([5.31e-4, 4.108e-3, 2.1087e-2, 8.1853e-2, 0.234817, 0.434401, 0.346129, -8.983e-3])

c_alpha_2S = np.array([8236.0, 1235.0, 280.8, 79.27, 25.59, 8.997, 3.319, 0.3643])
c_coeffs_2S = np.array([-1.13e-4, -8.78e-4, -4.54e-3, -1.8133e-2, -5.576e-2, -0.126895, -0.170352, 0.598684])

c_alpha_3S = np.array([0.9059])
c_coeffs_3S = np.array([1.0])

c_alpha_4S = np.array([0.1285])
c_coeffs_4S = np.array([1.0])

c_alpha_5S = np.array([0.04402])
c_coeffs_5S = np.array([1.0])

# --- C P shells ---
c_alpha_1P = np.array([18.71, 4.133, 1.2])
c_coeffs_1P = np.array([0.014031, 0.086866, 0.290216])

c_alpha_2P = np.array([0.3827])
c_coeffs_2P = np.array([1.0])

c_alpha_3P = np.array([0.1209])
c_coeffs_3P = np.array([1.0])

c_alpha_4P = np.array([0.03569])
c_coeffs_4P = np.array([1.0])

# --- C D shells ---
c_alpha_1D = np.array([1.097])
c_coeffs_1D = np.array([1.0])

c_alpha_2D = np.array([0.318])
c_coeffs_2D = np.array([1.0])

c_alpha_3D = np.array([0.1])
c_coeffs_3D = np.array([1.0])

# --- C F shells ---
c_alpha_1F = np.array([0.761])
c_coeffs_1F = np.array([1.0])

c_alpha_2F = np.array([0.268])
c_coeffs_2F = np.array([1.0])

# --- N S shells ---
n_alpha_1S = np.array([11420.0, 1712.0, 389.3, 110.0, 35.57, 12.54, 4.644, 0.5118])
n_coeffs_1S = np.array([5.23e-4, 4.045e-3, 2.0775e-2, 8.0727e-2, 0.233074, 0.433501, 0.347472, -8.508e-3])

n_alpha_2S = np.array([11420.0, 1712.0, 389.3, 110.0, 35.57, 12.54, 4.644, 0.5118])
n_coeffs_2S = np.array([-1.15e-4, -8.95e-4, -4.624e-3, -1.8528e-2, -5.7339e-2, -0.132076, -0.17251, 0.599944])

n_alpha_3S = np.array([1.293])
n_coeffs_3S = np.array([1.0])

n_alpha_4S = np.array([0.1787])
n_coeffs_4S = np.array([1.0])

n_alpha_5S = np.array([0.0576])
n_coeffs_5S = np.array([1.0])

# --- N P shells ---
n_alpha_1P = np.array([26.63, 5.948, 1.742])
n_coeffs_1P = np.array([0.01467, 0.091764, 0.298683])

n_alpha_2P = np.array([0.555])
n_coeffs_2P = np.array([1.0])

n_alpha_3P = np.array([0.1725])
n_coeffs_3P = np.array([1.0])

n_alpha_4P = np.array([0.0491])
n_coeffs_4P = np.array([1.0])

# --- N D shells ---
n_alpha_1D = np.array([1.654])
n_coeffs_1D = np.array([1.0])

n_alpha_2D = np.array([0.469])
n_coeffs_2D = np.array([1.0])

n_alpha_3D = np.array([0.151])
n_coeffs_3D = np.array([1.0])

# --- N F shells ---
n_alpha_1F = np.array([1.093])
n_coeffs_1F = np.array([1.0])

n_alpha_2F = np.array([0.364])
n_coeffs_2F = np.array([1.0])

DO_coeffs_L = np.array([
     3.16457395e-04,
    -1.89930358e-01,
    -3.50233243e-02,
    -3.58713314e-01,
    -2.49990689e-01,
    -3.03663240e-15,
     6.41779414e-15,
    -2.42500851e-01,
    -4.75820011e-15,
     4.41906834e-14,
    -2.79193372e-01,
     8.81137522e-16,
     8.64041334e-14,
    -1.47366462e-01,
    -1.05372968e-15,
    -2.14744508e-14,
    -2.78468517e-02,
    -6.98568243e-17,
     1.24730248e-15,
     1.85075949e-02,
    -4.28511758e-17,
     1.54730177e-12,
    -2.34157912e-16,
    -5.63944497e-14,
     5.65877551e-03,
     9.59130747e-16,
    -3.99306125e-12,
     2.56604399e-16,
     2.93046446e-14,
    -1.34127809e-02,
    -1.85117139e-15,
     4.32062498e-11,
     4.88826008e-18,
    -3.01323087e-17,
     6.25882949e-15,
    -4.30515578e-03,
     1.49122565e-16,
     2.83627789e-14,
    -1.26459852e-17,
     9.05400608e-17,
     1.26815467e-16,
     1.12414448e-14,
     7.69190266e-03,
    -1.02087070e-16,
     1.30995702e-12,
     1.24571426e-16,
     3.78215640e-04,
    -9.52229509e-04,
     5.37506886e-03,
    -6.36037424e-02,
    -1.11538630e-02,
    -4.00728310e-15,
     6.20403387e-15,
     1.90263932e-01,
    -4.07635585e-15,
     6.43897039e-15,
     2.19452805e-01,
    -8.49060431e-15,
    -1.20313677e-13,
     1.51555287e-01,
    -1.04765063e-15,
     2.04918419e-15,
     5.50777242e-02,
    -7.70158490e-17,
     8.47518635e-16,
     8.88876540e-03,
    -2.19315244e-16,
     9.03595676e-14,
     6.05136703e-16,
    -2.51468644e-14,
     7.68754620e-03,
    -6.47600463e-16,
     1.35543522e-12,
    -1.32311267e-15,
    -4.79881527e-14,
    -4.87346023e-03,
    -2.55342304e-15,
    -4.03394201e-12,
     2.18537066e-16,
    -5.61713048e-17,
    -1.56203265e-15,
     1.76705276e-03,
    -1.67033786e-16,
     2.76310503e-13,
     4.38463794e-17,
    -1.56201810e-16,
     2.83876725e-16,
    -1.45705647e-14,
     1.44059310e-03,
    -4.71456323e-16,
    -1.72661546e-12,
    -9.11873038e-17
])

norm_L = 0.9492

DO_coeffs_R = np.array([
    -7.44795928e-06,
    -1.89986581e-01,
    -3.52877156e-02,
    -3.58001633e-01,
    -2.40966948e-01,
    -2.86794246e-15,
    5.42610074e-15,
    -2.40969120e-01,
    -5.71568106e-15,
    2.10538992e-14,
    -2.80151343e-01,
    -1.52887894e-15,
    1.24146531e-14,
    -1.49583776e-01,
    -1.33355472e-15,
    -2.66433710e-14,
    -2.81279911e-02,
    -7.00415150e-17,
    9.48614459e-16,
    1.90613388e-02,
    6.74965305e-17,
    1.51450585e-12,
    -3.71296890e-16,
    -3.08967367e-14,
    7.07509010e-03,
    1.46565686e-15,
    -3.76905492e-12,
    2.88010857e-16,
    5.72229519e-14,
    -1.27558660e-02,
    -6.89664831e-16,
    4.25353085e-11,
    -1.26981412e-17,
    -3.96918101e-18,
    4.77450703e-15,
    -4.80944406e-03,
    -3.74531902e-18,
    3.77361194e-14,
    -2.03180916e-17,
    7.85595898e-17,
    2.41679452e-16,
    -2.57627668e-15,
    7.42399997e-03,
    -1.43777137e-16,
    1.30359998e-12,
    1.24398636e-16,
    5.04272426e-05,
    -3.18006425e-03,
    4.00801753e-03,
    -6.70167109e-02,
    -1.38558028e-02,
    -4.16854704e-15,
    2.30678269e-15,
    1.90949090e-01,
    -4.54672669e-15,
    6.27792961e-17,
    2.22751816e-01,
    -5.22868482e-15,
    -2.02409163e-14,
    1.52259773e-01,
    -3.87792537e-16,
    1.66232483e-14,
    5.18366878e-02,
    -7.00651963e-17,
    2.51514408e-17,
    9.15228869e-03,
    -2.44194802e-16,
    9.16220356e-14,
    6.97973958e-16,
    -1.15258857e-14,
    8.02230448e-03,
    -4.79019912e-16,
    1.38830659e-12,
    -1.24848608e-15,
    -1.49016592e-14,
    -5.38455344e-03,
    -1.05611954e-15,
    -4.11764600e-12,
    2.29500182e-16,
    -6.39377552e-17,
    -9.92413523e-16,
    1.87635685e-03,
    -9.78234408e-17,
    2.91068487e-13,
    4.34705307e-17,
    -1.62594082e-16,
    4.07010786e-16,
    -1.11453300e-14,
    1.26176589e-03,
    -9.00676731e-17,
    -1.74207006e-12,
    -8.61067160e-17
])

norm_R = 0.9365
def build_DO(DO_coeffs, x, y, z, dV, recenter=True, print_geom=False):
    global R_C, R_N
    basis_info = [
    # ---- Carbon S orbitals ----
        (c_alpha_1S, c_coeffs_1S, 0, 0, 0, R_C, 0),
        (c_alpha_2S, c_coeffs_2S, 0, 0, 0, R_C, 1),
        (c_alpha_3S, c_coeffs_3S, 0, 0, 0, R_C, 2),
        (c_alpha_4S, c_coeffs_4S, 0, 0, 0, R_C, 3),
        (c_alpha_5S, c_coeffs_5S, 0, 0, 0, R_C, 4),

    # ---- Carbon P orbitals ----
        (c_alpha_1P, c_coeffs_1P, 1, 0, 0, R_C, 5),
        (c_alpha_1P, c_coeffs_1P, 0, 1, 0, R_C, 6),
        (c_alpha_1P, c_coeffs_1P, 0, 0, 1, R_C, 7),

        (c_alpha_2P, c_coeffs_2P, 1, 0, 0, R_C, 8),
        (c_alpha_2P, c_coeffs_2P, 0, 1, 0, R_C, 9),
        (c_alpha_2P, c_coeffs_2P, 0, 0, 1, R_C, 10),

        (c_alpha_3P, c_coeffs_3P, 1, 0, 0, R_C, 11),
        (c_alpha_3P, c_coeffs_3P, 0, 1, 0, R_C, 12),
        (c_alpha_3P, c_coeffs_3P, 0, 0, 1, R_C, 13),

        (c_alpha_4P, c_coeffs_4P, 1, 0, 0, R_C, 14),
        (c_alpha_4P, c_coeffs_4P, 0, 1, 0, R_C, 15),
        (c_alpha_4P, c_coeffs_4P, 0, 0, 1, R_C, 16),

    # ---- Carbon D orbitals ----
        (c_alpha_1D, c_coeffs_1D, 1, 1, 0, R_C, 17),
        (c_alpha_1D, c_coeffs_1D, 0, 1, 1, R_C, 18),
    # d_z^2 (19) done manually
        (c_alpha_1D, c_coeffs_1D, 1, 0, 1, R_C, 20),
    # d_x2-y2 (21) done manually

        (c_alpha_2D, c_coeffs_2D, 1, 1, 0, R_C, 22),
        (c_alpha_2D, c_coeffs_2D, 0, 1, 1, R_C, 23),
    # d_z^2 (24) done manually
        (c_alpha_2D, c_coeffs_2D, 1, 0, 1, R_C, 25),
    # d_x2-y2 (26) done manually

        (c_alpha_3D, c_coeffs_3D, 1, 1, 0, R_C, 27),
        (c_alpha_3D, c_coeffs_3D, 0, 1, 1, R_C, 28),
    # d_z^2 (29) done manually
        (c_alpha_3D, c_coeffs_3D, 1, 0, 1, R_C, 30),
    # d_x2-y2 (31) done manually

    # ---- Carbon F orbitals ----
    # f_y3x2 (32) done manually
        (c_alpha_1F, c_coeffs_1F, 1, 1, 1, R_C, 33),
    # f_yz2 (34) done manually
    # f_z3 (35) done manually
    # f_xz2 (36) done manually
    # f_x2y2z (37) done manually
    # f_x3y2 (38) done manually

    # f_y3x2 (39) done manually
        (c_alpha_2F, c_coeffs_2F, 1, 1, 1, R_C, 40),
    # f_yz2 (41) done manually
    # f_z3 (42) done manually
    # f_xz2 (43) done manually
    # f_x2y2z (44) done manually
    # f_x3y2 (45) done manually

    # ---- Nitrogen S orbitals ----
        (n_alpha_1S, n_coeffs_1S, 0, 0, 0, R_N, 46),
        (n_alpha_2S, n_coeffs_2S, 0, 0, 0, R_N, 47),
        (n_alpha_3S, n_coeffs_3S, 0, 0, 0, R_N, 48),
        (n_alpha_4S, n_coeffs_4S, 0, 0, 0, R_N, 49),
        (n_alpha_5S, n_coeffs_5S, 0, 0, 0, R_N, 50),

    # ---- Nitrogen P orbitals ----
        (n_alpha_1P, n_coeffs_1P, 1, 0, 0, R_N, 51),
        (n_alpha_1P, n_coeffs_1P, 0, 1, 0, R_N, 52),
        (n_alpha_1P, n_coeffs_1P, 0, 0, 1, R_N, 53),

        (n_alpha_2P, n_coeffs_2P, 1, 0, 0, R_N, 54),
        (n_alpha_2P, n_coeffs_2P, 0, 1, 0, R_N, 55),
        (n_alpha_2P, n_coeffs_2P, 0, 0, 1, R_N, 56),

        (n_alpha_3P, n_coeffs_3P, 1, 0, 0, R_N, 57),
        (n_alpha_3P, n_coeffs_3P, 0, 1, 0, R_N, 58),
        (n_alpha_3P, n_coeffs_3P, 0, 0, 1, R_N, 59),

        (n_alpha_4P, n_coeffs_4P, 1, 0, 0, R_N, 60),
        (n_alpha_4P, n_coeffs_4P, 0, 1, 0, R_N, 61),
        (n_alpha_4P, n_coeffs_4P, 0, 0, 1, R_N, 62),

    # ---- Nitrogen D orbitals ----
        (n_alpha_1D, n_coeffs_1D, 1, 1, 0, R_N, 63),
        (n_alpha_1D, n_coeffs_1D, 0, 1, 1, R_N, 64),
    # d_z^2 (65) done manually
        (n_alpha_1D, n_coeffs_1D, 1, 0, 1, R_N, 66),
    # d_x2-y2 (67) done manually

        (n_alpha_2D, n_coeffs_2D, 1, 1, 0, R_N, 68),
        (n_alpha_2D, n_coeffs_2D, 0, 1, 1, R_N, 69),
    # d_z^2 (70) done manually
        (n_alpha_2D, n_coeffs_2D, 1, 0, 1, R_N, 71),
    # d_x2-y2 (72) done manually

        (n_alpha_3D, n_coeffs_3D, 1, 1, 0, R_N, 73),
        (n_alpha_3D, n_coeffs_3D, 0, 1, 1, R_N, 74),
    # d_z^2 (75) done manually
        (n_alpha_3D, n_coeffs_3D, 1, 0, 1, R_N, 76),
    # d_x2-y2 (77) done manually

    # ---- Nitrogen F orbitals ----
    # f_y3x2 (78) done manually
        (n_alpha_1F, n_coeffs_1F, 1, 1, 1, R_N, 79),
    # f_yz2 (80) done manually
    # f_z3 (81) done manually
    # f_xz2 (82) done manually
    # f_x2y2z (83) done manually
    # f_x3y2 (84) done manually

    # f_y3x2 (85) done manually
        (n_alpha_2F, n_coeffs_2F, 1, 1, 1, R_N, 86),
    # f_yz2 (87) done manually
    # f_z3 (88) done manually
    # f_xz2 (89) done manually
    # f_x2y2z (90) done manually
    # f_x3y2 (91) done manually
]
 
    DO = np.zeros_like(x, dtype=np.float64)
    threshold = 1e-3
    for alphas, coeffs, a, b, c, center, i in basis_info:
        coeff = DO_coeffs[i]
        if abs(coeff) >= threshold:
            ao = AO(alphas, coeffs, a, b, c, center, x, y, z, dV)
            DO += coeff * ao

    def norm(orb):
        return 1 / np.sqrt(integrate_3d(abs(orb)**2, dV))
    
    # Handle weird orbitals manually
    if abs(DO_coeffs[19]) >= threshold:
        ao_zz = AO(c_alpha_1D, c_coeffs_1D, 0, 0, 2, R_C, x, y, z, dV)
        ao_xx = AO(c_alpha_1D, c_coeffs_1D, 2, 0, 0, R_C, x, y, z, dV)
        ao_yy = AO(c_alpha_1D, c_coeffs_1D, 0, 2, 0, R_C, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[19] * normed_AO

    if abs(DO_coeffs[21]) >= threshold:
        ao_xx = AO(c_alpha_1D, c_coeffs_1D, 2, 0, 0, R_C, x, y, z, dV)
        ao_yy = AO(c_alpha_1D, c_coeffs_1D, 0, 2, 0, R_C, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[21] * normed_AO

    if abs(DO_coeffs[24]) >= threshold:
        ao_zz = AO(c_alpha_2D, c_coeffs_2D, 0, 0, 2, R_C, x, y, z, dV)
        ao_xx = AO(c_alpha_2D, c_coeffs_2D, 2, 0, 0, R_C, x, y, z, dV)
        ao_yy = AO(c_alpha_2D, c_coeffs_2D, 0, 2, 0, R_C, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[24] * normed_AO

    if abs(DO_coeffs[26]) >= threshold:
        ao_xx = AO(c_alpha_2D, c_coeffs_2D, 2, 0, 0, R_C, x, y, z, dV)
        ao_yy = AO(c_alpha_2D, c_coeffs_2D, 0, 2, 0, R_C, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[26] * normed_AO

    if abs(DO_coeffs[29]) >= threshold:
        ao_zz = AO(c_alpha_3D, c_coeffs_3D, 0, 0, 2, R_C, x, y, z, dV)
        ao_xx = AO(c_alpha_3D, c_coeffs_3D, 2, 0, 0, R_C, x, y, z, dV)
        ao_yy = AO(c_alpha_3D, c_coeffs_3D, 0, 2, 0, R_C, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[29] * normed_AO

    if abs(DO_coeffs[31]) >= threshold:
        ao_xx = AO(c_alpha_3D, c_coeffs_3D, 2, 0, 0, R_C, x, y, z, dV)
        ao_yy = AO(c_alpha_3D, c_coeffs_3D, 0, 2, 0, R_C, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[31] * normed_AO

    if abs(DO_coeffs[32]) >= threshold:
        ao_yx2 = AO(c_alpha_1F, c_coeffs_1F, 2, 1, 0, R_C, x, y, z, dV)
        ao_y3  = AO(c_alpha_1F, c_coeffs_1F, 0, 3, 0, R_C, x, y, z, dV)
        orb = np.sqrt(5/8) * (3 * ao_yx2 - ao_y3)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[32] * normed_AO

    if abs(DO_coeffs[34]) >= threshold:
        ao_yz2 = AO(c_alpha_1F, c_coeffs_1F, 0, 1, 2, R_C, x, y, z, dV)
        ao_y3  = AO(c_alpha_1F, c_coeffs_1F, 0, 3, 0, R_C, x, y, z, dV)
        ao_yx2 = AO(c_alpha_1F, c_coeffs_1F, 2, 1, 0, R_C, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_yz2 - ao_y3 - ao_yx2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[34] * normed_AO

    if abs(DO_coeffs[35]) >= threshold:
        ao_z3   = AO(c_alpha_1F, c_coeffs_1F, 0, 0, 3, R_C, x, y, z, dV)
        ao_x2z  = AO(c_alpha_1F, c_coeffs_1F, 2, 0, 1, R_C, x, y, z, dV)
        ao_y2z  = AO(c_alpha_1F, c_coeffs_1F, 0, 2, 1, R_C, x, y, z, dV)
        orb = 0.5 * (2 * ao_z3 - 3 * ao_x2z - 3 * ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[35] * normed_AO

    if abs(DO_coeffs[36]) >= threshold:
        ao_xz2 = AO(c_alpha_1F, c_coeffs_1F, 1, 0, 2, R_C, x, y, z, dV)
        ao_x3   = AO(c_alpha_1F, c_coeffs_1F, 3, 0, 0, R_C, x, y, z, dV)
        ao_xy2  = AO(c_alpha_1F, c_coeffs_1F, 1, 2, 0, R_C, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_xz2 - ao_x3 - ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[36] * normed_AO

    if abs(DO_coeffs[37]) >= threshold:
        ao_x2z = AO(c_alpha_1F, c_coeffs_1F, 2, 0, 1, R_C, x, y, z, dV)
        ao_y2z = AO(c_alpha_1F, c_coeffs_1F, 0, 2, 1, R_C, x, y, z, dV)
        orb = (np.sqrt(15) / 2) * (ao_x2z - ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[37] * normed_AO

    if abs(DO_coeffs[38]) >= threshold:
        ao_x3 = AO(c_alpha_1F, c_coeffs_1F, 3, 0, 0, R_C, x, y, z, dV)
        ao_xy2 = AO(c_alpha_1F, c_coeffs_1F, 1, 2, 0, R_C, x, y, z, dV)
        orb = np.sqrt(5/8) * (ao_x3 - 3 * ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[38] * normed_AO

    if abs(DO_coeffs[39]) >= threshold:
        ao_yx2 = AO(c_alpha_2F, c_coeffs_2F, 2, 1, 0, R_C, x, y, z, dV)
        ao_y3  = AO(c_alpha_2F, c_coeffs_2F, 0, 3, 0, R_C, x, y, z, dV)
        orb = np.sqrt(5/8) * (3 * ao_yx2 - ao_y3)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[39] * normed_AO

    if abs(DO_coeffs[41]) >= threshold:
        ao_yz2 = AO(c_alpha_2F, c_coeffs_2F, 0, 1, 2, R_C, x, y, z, dV)
        ao_y3  = AO(c_alpha_2F, c_coeffs_2F, 0, 3, 0, R_C, x, y, z, dV)
        ao_yx2 = AO(c_alpha_2F, c_coeffs_2F, 2, 1, 0, R_C, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_yz2 - ao_y3 - ao_yx2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[41] * normed_AO

    if abs(DO_coeffs[42]) >= threshold:
        ao_z3   = AO(c_alpha_2F, c_coeffs_2F, 0, 0, 3, R_C, x, y, z, dV)
        ao_x2z  = AO(c_alpha_2F, c_coeffs_2F, 2, 0, 1, R_C, x, y, z, dV)
        ao_y2z  = AO(c_alpha_2F, c_coeffs_2F, 0, 2, 1, R_C, x, y, z, dV)
        orb = 0.5 * (2 * ao_z3 - 3 * ao_x2z - 3 * ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[42] * normed_AO

    if abs(DO_coeffs[43]) >= threshold:
        ao_xz2 = AO(c_alpha_2F, c_coeffs_2F, 1, 0, 2, R_C, x, y, z, dV)
        ao_x3   = AO(c_alpha_2F, c_coeffs_2F, 3, 0, 0, R_C, x, y, z, dV)
        ao_xy2  = AO(c_alpha_2F, c_coeffs_2F, 1, 2, 0, R_C, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_xz2 - ao_x3 - ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[43] * normed_AO

    if abs(DO_coeffs[44]) >= threshold:
        ao_x2z = AO(c_alpha_2F, c_coeffs_2F, 2, 0, 1, R_C, x, y, z, dV)
        ao_y2z = AO(c_alpha_2F, c_coeffs_2F, 0, 2, 1, R_C, x, y, z, dV)
        orb = (np.sqrt(15) / 2) * (ao_x2z - ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[44] * normed_AO

    if abs(DO_coeffs[45]) >= threshold:
        ao_x3 = AO(c_alpha_2F, c_coeffs_2F, 3, 0, 0, R_C, x, y, z, dV)
        ao_xy2 = AO(c_alpha_2F, c_coeffs_2F, 1, 2, 0, R_C, x, y, z, dV)
        orb = np.sqrt(5/8) * (ao_x3 - 3 * ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[45] * normed_AO

    if abs(DO_coeffs[65]) >= threshold:
        ao_zz = AO(n_alpha_1D, n_coeffs_1D, 0, 0, 2, R_N, x, y, z, dV)
        ao_xx = AO(n_alpha_1D, n_coeffs_1D, 2, 0, 0, R_N, x, y, z, dV)
        ao_yy = AO(n_alpha_1D, n_coeffs_1D, 0, 2, 0, R_N, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[65] * normed_AO

    if abs(DO_coeffs[67]) >= threshold:
        ao_xx = AO(n_alpha_1D, n_coeffs_1D, 2, 0, 0, R_N, x, y, z, dV)
        ao_yy = AO(n_alpha_1D, n_coeffs_1D, 0, 2, 0, R_N, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[67] * normed_AO

    if abs(DO_coeffs[70]) >= threshold:
        ao_zz = AO(n_alpha_2D, n_coeffs_2D, 0, 0, 2, R_N, x, y, z, dV)
        ao_xx = AO(n_alpha_2D, n_coeffs_2D, 2, 0, 0, R_N, x, y, z, dV)
        ao_yy = AO(n_alpha_2D, n_coeffs_2D, 0, 2, 0, R_N, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[70] * normed_AO

    if abs(DO_coeffs[72]) >= threshold:
        ao_xx = AO(n_alpha_2D, n_coeffs_2D, 2, 0, 0, R_N, x, y, z, dV)
        ao_yy = AO(n_alpha_2D, n_coeffs_2D, 0, 2, 0, R_N, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[72] * normed_AO

    if abs(DO_coeffs[75]) >= threshold:
        ao_zz = AO(n_alpha_3D, n_coeffs_3D, 0, 0, 2, R_N, x, y, z, dV)
        ao_xx = AO(n_alpha_3D, n_coeffs_3D, 2, 0, 0, R_N, x, y, z, dV)
        ao_yy = AO(n_alpha_3D, n_coeffs_3D, 0, 2, 0, R_N, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[75] * normed_AO

    if abs(DO_coeffs[77]) >= threshold:
        ao_xx = AO(n_alpha_3D, n_coeffs_3D, 2, 0, 0, R_N, x, y, z, dV)
        ao_yy = AO(n_alpha_3D, n_coeffs_3D, 0, 2, 0, R_N, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[77] * normed_AO

    if abs(DO_coeffs[78]) >= threshold:
        ao_yx2 = AO(n_alpha_1F, n_coeffs_1F, 2, 1, 0, R_N, x, y, z, dV)
        ao_y3  = AO(n_alpha_1F, n_coeffs_1F, 0, 3, 0, R_N, x, y, z, dV)
        orb = np.sqrt(5/8) * (3 * ao_yx2 - ao_y3)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[78] * normed_AO

    if abs(DO_coeffs[80]) >= threshold:
        ao_yz2 = AO(n_alpha_1F, n_coeffs_1F, 0, 1, 2, R_N, x, y, z, dV)
        ao_y3  = AO(n_alpha_1F, n_coeffs_1F, 0, 3, 0, R_N, x, y, z, dV)
        ao_yx2 = AO(n_alpha_1F, n_coeffs_1F, 2, 1, 0, R_N, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_yz2 - ao_y3 - ao_yx2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[80] * normed_AO

    if abs(DO_coeffs[81]) >= threshold:
        ao_z3   = AO(n_alpha_1F, n_coeffs_1F, 0, 0, 3, R_N, x, y, z, dV)
        ao_x2z  = AO(n_alpha_1F, n_coeffs_1F, 2, 0, 1, R_N, x, y, z, dV)
        ao_y2z  = AO(n_alpha_1F, n_coeffs_1F, 0, 2, 1, R_N, x, y, z, dV)
        orb = 0.5 * (2 * ao_z3 - 3 * ao_x2z - 3 * ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[81] * normed_AO

    if abs(DO_coeffs[82]) >= threshold:
        ao_xz2 = AO(n_alpha_1F, n_coeffs_1F, 1, 0, 2, R_N, x, y, z, dV)
        ao_x3   = AO(n_alpha_1F, n_coeffs_1F, 3, 0, 0, R_N, x, y, z, dV)
        ao_xy2  = AO(n_alpha_1F, n_coeffs_1F, 1, 2, 0, R_N, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_xz2 - ao_x3 - ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[82] * normed_AO

    if abs(DO_coeffs[83]) >= threshold:
        ao_x2z = AO(n_alpha_1F, n_coeffs_1F, 2, 0, 1, R_N, x, y, z, dV)
        ao_y2z = AO(n_alpha_1F, n_coeffs_1F, 0, 2, 1, R_N, x, y, z, dV)
        orb = (np.sqrt(15) / 2) * (ao_x2z - ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[83] * normed_AO

    if abs(DO_coeffs[84]) >= threshold:
        ao_x3 = AO(n_alpha_1F, n_coeffs_1F, 3, 0, 0, R_N, x, y, z, dV)
        ao_xy2 = AO(n_alpha_1F, n_coeffs_1F, 1, 2, 0, R_N, x, y, z, dV)
        orb = np.sqrt(5/8) * (ao_x3 - 3 * ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[84] * normed_AO

    if abs(DO_coeffs[85]) >= threshold:
        ao_yx2 = AO(n_alpha_2F, n_coeffs_2F, 2, 1, 0, R_N, x, y, z, dV)
        ao_y3  = AO(n_alpha_2F, n_coeffs_2F, 0, 3, 0, R_N, x, y, z, dV)
        orb = np.sqrt(5/8) * (3 * ao_yx2 - ao_y3)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[85] * normed_AO

    if abs(DO_coeffs[87]) >= threshold:
        ao_yz2 = AO(n_alpha_2F, n_coeffs_2F, 0, 1, 2, R_N, x, y, z, dV)
        ao_y3  = AO(n_alpha_2F, n_coeffs_2F, 0, 3, 0, R_N, x, y, z, dV)
        ao_yx2 = AO(n_alpha_2F, n_coeffs_2F, 2, 1, 0, R_N, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_yz2 - ao_y3 - ao_yx2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[87] * normed_AO

    if abs(DO_coeffs[88]) >= threshold:
        ao_z3   = AO(n_alpha_2F, n_coeffs_2F, 0, 0, 3, R_N, x, y, z, dV)
        ao_x2z  = AO(n_alpha_2F, n_coeffs_2F, 2, 0, 1, R_N, x, y, z, dV)
        ao_y2z  = AO(n_alpha_2F, n_coeffs_2F, 0, 2, 1, R_N, x, y, z, dV)
        orb = 0.5 * (2 * ao_z3 - 3 * ao_x2z - 3 * ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[88] * normed_AO

    if abs(DO_coeffs[89]) >= threshold:
        ao_xz2 = AO(n_alpha_2F, n_coeffs_2F, 1, 0, 2, R_N, x, y, z, dV)
        ao_x3   = AO(n_alpha_2F, n_coeffs_2F, 3, 0, 0, R_N, x, y, z, dV)
        ao_xy2  = AO(n_alpha_2F, n_coeffs_2F, 1, 2, 0, R_N, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_xz2 - ao_x3 - ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[89] * normed_AO

    if abs(DO_coeffs[90]) >= threshold:
        ao_x2z = AO(n_alpha_2F, n_coeffs_2F, 2, 0, 1, R_N, x, y, z, dV)
        ao_y2z = AO(n_alpha_2F, n_coeffs_2F, 0, 2, 1, R_N, x, y, z, dV)
        orb = (np.sqrt(15) / 2) * (ao_x2z - ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[90] * normed_AO

    if abs(DO_coeffs[91]) >= threshold:
        ao_x3 = AO(n_alpha_2F, n_coeffs_2F, 3, 0, 0, R_N, x, y, z, dV)
        ao_xy2 = AO(n_alpha_2F, n_coeffs_2F, 1, 2, 0, R_N, x, y, z, dV)
        orb = np.sqrt(5/8) * (ao_x3 - 3 * ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[91] * normed_AO

    Norm = norm(DO)
    DO /= Norm

    if recenter:
        centroid = recenter_DO(DO, x, y, z, dV)
        
        # Shift grid
        x_new = x - centroid[0]
        y_new = y - centroid[1]
        z_new = z - centroid[2]

        # Shift molecular geometry
        R_C = R_C - centroid
        R_N = R_N - centroid

        if print_geom:
            bohr_to_angstrom = 0.529177
            R_C_ang = R_C * bohr_to_angstrom
            R_N_ang = R_N * bohr_to_angstrom

            print("\nShifting molecular geometry to the new center")
            print("New molecular geometry is (in Ångströms):")
            print(" atom         X             Y             Z")
            print("   1      {: .6f}      {: .6f}      {: .6f}".format(*R_C_ang))
            print("   2      {: .6f}      {: .6f}      {: .6f}".format(*R_N_ang))

        # Rebuild on shifted grid
        return build_DO(DO_coeffs, x_new, y_new, z_new, dV, recenter=False)

    return DO

# ---------------------------------------------
# Plane-wave continuum
# ---------------------------------------------
def plane_wave(k_vec, x, y, z):
    k_dot_r = k_vec[0]*x + k_vec[1]*y + k_vec[2]*z
    return np.exp(1j * k_dot_r)

# ---------------------------------------------
# Build rotation matrix from Euler angles (ZYZ convention)
# ---------------------------------------------
def rotation_matrix(alpha, beta, gamma):
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta),  np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)
    Rz1 = np.array([[ ca, -sa, 0], [ sa,  ca, 0], [  0,   0, 1]])
    Ry  = np.array([[ cb,   0, sb], [  0,   1,  0], [-sb,  0, cb]])
    Rz2 = np.array([[ cg, -sg, 0], [ sg,  cg, 0], [  0,   0, 1]])
    return Rz1 @ Ry @ Rz2

# ---------------------------------------------
# Orientational averaging routines
# ---------------------------------------------
def compute_average_amplitudes(X, Y, Z, DO_coeffs_L, DO_coeffs_R, ang_grid, k_lab, pol_lab, mode='parallel'):
    DO_L = build_DO(DO_coeffs_L, X, Y, Z, dV)
    DO_R = build_DO(DO_coeffs_R, X, Y, Z, dV)
    N = ang_grid.shape[1]
    amps = np.zeros(N, dtype=complex)
    for i in range(N):
        a, b, g = ang_grid[:, i]
        R = rotation_matrix(a, b, g)
        # rotate polarization and k into molecular frame
        eps_mol = R.T @ pol_lab
        k_mol   = R.T @ k_lab
        # dipole operator mu·eps = eps_x x + eps_y y + eps_z z
        mu = eps_mol[0]*X + eps_mol[1]*Y + eps_mol[2]*Z
        # plane wave in molecular frame
        psi_f = plane_wave(k_mol, X, Y, Z)
        # transition amplitude A = ∫ psi_f * mu * DO d^3r
        integrand_L = np.conj(psi_f) * mu * DO_L
        integrand_R = psi_f * mu * np.conj(DO_R)
        A_L = integrate_3d(integrand_L, dV)
        A_R = integrate_3d(integrand_R, dV)
        amps[i] = A_L * A_R
    return np.real(amps)

# ---------------------------------------------
# Main script
# ---------------------------------------------
if __name__ == '__main__':

    # lab-frame vectors
    pol_lab = np.array([0.0, 0.0, 1.0])  # z-polarization
    k_par_lab = np.array([0.0, 0.0, 1.0])
    k_perp1_lab = np.array([1.0, 0.0, 0.0])
    k_perp2_lab = np.array([0.0, 1.0, 0.0])

    # Euler angle grid (3, N)
    ang_grid = np.array([
        [1.07393, 3.31052, 1.74962, 1.64604, 6.10750, 1.75961, 5.53027, 5.52842, 2.28479, 3.58752, 1.66523, 0.84226, 3.29536, 5.57673, 2.47174, 0.92423, 3.93291, 3.31288, 1.90607, 5.84085, 0.38988, 3.04104, 3.73432, 0.06801, 4.53341, 3.95200, 4.40989, 3.00258, 0.72053, 3.05150,  3.43130, 4.27813, 0.33053, 0.21824, 1.26162, 2.18568, 4.36603, 0.72968, 0.09293, 1.39264, 0.61877, 3.83971, 1.33503, 2.72987, 5.64007, 0.64928, 4.61242, 3.37637, 5.84322, 3.88503, 3.30267, 2.80485, 5.96564, 3.63200, 3.21888, 5.15715, 6.16719, 5.35457, 0.54376, 4.22616, 4.85400, 3.31674, 1.37092, 4.28857, 1.13002, 0.71738, 0.54974, 0.55092, 1.24132, 3.64290, 1.64918, 0.21685, 3.08618, 0.92048, 2.12227, 2.48303, 5.11643, 4.98922, 5.85691, 5.20647, 0.55287, 4.78971, 0.20077, 1.62743, 1.16900, 5.35471, 3.02272, 6.06776, 2.18861, 4.82467, 5.11437, 2.76769, 1.87105, 4.61828, 4.49672, 3.57552, 2.38265, 2.76245, 4.72005, 4.04847, 4.67308, 0.03871, 5.47292, 0.39803, 5.04133, 1.91134, 2.47181, 2.18850, 0.02788, 1.64398, 2.99010, 5.39785, 2.22545, 0.13728, 4.95468, 3.95088, 4.48075, 3.89107, 5.28079, 1.04117, 1.52149, 6.18079, 5.73604, 1.91169, 5.92459, 5.51662, 3.78850, 2.46795, 4.85331, 2.43693, 0.84941, 4.11191, 1.05014, 3.60250, 2.66381, 5.14617, 1.60611, 5.85919, 6.20983, 3.03865, 5.79858, 4.28145, 2.75753, 4.16086, 5.81163, 4.79705, 1.95620, 2.56970, 1.91953, 1.30330],
        [2.37487, 2.01802, 0.19298, 1.42732, 0.94148, 2.87768, 0.43346, 2.23837, 2.38175, 1.86742, 0.48430, 1.47576, 2.33447, 1.61340, 1.85509, 2.01268, 0.75689, 1.06096, 1.28682, 2.45608, 1.22613, 0.85503, 2.48750, 2.59037, 1.12758, 1.08387, 1.71823, 2.89856, 1.19880, 1.82581, 0.72273, 2.03217, 2.31736, 1.47451, 2.64602, 1.75559, 0.58855, 1.75632, 1.74660, 1.57420, 0.64023, 1.67455, 1.85446, 1.02968, 1.28404, 2.56314, 1.94129, 0.21312, 0.13739, 1.37820, 1.37309, 2.42800, 1.21872, 2.18767, 2.62084, 2.30207, 2.00292, 1.40952, 1.48686, 1.22480, 2.13311, 1.69066, 1.25267, 2.52738, 1.45655, 2.25519, 2.00568, 2.85707, 2.12000, 0.96483, 1.75524, 2.03425, 0.51224, 0.89061, 0.67305, 1.24734, 0.87917, 0.58117, 1.50352, 2.01029, 0.93519, 1.36289, 0.97133, 0.77734, 0.64195, 2.56591, 1.19311, 1.72994, 1.16124, 2.76802, 1.23525, 1.94278, 2.56288, 0.31374, 2.25460, 1.54648, 0.91810, 1.35613, 1.64956, 2.27038, 0.80996, 0.41045, 0.74431, 1.74535, 1.53611, 1.60435, 1.55181, 1.45794, 1.21354, 1.09534, 2.15902, 1.07567, 2.06928, 0.69204, 1.83188, 2.77746, 1.42929, 1.99852, 1.71301, 1.17178, 2.35523, 1.46674, 0.97673, 2.24095, 0.67027, 1.93475, 0.47169, 2.66075, 1.07081, 0.42167, 0.36992, 1.78975, 1.74413, 1.25456, 0.70579, 3.06199, 2.05766, 2.12962, 2.29422, 1.51664, 1.83453, 0.90821, 1.65087, 1.50476, 2.78176, 2.47069, 0.94743, 2.17897, 1.92705, 0.94949],
        [0.00649606, 0.00677567, 0.00640311, 0.00666473, 0.00669922, 0.00663394, 0.00674343, 0.00667127, 0.00668776, 0.00676721, 0.00665191, 0.00663109, 0.00663458, 0.00674962, 0.00667875, 0.00676506, 0.00663678, 0.00664742, 0.00663493, 0.00677988, 0.00666949, 0.00678447, 0.00678442, 0.00674655, 0.00660888, 0.00662420, 0.00663710, 0.00638865, 0.00674299, 0.00667288, 0.00663664, 0.00663669, 0.00669460, 0.00674184, 0.00666331, 0.00668894, 0.00668811, 0.00678161, 0.00666844, 0.00649848, 0.00672334, 0.00663329, 0.00657182, 0.00678480, 0.00678987, 0.00668689, 0.00668108, 0.00665181, 0.00659894, 0.00675097, 0.00676891, 0.00677254, 0.00678260, 0.00672571, 0.00664667, 0.00640289, 0.00656302, 0.00671165, 0.00670240, 0.00655946, 0.00665147, 0.00671663, 0.00668897, 0.00670018, 0.00646668, 0.00646788, 0.00663172, 0.00672750, 0.00657072, 0.00637692, 0.00672365, 0.00674388, 0.00667000, 0.00666926, 0.00675168, 0.00664912, 0.00663779, 0.00677741, 0.00670043, 0.00660064, 0.00656175, 0.00666901, 0.00638311, 0.00668824, 0.00678113, 0.00665192, 0.00672800, 0.00638353, 0.00639049, 0.00675165, 0.00663724, 0.00639167, 0.00676403, 0.00675114, 0.00666951, 0.00662358, 0.00666779, 0.00663807, 0.00668906, 0.00678518, 0.00666891, 0.00666294, 0.00671233, 0.00670227, 0.00677799, 0.00676458, 0.00677420, 0.00666742, 0.00659765, 0.00672893, 0.00665960, 0.00677759, 0.00676543, 0.00657994, 0.00675244, 0.00666695, 0.00672720, 0.00664683, 0.00674483, 0.00669447, 0.00672329, 0.00659747, 0.00678974, 0.00678624, 0.00674990, 0.00666400, 0.00668276, 0.00666590, 0.00637544, 0.00664386, 0.00667255, 0.00637753, 0.00676400, 0.00663239, 0.00669928, 0.00664242, 0.00676106, 0.00672395, 0.00666835, 0.00677843, 0.00658002, 0.00672617, 0.00666076, 0.00662597, 0.00668821, 0.00664327, 0.00664268, 0.00667836, 0.00678659, 0.00674897]
    ])

    N_ang = ang_grid.shape[1]

    # photoelectron energies (eV)
    E_eV = np.linspace(0.001, 0.501, 11)
    hartree = 27.2114
    beta_vals = []

    for E in E_eV:
        E_au = E / hartree
        k_mag = np.sqrt(2 * E_au)
        # scale k_lab directions
        k_par_lab_scaled = k_par_lab * k_mag
        k_perp1_lab_scaled = k_perp1_lab * k_mag
        k_perp2_lab_scaled = k_perp2_lab * k_mag

        # compute amplitudes
        A_par = compute_average_amplitudes(X, Y, Z, DO_coeffs_L, DO_coeffs_R, ang_grid, k_par_lab_scaled, pol_lab, mode='parallel')
        A_perp1 = compute_average_amplitudes(X, Y, Z, DO_coeffs_L, DO_coeffs_R, ang_grid, k_perp1_lab_scaled, pol_lab, mode='perpendicular')
        A_perp2 = compute_average_amplitudes(X, Y, Z, DO_coeffs_L, DO_coeffs_R, ang_grid, k_perp2_lab_scaled, pol_lab, mode='perpendicular')

        # intensities
        sigma_par = np.mean(A_par)
        sigma_perp = 0.5 * (np.mean(A_perp1) + np.mean(A_perp2))

        # anisotropy parameter
        beta = 2 * (sigma_par - sigma_perp) / (sigma_par + 2 * sigma_perp)
        beta_vals.append(beta)

    # plot
    print(beta_vals)
    plt.plot(E_eV, beta_vals, marker='o')
    plt.xlabel('Photoelectron KE (eV)')
    plt.ylabel(r'$\beta$')
    plt.ylim(-1, 2)
    plt.grid(True)
    plt.show()

