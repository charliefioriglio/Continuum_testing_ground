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
bl_A = 1.67  # Bond length in angstroms
bl_au = bl_A / 52.9e-2  # bond length in au
a = bl_au / 2
R_A = np.array([0, 0, -a])
R_B = np.array([0, 0, a])

def recenter_DO(DO, x, y, z, dV):
    density = np.abs(DO)**2
    x_c = np.sum(density * x) * dV
    y_c = np.sum(density * y) * dV
    z_c = np.sum(density * z) * dV
    return np.array([x_c, y_c, z_c])

# --- Cu S shells ---
cu_alpha_1S = np.array([560.088, 56.6486, 35.4258, 11.0546, 2.30682, 0.951429, 0.145184, 0.050325])
cu_coeffs_1S = np.array([6.37e-4, -9.735e-3, 6.5793e-2, -0.415035, 0.746611, 0.462173, 0.015983, -0.002767])

cu_alpha_2S = np.array([560.088, 56.6486, 35.4258, 11.0546, 2.30682, 0.951429, 0.145184, 0.050325])
cu_coeffs_2S = np.array([-1.36e-4, 0.001401, -0.013174, 0.095695, -0.211874, -0.235944, 0.508115, 0.621519])

cu_alpha_3S = np.array([0.050325])
cu_coeffs_3S = np.array([1.0])

cu_alpha_4S = np.array([560.088, 56.6486, 35.4258, 11.0546, 2.30682, 0.951429, 0.145184, 0.050325])
cu_coeffs_4S = np.array([-3.33e-4, 0.00593, -0.032549, 0.211071, -0.730556, 0.177242, 1.714873, -1.604573])

cu_alpha_5S = np.array([0.0174])
cu_coeffs_5S = np.array([1.0])

# --- Cu P shells ---
cu_alpha_1P = np.array([70.9739, 17.851, 4.24679, 1.8776, 0.793335, 0.193476, 0.057393])
cu_coeffs_1P = np.array([0.003682, -0.082128, 0.375379, 0.508409, 0.239095, 0.01585, -0.002527])

cu_alpha_2P = np.array([70.9739, 17.851, 4.24679, 1.8776, 0.793335, 0.193476, 0.057393])
cu_coeffs_2P = np.array([-9.1e-4, 0.029914, -0.166094, -0.278, 0.157584, 0.830488, 0.165869])

cu_alpha_3P = np.array([70.9739, 17.851, 4.24679, 1.8776, 0.793335, 0.193476, 0.057393])
cu_coeffs_3P = np.array([-6.28e-4, 0.016563, -0.084572, -0.141283, -0.003571, 0.519005, 0.592475])

cu_alpha_4P = np.array([0.057393])
cu_coeffs_4P = np.array([1.0])

cu_alpha_5P = np.array([0.017])
cu_coeffs_5P = np.array([1.0])

# --- Cu D shells ---
cu_alpha_1D = np.array([60.3804, 19.1121, 6.95288, 2.60994, 0.922567, 0.283642])
cu_coeffs_1D = np.array([0.017564, 0.099134, 0.271171, 0.40618, 0.381427, 0.200626])

cu_alpha_2D = np.array([60.3804, 19.1121, 6.95288, 2.60994, 0.922567, 0.283642])
cu_coeffs_2D = np.array([-0.022286, -0.128274, -0.362797, -0.325722, 0.327087, 0.656809])

cu_alpha_3D = np.array([0.283642])
cu_coeffs_3D = np.array([1.0])

cu_alpha_4D = np.array([0.0872])
cu_coeffs_4D = np.array([1.0])

# --- Cu F shells ---
cu_alpha_1F = np.array([2.7482])
cu_coeffs_1F = np.array([1.0])

cu_alpha_2F = np.array([0.9066])
cu_coeffs_2F = np.array([1.0])

# --- O S shells ---
o_alpha_1S = np.array([1.172e4, 1.759e3, 400.8, 113.7, 37.03, 13.27, 5.025, 1.013, 0.3023])
o_coeffs_1S = np.array([0.00071, 0.00547, 0.027837, 0.1048, 0.283062, 0.448719, 0.270952, 0.015458, -0.002585])

o_alpha_2S = np.array([1.172e4, 1.759e3, 400.8, 113.7, 37.03, 13.27, 5.025, 1.013, 0.3023])
o_coeffs_2S = np.array([-0.00016, -0.001263, -0.006267, -0.025716, -0.070924, -0.165411, -0.116955, 0.557368, 0.572759])

o_alpha_3S = np.array([0.3023])
o_coeffs_3S = np.array([1.0])

o_alpha_4S = np.array([0.07896])
o_coeffs_4S = np.array([1.0])

# --- O P shells ---
o_alpha_1P = np.array([17.7, 3.854, 1.046, 0.2753])
o_coeffs_1P = np.array([0.043018, 0.228913, 0.508728, 0.460531])

o_alpha_2P = np.array([0.2753])
o_coeffs_2P = np.array([1.0])

o_alpha_3P = np.array([0.06856])
o_coeffs_3P = np.array([1.0])

# --- O D shells ---
o_alpha_1D = np.array([1.185])
o_coeffs_1D = np.array([1.0])

o_alpha_2D = np.array([0.332])
o_coeffs_2D = np.array([1.0])

DO_coeffs_b2_L = np.array([
    -6.94442091e-16,
    1.54322894e-14,
    -8.16556534e-15,
    -6.85466315e-15,
    6.89892231e-16,
    -2.00622426e-16,
    3.92677557e-02,
    2.24687311e-15,
    -4.36093808e-15,
    6.35512242e-03,
    -4.95935142e-15,
    8.22533601e-15,
    -2.24451787e-01,
    2.46997604e-14,
    -3.21570989e-15,
    1.67383041e-02,
    -1.11146717e-14,
    -1.76453278e-15,
    -5.77952387e-02,
    4.37137048e-15,
    8.80681088e-16,
    3.74719945e-01,
    -3.69808251e-15,
    -1.51255806e-15,
    -1.18363905e-15,
    -6.60329875e-16,
    8.61318631e-04,
    -3.95923740e-16,
    1.25755553e-16,
    -6.15662073e-17,
    5.12234002e-17,
    7.77636248e-03,
    4.10606528e-15,
    -1.11220031e-15,
    -7.24981056e-17,
    3.56537965e-15,
    3.46011275e-02,
    -6.58705002e-17,
    6.74935290e-16,
    1.37075857e-16,
    7.21917942e-10,
    -9.58442026e-17,
    -5.44486028e-04,
    -1.40521424e-16,
    -4.12238088e-18,
    9.05886080e-18,
    3.18989045e-17,
    -1.34958316e-09,
    1.13568400e-16,
    -9.65162476e-03,
    9.67601693e-16,
    -1.49803523e-16,
    2.91498133e-17,
    2.18169470e-16,
    8.46891069e-16,
    -1.59381904e-16,
    -1.15497001e-14,
    1.52600728e-15,
    9.82951747e-15,
    -7.16525899e-01,
    -5.66361030e-15,
    -1.74915127e-15,
    -6.38913608e-02,
    8.08100817e-15,
    3.70655324e-15,
    -2.00796324e-01,
    -1.94113959e-14,
    1.43967332e-16,
    4.18305726e-03,
    2.85210011e-16,
    -3.39097891e-17,
    -2.75850999e-17,
    -6.90926008e-17,
    1.83280178e-02,
    -1.88694495e-15,
    -3.68335650e-16,
    6.41325405e-17
])

DO_coeffs_b2_R = np.array([
    -5.67517534e-16,
    1.37168744e-14,
    -2.61636403e-15,
    -3.83592093e-15,
    -1.23809657e-15,
    -1.96478612e-16,
    3.61566216e-02,
    7.36131034e-16,
    -3.01992063e-15,
    -1.86879450e-03,
    -6.45379424e-16,
    6.14249301e-15,
    -1.81901093e-01,
    1.05292060e-14,
    -3.64010491e-15,
    3.48739665e-02,
    -5.97244272e-15,
    -9.95646768e-16,
    -2.51496640e-02,
    1.78901003e-15,
    1.20846480e-15,
    3.85856779e-01,
    -5.51793273e-15,
    -1.75974668e-15,
    -1.12521797e-15,
    -2.91605802e-16,
    -1.36540925e-03,
    -5.05743184e-17,
    1.05295127e-16,
    -3.53994058e-17,
    9.30025278e-17,
    5.25478068e-04,
    1.74917271e-15,
    -7.09962143e-16,
    -4.01156839e-17,
    1.70629576e-15,
    2.30089831e-02,
    7.71045511e-16,
    -1.38766219e-16,
    -7.78683915e-18,
    8.60088397e-10,
    -9.17465035e-17,
    -4.76234032e-04,
    -7.64734718e-17,
    -4.10458015e-18,
    9.42446733e-18,
    2.19693069e-17,
    -2.18539923e-09,
    9.33836303e-17,
    -9.34782368e-03,
    7.87127947e-16,
    -2.49434071e-16,
    2.37541421e-17,
    1.59799919e-16,
    5.59960385e-16,
    -1.17852092e-16,
    -6.98058506e-15,
    -3.80705740e-16,
    1.00430907e-14,
    -7.62587387e-01,
    -7.94540820e-15,
    -9.18810199e-16,
    -6.40185831e-02,
    3.20769696e-15,
    3.19677164e-15,
    -1.81053453e-01,
    -8.60204916e-15,
    1.36939993e-16,
    4.50762642e-03,
    2.89632175e-16,
    -7.90791320e-18,
    -2.68936331e-17,
    -3.93183828e-17,
    1.76001171e-02,
    -1.03227060e-15,
    -6.25103850e-16,
    1.33959388e-17
])

DO_coeffs_b1_L = np.array([
    3.28874313e-16,
    -3.33496034e-15,
    1.64696964e-13,
    2.08116690e-15,
    -1.10721182e-14,
    -3.92677534e-02,
    8.45616805e-17,
    -3.12285122e-15,
    -6.35514044e-03,
    -7.12018444e-15,
    8.05790473e-15,
    2.24451844e-01,
    1.56007133e-14,
    -3.97047781e-14,
    -1.67383195e-02,
    -7.97721882e-15,
    1.49560514e-13,
    5.77952419e-02,
    1.08453927e-15,
    -1.16047452e-14,
    1.91120792e-15,
    -7.73733204e-16,
    -1.53127795e-15,
    -3.74720372e-01,
    -1.36669792e-15,
    2.60535889e-17,
    1.47750931e-16,
    5.48962605e-16,
    -8.61324515e-04,
    -1.95388245e-16,
    2.31000504e-15,
    -5.79539922e-16,
    -7.39336323e-15,
    -7.77638753e-03,
    -8.55000428e-17,
    -1.49458444e-14,
    -1.05359434e-15,
    6.34606496e-14,
    -3.46011603e-02,
    -1.30571521e-15,
    -1.52468458e-17,
    6.83568851e-17,
    -5.64412754e-17,
    4.99810423e-17,
    5.44486930e-04,
    1.72315692e-17,
    6.38270907e-10,
    1.77925947e-17,
    -3.35611776e-17,
    1.30666930e-16,
    -2.15425433e-16,
    9.65162970e-03,
    1.93506973e-17,
    -1.17833323e-09,
    -1.15958665e-15,
    -4.17720741e-16,
    1.42901594e-14,
    -2.09675347e-13,
    7.16525741e-01,
    1.00441894e-14,
    5.02055686e-15,
    6.38913666e-02,
    -4.79125998e-16,
    -5.21147492e-15,
    2.00796320e-01,
    4.03198203e-15,
    3.49887218e-14,
    -1.29187528e-16,
    3.92693182e-17,
    -2.09659276e-16,
    -4.18305449e-03,
    -3.03225438e-17,
    7.74488202e-16,
    -4.91249516e-17,
    4.07730923e-15,
    -1.83280211e-02,
    5.20602657e-16
])

DO_coeffs_b1_R = np.array([
    1.03592157e-16,
    -5.81972234e-15,
    8.18574859e-14,
    1.27371952e-15,
    -7.36789917e-15,
    -3.61566170e-02,
    8.65745417e-17,
    -1.55613578e-15,
    1.86874955e-03,
    -6.23925440e-15,
    5.47235750e-15,
    1.81901127e-01,
    1.41098088e-14,
    -2.40922044e-14,
    -3.48739798e-02,
    -7.81578655e-15,
    7.38684569e-14,
    2.51496624e-02,
    1.71467137e-15,
    -6.15295520e-15,
    8.56311864e-16,
    -1.10731589e-15,
    1.29823900e-15,
    -3.85857047e-01,
    -1.49725417e-15,
    2.32222208e-18,
    1.75187911e-16,
    3.64988633e-16,
    1.36540431e-03,
    -8.94426195e-17,
    1.12924353e-15,
    -4.25654993e-16,
    -3.60182752e-15,
    -5.25477146e-04,
    -1.14300558e-16,
    -6.80655073e-15,
    -6.29621637e-16,
    2.91940791e-14,
    -2.30089870e-02,
    -6.14531540e-16,
    -2.34467655e-17,
    7.97351891e-17,
    -5.08487529e-17,
    8.79846310e-18,
    4.76232904e-04,
    1.85161207e-17,
    8.89589895e-10,
    -7.57192691e-18,
    -1.60365605e-17,
    3.87429584e-17,
    -1.62113259e-16,
    9.34781514e-03,
    1.03116877e-17,
    -2.34098590e-09,
    -7.73039003e-16,
    -6.54630781e-16,
    8.86366927e-15,
    -9.91973544e-14,
    7.62587340e-01,
    1.03667931e-14,
    6.05232635e-15,
    6.40185595e-02,
    1.64147935e-16,
    -2.93380795e-15,
    1.81053434e-01,
    2.76965944e-15,
    1.88953511e-14,
    -6.89835221e-17,
    7.66481233e-17,
    -1.79744907e-16,
    -4.50762041e-03,
    -1.25739786e-17,
    5.77775919e-16,
    -2.89940224e-16,
    2.01440663e-15,
    -1.76001042e-02,
    2.89323614e-16
])

DO_coeffs_a1_L = np.array([
    2.81683601e-03,
    -8.44160150e-01,
    7.62010796e-02,
    1.15442635e-03,
    -1.03450947e-01,
    2.61741172e-15,
    4.87390669e-16,
    1.90430124e-02,
    -2.35859781e-14,
    1.70848261e-14,
    -4.49532781e-02,
    4.72449955e-14,
    -2.13543686e-14,
    8.36156087e-02,
    -2.42230803e-14,
    1.42599106e-14,
    2.36091911e-02,
    1.73482336e-15,
    1.85393272e-15,
    4.60466570e-02,
    2.15498207e-15,
    4.98394658e-15,
    3.66901353e-01,
    3.81808720e-15,
    4.54909975e-09,
    5.45362087e-17,
    -6.77213525e-17,
    1.62088185e-04,
    3.48418443e-16,
    2.28514406e-10,
    7.07622212e-17,
    -5.95284635e-16,
    7.06102235e-03,
    -9.06061626e-16,
    -1.12769700e-09,
    -1.16201636e-15,
    7.98219930e-15,
    2.75880516e-02,
    4.49783208e-15,
    -6.16115062e-10,
    -1.88040530e-17,
    5.17605553e-17,
    -9.11936395e-17,
    -6.61522050e-04,
    1.09281555e-16,
    3.16929905e-10,
    3.90654964e-17,
    1.36367503e-16,
    6.59380408e-17,
    9.40617081e-16,
    -8.30480590e-03,
    -9.51504151e-16,
    4.95822385e-10,
    -4.45526181e-17,
    3.19882161e-04,
    9.71662396e-02,
    2.06903207e-02,
    1.63212688e-01,
    -6.89958577e-15,
    -1.17442472e-14,
    3.79712481e-01,
    1.10681413e-16,
    2.72231209e-16,
    1.60932727e-03,
    -1.06865662e-14,
    -2.05889005e-14,
    3.14262200e-02,
    -2.50670278e-16,
    -2.63754278e-16,
    -1.54627744e-03,
    7.12133238e-18,
    -8.56554038e-11,
    -5.71663957e-17,
    1.03131538e-15,
    4.10371806e-03,
    3.12118885e-16,
    -5.73112255e-10
])

DO_coeffs_a1_R = np.array([
    3.85903190e-03,
    -8.03292957e-01,
    1.27137139e-01,
    -6.56816389e-03,
    -5.03541798e-02,
    2.01073403e-15,
    3.22924388e-16,
    3.62225932e-02,
    -2.06797702e-14,
    1.35898768e-14,
    -4.28996322e-02,
    3.91054342e-14,
    -1.92504912e-14,
    1.35535735e-02,
    -1.56850229e-14,
    7.89387428e-15,
    7.80106253e-03,
    -4.97088219e-17,
    -4.35931773e-16,
    2.44106732e-02,
    1.42871048e-15,
    5.36218180e-15,
    4.03268746e-01,
    3.62721432e-15,
    -8.55295921e-09,
    8.53241153e-17,
    -2.24835547e-16,
    -1.77081955e-03,
    6.97427601e-16,
    6.69708379e-10,
    -1.05581758e-16,
    1.95135221e-17,
    6.06801371e-03,
    -2.11404415e-15,
    4.04485655e-10,
    -3.57677695e-17,
    2.52407790e-15,
    2.92718852e-02,
    1.02100472e-14,
    3.44397485e-12,
    -2.76032596e-17,
    5.10297386e-17,
    -1.22438109e-16,
    -8.18887662e-04,
    1.07199893e-16,
    -8.65200396e-11,
    5.28525364e-17,
    1.22875591e-16,
    2.73323550e-17,
    8.92819161e-16,
    -1.07472456e-02,
    -8.98237230e-16,
    1.88521407e-10,
    -5.71429240e-17,
    8.91680077e-04,
    1.10717832e-01,
    1.96320528e-02,
    1.66214581e-01,
    -9.33187995e-15,
    -1.41030156e-14,
    5.05206548e-01,
    -1.83554143e-15,
    -1.59892996e-15,
    1.31918227e-02,
    -9.32621005e-15,
    -4.92488327e-15,
    5.27977109e-02,
    -2.33370353e-16,
    -9.80229436e-17,
    -2.92537503e-03,
    7.88447352e-17,
    8.85740504e-11,
    -6.29617127e-17,
    8.70130642e-16,
    1.36584197e-03,
    -3.88290595e-17,
    5.94997309e-10
])

b2_L_norm = 0.9250
b2_R_norm = 0.9280
b1_L_norm = 0.9250
b1_R_norm = 0.9280
a1_L_norm = 0.9734
a1_R_norm = 0.9296

def build_DO(DO_coeffs, x, y, z, dV, recenter=True, print_geom=False):
    global R_A, R_B

    basis_info = [
    # ---- Cu S orbitals ----
        (cu_alpha_1S, cu_coeffs_1S, 0, 0, 0, R_A, 0),
        (cu_alpha_2S, cu_coeffs_2S, 0, 0, 0, R_A, 1),
        (cu_alpha_3S, cu_coeffs_3S, 0, 0, 0, R_A, 2),
        (cu_alpha_4S, cu_coeffs_4S, 0, 0, 0, R_A, 3),
        (cu_alpha_5S, cu_coeffs_5S, 0, 0, 0, R_A, 4),

    # ---- Cu P orbitals ----
        (cu_alpha_1P, cu_coeffs_1P, 1, 0, 0, R_A, 5),
        (cu_alpha_1P, cu_coeffs_1P, 0, 1, 0, R_A, 6),
        (cu_alpha_1P, cu_coeffs_1P, 0, 0, 1, R_A, 7),

        (cu_alpha_2P, cu_coeffs_2P, 1, 0, 0, R_A, 8),
        (cu_alpha_2P, cu_coeffs_2P, 0, 1, 0, R_A, 9),
        (cu_alpha_2P, cu_coeffs_2P, 0, 0, 1, R_A, 10),

        (cu_alpha_3P, cu_coeffs_3P, 1, 0, 0, R_A, 11),
        (cu_alpha_3P, cu_coeffs_3P, 0, 1, 0, R_A, 12),
        (cu_alpha_3P, cu_coeffs_3P, 0, 0, 1, R_A, 13),

        (cu_alpha_4P, cu_coeffs_4P, 1, 0, 0, R_A, 14),
        (cu_alpha_4P, cu_coeffs_4P, 0, 1, 0, R_A, 15),
        (cu_alpha_4P, cu_coeffs_4P, 0, 0, 1, R_A, 16),
    
        (cu_alpha_5P, cu_coeffs_5P, 1, 0, 0, R_A, 17),
        (cu_alpha_5P, cu_coeffs_5P, 0, 1, 0, R_A, 18),
        (cu_alpha_5P, cu_coeffs_5P, 0, 0, 1, R_A, 19),

    # ---- Cu D orbitals ----
        (cu_alpha_1D, cu_coeffs_1D, 1, 1, 0, R_A, 20),
        (cu_alpha_1D, cu_coeffs_1D, 0, 1, 1, R_A, 21),
    # d_z^2 (22) done manually
        (cu_alpha_1D, cu_coeffs_1D, 1, 0, 1, R_A, 23),
    # d_x2-y2 (24) done manually

        (cu_alpha_2D, cu_coeffs_2D, 1, 1, 0, R_A, 25),
        (cu_alpha_2D, cu_coeffs_2D, 0, 1, 1, R_A, 26),
    # d_z^2 (27) done manually
        (cu_alpha_2D, cu_coeffs_2D, 1, 0, 1, R_A, 28),
    # d_x2-y2 (29) done manually

        (cu_alpha_3D, cu_coeffs_3D, 1, 1, 0, R_A, 30),
        (cu_alpha_3D, cu_coeffs_3D, 0, 1, 1, R_A, 31),
    # d_z^2 (32) done manually
        (cu_alpha_3D, cu_coeffs_3D, 1, 0, 1, R_A, 33),
    # d_x2-y2 (34) done manually

        (cu_alpha_4D, cu_coeffs_4D, 1, 1, 0, R_A, 35),
        (cu_alpha_4D, cu_coeffs_4D, 0, 1, 1, R_A, 36),
    # d_z^2 (37) done manually
        (cu_alpha_4D, cu_coeffs_4D, 1, 0, 1, R_A, 38),
    # d_x2-y2 (39) done manually

    # f_y3x2 (40) done manually
        (cu_alpha_1F, cu_coeffs_1F, 1, 1, 1, R_A, 41),
    # f_yz2 (42) done manually
    # f_z3 (43) done manually
    # f_xz2 (44) done manually
    # f_x2y2z (45) done manually
    # f_x3y2 (46) done manually

    # f_y3x2 (47) done manually
        (cu_alpha_2F, cu_coeffs_2F, 1, 1, 1, R_A, 48),
    # f_yz2 (49) done manually
    # f_z3 (50) done manually
    # f_xz2 (51) done manually
    # f_x2y2z (52) done manually
    # f_x3y2 (53) done manually

    # ---- O S orbitals ----
        (o_alpha_1S, o_coeffs_1S, 0, 0, 0, R_B, 54),
        (o_alpha_2S, o_coeffs_2S, 0, 0, 0, R_B, 55),
        (o_alpha_3S, o_coeffs_3S, 0, 0, 0, R_B, 56),
        (o_alpha_4S, o_coeffs_4S, 0, 0, 0, R_B, 57),

    # ---- O P orbitals ----
        (o_alpha_1P, o_coeffs_1P, 1, 0, 0, R_B, 58),
        (o_alpha_1P, o_coeffs_1P, 0, 1, 0, R_B, 59),
        (o_alpha_1P, o_coeffs_1P, 0, 0, 1, R_B, 60),

        (o_alpha_2P, o_coeffs_2P, 1, 0, 0, R_B, 61),
        (o_alpha_2P, o_coeffs_2P, 0, 1, 0, R_B, 62),
        (o_alpha_2P, o_coeffs_2P, 0, 0, 1, R_B, 63),

        (o_alpha_3P, o_coeffs_3P, 1, 0, 0, R_B, 64),
        (o_alpha_3P, o_coeffs_3P, 0, 1, 0, R_B, 65),
        (o_alpha_3P, o_coeffs_3P, 0, 0, 1, R_B, 66),

    # ---- O D orbitals ----
        (o_alpha_1D, o_coeffs_1D, 1, 1, 0, R_B, 67),
        (o_alpha_1D, o_coeffs_1D, 0, 1, 1, R_B, 68),
    # d_z^2 (69) done manually
        (o_alpha_1D, o_coeffs_1D, 1, 0, 1, R_B, 70),
    # d_x2-y2 (71) done manually

        (o_alpha_2D, o_coeffs_2D, 1, 1, 0, R_B, 72),
        (o_alpha_2D, o_coeffs_2D, 0, 1, 1, R_B, 73),
    # d_z^2 (74) done manually
        (o_alpha_2D, o_coeffs_2D, 1, 0, 1, R_B, 75),
    # d_x2-y2 (76) done manually
]
    
    DO = np.zeros_like(x, dtype=np.float64)
    threshold = 1e-5
    for alphas, coeffs, a, b, c, center, i in basis_info:
        coeff = DO_coeffs[i]
        if abs(coeff) >= threshold:
            ao = AO(alphas, coeffs, a, b, c, center, x, y, z, dV)
            DO += coeff * ao

    def norm(orb):
        return 1 / np.sqrt(integrate_3d(abs(orb)**2, dV))
    
    # Handle weird orbitals manually
    if abs(DO_coeffs[22]) >= threshold:
        ao_zz = AO(cu_alpha_1D, cu_coeffs_1D, 0, 0, 2, R_A, x, y, z, dV)
        ao_xx = AO(cu_alpha_1D, cu_coeffs_1D, 2, 0, 0, R_A, x, y, z, dV)
        ao_yy = AO(cu_alpha_1D, cu_coeffs_1D, 0, 2, 0, R_A, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[22] * normed_AO

    if abs(DO_coeffs[24]) >= threshold:
        ao_xx = AO(cu_alpha_1D, cu_coeffs_1D, 2, 0, 0, R_A, x, y, z, dV)
        ao_yy = AO(cu_alpha_1D, cu_coeffs_1D, 0, 2, 0, R_A, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[24] * normed_AO

    if abs(DO_coeffs[27]) >= threshold:
        ao_zz = AO(cu_alpha_2D, cu_coeffs_2D, 0, 0, 2, R_A, x, y, z, dV)
        ao_xx = AO(cu_alpha_2D, cu_coeffs_2D, 2, 0, 0, R_A, x, y, z, dV)
        ao_yy = AO(cu_alpha_2D, cu_coeffs_2D, 0, 2, 0, R_A, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[27] * normed_AO

    if abs(DO_coeffs[29]) >= threshold:
        ao_xx = AO(cu_alpha_2D, cu_coeffs_2D, 2, 0, 0, R_A, x, y, z, dV)
        ao_yy = AO(cu_alpha_2D, cu_coeffs_2D, 0, 2, 0, R_A, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[29] * normed_AO

    if abs(DO_coeffs[32]) >= threshold:
        ao_zz = AO(cu_alpha_3D, cu_coeffs_3D, 0, 0, 2, R_A, x, y, z, dV)
        ao_xx = AO(cu_alpha_3D, cu_coeffs_3D, 2, 0, 0, R_A, x, y, z, dV)
        ao_yy = AO(cu_alpha_3D, cu_coeffs_3D, 0, 2, 0, R_A, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[32] * normed_AO

    if abs(DO_coeffs[34]) >= threshold:
        ao_xx = AO(cu_alpha_3D, cu_coeffs_3D, 2, 0, 0, R_A, x, y, z, dV)
        ao_yy = AO(cu_alpha_3D, cu_coeffs_3D, 0, 2, 0, R_A, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[34] * normed_AO

    if abs(DO_coeffs[37]) >= threshold:
        ao_zz = AO(cu_alpha_4D, cu_coeffs_4D, 0, 0, 2, R_A, x, y, z, dV)
        ao_xx = AO(cu_alpha_4D, cu_coeffs_4D, 2, 0, 0, R_A, x, y, z, dV)
        ao_yy = AO(cu_alpha_4D, cu_coeffs_4D, 0, 2, 0, R_A, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[37] * normed_AO

    if abs(DO_coeffs[39]) >= threshold:
        ao_xx = AO(cu_alpha_4D, cu_coeffs_4D, 2, 0, 0, R_A, x, y, z, dV)
        ao_yy = AO(cu_alpha_4D, cu_coeffs_4D, 0, 2, 0, R_A, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[39] * normed_AO

    if abs(DO_coeffs[40]) >= threshold:
        ao_yx2 = AO(cu_alpha_1F, cu_coeffs_1F, 2, 1, 0, R_A, x, y, z, dV)
        ao_y3  = AO(cu_alpha_1F, cu_coeffs_1F, 0, 3, 0, R_A, x, y, z, dV)
        orb = np.sqrt(5/8) * (3 * ao_yx2 - ao_y3)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[40] * normed_AO

    if abs(DO_coeffs[42]) >= threshold:
        ao_yz2 = AO(cu_alpha_1F, cu_coeffs_1F, 0, 1, 2, R_A, x, y, z, dV)
        ao_y3  = AO(cu_alpha_1F, cu_coeffs_1F, 0, 3, 0, R_A, x, y, z, dV)
        ao_yx2 = AO(cu_alpha_1F, cu_coeffs_1F, 2, 1, 0, R_A, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_yz2 - ao_y3 - ao_yx2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[42] * normed_AO

    if abs(DO_coeffs[43]) >= threshold:
        ao_z3   = AO(cu_alpha_1F, cu_coeffs_1F, 0, 0, 3, R_A, x, y, z, dV)
        ao_x2z  = AO(cu_alpha_1F, cu_coeffs_1F, 2, 0, 1, R_A, x, y, z, dV)
        ao_y2z  = AO(cu_alpha_1F, cu_coeffs_1F, 0, 2, 1, R_A, x, y, z, dV)
        orb = 0.5 * (2 * ao_z3 - 3 * ao_x2z - 3 * ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[43] * normed_AO

    if abs(DO_coeffs[44]) >= threshold:
        ao_xz2 = AO(cu_alpha_1F, cu_coeffs_1F, 1, 0, 2, R_A, x, y, z, dV)
        ao_x3   = AO(cu_alpha_1F, cu_coeffs_1F, 3, 0, 0, R_A, x, y, z, dV)
        ao_xy2  = AO(cu_alpha_1F, cu_coeffs_1F, 1, 2, 0, R_A, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_xz2 - ao_x3 - ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[44] * normed_AO

    if abs(DO_coeffs[45]) >= threshold:
        ao_x2z = AO(cu_alpha_1F, cu_coeffs_1F, 2, 0, 1, R_A, x, y, z, dV)
        ao_y2z = AO(cu_alpha_1F, cu_coeffs_1F, 0, 2, 1, R_A, x, y, z, dV)
        orb = (np.sqrt(15) / 2) * (ao_x2z - ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[45] * normed_AO

    if abs(DO_coeffs[46]) >= threshold:
        ao_x3 = AO(cu_alpha_1F, cu_coeffs_1F, 3, 0, 0, R_A, x, y, z, dV)
        ao_xy2 = AO(cu_alpha_1F, cu_coeffs_1F, 1, 2, 0, R_A, x, y, z, dV)
        orb = np.sqrt(5/8) * (ao_x3 - 3 * ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[46] * normed_AO

    if abs(DO_coeffs[47]) >= threshold:
        ao_yx2 = AO(cu_alpha_2F, cu_coeffs_2F, 2, 1, 0, R_A, x, y, z, dV)
        ao_y3  = AO(cu_alpha_2F, cu_coeffs_2F, 0, 3, 0, R_A, x, y, z, dV)
        orb = np.sqrt(5/8) * (3 * ao_yx2 - ao_y3)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[47] * normed_AO

    if abs(DO_coeffs[49]) >= threshold:
        ao_yz2 = AO(cu_alpha_2F, cu_coeffs_2F, 0, 1, 2, R_A, x, y, z, dV)
        ao_y3  = AO(cu_alpha_2F, cu_coeffs_2F, 0, 3, 0, R_A, x, y, z, dV)
        ao_yx2 = AO(cu_alpha_2F, cu_coeffs_2F, 2, 1, 0, R_A, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_yz2 - ao_y3 - ao_yx2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[49] * normed_AO

    if abs(DO_coeffs[50]) >= threshold:
        ao_z3   = AO(cu_alpha_2F, cu_coeffs_2F, 0, 0, 3, R_A, x, y, z, dV)
        ao_x2z  = AO(cu_alpha_2F, cu_coeffs_2F, 2, 0, 1, R_A, x, y, z, dV)
        ao_y2z  = AO(cu_alpha_2F, cu_coeffs_2F, 0, 2, 1, R_A, x, y, z, dV)
        orb = (2 * ao_z3 - 3 * ao_x2z - 3 * ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[50] * 0.5 * normed_AO

    if abs(DO_coeffs[51]) >= threshold:
        ao_xz2 = AO(cu_alpha_2F, cu_coeffs_2F, 1, 0, 2, R_A, x, y, z, dV)
        ao_x3   = AO(cu_alpha_2F, cu_coeffs_2F, 3, 0, 0, R_A, x, y, z, dV)
        ao_xy2  = AO(cu_alpha_2F, cu_coeffs_2F, 1, 2, 0, R_A, x, y, z, dV)
        orb = np.sqrt(3/8) * (4 * ao_xz2 - ao_x3 - ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[51] * normed_AO

    if abs(DO_coeffs[52]) >= threshold:
        ao_x2z = AO(cu_alpha_2F, cu_coeffs_2F, 2, 0, 1, R_A, x, y, z, dV)
        ao_y2z = AO(cu_alpha_2F, cu_coeffs_2F, 0, 2, 1, R_A, x, y, z, dV)
        orb = (np.sqrt(15) / 2) * (ao_x2z - ao_y2z)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[52] * normed_AO

    if abs(DO_coeffs[53]) >= threshold:
        ao_x3 = AO(cu_alpha_2F, cu_coeffs_2F, 3, 0, 0, R_A, x, y, z, dV)
        ao_xy2 = AO(cu_alpha_2F, cu_coeffs_2F, 1, 2, 0, R_A, x, y, z, dV)
        orb = np.sqrt(5/8) * (ao_x3 - 3 * ao_xy2)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[53] * normed_AO

    if abs(DO_coeffs[69]) >= threshold:
        ao_zz = AO(o_alpha_1D, o_coeffs_1D, 0, 0, 2, R_B, x, y, z, dV)
        ao_xx = AO(o_alpha_1D, o_coeffs_1D, 2, 0, 0, R_B, x, y, z, dV)
        ao_yy = AO(o_alpha_1D, o_coeffs_1D, 0, 2, 0, R_B, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[69] * normed_AO

    if abs(DO_coeffs[71]) >= threshold:
        ao_xx = AO(o_alpha_1D, o_coeffs_1D, 2, 0, 0, R_B, x, y, z, dV)
        ao_yy = AO(o_alpha_1D, o_coeffs_1D, 0, 2, 0, R_B, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[71] * normed_AO

    if abs(DO_coeffs[74]) >= threshold:
        ao_zz = AO(o_alpha_2D, o_coeffs_2D, 0, 0, 2, R_B, x, y, z, dV)
        ao_xx = AO(o_alpha_2D, o_coeffs_2D, 2, 0, 0, R_B, x, y, z, dV)
        ao_yy = AO(o_alpha_2D, o_coeffs_2D, 0, 2, 0, R_B, x, y, z, dV)
        orb = 0.5 * (2 * ao_zz - ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb
        DO += DO_coeffs[74] * normed_AO

    if abs(DO_coeffs[76]) >= threshold:
        ao_xx = AO(o_alpha_2D, o_coeffs_2D, 2, 0, 0, R_B, x, y, z, dV)
        ao_yy = AO(o_alpha_2D, o_coeffs_2D, 0, 2, 0, R_B, x, y, z, dV)
        orb = (np.sqrt(3)/2) * (ao_xx - ao_yy)
        N = norm(orb)
        normed_AO = N * orb   
        DO += DO_coeffs[76] * normed_AO

    Norm = norm(DO)
    DO /= Norm

    if recenter:
        centroid = recenter_DO(DO, x, y, z, dV)
        
        # Shift grid
        x_new = x - centroid[0]
        y_new = y - centroid[1]
        z_new = z - centroid[2]

        # Shift molecular geometry
        R_A = R_A - centroid
        R_B = R_B - centroid

        if print_geom:
            bohr_to_angstrom = 0.529177
            R_A_ang = R_A * bohr_to_angstrom
            R_B_ang = R_B * bohr_to_angstrom

            print("\nShifting molecular geometry to the new center")
            print("New molecular geometry is (in Ångströms):")
            print(" atom         X             Y             Z")
            print("   1      {: .6f}      {: .6f}      {: .6f}".format(*R_A_ang))
            print("   2      {: .6f}      {: .6f}      {: .6f}".format(*R_B_ang))

        # Rebuild on shifted grid
        return build_DO(DO_coeffs, x_new, y_new, z_new, dV, recenter=False)

    return DO
