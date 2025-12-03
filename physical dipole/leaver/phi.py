import numpy as np

def compute_Phi_m(m, phi):
    return (1 / np.sqrt(2 * np.pi)) * np.exp(1j * m * phi)