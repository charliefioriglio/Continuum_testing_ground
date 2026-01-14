import numpy as np
from CuO_do import (
    build_DO,
    DO_coeffs_b1_L
)

# Test point
x = np.array([1.0])
y = np.array([1.0])
z = np.array([1.0])
dV = 1.0 # Dummy

val = build_DO(DO_coeffs_b1_L, x, y, z, dV)
print(f"Python DO(b1_L) at (1,1,1): {val[0]:.8f}")
