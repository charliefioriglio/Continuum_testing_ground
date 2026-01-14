import numpy as np
from CuO_do import (
    build_DO,
    DO_coeffs_b1_L
)

# Global grid definition matching calculate.py
L = 30
n_pts = 100
x = np.linspace(-L, L, n_pts)
y = np.linspace(-L, L, n_pts)
z = np.linspace(-L, L, n_pts)
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]
dV = dx * dy * dz

# Check indices (51, 51, 51)
idx = 51
print(f"Checking Grid Point: ({x[idx]:.4f}, {y[idx]:.4f}, {z[idx]:.4f})")

# We only need to pass the point, BUT build_DO calculates norm over the arrays passed.
# So we must pass the FULL arrays x, y, z to ensure correct normalization.
# build_DO assumes 3D arrays for eval?
# CuO_do.py: AO(..., x, y, z, dV).
# primitives = ... exp(-alpha * r2). 
# If x,y,z are 1D arrays, numpy broadcasting might behave differently than meshgrid?
# calculate.py passes X, Y, Z from meshgrid.
# Let's pass meshgrid to be safe and match calculate.py.

X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

val_cube = build_DO(DO_coeffs_b1_L, X, Y, Z, dV)
val = val_cube[idx, idx, idx]

print(f"Python DO(b1_L) at index ({idx},{idx},{idx}): {val:.8f}")
