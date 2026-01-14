import numpy as np

# Grid
L = 30
n_pts = 100
x = np.linspace(-L, L, n_pts)
y = np.linspace(-L, L, n_pts)
z = np.linspace(-L, L, n_pts)
dx = x[1] - x[0]
dV = dx**3

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Primitive dxy: x*y * exp(-alpha * r^2)
alpha = 0.5 # Test alpha
r2 = X**2 + Y**2 + Z**2
prim = X * Y * np.exp(-alpha * r2)

norm_sq = np.sum(prim**2) * dV
N_num = 1.0 / np.sqrt(norm_sq)

# Evaluate at index 51 (0.909...)
idx = 51
val = prim[idx, idx, idx]
print(f"Index 51 Coord: {x[idx]:.4f}")
print(f"Primitive Value: {val:.6f}")
print(f"Numerical Norm Factor: {N_num:.6f}")
print(f"Normalized Value: {val * N_num:.6f}")

# Analytical Norm for dxy (1,1,0)
# Formula: N = (2a/pi)^0.75 * sqrt( (4a)^L / (2lx-1)!! ... )
# L=2. lx=1, ly=1, lz=0.
# (2*1-1)!! = 1.
# (4a)^2 = 16a^2.
# N = (2a/pi)^0.75 * sqrt(16 a^2) = (2a/pi)^0.75 * 4a.
# Wait. Formula check.
# Integrated x^2 y^2 exp(-2a r^2).
# = Int x^2 exp(-2ax^2) * Int y^2 exp(-2ay^2) * Int exp(-2az^2).
# Int exp(-2au^2) = sqrt(pi/2a).
# Int u^2 exp(-2au^2) = 1/4a * sqrt(pi/2a).
# Norm^2 = (1/4a sqrt)^2 * sqrt.
# = 1/16a^2 * (pi/2a)^(3/2).
# N = 4a * (2a/pi)^0.75.
# Let's see if Python matches.
pass
