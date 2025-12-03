"""
main.py
=======

Main driver script for the point-dipole photodetachment calculation,
translated from MATLAB extendingdpmodel.m.

Running this script reproduces the workflow:
1. Load / compute Clebsch-Gordan coefficient arrays.
2. Build complex spherical harmonic grids.
3. Compute parallel/perpendicular harmonic evaluation points.
4. Load Dyson orbitals (left & right).
5. Loop over dipole moments and kinetic energies:
   a. Solve the point-dipole eigenproblem.
   b. Build continuum functions on the grid.
   c. Compute overlap integrals.
   d. Calculate β via analytical CG summation.

Usage
-----
Edit the ``USER DEFINED PARAMETERS`` section, then run::

    python main.py
"""

from __future__ import annotations

import numpy as np

from .beta import beta2
from .clebsch_array import clebarray1
from .continuum import build_S_grid, continuum
from .dipole_eigensystem import pointdpmatrixc
from .dyson import dyson
from .harmonics import Ylmc, thph
from .integration import simp3D
from .parper import parper

# ---------------------------------------------------------------------------
# USER DEFINED PARAMETERS
# ---------------------------------------------------------------------------
# Molecule data function (returns geometry, basis, Dyson coefficients)
from .CuO_NTOs import CuO_NTOs as qcoutput

LMAX = 6                          # highest l partial wave
E_RANGE = (0.1, 1.0)              # first and last eKE values (eV)
E_STEP = 0.2                      # step size (eV) for β calculation
D_START = 0.639                   # first dipole moment (a.u.)
D_END = 0.639                     # last dipole moment (a.u.)
D_STEP = 0.1                      # dipole moment increment
GRIDMAX = 18.897260               # grid dimension +/- in a0
PS = 100                          # number of points on the integral grid
DYSTOL = 1e-8                     # threshold for Dyson orbital coefficients
ISOVAL = 0.001                    # isovalue for plotting (not used here)
TOL = 10                          # tolerance for rounding integrals to 0

# ---------------------------------------------------------------------------
# CONSTANTS AND DERIVED VARIABLES
# ---------------------------------------------------------------------------
A0 = 0.5291772109                 # Bohr radius (Å)
NPTS = PS + 1                     # number of grid points per dimension
STEP = 2 * GRIDMAX / PS           # grid step size

# Energy / momentum arrays
N_BETA = int(round((E_RANGE[1] - E_RANGE[0]) / E_STEP)) + 1
E_VAL = np.array([E_RANGE[0] + i * E_STEP for i in range(N_BETA)])
Q_VAL = np.sqrt(2 * E_VAL / 27.211)  # momentum in a.u.

N_DIPOLES = int(round((D_END - D_START) / D_STEP)) + 1 if D_STEP > 0 else 1


def main() -> None:
    """Run the full point-dipole β calculation."""

    # -------------------- Build coordinate grid --------------------
    x1 = np.linspace(-GRIDMAX, GRIDMAX, NPTS)
    X, Y, Z = np.meshgrid(x1, x1, x1, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)

    # -------------------- Clebsch-Gordan coefficients --------------------
    print(f"Building Clebsch-Gordan coefficient arrays (lmax = {LMAX})")
    CG: dict = {}
    for l in range(LMAX + 1):
        for L in range(abs(l - 1), l + 2):
            if L > LMAX + 1:
                continue
            CG[(l, L)] = clebarray1(l, L, LMAX)

    # -------------------- Complex spherical harmonics --------------------
    print(f"Loading complex spherical harmonics (up to l = {LMAX})")
    _, theta, phi = thph(X, Y, Z)
    n_m = 2 * LMAX + 1
    Sc = np.zeros((NPTS, NPTS, NPTS, LMAX + 1, n_m), dtype=np.complex128)
    for l in range(LMAX + 1):
        for m in range(-l, l + 1):
            Sc[..., l, m + LMAX] = Ylmc(l, m, theta, phi)

    # -------------------- Parallel / perpendicular values --------------------
    print(f"Generating parallel/perpendicular values (up to l = {LMAX})")
    ppar, pper = parper(LMAX)
    # MATLAB SPV indexed as SPV(l+1, m+lmax+2, pol) where pol=1,2
    # We store as SPV[l, m+lmax+1, pol] with pol=0 (par), 1 (per)
    SPV = np.zeros((LMAX + 1, 2 * LMAX + 3, 2), dtype=np.float64)
    SPV[..., 0] = ppar
    SPV[..., 1] = pper

    # -------------------- Load Dyson orbitals --------------------
    print(f"Loading Dyson orbitals (Dystol = {DYSTOL:.2e})")
    orbL, orbR, llim = dyson(qcoutput, PS, STEP, X, Y, Z, TOL, DYSTOL)

    # -------------------- Main loop over dipole moments --------------------
    betas = np.zeros((N_BETA, N_DIPOLES + 1), dtype=np.float64)
    betas[:, 0] = E_VAL  # first column holds kinetic energies

    dind = 0
    dp_values = (
        np.arange(D_START, D_END + D_STEP / 2, D_STEP)
        if D_STEP > 0
        else np.array([D_START])
    )

    for dp in dp_values:
        dind += 1
        print(f"\n===== Dipole moment = {dp:.3f} a.u. =====")

        # ---------- Solve eigenproblem for each λ ----------
        print("Setting up and normalizing eigenfunctions/eigenvalues")
        eign = np.zeros((LMAX + 1, 2 * LMAX + 1), dtype=np.float64)
        eigv = np.zeros((LMAX + 1, LMAX + 1, 2 * LMAX + 1), dtype=np.float64)

        for lam in range(-LMAX, LMAX + 1):
            eign[:, lam + LMAX] = pointdpmatrixc(LMAX, dp, lam, 1)
            eigv[:, :, lam + LMAX] = pointdpmatrixc(LMAX, dp, lam, 2)

        # Normalize eigenvectors (sum over l of eigv^2 = 1)
        for lam in range(-LMAX, LMAX + 1):
            for N in range(abs(lam), LMAX + 1):
                norm2 = 0.0
                for l in range(abs(lam), LMAX + 1):
                    norm2 += eigv[N, l, lam + LMAX] ** 2
                if norm2 > 0:
                    eigv[N, :, lam + LMAX] /= np.sqrt(norm2)

        # ---------- Loop over kinetic energies ----------
        for i, (q, E_eV) in enumerate(zip(Q_VAL, E_VAL)):
            print(f"  k = {q:.4f} a.u., E = {E_eV:.2f} eV  ", end="")

            # Build continuum grids F and Fd
            F = np.zeros((NPTS, NPTS, NPTS, LMAX + 1, 2 * LMAX + 1), dtype=np.float64)
            Fd = np.zeros_like(F)

            for N in range(LMAX + 1):
                for lam in range(-N, N + 1):
                    F[..., N, lam + LMAX] = continuum(
                        N, lam, LMAX, eign, eigv, q, r, Sc, NPTS
                    )
                    # Fd uses (-1)^λ * continuum(N, -λ, ...)
                    Fd[..., N, lam + LMAX] = ((-1) ** lam) * continuum(
                        N, -lam, LMAX, eign, eigv, q, r, Sc, NPTS
                    )

            # ---------- Compute integrals ----------
            Integrals = np.zeros((LMAX + 1, 2 * LMAX + 1, 3), dtype=np.float64)
            Integralsd = np.zeros_like(Integrals)

            for lam in range(-LMAX, LMAX + 1):
                for N in range(abs(lam), LMAX + 1):
                    # m2 = -1, 0, +1 → Sc indices lmax-1, lmax, lmax+1 (for l=1 harmonics)
                    # Sc[..., 1, lmax + m] for Y_1^m
                    Y1m1 = Sc[..., 1, LMAX - 1]  # m = -1
                    Y10 = Sc[..., 1, LMAX]       # m = 0
                    Y1p1 = Sc[..., 1, LMAX + 1]  # m = +1

                    # Left integrals
                    Integrands = [
                        Y1m1 * orbL * F[..., N, lam + LMAX],
                        Y10 * orbL * F[..., N, lam + LMAX],
                        Y1p1 * orbL * F[..., N, lam + LMAX],
                    ]
                    # Right integrals (note sign pattern from MATLAB)
                    Integrandsd = [
                        -Y1p1 * orbR * Fd[..., N, lam + LMAX],
                        +Y10 * orbR * Fd[..., N, lam + LMAX],
                        -Y1m1 * orbR * Fd[..., N, lam + LMAX],
                    ]

                    for idx, (intL, intR) in enumerate(zip(Integrands, Integrandsd)):
                        val_L = simp3D(intL.real, STEP, STEP, STEP)
                        val_R = simp3D(intR.real, STEP, STEP, STEP)
                        Integrals[N, lam + LMAX, idx] = round(val_L, TOL)
                        Integralsd[N, lam + LMAX, idx] = round(val_R, TOL)

            # ---------- β calculation ----------
            beta_val = beta2(
                CG, Integrals, Integralsd, eigv, eign, SPV, LMAX
            )
            betas[i, dind] = beta_val
            print(f"→ β = {beta_val:.4f}")

    # -------------------- Summary output --------------------
    print("\n===== Results Summary =====")
    header = "E (eV)" + "".join(f"  D={d:.2f}" for d in dp_values)
    print(header)
    for row in betas:
        print("  ".join(f"{v:7.4f}" for v in row))

    return betas


if __name__ == "__main__":
    main()
