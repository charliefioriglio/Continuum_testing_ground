# Point Dipole Program (Python)

A faithful Python translation of the MATLAB point dipole photodetachment
program. This package computes photoelectron angular distributions (β
parameters) for anion photodetachment in the presence of a permanent
molecular dipole using the analytical Clebsch-Gordan coupling approach.

## Installation

```bash
# From the code directory:
pip install numpy scipy matplotlib scikit-image
```

## Package Structure

```
point dipole program python/
├── __init__.py              # Package metadata
├── utilities.py             # dfac, tricoeff, x (low-level helpers)
├── wigner.py                # Wigner 3-j symbols, Clebsch-Gordan coefficients
├── clebsch_array.py         # Pre-compute CG arrays for fast β evaluation
├── harmonics.py             # Spherical harmonics (real and complex)
├── parper.py                # Parallel/perpendicular harmonic values
├── integration.py           # 3-D Simpson's rule integration
├── dipole_eigensystem.py    # Point-dipole Hamiltonian eigenproblem
├── continuum.py             # Continuum wavefunction construction
├── dyson.py                 # Dyson orbital builder from Gaussian basis
├── beta.py                  # Analytical β parameter calculation
├── main.py                  # Main driver script
├── plotting.py              # Isosurface visualization
├── CuO_NTOs.py              # Example molecule data (CuO)
└── README.md                # This file
```

## MATLAB → Python File Mapping

| MATLAB File         | Python Module            | Description                                |
|---------------------|-------------------------|--------------------------------------------|
| `dfac.m`            | `utilities.dfac`        | Double factorial                           |
| `tricoeff.m`        | `utilities.tricoeff`    | Triangle coefficient for Wigner symbols    |
| `x.m`               | `utilities.x`           | Summation helper for Wigner 3-j            |
| `wigner3j.m`        | `wigner.wigner3j`       | Wigner 3-j symbol                          |
| `clebschgordan.m`   | `wigner.clebschgordan`  | Clebsch-Gordan coefficient                 |
| `clebarray1.m`      | `clebsch_array.clebarray1` | Pre-computed CG array                   |
| `prefactor.m`       | `harmonics.prefactor`   | Normalization prefactor                    |
| `thph.m`            | `harmonics.thph`        | Cartesian → spherical coordinate conversion|
| `Ylm.m`             | `harmonics.Ylm`         | Real (cubic) spherical harmonics           |
| `Ylmc.m`            | `harmonics.Ylmc`        | Complex spherical harmonics                |
| `parper.m`          | `parper.parper`         | P_l^m(±1) for parallel/perpendicular       |
| `simp3D.m`          | `integration.simp3D`    | 3-D Simpson integration                    |
| `pointdpmatrixc.m`  | `dipole_eigensystem.pointdpmatrixc` | Eigenproblem solver      |
| `continuum.m`       | `continuum.continuum`   | Continuum wavefunction                     |
| `Dyson.m`           | `dyson.dyson`           | Dyson orbital builder                      |
| `beta2.m`           | `beta.beta2`            | β parameter calculation                    |
| `extendingdpmodel.m`| `main.main`             | Main driver                                |
| `figplot.m`         | `plotting.figplot`      | Isosurface plotting                        |
| `CuO_NTOs.m`        | `CuO_NTOs.CuO_NTOs`     | CuO molecule data                          |

## Quick Start

```python
import numpy as np
from point_dipole_program_python import main

# Run the full calculation
betas = main.main()
```

Or modify `main.py` to change parameters:

```python
# Edit these in main.py:
LMAX = 6                  # highest partial wave
E_RANGE = (0.1, 1.0)      # kinetic energy range (eV)
E_STEP = 0.2              # energy step
D_START = 0.639           # dipole moment (a.u.)
GRIDMAX = 18.897260       # grid extent (a₀)
PS = 100                  # grid points
```

## Key Algorithms

### Point-Dipole Eigenproblem

The angular eigenproblem is solved by diagonalizing the matrix:

```
H_{ll'} = δ_{ll'} l(l+1) - 2D √(4π/3) √((2l+1)·3·(2l'+1)/(4π))
          × (l 1 l' | 0 0 0) × (l 1 l' | -λ 0 λ)
```

where D is the dipole moment and λ is the projection quantum number.

### Continuum Wavefunction

The point-dipole continuum uses Gallup normalization (PRA 23, 632, 1981):

```
ψ_{N,λ}(r) = √r · √(π/2k) · J_{ν_N+1/2}(kr) · Ω_{N,λ}(θ,φ)
```

where ν_N is the effective angular momentum from the eigenproblem.

### β Parameter

Computed analytically via nested Clebsch-Gordan sums:

```
β = 2(σ_∥ - σ_⊥) / (σ_∥ + 2σ_⊥)
```

with parallel/perpendicular cross sections from CG products and phase factors.

## Adding New Molecules

1. Create a new file (e.g., `MyMolecule.py`) with a function returning:
   - `A`: atom coordinates (Å)
   - `N`: Dyson norms
   - `CZ`: list of basis set arrays
   - `wL`, `wR`: left/right Dyson coefficients
   - `llim`: max angular momentum

2. Import it in `main.py`:
   ```python
   from .MyMolecule import MyMolecule as qcoutput
   ```

## References

- G.A. Gallup, Phys. Rev. A **23**, 632 (1981) — Point-dipole continuum
- E.P. Wigner, *Group Theory* — Wigner 3-j symbols
- Sanov & Mabbs, Int. Rev. Phys. Chem. **27**, 53 (2008) — β parameters

## License

Research code — use at your own risk.
