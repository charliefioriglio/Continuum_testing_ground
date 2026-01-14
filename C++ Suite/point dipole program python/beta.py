"""
beta.py
=======

Analytical β parameter calculation via Clebsch-Gordan coupling, translated
from MATLAB beta2.m.

Functions
---------
beta2(CG, B, Bd, eigv, eign, SPV, lmax, kr, r)
    Compute the asymmetry parameter β using nested sums over (N, λ, l, m, L).
"""

from __future__ import annotations

import numpy as np

__all__ = ["beta2"]


def beta2(
    CG: dict,
    B: np.ndarray,
    Bd: np.ndarray,
    eigv: np.ndarray,
    eign: np.ndarray,
    SPV: np.ndarray,
    lmax: int,
    kr: float = 0.0,  # not used directly, kept for signature parity
    r: float = 0.0,   # not used directly, kept for signature parity
) -> float:
    r"""
    Compute the asymmetry parameter β using the analytical Clebsch-Gordan approach.

    Parameters
    ----------
    CG : dict
        Pre-computed Clebsch-Gordan arrays indexed by ``(l, L)`` keys, each of
        shape ``(3, 2*lmax+1, 2*lmax+3)`` indexed as ``CG[l,L][m2+1, m1+lmax, m3+lmax+1]``.
    B : ndarray, shape (lmax+1, 2*lmax+1, 3)
        Left Dyson overlap integrals indexed as ``B[N, lam+lmax, m2+1]``.
    Bd : ndarray, shape (lmax+1, 2*lmax+1, 3)
        Right Dyson overlap integrals indexed as ``Bd[Nd, lamd+lmax, m2d+1]``.
    eigv : ndarray, shape (lmax+1, lmax+1, 2*lmax+1)
        Eigenvector array indexed as ``eigv[N, l, lam+lmax]``.
    eign : ndarray, shape (lmax+1, 2*lmax+1)
        Eigenvalue array (ν_N) indexed as ``eign[N, lam+lmax]``.
    SPV : ndarray, shape (lmax+1, 2*lmax+3, 2)
        Parallel/perpendicular harmonic values indexed as
        ``SPV[l, m+lmax+1, pol]`` where pol=0 is parallel, pol=1 is perpendicular.
    lmax : int
        Maximum angular momentum.
    kr : float
        Unused (present for signature parity with MATLAB).
    r : float
        Unused (present for signature parity with MATLAB).

    Returns
    -------
    float
        The asymmetry parameter β, rounded to 4 decimal places.

    Notes
    -----
    Translates MATLAB beta2.m which performs a nested summation:

    .. math::
        \sigma_\parallel &= \sum_{N,N',\lambda,\lambda',m_2,m_2',l,l',m,L}
            (\text{eigvectors}) \times (\text{integrals}) \times
            e^{i\pi(\nu_N - \nu_{N'})/2} \times (\text{CG products}) \times (\text{SPV}_\parallel) \\
        \sigma_\perp &= \text{(same with perpendicular SPV)} \\
        \beta &= \frac{2(\sigma_\parallel - \sigma_\perp)}{\sigma_\parallel + 2\sigma_\perp}

    Selection rules enforced:
    - λ + m₂ == λ' + m₂'
    - m == m'
    - L == L'
    """
    sumpar = 0.0
    sumper = 0.0

    for N in range(lmax + 1):
        for Nd in range(lmax + 1):
            for lam in range(-N, N + 1):
                for lamd in range(-Nd, Nd + 1):
                    for m2 in range(-1, 2):  # -1, 0, 1
                        for m2d in range(-1, 2):
                            # Integral product
                            In = (
                                B[N, lam + lmax, m2 + 1]
                                * Bd[Nd, lamd + lmax, m2d + 1]
                            )
                            if In == 0:
                                continue
                            if (lam + m2) != (lamd + m2d):
                                continue

                            for l in range(abs(lam), lmax + 1):
                                for ld in range(abs(lamd), lmax + 1):
                                    vv = (
                                        eigv[N, l, lam + lmax]
                                        * eigv[Nd, ld, lamd + lmax]
                                    )
                                    if vv == 0:
                                        continue

                                    for m in range(-l, l + 1):
                                        md = m  # selection rule m == md
                                        if abs(md) > ld:
                                            continue

                                        for L in range(abs(l - 1), l + 2):
                                            Ld = L  # selection rule L == Ld
                                            if Ld < abs(ld - 1) or Ld > ld + 1:
                                                continue

                                            # CG coefficients
                                            # Index: CG[l,L][m2+1, m1+lmax, m3+lmax+1]
                                            # MATLAB: CG(l+1,L+1,m2+2,lam+lmax+1,(m2+lam)+lmax+2)
                                            key = (l, L)
                                            if key not in CG:
                                                continue
                                            cg_arr = CG[key]

                                            c1 = cg_arr[m2 + 1, lam + lmax, (m2 + lam) + lmax + 1]
                                            c2 = cg_arr[1, m + lmax, m + lmax + 1]  # m2=0
                                            key_d = (ld, Ld)
                                            if key_d not in CG:
                                                continue
                                            cg_arr_d = CG[key_d]
                                            c3 = cg_arr_d[m2d + 1, lamd + lmax, (m2d + lamd) + lmax + 1]
                                            c4 = cg_arr_d[1, md + lmax, md + lmax + 1]

                                            cleb = c1 * c2 * c3 * c4
                                            if cleb == 0:
                                                continue

                                            # Phase from eigenvalues
                                            phase = np.real(
                                                np.exp(
                                                    0.5j
                                                    * np.pi
                                                    * (
                                                        eign[N, lam + lmax]
                                                        - eign[Nd, lamd + lmax]
                                                    )
                                                )
                                            )

                                            # Prefactors
                                            prefac = (
                                                8.0
                                                / 3.0
                                                * 8.0
                                                * np.pi**2
                                                / (2 * L + 1)
                                            )
                                            signm = 2 * ((-1) ** m)

                                            # SPV indexed as SPV[l, m+lmax+1, pol]
                                            sp_l = SPV[l, -m + lmax + 1, 0]
                                            sp_ld = SPV[ld, md + lmax + 1, 0]
                                            sp_l_per = SPV[l, -m + lmax + 1, 1]
                                            sp_ld_per = SPV[ld, md + lmax + 1, 1]

                                            common = vv * In * phase * prefac * cleb * signm
                                            sumpar += common * sp_l * sp_ld
                                            sumper += common * sp_l_per * sp_ld_per

    denom = sumpar + 2 * sumper
    if denom == 0:
        return 0.0
    beta_val = 2 * (sumpar - sumper) / denom
    return round(beta_val, 4)
