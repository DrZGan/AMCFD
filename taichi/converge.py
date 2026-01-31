"""
AM-CFD Taichi Implementation - Convergence Acceleration

Python port of Fortran enhance_converge_speed (mod_converge.f90).

Performs an x-direction line correction (TDMA-like sweep) to accelerate
convergence of the enthalpy equation. Aggregates coefficients over j,k
slices to form a 1D tri-diagonal system in i, computes the correction
u0394h(i), and applies it in-place to `state.enthalpy`.

Inputs are Taichi data structures; computation is done via NumPy views
for simplicity, then written back to Taichi fields.
"""

from __future__ import annotations
from typing import Optional
import numpy as np

# Taichi structures
from data_structures import State, DiscretCoeffs, GridParams

def enhance_converge_speed(state: State, coeffs: DiscretCoeffs, grid: Optional[GridParams] = None) -> State:
    """Apply slice-wise enthalpy correction (x-direction TDMA).

    Mirrors logic from Fortran `enhance_converge_speed`:
    - Build aggregated arrays `bl`, `blp`, `blm`, `blc` across j,k for each i.
    - Forward sweep to compute `pib`, `qib`.
    - Backward substitution to get `delh(i)`.
    - Update interior enthalpy: `enthalpy(i,j,k) += delh(i)`.

    Args:
        state: Flow `State` containing `enthalpy` field.
        coeffs: Discretization coefficients (`ap, ae, aw, an, as_, at, ab, su`).
        grid: Optional `GridParams` for dimension hints; not required.

    Returns:
        The same `state` object with `enthalpy` modified in place.
    """
    # Dimensions (prefer explicit grid if provided, otherwise state)
    ni = grid.ni if grid is not None else state.ni
    nj = grid.nj if grid is not None else state.nj
    nk = grid.nk if grid is not None else state.nk

    # Interior index bounds (match Fortran: 2..nim1 etc.) using 0-based Python indices
    ist, ien = 1, ni - 2
    jst, jen = 1, nj - 2
    kst, ken = 1, nk - 2

    if ien < ist or jen < jst or ken < kst:
        # Too small grid to apply correction
        return state

    # Pull NumPy views
    ap = coeffs.ap.to_numpy()
    ae = coeffs.ae.to_numpy()
    aw = coeffs.aw.to_numpy()
    an = coeffs.an.to_numpy()
    as_ = coeffs.as_.to_numpy()
    at = coeffs.at.to_numpy()
    ab = coeffs.ab.to_numpy()
    su = coeffs.su.to_numpy()
    h = state.enthalpy.to_numpy()

    # Aggregated line arrays (full ni length to keep indices consistent)
    bl = np.zeros(ni, dtype=np.float64)
    blp = np.zeros(ni, dtype=np.float64)
    blm = np.zeros(ni, dtype=np.float64)
    blc = np.zeros(ni, dtype=np.float64)

    # Accumulate across interior j,k for each i
    # Vectorized sum over j,k for speed where possible
    # Note: we compute per-i terms; slicing ensures bounds correctness.
    for i in range(ist, ien + 1):
        # Slices
        js = slice(jst, jen + 1)
        ks = slice(kst, ken + 1)

        # Extract 2D interior planes at i
        ap_i = ap[i, js, ks]
        ae_i = ae[i, js, ks]
        aw_i = aw[i, js, ks]
        an_i = an[i, js, ks]
        as_i = as_[i, js, ks]
        at_i = at[i, js, ks]
        ab_i = ab[i, js, ks]
        su_i = su[i, js, ks]

        # Neighbor enthalpies
        h_e = h[i + 1, js, ks]
        h_w = h[i - 1, js, ks]
        h_n = h[i, jst + 1:jen + 2, ks]  # temp slice; will fix below
        h_s = h[i, jst - 1:jen, ks]      # temp slice; will fix below
        h_t = h[i, js, kst + 1:ken + 2]
        h_b = h[i, js, kst - 1:ken]

        # Correct neighbor slices for exact alignment
        # `an` multiplies h(i,j+1,k), `as` multiplies h(i,j-1,k)
        h_n = h[i, slice(jst + 1, jen + 2), ks]
        h_s = h[i, slice(jst - 1, jen), ks]
        # `at` multiplies h(i,j,k+1), `ab` multiplies h(i,j,k-1)
        h_t = h[i, js, slice(kst + 1, ken + 2)]
        h_b = h[i, js, slice(kst - 1, ken)]

        # Sums
        bl[i] = np.sum(ap_i - an_i - as_i - at_i - ab_i)
        blp[i] = np.sum(ae_i)
        blm[i] = np.sum(aw_i)
        blc[i] = (
            np.sum(ae_i * h_e) +
            np.sum(aw_i * h_w) +
            np.sum(an_i * h_n) +
            np.sum(as_i * h_s) +
            np.sum(at_i * h_t) +
            np.sum(ab_i * h_b) +
            np.sum(su_i) -
            np.sum(ap_i * h[i, js, ks])
        )

    # TDMA forward/backward sweep
    eps = 1e-20
    pib = np.zeros(ni, dtype=np.float64)
    qib = np.zeros(ni, dtype=np.float64)
    delh = np.zeros(ni, dtype=np.float64)

    # Seed at i=2 (Fortran indexing) -> Python index 2
    i_seed = 2
    if bl[i_seed] != 0.0:
        pib[i_seed] = blp[i_seed] / (bl[i_seed] + eps)
        qib[i_seed] = blc[i_seed] / (bl[i_seed] + eps)

    # Forward sweep: i=3..nim1 -> Python 3..(ni-2)
    for i in range(3, ni - 1):
        denom = bl[i] - blm[i] * pib[i - 1]
        denom = denom if abs(denom) > eps else eps
        pib[i] = blp[i] / denom
        qib[i] = (blc[i] + blm[i] * qib[i - 1]) / denom

    # Back substitution
    delh[ni - 2] = qib[ni - 2]
    for i in range(ni - 3, 1, -1):
        delh[i] = pib[i] * delh[i + 1] + qib[i]

    # Apply correction to interior cells
    for i in range(ist, ien + 1):
        h[i, jst:jen + 1, kst:ken + 1] += delh[i]

    # Write back to Taichi field
    state.enthalpy.from_numpy(h)
    return state


__all__ = ["enhance_converge_speed"]
