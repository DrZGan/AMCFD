import numpy as np
import taichi as ti
from typing import Optional

from data_structures import GridParams, State, DiscretCoeffs, ConvergenceState, PoolDimensions

ti.init(arch=ti.cpu)

def _to_np(field):
    """Safe convert Taichi field to numpy array."""
    return field.to_numpy() if hasattr(field, "to_numpy") else np.asarray(field)

def compute_residuals(
    ivar: int,
    gp: GridParams,
    st: State,
    coeff: DiscretCoeffs,
    denl: float,
    pi: float,
    small: float,
    pool: PoolDimensions,
    conv: ConvergenceState,
) -> ConvergenceState:
    """
    Python/Taichi conversion of Fortran mod_resid.f90 `residual`.

    Updates residuals in `conv`:
    - residual_u, residual_v, residual_w, residual_h
    - residual_p (normalized mass residual)
    - max_residual (max of u,v,w,p)

    Assumptions:
    - Interior loops use indices [1..nim1-1], [1..njm1-1], [1..nkm1-1] to avoid boundaries.
    - `conv.residual_p` holds the pre-normalized mass residual prior to this call.
    """
    # Dimensions
    ni, nj, nk = gp.ni, gp.nj, gp.nk
    nim1, njm1, nkm1 = ni - 1, nj - 1, nk - 1

    # Fields (numpy views)
    u = _to_np(st.uVel)
    v = _to_np(st.vVel)
    w = _to_np(st.wVel)
    H = _to_np(st.enthalpy)

    an = _to_np(coeff.an)
    as_ = _to_np(coeff.as_)  # Fortran `as` -> Python `as_`
    ae = _to_np(coeff.ae)
    aw = _to_np(coeff.aw)
    at = _to_np(coeff.at)
    ab = _to_np(coeff.ab)
    ap = _to_np(coeff.ap)
    su = _to_np(coeff.su)

    areaij = _to_np(gp.areaij)
    areaik = _to_np(gp.areaik)
    areajk = _to_np(gp.areajk)

    # Helper: reference momentum (uses last k-plane velocity max)
    def _ref_momentum(vel_last_plane: np.ndarray) -> float:
        umaxt = np.max(np.abs(vel_last_plane))
        diameter = min(pool.width, pool.length, pool.depth)
        return 0.25 * pi * (diameter ** 2) * denl * (umaxt ** 2)

    # Interior index ranges (avoid boundaries for +/- neighbors)
    i_rng = range(1, nim1)
    j_rng = range(1, njm1)
    k_rng = range(1, nkm1)

    # Fortran computed different ranges for some cases; adopt interior consistently here.
    if ivar == 5:
        # Energy equation residual (case 500 in Fortran)
        sumd = 0.0
        sumh = 0.0
        for k in k_rng:
            for j in j_rng:
                for i in i_rng:
                    resor = (
                        an[i, j, k] * H[i, j + 1, k]
                        + as_[i, j, k] * H[i, j - 1, k]
                        + ae[i, j, k] * H[i + 1, j, k]
                        + aw[i, j, k] * H[i - 1, j, k]
                        + at[i, j, k] * H[i, j, k + 1]
                        + ab[i, j, k] * H[i, j, k - 1]
                        + su[i, j, k]
                    ) / ap[i, j, k] - H[i, j, k]
                    sumd += abs(resor)
                    sumh += abs(H[i, j, k])
        conv.residual_h = sumd / (sumh + small)
    elif ivar == 1:
        # u-velocity residual (case 100)
        sumd = 0.0
        for k in k_rng:
            for j in j_rng:
                for i in i_rng:
                    resor = (
                        an[i, j, k] * u[i, j + 1, k]
                        + as_[i, j, k] * u[i, j - 1, k]
                        + ae[i, j, k] * u[i + 1, j, k]
                        + aw[i, j, k] * u[i - 1, j, k]
                        + at[i, j, k] * u[i, j, k + 1]
                        + ab[i, j, k] * u[i, j, k - 1]
                        + su[i, j, k]
                        - ap[i, j, k] * u[i, j, k]
                    )
                    sumd += abs(resor)
        refmom = _ref_momentum(u[1:nim1, 1:njm1, nk - 1])
        conv.residual_u = sumd / (refmom if refmom > 0 else small)
    elif ivar == 2:
        # v-velocity residual (case 200)
        sumd = 0.0
        for k in k_rng:
            for j in j_rng:
                for i in i_rng:
                    resor = (
                        an[i, j, k] * v[i, j + 1, k]
                        + as_[i, j, k] * v[i, j - 1, k]
                        + ae[i, j, k] * v[i + 1, j, k]
                        + aw[i, j, k] * v[i - 1, j, k]
                        + at[i, j, k] * v[i, j, k + 1]
                        + ab[i, j, k] * v[i, j, k - 1]
                        + su[i, j, k]
                        - ap[i, j, k] * v[i, j, k]
                    )
                    sumd += abs(resor)
        refmom = _ref_momentum(v[1:nim1, 1:njm1, nk - 1])
        conv.residual_v = sumd / (refmom if refmom > 0 else small)
    elif ivar == 3:
        # w-velocity residual (case 300)
        sumd = 0.0
        for k in k_rng:
            for j in j_rng:
                for i in i_rng:
                    resor = (
                        an[i, j, k] * w[i, j + 1, k]
                        + as_[i, j, k] * w[i, j - 1, k]
                        + ae[i, j, k] * w[i + 1, j, k]
                        + aw[i, j, k] * w[i - 1, j, k]
                        + at[i, j, k] * w[i, j, k + 1]
                        + ab[i, j, k] * w[i, j, k - 1]
                        + su[i, j, k]
                        - ap[i, j, k] * w[i, j, k]
                    )
                    sumd += abs(resor)
        refmom = _ref_momentum(w[1:nim1, 1:njm1, nk - 1])
        conv.residual_w = sumd / (refmom if refmom > 0 else small)
    elif ivar == 4:
        # Normalized mass residual (case 400)
        denom = 0.0
        for k in k_rng:
            for j in j_rng:
                for i in i_rng:
                    dtpvar = (
                        (abs(u[i, j, k]) + abs(u[i + 1, j, k])) * areajk[j, k]
                        + (abs(v[i, j, k]) + abs(v[i, j + 1, k])) * areaik[i, k]
                        + (abs(w[i, j, k]) + abs(w[i, j, k + 1])) * areaij[i, j]
                    )
                    denom += 0.5 * abs(dtpvar)
        denom *= denl
        # Normalize existing residual_p with denom
        base = conv.residual_p
        conv.residual_p = base / (denom + small)
    else:
        # No-op for other ivar values
        pass

    # Update max residual (consistent with Fortran main logic)
    conv.max_residual = max(
        abs(conv.residual_u),
        abs(conv.residual_v),
        abs(conv.residual_w),
        abs(conv.residual_p),
    )

    return conv


__all__ = ["compute_residuals"]
