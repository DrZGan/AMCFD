"""
AM-CFD Python Implementation - Solvers (TDMA line-by-line, NumPy)

Converted from Fortran mod_solve.f90 without Taichi kernels.

Functions:
- solution_uvw: momentum equations line TDMA along x
- solution_enthalpy: enthalpy equation line TDMA along x
- clean_uvw: zero velocities in solid or boiling regions

Inputs:
- state: State object (fields updated in place)
- coeffs: DiscretCoeffs (ap/ae/aw/an/as_/at/ab/su)
- physics: PhysicsParams (tsolid/tvapor from input_param.yaml) — only for clean_uvw

Outputs:
- Returns the same State object with modified fields.
"""

from __future__ import annotations
import numpy as np

from data_structures import State, DiscretCoeffs, PhysicsParams


def _tdma_line_update_scalar(field: np.ndarray, ap: np.ndarray, ae: np.ndarray, aw: np.ndarray,
                             an: np.ndarray, as_: np.ndarray, at: np.ndarray, ab: np.ndarray,
                             su: np.ndarray, ni: int, nj: int, nk: int) -> None:
    """Line-by-line TDMA update along i for a scalar field (enthalpy)."""
    for _ in range(2):
        for k in range(nk - 2, 0, -1):
            for _ in range(2):
                for j in range(1, nj - 1):
                    pr = np.zeros(ni)
                    qr = np.zeros(ni)
                    pr[0] = 0.0
                    qr[0] = field[0, j, k]
                    for i in range(1, ni - 1):
                        d = (at[i, j, k] * field[i, j, k + 1] +
                             ab[i, j, k] * field[i, j, k - 1] +
                             an[i, j, k] * field[i, j + 1, k] +
                             as_[i, j, k] * field[i, j - 1, k] +
                             su[i, j, k])
                        denom = ap[i, j, k] - aw[i, j, k] * pr[i - 1]
                        if denom <= 1e-12 and denom >= 0.0:
                            denom = denom + 1e-13
                        if denom >= -1e-12 and denom < 0.0:
                            denom = denom - 1e-13
                        pr[i] = ae[i, j, k] / denom
                        qr[i] = (d + aw[i, j, k] * qr[i - 1]) / denom
                    for i in range(ni - 2, 0, -1):
                        field[i, j, k] = pr[i] * field[i + 1, j, k] + qr[i]


def _tdma_line_update_vector(field: np.ndarray, ap: np.ndarray, ae: np.ndarray, aw: np.ndarray,
                             an: np.ndarray, as_: np.ndarray, at: np.ndarray, ab: np.ndarray,
                             su: np.ndarray, ni: int, nj: int, nk: int) -> None:
    """Line-by-line TDMA update along i for a vector component (u/v/w)."""
    for _ in range(2):
        for k in range(nk - 2, 0, -1):
            for _ in range(2):
                for j in range(1, nj - 1):
                    pr = np.zeros(ni)
                    qr = np.zeros(ni)
                    pr[0] = 0.0
                    qr[0] = field[0, j, k]
                    for i in range(1, ni - 1):
                        d = (at[i, j, k] * field[i, j, k + 1] +
                             ab[i, j, k] * field[i, j, k - 1] +
                             an[i, j, k] * field[i, j + 1, k] +
                             as_[i, j, k] * field[i, j - 1, k] +
                             su[i, j, k])
                        denom = ap[i, j, k] - aw[i, j, k] * pr[i - 1]
                        if denom <= 1e-12 and denom >= 0.0:
                            denom = denom + 1e-13
                        if denom >= -1e-12 and denom < 0.0:
                            denom = denom - 1e-13
                        pr[i] = ae[i, j, k] / denom
                        qr[i] = (d + aw[i, j, k] * qr[i - 1]) / denom
                    for i in range(ni - 2, 0, -1):
                        field[i, j, k] = pr[i] * field[i + 1, j, k] + qr[i]


def _tdma_scalar_with_coeffs(field_ti, coeffs: DiscretCoeffs, ni: int, nj: int, nk: int) -> None:
    """Wrapper: accept `coeffs` object, unpack arrays once, and update scalar field."""
    field = field_ti.to_numpy()
    ap = coeffs.ap.to_numpy()
    ae = coeffs.ae.to_numpy()
    aw = coeffs.aw.to_numpy()
    an = coeffs.an.to_numpy()
    as_ = coeffs.as_.to_numpy()
    at = coeffs.at.to_numpy()
    ab = coeffs.ab.to_numpy()
    su = coeffs.su.to_numpy()
    _tdma_line_update_scalar(field, ap, ae, aw, an, as_, at, ab, su, ni, nj, nk)
    field_ti.from_numpy(field)


def _tdma_vector_with_coeffs(field_ti, coeffs: DiscretCoeffs, ni: int, nj: int, nk: int) -> None:
    """Wrapper: accept `coeffs` object, unpack arrays once, and update vector component."""
    field = field_ti.to_numpy()
    ap = coeffs.ap.to_numpy()
    ae = coeffs.ae.to_numpy()
    aw = coeffs.aw.to_numpy()
    an = coeffs.an.to_numpy()
    as_ = coeffs.as_.to_numpy()
    at = coeffs.at.to_numpy()
    ab = coeffs.ab.to_numpy()
    su = coeffs.su.to_numpy()
    _tdma_line_update_vector(field, ap, ae, aw, an, as_, at, ab, su, ni, nj, nk)
    field_ti.from_numpy(field)


def solution_enthalpy(state: State, coeffs: DiscretCoeffs) -> State:
    """Solve enthalpy equation using TDMA wrapper that accepts `coeffs`."""
    _tdma_scalar_with_coeffs(state.enthalpy, coeffs, state.ni, state.nj, state.nk)
    return state


def solution_uvw(ivar: int, state: State, coeffs: DiscretCoeffs) -> State:
    """Solve momentum equations for ivar in {1,2,3} using TDMA wrapper that accepts `coeffs`."""
    ni, nj, nk = state.ni, state.nj, state.nk
    if ivar == 1:
        _tdma_vector_with_coeffs(state.uVel, coeffs, ni, nj, nk)
    elif ivar == 2:
        _tdma_vector_with_coeffs(state.vVel, coeffs, ni, nj, nk)
    elif ivar == 3:
        _tdma_vector_with_coeffs(state.wVel, coeffs, ni, nj, nk)
    return state


def clean_uvw(state: State, physics: PhysicsParams) -> State:
    """Zero velocity components in solid or boiling regions using thresholds from physics."""
    ni, nj, nk = state.ni, state.nj, state.nk
    temp = state.temp.to_numpy()
    u = state.uVel.to_numpy()
    v = state.vVel.to_numpy()
    w = state.wVel.to_numpy()
    for k in range(0, nk - 1):
        for j in range(0, nj - 1):
            for i in range(0, ni - 1):
                tulc = min(temp[i, j, k], temp[i + 1, j, k])
                tvlc = min(temp[i, j, k], temp[i, j + 1, k])
                twlc = min(temp[i, j, k], temp[i, j, k + 1])
                if tulc <= physics.tsolid:
                    u[i + 1, j, k] = 0.0
                if tvlc <= physics.tsolid:
                    v[i, j + 1, k] = 0.0
                if twlc <= physics.tsolid:
                    w[i, j, k + 1] = 0.0
                if temp[i, j, nk - 1] >= physics.tvapor:
                    u[i, j, k] = 0.0
                    v[i, j, k] = 0.0
                    w[i, j, k] = 0.0
    state.uVel.from_numpy(u)
    state.vVel.from_numpy(v)
    state.wVel.from_numpy(w)
    return state


__all__ = ["solution_enthalpy", "solution_uvw", "clean_uvw"]
