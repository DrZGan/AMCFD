"""
Unit tests for solve: solution_enthalpy, solution_uvw, and clean_uvw.

Uses small uniform grids where TDMA results are analytically reproducible.
"""

import taichi as ti
import numpy as np

from data_structures import State, DiscretCoeffs, PhysicsParams
from solve import solution_enthalpy, solution_uvw, clean_uvw


def _build_uniform_coeffs(ni, nj, nk, a=2.0, b=0.5, c=0.25, s=1.0):
    coeffs = DiscretCoeffs(ni, nj, nk)
    ap = np.full((ni, nj, nk), a, dtype=np.float64)
    ae = np.full((ni, nj, nk), b, dtype=np.float64)
    aw = np.full((ni, nj, nk), c, dtype=np.float64)
    an = np.zeros((ni, nj, nk), dtype=np.float64)
    as_ = np.zeros((ni, nj, nk), dtype=np.float64)
    at = np.zeros((ni, nj, nk), dtype=np.float64)
    ab = np.zeros((ni, nj, nk), dtype=np.float64)
    su = np.full((ni, nj, nk), s, dtype=np.float64)
    coeffs.ap.from_numpy(ap)
    coeffs.ae.from_numpy(ae)
    coeffs.aw.from_numpy(aw)
    coeffs.an.from_numpy(an)
    coeffs.as_.from_numpy(as_)
    coeffs.at.from_numpy(at)
    coeffs.ab.from_numpy(ab)
    coeffs.su.from_numpy(su)
    return coeffs


def _tdma_line_reference(ni, a, b, c, s):
    """1D TDMA along i for uniform coefficients with zero neighbors and zero initial field."""
    pr = np.zeros(ni, dtype=np.float64)
    qr = np.zeros(ni, dtype=np.float64)
    phi = np.zeros(ni, dtype=np.float64)
    pr[0] = 0.0
    qr[0] = 0.0
    for i in range(1, ni - 1):
        d = s
        denom = a - c * pr[i - 1]
        if denom <= 1e-12 and denom >= 0.0:
            denom = denom + 1e-13
        if denom >= -1e-12 and denom < 0.0:
            denom = denom - 1e-13
        pr[i] = b / denom
        qr[i] = (d + c * qr[i - 1]) / denom
    for i in range(ni - 2, 0, -1):
        phi[i] = pr[i] * phi[i + 1] + qr[i]
    return phi


def test_solution_enthalpy_uniform():
    ti.init(arch=ti.cpu)
    ni, nj, nk = 6, 5, 4
    a, b, c, s = 2.0, 0.5, 0.25, 1.0
    state = State(ni, nj, nk)
    coeffs = _build_uniform_coeffs(ni, nj, nk, a=a, b=b, c=c, s=s)
    state.enthalpy.from_numpy(np.zeros((ni, nj, nk)))
    solution_enthalpy(state, coeffs)
    # Expected: same line solution applied to all interior j,k
    line = _tdma_line_reference(ni, a, b, c, s)
    expected = np.zeros((ni, nj, nk), dtype=np.float64)
    for k in range(1, nk - 1):
        for j in range(1, nj - 1):
            for i in range(1, ni - 1):
                expected[i, j, k] = line[i]
    h_result = state.enthalpy.to_numpy()
    assert np.allclose(h_result, expected, rtol=1e-12, atol=1e-12)


def test_solution_uvw_uniform():
    ti.init(arch=ti.cpu)
    ni, nj, nk = 6, 5, 4
    a, b, c, s = 2.0, 0.5, 0.25, 1.0
    state = State(ni, nj, nk)
    coeffs = _build_uniform_coeffs(ni, nj, nk, a=a, b=b, c=c, s=s)
    # Zero velocities initially
    zeros = np.zeros((ni, nj, nk), dtype=np.float64)
    state.uVel.from_numpy(zeros)
    state.vVel.from_numpy(zeros)
    state.wVel.from_numpy(zeros)
    # Solve for u, then v
    solution_uvw(1, state, coeffs)
    solution_uvw(2, state, coeffs)
    line = _tdma_line_reference(ni, a, b, c, s)
    expected = np.zeros((ni, nj, nk), dtype=np.float64)
    for k in range(1, nk - 1):
        for j in range(1, nj - 1):
            for i in range(1, ni - 1):
                expected[i, j, k] = line[i]
    u_res = state.uVel.to_numpy()
    v_res = state.vVel.to_numpy()
    w_res = state.wVel.to_numpy()
    assert np.allclose(u_res, expected, rtol=1e-12, atol=1e-12)
    assert np.allclose(v_res, expected, rtol=1e-12, atol=1e-12)
    # w remains zero (never solved for ivar=3 here)
    assert np.allclose(w_res, zeros)


def test_clean_uvw_thresholds():
    ti.init(arch=ti.cpu)
    ni, nj, nk = 4, 4, 3
    state = State(ni, nj, nk)
    physics = PhysicsParams(tsolid=1200.0, tvapor=3000.0)
    # Initialize velocities to ones
    ones = np.ones((ni, nj, nk), dtype=np.float64)
    state.uVel.from_numpy(ones.copy())
    state.vVel.from_numpy(ones.copy())
    state.wVel.from_numpy(ones.copy())
    # Temperature field: default above solidus, below vapor
    temp = np.full((ni, nj, nk), 1600.0, dtype=np.float64)
    # Solid threshold triggers
    temp[1, 1, 1] = 1000.0
    temp[2, 1, 1] = 900.0   # min(i,i+1) -> u at [2,1,1] zero
    temp[1, 2, 1] = 900.0   # min(j,j+1) -> v at [1,2,1] zero
    temp[1, 1, 2] = 900.0   # min(k,k+1) -> w at [1,1,2] zero
    # Vapor threshold: zero entire column at (0,0,:)
    temp[0, 0, nk - 1] = 3500.0
    state.temp.from_numpy(temp)
    clean_uvw(state, physics)
    u = state.uVel.to_numpy()
    v = state.vVel.to_numpy()
    w = state.wVel.to_numpy()
    # Solid checks
    assert u[2, 1, 1] == 0.0
    assert v[1, 2, 1] == 0.0
    assert w[1, 1, 2] == 0.0
    # Vapor checks: entire column (0,0,:) zero
    assert u[0, 0, 0] == 0.0 and u[0, 0, 1] == 0.0
    assert v[0, 0, 0] == 0.0 and v[0, 0, 1] == 0.0
    assert w[0, 0, 0] == 0.0 and w[0, 0, 1] == 0.0


if __name__ == "__main__":
    test_solution_enthalpy_uniform()
    test_solution_uvw_uniform()
    test_clean_uvw_thresholds()
    print("All solve tests: PASSED")
