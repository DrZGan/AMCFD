"""
Test residuals: Verify residual.py matches Fortran mod_resid.f90 logic

Constructs synthetic fields and coefficients, computes expected residuals
per the Fortran equations, and compares with compute_residuals().
"""

import sys
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

import math
import numpy as np
import taichi as ti

ti.init(arch=ti.cpu)

from data_structures import GridParams, State, DiscretCoeffs, ConvergenceState, PoolDimensions
from residual import compute_residuals


def _make_uniform_grid(ni, nj, nk):
    gp = GridParams(ni, nj, nk)
    # Fill faces uniformly (not used for residuals directly)
    gp.xu.from_numpy(np.linspace(0.0, float(ni - 1), ni))
    gp.yv.from_numpy(np.linspace(0.0, float(nj - 1), nj))
    gp.zw.from_numpy(np.linspace(0.0, float(nk - 1), nk))
    # Centers as midpoints
    x = np.zeros(ni); x[:-1] = 0.5 * (gp.xu.to_numpy()[1:] + gp.xu.to_numpy()[:-1]); x[-1] = gp.xu.to_numpy()[-1]
    y = np.zeros(nj); y[:-1] = 0.5 * (gp.yv.to_numpy()[1:] + gp.yv.to_numpy()[:-1]); y[-1] = gp.yv.to_numpy()[-1]
    z = np.zeros(nk); z[:-1] = 0.5 * (gp.zw.to_numpy()[1:] + gp.zw.to_numpy()[:-1]); z[-1] = gp.zw.to_numpy()[-1]
    gp.x.from_numpy(x)
    gp.y.from_numpy(y)
    gp.z.from_numpy(z)
    # Areas for mass residual calculation
    gp.areaij.from_numpy(np.ones((ni, nj)))
    gp.areaik.from_numpy(np.ones((ni, nk)))
    gp.areajk.from_numpy(np.ones((nj, nk)))
    return gp


def _fill_state(ni, nj, nk):
    st = State(ni, nj, nk)
    # Deterministic pattern for velocities/enthalpy
    u = np.zeros((ni, nj, nk)); v = np.zeros_like(u); w = np.zeros_like(u); H = np.zeros_like(u)
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                u[i, j, k] = i + j + k
                v[i, j, k] = 2 * i - j + k
                w[i, j, k] = -i + 3 * j - 2 * k
                H[i, j, k] = 0.5 * (i + j + k)
    st.uVel.from_numpy(u)
    st.vVel.from_numpy(v)
    st.wVel.from_numpy(w)
    st.enthalpy.from_numpy(H)
    return st


def _fill_coeffs(ni, nj, nk):
    coeff = DiscretCoeffs(ni, nj, nk)
    # Use simple coefficients to avoid degeneracy
    an = np.ones((ni, nj, nk))
    as_ = np.ones((ni, nj, nk))
    ae = np.ones((ni, nj, nk))
    aw = np.ones((ni, nj, nk))
    at = np.ones((ni, nj, nk))
    ab = np.ones((ni, nj, nk))
    ap = np.full((ni, nj, nk), 6.0)
    su = np.zeros((ni, nj, nk))

    coeff.an.from_numpy(an)
    coeff.as_.from_numpy(as_)
    coeff.ae.from_numpy(ae)
    coeff.aw.from_numpy(aw)
    coeff.at.from_numpy(at)
    coeff.ab.from_numpy(ab)
    coeff.ap.from_numpy(ap)
    coeff.su.from_numpy(su)
    return coeff


def _expected_refmom(denl, pi, pool, umaxt):
    diameter = min(pool.width, pool.length, pool.depth)
    return 0.25 * pi * (diameter ** 2) * denl * (umaxt ** 2)


def _expected_sumd_velocity(st_arr, coeff, i_rng, j_rng, k_rng):
    an, as_, ae, aw, at, ab, ap, su = [
        coeff.an.to_numpy(), coeff.as_.to_numpy(), coeff.ae.to_numpy(), coeff.aw.to_numpy(),
        coeff.at.to_numpy(), coeff.ab.to_numpy(), coeff.ap.to_numpy(), coeff.su.to_numpy()
    ]
    sumd = 0.0
    for k in k_rng:
        for j in j_rng:
            for i in i_rng:
                resor = (
                    an[i, j, k] * st_arr[i, j + 1, k]
                    + as_[i, j, k] * st_arr[i, j - 1, k]
                    + ae[i, j, k] * st_arr[i + 1, j, k]
                    + aw[i, j, k] * st_arr[i - 1, j, k]
                    + at[i, j, k] * st_arr[i, j, k + 1]
                    + ab[i, j, k] * st_arr[i, j, k - 1]
                    + su[i, j, k]
                    - ap[i, j, k] * st_arr[i, j, k]
                )
                sumd += abs(resor)
    return sumd


def _expected_sumd_h(H, coeff, i_rng, j_rng, k_rng, small):
    an, as_, ae, aw, at, ab, ap, su = [
        coeff.an.to_numpy(), coeff.as_.to_numpy(), coeff.ae.to_numpy(), coeff.aw.to_numpy(),
        coeff.at.to_numpy(), coeff.ab.to_numpy(), coeff.ap.to_numpy(), coeff.su.to_numpy()
    ]
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
    return sumd / (sumh + small)


def _expected_mass_norm(gp, st, denl, small):
    u = st.uVel.to_numpy(); v = st.vVel.to_numpy(); w = st.wVel.to_numpy()
    areaij = gp.areaij.to_numpy(); areaik = gp.areaik.to_numpy(); areajk = gp.areajk.to_numpy()
    ni, nj, nk = u.shape
    i_rng = range(1, ni - 1); j_rng = range(1, nj - 1); k_rng = range(1, nk - 1)
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
    return denom + small


def test_velocity_residuals_match_fortran():
    ni, nj, nk = 6, 6, 5
    gp = _make_uniform_grid(ni, nj, nk)
    st = _fill_state(ni, nj, nk)
    coeff = _fill_coeffs(ni, nj, nk)
    pool = PoolDimensions(length=1.0, width=1.0, depth=1.0)
    conv = ConvergenceState()

    denl = 1.0
    pi = math.pi
    small = 1.0e-6

    # Interior ranges equivalent to Fortran defaults
    i_rng = range(1, ni - 1); j_rng = range(1, nj - 1); k_rng = range(1, nk - 1)

    # u-residual
    sumd_u = _expected_sumd_velocity(st.uVel.to_numpy(), coeff, i_rng, j_rng, k_rng)
    umaxt_u = np.max(np.abs(st.uVel.to_numpy()[1:ni-1, 1:nj-1, nk - 1]))
    refmom_u = _expected_refmom(denl, pi, pool, umaxt_u)
    expected_u = sumd_u / refmom_u
    conv = compute_residuals(1, gp, st, coeff, denl, pi, small, pool, conv)
    assert abs(conv.residual_u - expected_u) / (expected_u + small) < 1e-12

    # v-residual
    sumd_v = _expected_sumd_velocity(st.vVel.to_numpy(), coeff, i_rng, j_rng, k_rng)
    umaxt_v = np.max(np.abs(st.vVel.to_numpy()[1:ni-1, 1:nj-1, nk - 1]))
    refmom_v = _expected_refmom(denl, pi, pool, umaxt_v)
    expected_v = sumd_v / refmom_v
    conv = compute_residuals(2, gp, st, coeff, denl, pi, small, pool, conv)
    assert abs(conv.residual_v - expected_v) / (expected_v + small) < 1e-12

    # w-residual
    sumd_w = _expected_sumd_velocity(st.wVel.to_numpy(), coeff, i_rng, j_rng, k_rng)
    umaxt_w = np.max(np.abs(st.wVel.to_numpy()[1:ni-1, 1:nj-1, nk - 1]))
    refmom_w = _expected_refmom(denl, pi, pool, umaxt_w)
    expected_w = sumd_w / refmom_w
    conv = compute_residuals(3, gp, st, coeff, denl, pi, small, pool, conv)
    assert abs(conv.residual_w - expected_w) / (expected_w + small) < 1e-12

    print("✓ Velocity residuals match Fortran logic")


def test_energy_residual_matches_fortran():
    ni, nj, nk = 6, 6, 5
    gp = _make_uniform_grid(ni, nj, nk)
    st = _fill_state(ni, nj, nk)
    coeff = _fill_coeffs(ni, nj, nk)
    pool = PoolDimensions(length=1.0, width=1.0, depth=1.0)
    conv = ConvergenceState()

    small = 1.0e-6
    i_rng = range(1, ni - 1); j_rng = range(1, nj - 1); k_rng = range(1, nk - 1)

    expected_h = _expected_sumd_h(st.enthalpy.to_numpy(), coeff, i_rng, j_rng, k_rng, small)
    conv = compute_residuals(5, gp, st, coeff, denl=1.0, pi=math.pi, small=small, pool=pool, conv=conv)
    assert abs(conv.residual_h - expected_h) / (expected_h + small) < 1e-12

    print("✓ Energy residual matches Fortran logic")


def test_mass_residual_normalization_matches_fortran():
    ni, nj, nk = 6, 6, 5
    gp = _make_uniform_grid(ni, nj, nk)
    st = _fill_state(ni, nj, nk)
    coeff = _fill_coeffs(ni, nj, nk)
    pool = PoolDimensions(length=1.0, width=1.0, depth=1.0)
    conv = ConvergenceState(residual_p=2.5)  # emulate existing resorm prior to normalization

    small = 1.0e-6
    denl = 1.23

    denom_plus_small = _expected_mass_norm(gp, st, denl, small)
    conv = compute_residuals(4, gp, st, coeff, denl=denl, pi=math.pi, small=small, pool=pool, conv=conv)
    expected_p = 2.5 / denom_plus_small
    assert abs(conv.residual_p - expected_p) / (expected_p + small) < 1e-12

    print("✓ Mass residual normalization matches Fortran logic")


def run_all_tests():
    print("\n" + "=" * 60)
    print("Testing residuals vs Fortran mod_resid.f90")
    print("=" * 60 + "\n")

    test_velocity_residuals_match_fortran()
    test_energy_residual_matches_fortran()
    test_mass_residual_normalization_matches_fortran()

    print("\n" + "=" * 60)
    print("All residual tests passed! ✓")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    run_all_tests()
