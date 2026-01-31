"""
AM-CFD Taichi Implementation - Source Terms

Converted from Fortran module: mod_sour.f90
"""

import taichi as ti

from data_structures import (
    State,
    StatePrev,
    MaterialProps,
    GridParams,
    DiscretCoeffs,
    SimulationParams,
    PhysicsParams,
    LaserParams,
    LaserState,
)


@ti.kernel
def _source_u(
    state: ti.template(),
    state_prev: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
    physics: ti.template(),
):
    i_start = ti.static(simu_params.istatp1)
    i_end = ti.static(min(simu_params.iendm1, simu_params.ni - 2))
    j_start = ti.static(simu_params.jstatp1)
    j_end = ti.static(min(simu_params.jendm1, simu_params.nj - 2))
    k_start = ti.static(simu_params.kstatp1)
    k_end = ti.static(min(simu_params.nkm1, simu_params.nk - 2))

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                fracx = grid.fracx[i - 1]
                fraclu = state.fracl[i, j, k] * (1.0 - fracx) + state.fracl[i - 1, j, k] * fracx
                if fraclu > 0.0:
                    term = 180.0 * physics.vis0 / (1.0e-5 ** 2) * (1.0 - fraclu) ** 2 / (fraclu + 1.0e-3)
                    coeffs.sp[i, j, k] = coeffs.sp[i, j, k] - term * grid.volume_u[i, j, k]

    # k = nkm1 (top interior)
    k_top_inner = ti.static(min(simu_params.nkm1 - 1, simu_params.nk - 2))
    k_top = ti.static(simu_params.nkm1)
    for j in range(j_start, j_end + 1):
        for i in range(i_start, i_end + 1):
            coeffs.su[i, j, k_top_inner] = coeffs.su[i, j, k_top_inner] + coeffs.at[i, j, k_top_inner] * state.uVel[i, j, k_top]
            coeffs.sp[i, j, k_top_inner] = coeffs.sp[i, j, k_top_inner] - coeffs.at[i, j, k_top_inner]
            coeffs.at[i, j, k_top_inner] = 0.0

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                coeffs.ap[i, j, k] = (
                    coeffs.an[i, j, k]
                    + coeffs.as_[i, j, k]
                    + coeffs.ae[i, j, k]
                    + coeffs.aw[i, j, k]
                    + coeffs.at[i, j, k]
                    + coeffs.ab[i, j, k]
                    + coeffs.apnot[i, j, k]
                    - coeffs.sp[i, j, k]
                )
                coeffs.dux[i, j, k] = grid.areajk[j, k] / coeffs.ap[i, j, k]

                # Under-relaxation
                coeffs.ap[i, j, k] = coeffs.ap[i, j, k] / simu_params.urf_vel
                coeffs.su[i, j, k] = coeffs.su[i, j, k] + (1.0 - simu_params.urf_vel) * coeffs.ap[i, j, k] * state.uVel[i, j, k]
                coeffs.dux[i, j, k] = coeffs.dux[i, j, k] * simu_params.urf_vel

                # Zero velocity in solid
                tulc = ti.min(state.temp[i, j, k], state.temp[i - 1, j, k])
                if tulc <= physics.tsolid:
                    coeffs.su[i, j, k] = 0.0
                    coeffs.an[i, j, k] = 0.0
                    coeffs.as_[i, j, k] = 0.0
                    coeffs.ae[i, j, k] = 0.0
                    coeffs.aw[i, j, k] = 0.0
                    coeffs.at[i, j, k] = 0.0
                    coeffs.ab[i, j, k] = 0.0
                    coeffs.ap[i, j, k] = 1.0e20


@ti.kernel
def _source_v(
    state: ti.template(),
    state_prev: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
    physics: ti.template(),
):
    i_start = ti.static(simu_params.istatp1)
    i_end = ti.static(min(simu_params.iendm1, simu_params.ni - 2))
    j_start = ti.static(simu_params.jstatp1)
    j_end = ti.static(min(simu_params.jendm1, simu_params.nj - 2))
    k_start = ti.static(simu_params.kstatp1)
    k_end = ti.static(min(simu_params.nkm1, simu_params.nk - 2))

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                fracy = grid.fracy[j - 1]
                fraclv = state.fracl[i, j, k] * (1.0 - fracy) + state.fracl[i, j - 1, k] * fracy
                if fraclv > 0.0:
                    term = 180.0 * physics.vis0 / (1.0e-5 ** 2) * (1.0 - fraclv) ** 2 / (fraclv + 1.0e-3)
                    coeffs.sp[i, j, k] = coeffs.sp[i, j, k] - term * grid.volume_v[i, j, k]

    # k = nkm1 (top interior)
    k_top_inner = ti.static(min(simu_params.nkm1 - 1, simu_params.nk - 2))
    k_top = ti.static(simu_params.nkm1)
    for j in range(j_start, j_end + 1):
        for i in range(i_start, i_end + 1):
            coeffs.su[i, j, k_top_inner] = coeffs.su[i, j, k_top_inner] + coeffs.at[i, j, k_top_inner] * state.vVel[i, j, k_top]
            coeffs.sp[i, j, k_top_inner] = coeffs.sp[i, j, k_top_inner] - coeffs.at[i, j, k_top_inner]
            coeffs.at[i, j, k_top_inner] = 0.0

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                coeffs.ap[i, j, k] = (
                    coeffs.an[i, j, k]
                    + coeffs.as_[i, j, k]
                    + coeffs.ae[i, j, k]
                    + coeffs.aw[i, j, k]
                    + coeffs.at[i, j, k]
                    + coeffs.ab[i, j, k]
                    + coeffs.apnot[i, j, k]
                    - coeffs.sp[i, j, k]
                )
                coeffs.dvy[i, j, k] = grid.areaik[i, k] / coeffs.ap[i, j, k]

                # Under-relaxation
                coeffs.ap[i, j, k] = coeffs.ap[i, j, k] / simu_params.urf_vel
                coeffs.su[i, j, k] = coeffs.su[i, j, k] + (1.0 - simu_params.urf_vel) * coeffs.ap[i, j, k] * state.vVel[i, j, k]
                coeffs.dvy[i, j, k] = coeffs.dvy[i, j, k] * simu_params.urf_vel

                # Zero velocity in solid
                tvlc = ti.min(state.temp[i, j, k], state.temp[i, j - 1, k])
                if tvlc <= physics.tsolid:
                    coeffs.su[i, j, k] = 0.0
                    coeffs.an[i, j, k] = 0.0
                    coeffs.as_[i, j, k] = 0.0
                    coeffs.ae[i, j, k] = 0.0
                    coeffs.aw[i, j, k] = 0.0
                    coeffs.at[i, j, k] = 0.0
                    coeffs.ab[i, j, k] = 0.0
                    coeffs.ap[i, j, k] = 1.0e20


@ti.kernel
def _source_w(
    state: ti.template(),
    state_prev: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
    physics: ti.template(),
):
    i_start = ti.static(simu_params.istatp1)
    i_end = ti.static(min(simu_params.iendm1, simu_params.ni - 2))
    j_start = ti.static(simu_params.jstatp1)
    j_end = ti.static(min(simu_params.jendm1, simu_params.nj - 2))
    k_start = ti.static(simu_params.kstatp1)
    k_end = ti.static(min(simu_params.nkm1, simu_params.nk - 2))

    boufac = physics.rholiq * physics.grav * physics.beta

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                fracz = grid.fracz[k - 1]
                fraclw = state.fracl[i, j, k] * (1.0 - fracz) + state.fracl[i, j, k - 1] * fracz
                if fraclw > 0.0:
                    term = 180.0 * physics.vis0 / (1.0e-5 ** 2) * (1.0 - fraclw) ** 2 / (fraclw + 1.0e-3)
                    coeffs.sp[i, j, k] = coeffs.sp[i, j, k] - term * grid.volume_w[i, j, k]

                    # Buoyancy
                    tw = state.temp[i, j, k] * (1.0 - fracz) + state.temp[i, j, k - 1] * fracz
                    coeffs.su[i, j, k] = coeffs.su[i, j, k] + boufac * grid.volume_w[i, j, k] * (tw - physics.tsolid)

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                coeffs.ap[i, j, k] = (
                    coeffs.an[i, j, k]
                    + coeffs.as_[i, j, k]
                    + coeffs.ae[i, j, k]
                    + coeffs.aw[i, j, k]
                    + coeffs.at[i, j, k]
                    + coeffs.ab[i, j, k]
                    + coeffs.apnot[i, j, k]
                    - coeffs.sp[i, j, k]
                )
                coeffs.dwz[i, j, k] = grid.areaij[i, j] / coeffs.ap[i, j, k]

                # Under-relaxation
                coeffs.ap[i, j, k] = coeffs.ap[i, j, k] / simu_params.urf_vel
                coeffs.su[i, j, k] = coeffs.su[i, j, k] + (1.0 - simu_params.urf_vel) * coeffs.ap[i, j, k] * state.wVel[i, j, k]
                coeffs.dwz[i, j, k] = coeffs.dwz[i, j, k] * simu_params.urf_vel

                # Zero velocity in solid
                twlc = ti.min(state.temp[i, j, k], state.temp[i, j, k - 1])
                if twlc <= physics.tsolid:
                    coeffs.su[i, j, k] = 0.0
                    coeffs.an[i, j, k] = 0.0
                    coeffs.as_[i, j, k] = 0.0
                    coeffs.ae[i, j, k] = 0.0
                    coeffs.aw[i, j, k] = 0.0
                    coeffs.at[i, j, k] = 0.0
                    coeffs.ab[i, j, k] = 0.0
                    coeffs.ap[i, j, k] = 1.0e20


@ti.kernel
def _source_p(
    state: ti.template(),
    mat_props: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
    physics: ti.template(),
):
    i_start = ti.static(simu_params.istatp1)
    i_end = ti.static(min(simu_params.iendm1, simu_params.ni - 2))
    j_start = ti.static(simu_params.jstatp1)
    j_end = ti.static(min(simu_params.jendm1, simu_params.nj - 2))
    k_start = ti.static(simu_params.kstatp1)
    k_end = ti.static(min(simu_params.nkm1, simu_params.nk - 2))

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                coeffs.ap[i, j, k] = (
                    coeffs.an[i, j, k]
                    + coeffs.as_[i, j, k]
                    + coeffs.ae[i, j, k]
                    + coeffs.aw[i, j, k]
                    + coeffs.at[i, j, k]
                    + coeffs.ab[i, j, k]
                    - coeffs.sp[i, j, k]
                )
                if state.temp[i, j, k] <= physics.tsolid:
                    coeffs.su[i, j, k] = 0.0
                    coeffs.ap[i, j, k] = 1.0e20
                    coeffs.an[i, j, k] = 0.0
                    coeffs.as_[i, j, k] = 0.0
                    coeffs.ae[i, j, k] = 0.0
                    coeffs.aw[i, j, k] = 0.0
                    coeffs.at[i, j, k] = 0.0
                    coeffs.ab[i, j, k] = 0.0


@ti.kernel
def _source_h(
    state: ti.template(),
    state_prev: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
    physics: ti.template(),
    laser_params: ti.template(),
    laser_state: ti.template(),
):
    i_start = ti.static(simu_params.istatp1)
    i_end = ti.static(min(simu_params.iendm1, simu_params.ni - 2))
    j_start = ti.static(simu_params.jstatp1)
    j_end = ti.static(min(simu_params.jendm1, simu_params.nj - 2))
    k_start = ti.static(simu_params.kstatp1)
    k_end = ti.static(min(simu_params.nkm1, simu_params.nk - 2))

    # RHF-modified volumetric parameters
    rhf = ti.cast(laser_state.rhf, ti.f64)
    source_depth_rhf = ti.cast(laser_params.source_depth, ti.f64) * rhf * rhf
    alasetavol_rhf = source_depth_rhf * ti.cast(3700.0, ti.f64)
    sourcerad_rhf = source_depth_rhf * ti.cast(0.37, ti.f64)
    pi = ti.cast(3.1415926, ti.f64)

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                sourceinput = ti.cast(0.0, ti.f64)
                if laser_state.laser_on:
                    if (grid.z[simu_params.nkm1] - grid.z[k]) <= source_depth_rhf and source_depth_rhf > 0.0:
                        dx = ti.cast(laser_state.beam_x, ti.f64) - grid.x[i]
                        dy = ti.cast(laser_state.beam_y, ti.f64) - grid.y[j]
                        power_vol = ti.cast(laser_params.power_vol, ti.f64)
                        factor = ti.cast(laser_params.factor, ti.f64)
                        sourceinput = (
                            power_vol
                            * factor
                            / pi
                            / (sourcerad_rhf ** 2)
                            / source_depth_rhf
                            * alasetavol_rhf
                            * ti.exp(-factor / (sourcerad_rhf ** 2) * (dx * dx + dy * dy))
                        )
                coeffs.su[i, j, k] = coeffs.su[i, j, k] + grid.vol[i, j, k] * sourceinput

    # Latent heat source terms
    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                volht = grid.vol[i, j, k] * physics.hlatent * mat_props.den[i, j, k] / simu_params.delt
                coeffs.su[i, j, k] = coeffs.su[i, j, k] - volht * (state.fracl[i, j, k] - state_prev.fraclnot[i, j, k])

                flew = grid.areajk[j, k] * (
                    ti.max(state.uVel[i, j, k], 0.0) * state.fracl[i - 1, j, k]
                    - ti.max(-state.uVel[i, j, k], 0.0) * state.fracl[i, j, k]
                    + ti.max(-state.uVel[i + 1, j, k], 0.0) * state.fracl[i + 1, j, k]
                    - ti.max(state.uVel[i + 1, j, k], 0.0) * state.fracl[i, j, k]
                )
                flns = grid.areaik[i, k] * (
                    ti.max(state.vVel[i, j, k], 0.0) * state.fracl[i, j - 1, k]
                    - ti.max(-state.vVel[i, j, k], 0.0) * state.fracl[i, j, k]
                    + ti.max(-state.vVel[i, j + 1, k], 0.0) * state.fracl[i, j + 1, k]
                    - ti.max(state.vVel[i, j + 1, k], 0.0) * state.fracl[i, j, k]
                )
                fltb = grid.areaij[i, j] * (
                    ti.max(state.wVel[i, j, k], 0.0) * state.fracl[i, j, k - 1]
                    - ti.max(-state.wVel[i, j, k], 0.0) * state.fracl[i, j, k]
                    + ti.max(-state.wVel[i, j, k + 1], 0.0) * state.fracl[i, j, k + 1]
                    - ti.max(state.wVel[i, j, k + 1], 0.0) * state.fracl[i, j, k]
                )
                coeffs.su[i, j, k] = coeffs.su[i, j, k] + mat_props.den[i, j, k] * physics.hlatent * (flew + flns + fltb)

    # k = 1 and k = nk boundaries (transfer to source term)
    k_bottom = ti.static(simu_params.kstatp1)
    k_top_inner = ti.static(min(simu_params.nkm1 - 1, simu_params.nk - 2))
    k_top = ti.static(simu_params.nkm1)
    for j in range(j_start, j_end + 1):
        for i in range(i_start, i_end + 1):
            coeffs.su[i, j, k_bottom] = coeffs.su[i, j, k_bottom] + coeffs.ab[i, j, k_bottom] * state.enthalpy[i, j, k_bottom - 1]
            coeffs.sp[i, j, k_bottom] = coeffs.sp[i, j, k_bottom] - coeffs.ab[i, j, k_bottom]
            coeffs.ab[i, j, k_bottom] = 0.0
            coeffs.su[i, j, k_top_inner] = coeffs.su[i, j, k_top_inner] + coeffs.at[i, j, k_top_inner] * state.enthalpy[i, j, k_top]
            coeffs.sp[i, j, k_top_inner] = coeffs.sp[i, j, k_top_inner] - coeffs.at[i, j, k_top_inner]
            coeffs.at[i, j, k_top_inner] = 0.0

    # j = 1 and j = nj boundaries
    j_bottom = ti.static(simu_params.jstatp1)
    j_top_inner = ti.static(min(simu_params.njm1 - 1, simu_params.nj - 2))
    j_top = ti.static(simu_params.njm1)
    for k in range(k_start, k_end + 1):
        for i in range(i_start, i_end + 1):
            coeffs.su[i, j_bottom, k] = coeffs.su[i, j_bottom, k] + coeffs.as_[i, j_bottom, k] * state.enthalpy[i, j_bottom - 1, k]
            coeffs.sp[i, j_bottom, k] = coeffs.sp[i, j_bottom, k] - coeffs.as_[i, j_bottom, k]
            coeffs.as_[i, j_bottom, k] = 0.0
            coeffs.su[i, j_top_inner, k] = coeffs.su[i, j_top_inner, k] + coeffs.an[i, j_top_inner, k] * state.enthalpy[i, j_top, k]
            coeffs.sp[i, j_top_inner, k] = coeffs.sp[i, j_top_inner, k] - coeffs.an[i, j_top_inner, k]
            coeffs.an[i, j_top_inner, k] = 0.0

    # i = 1 and i = ni boundaries
    i_bottom = ti.static(simu_params.istatp1)
    i_top_inner = ti.static(min(simu_params.nim1 - 1, simu_params.ni - 2))
    i_top = ti.static(simu_params.nim1)
    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            coeffs.su[i_bottom, j, k] = coeffs.su[i_bottom, j, k] + coeffs.aw[i_bottom, j, k] * state.enthalpy[i_bottom - 1, j, k]
            coeffs.sp[i_bottom, j, k] = coeffs.sp[i_bottom, j, k] - coeffs.aw[i_bottom, j, k]
            coeffs.aw[i_bottom, j, k] = 0.0
            coeffs.su[i_top_inner, j, k] = coeffs.su[i_top_inner, j, k] + coeffs.ae[i_top_inner, j, k] * state.enthalpy[i_top, j, k]
            coeffs.sp[i_top_inner, j, k] = coeffs.sp[i_top_inner, j, k] - coeffs.ae[i_top_inner, j, k]
            coeffs.ae[i_top_inner, j, k] = 0.0

    # Final ap and under-relaxation for enthalpy
    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                coeffs.ap[i, j, k] = (
                    coeffs.an[i, j, k]
                    + coeffs.as_[i, j, k]
                    + coeffs.ae[i, j, k]
                    + coeffs.aw[i, j, k]
                    + coeffs.at[i, j, k]
                    + coeffs.ab[i, j, k]
                    + coeffs.apnot[i, j, k]
                    - coeffs.sp[i, j, k]
                )
                coeffs.ap[i, j, k] = coeffs.ap[i, j, k] / simu_params.urf_h
                coeffs.su[i, j, k] = coeffs.su[i, j, k] + (1.0 - simu_params.urf_h) * coeffs.ap[i, j, k] * state.enthalpy[i, j, k]


def source_term(
    ivar: int,
    state: State,
    state_prev: StatePrev,
    grid: GridParams,
    coeffs: DiscretCoeffs,
    mat_props: MaterialProps,
    laser_state: LaserState,
    physics: PhysicsParams,
    sim: SimulationParams,
    laser_params: LaserParams,
) -> DiscretCoeffs:
    """Add source terms to discretized equations. (From mod_sour.f90)"""
    if ivar == 1:
        _source_u(state, state_prev, mat_props, grid, coeffs, sim, physics)
    elif ivar == 2:
        _source_v(state, state_prev, mat_props, grid, coeffs, sim, physics)
    elif ivar == 3:
        _source_w(state, state_prev, mat_props, grid, coeffs, sim, physics)
    elif ivar == 4:
        _source_p(state, mat_props, coeffs, sim, physics)
    elif ivar == 5:
        _source_h(state, state_prev, mat_props, grid, coeffs, sim, physics, laser_params, laser_state)
    else:
        raise ValueError(f"Invalid ivar={ivar}. Must be 1-5.")
    return coeffs

