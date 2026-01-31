"""
AM-CFD Taichi Implementation - Discretization

Converted from Fortran module: mod_discret.f90

Builds finite-volume discretization coefficients for momentum, pressure,
and enthalpy equations using a power-law scheme.
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
)


@ti.func
def _power_law_coeff(d: ti.f64, f: ti.f64) -> ti.f64:
    val = ti.cast(0.0, ti.f64)
    if d > ti.cast(0.0, ti.f64):
        one = ti.cast(1.0, ti.f64)
        coef = ti.cast(0.1, ti.f64)
        zero = ti.cast(0.0, ti.f64)
        val = d * ti.max(zero, (one - coef * ti.abs(f / d)) ** 5)
    return val


@ti.kernel
def _discretize_u(
    state: ti.template(),
    state_prev: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
):
    ni = ti.static(simu_params.ni)
    nj = ti.static(simu_params.nj)
    nk = ti.static(simu_params.nk)
    i_start = ti.static(simu_params.istatp1)
    i_end = ti.static(min(simu_params.iendm1, simu_params.ni - 2))
    j_start = ti.static(simu_params.jstatp1)
    j_end = ti.static(min(simu_params.jendm1, simu_params.nj - 2))
    k_start = ti.static(simu_params.kstatp1)
    k_end = ti.static(min(simu_params.nkm1, simu_params.nk - 2))

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                # Velocities at CV faces
                fracx = grid.fracx[i - 1]
                vn = state.vVel[i, j + 1, k] * (1.0 - fracx) + state.vVel[i - 1, j + 1, k] * fracx
                vs = state.vVel[i, j, k] * (1.0 - fracx) + state.vVel[i - 1, j, k] * fracx
                ue = 0.5 * (state.uVel[i + 1, j, k] + state.uVel[i, j, k])
                uw = 0.5 * (state.uVel[i - 1, j, k] + state.uVel[i, j, k])
                wt = state.wVel[i, j, k + 1] * (1.0 - fracx) + state.wVel[i - 1, j, k + 1] * fracx
                wb = state.wVel[i, j, k] * (1.0 - fracx) + state.wVel[i - 1, j, k] * fracx

                # Convection coefficients
                fn = vn * mat_props.den[i, j, k] * grid.areauik[i, k]
                fs = vs * mat_props.den[i, j, k] * grid.areauik[i, k]
                fe = ue * mat_props.den[i, j, k] * grid.areajk[j, k]
                fw = uw * mat_props.den[i, j, k] * grid.areajk[j, k]
                ft = wt * mat_props.den[i, j, k] * grid.areauij[i, j]
                fb = wb * mat_props.den[i, j, k] * grid.areauij[i, j]

                # Viscosity at CV faces
                visu = mat_props.vis[i, j, k] * mat_props.vis[i - 1, j, k] / (
                    fracx * mat_props.vis[i - 1, j, k] + (1.0 - fracx) * mat_props.vis[i, j, k]
                )
                visun = mat_props.vis[i, j + 1, k] * mat_props.vis[i - 1, j + 1, k] / (
                    fracx * mat_props.vis[i - 1, j + 1, k] + (1.0 - fracx) * mat_props.vis[i, j + 1, k]
                )
                visus = mat_props.vis[i, j - 1, k] * mat_props.vis[i - 1, j - 1, k] / (
                    fracx * mat_props.vis[i - 1, j - 1, k] + (1.0 - fracx) * mat_props.vis[i, j - 1, k]
                )
                visn = visu
                viss = visu
                if j != j_end:
                    visn = visu * visun / ((1.0 - grid.fracy[j]) * visun + grid.fracy[j] * visu)
                if j != j_start:
                    viss = visu * visus / (grid.fracy[j - 1] * visus + (1.0 - grid.fracy[j - 1]) * visu)

                vise = mat_props.vis[i, j, k]
                visw = mat_props.vis[i - 1, j, k]
                visut = mat_props.vis[i, j, k + 1] * mat_props.vis[i - 1, j, k + 1] / (
                    fracx * mat_props.vis[i - 1, j, k + 1] + (1.0 - fracx) * mat_props.vis[i, j, k + 1]
                )
                visub = mat_props.vis[i, j, k - 1] * mat_props.vis[i - 1, j, k - 1] / (
                    fracx * mat_props.vis[i - 1, j, k - 1] + (1.0 - fracx) * mat_props.vis[i, j, k - 1]
                )
                vist = visu
                visb = visu
                if k != k_end:
                    vist = visut * visu / (grid.fracz[k] * visu + (1.0 - grid.fracz[k]) * visut)
                if k != k_start:
                    visb = visub * visu / ((1.0 - grid.fracz[k - 1]) * visu + grid.fracz[k - 1] * visub)

                # Diffusion coefficients
                dn = visn * grid.areauik[i, k] * grid.dypsinv[j + 1]
                ds = viss * grid.areauik[i, k] * grid.dypsinv[j]
                de = vise * grid.areajk[j, k] / (grid.xu[i + 1] - grid.xu[i])
                dw = visw * grid.areajk[j, k] / (grid.xu[i] - grid.xu[i - 1])
                dt = vist * grid.areauij[i, j] * grid.dzpbinv[k + 1]
                db = visb * grid.areauij[i, j] * grid.dzpbinv[k]

                coeffs.an[i, j, k] = _power_law_coeff(dn, fn) + ti.max(0.0, -fn)
                coeffs.as_[i, j, k] = _power_law_coeff(ds, fs) + ti.max(0.0, fs)
                coeffs.ae[i, j, k] = _power_law_coeff(de, fe) + ti.max(0.0, -fe)
                coeffs.aw[i, j, k] = _power_law_coeff(dw, fw) + ti.max(0.0, fw)
                coeffs.at[i, j, k] = _power_law_coeff(dt, ft) + ti.max(0.0, -ft)
                coeffs.ab[i, j, k] = _power_law_coeff(db, fb) + ti.max(0.0, fb)
                coeffs.apnot[i, j, k] = mat_props.den[i, j, k] * grid.volume_u[i, j, k] / simu_params.delt

                delf = fn - fs + fe - fw + ft - fb
                cp0 = ti.max(0.0, delf)
                cp1 = ti.min(0.0, delf)
                coeffs.su[i, j, k] = -cp1 * state.uVel[i, j, k]
                coeffs.su[i, j, k] += grid.areajk[j, k] * (
                    state.pressure[i - 1, j, k] - state.pressure[i, j, k]
                )
                coeffs.sp[i, j, k] = -cp0
                coeffs.su[i, j, k] += coeffs.apnot[i, j, k] * state_prev.unot[i, j, k]

                dudxp = (state.uVel[i + 1, j, k] - state.uVel[i, j, k]) / (grid.xu[i + 1] - grid.xu[i])
                dudxm = (state.uVel[i, j, k] - state.uVel[i - 1, j, k]) / (grid.xu[i] - grid.xu[i - 1])
                coeffs.su[i, j, k] += (vise * dudxp - visw * dudxm) * grid.areajk[j, k]

                dvdxp = (state.vVel[i, j + 1, k] - state.vVel[i - 1, j + 1, k]) * grid.dxpwinv[i]
                dvdxm = (state.vVel[i, j, k] - state.vVel[i - 1, j, k]) * grid.dxpwinv[i]
                coeffs.su[i, j, k] += (visn * dvdxp - viss * dvdxm) * grid.areauik[i, k]

                dwdxp = (state.wVel[i, j, k + 1] - state.wVel[i - 1, j, k + 1]) * grid.dxpwinv[i]
                dwdxm = (state.wVel[i, j, k] - state.wVel[i - 1, j, k]) * grid.dxpwinv[i]
                coeffs.su[i, j, k] += (vist * dwdxp - visb * dwdxm) * grid.areauij[i, j]


@ti.kernel
def _discretize_v(
    state: ti.template(),
    state_prev: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
):
    ni = ti.static(simu_params.ni)
    nj = ti.static(simu_params.nj)
    nk = ti.static(simu_params.nk)
    i_start = ti.static(simu_params.istatp1)
    i_end = ti.static(min(simu_params.iendm1, simu_params.ni - 2))
    j_start = ti.static(simu_params.jstatp1)
    j_end = ti.static(min(simu_params.jendm1, simu_params.nj - 2))
    k_start = ti.static(simu_params.kstatp1)
    k_end = ti.static(min(simu_params.nkm1, simu_params.nk - 2))

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                vn = 0.5 * (state.vVel[i, j, k] + state.vVel[i, j + 1, k])
                vs = 0.5 * (state.vVel[i, j, k] + state.vVel[i, j - 1, k])
                fracy = grid.fracy[j - 1]
                ue = state.uVel[i + 1, j, k] * (1.0 - fracy) + state.uVel[i + 1, j - 1, k] * fracy
                uw = state.uVel[i, j, k] * (1.0 - fracy) + state.uVel[i, j - 1, k] * fracy
                wt = state.wVel[i, j, k + 1] * (1.0 - fracy) + state.wVel[i, j - 1, k + 1] * fracy
                wb = state.wVel[i, j, k] * (1.0 - fracy) + state.wVel[i, j - 1, k] * fracy

                fn = vn * mat_props.den[i, j, k] * grid.areaik[i, k]
                fs = vs * mat_props.den[i, j, k] * grid.areaik[i, k]
                fe = ue * mat_props.den[i, j, k] * grid.areavjk[j, k]
                fw = uw * mat_props.den[i, j, k] * grid.areavjk[j, k]
                ft = wt * mat_props.den[i, j, k] * grid.areavij[i, j]
                fb = wb * mat_props.den[i, j, k] * grid.areavij[i, j]

                visn = mat_props.vis[i, j, k]
                viss = mat_props.vis[i, j - 1, k]
                visv = mat_props.vis[i, j, k] * mat_props.vis[i, j - 1, k] / (
                    fracy * mat_props.vis[i, j - 1, k] + (1.0 - fracy) * mat_props.vis[i, j, k]
                )
                visve = mat_props.vis[i + 1, j, k] * mat_props.vis[i + 1, j - 1, k] / (
                    fracy * mat_props.vis[i + 1, j - 1, k] + (1.0 - fracy) * mat_props.vis[i + 1, j, k]
                )
                visvw = mat_props.vis[i - 1, j, k] * mat_props.vis[i - 1, j - 1, k] / (
                    fracy * mat_props.vis[i - 1, j - 1, k] + (1.0 - fracy) * mat_props.vis[i - 1, j, k]
                )
                vise = visv
                visw = visv
                if i != i_end:
                    vise = visv * visve / ((1.0 - grid.fracx[i]) * visve + grid.fracx[i] * visv)
                if i != i_start:
                    visw = visv * visvw / (grid.fracx[i - 1] * visvw + (1.0 - grid.fracx[i - 1]) * visv)

                visvt = mat_props.vis[i, j, k + 1] * mat_props.vis[i, j - 1, k + 1] / (
                    fracy * mat_props.vis[i, j - 1, k + 1] + (1.0 - fracy) * mat_props.vis[i, j, k + 1]
                )
                visvb = mat_props.vis[i, j, k - 1] * mat_props.vis[i, j - 1, k - 1] / (
                    fracy * mat_props.vis[i, j - 1, k - 1] + (1.0 - fracy) * mat_props.vis[i, j, k - 1]
                )
                vist = visv
                visb = visv
                if k != k_end:
                    vist = visv * visvt / ((1.0 - grid.fracz[k]) * visvt + grid.fracz[k] * visv)
                if k != k_start:
                    visb = visv * visvb / (grid.fracz[k - 1] * visvb + (1.0 - grid.fracz[k - 1]) * visv)

                dn = visn * grid.areaik[i, k] / (grid.yv[j + 1] - grid.yv[j])
                ds = viss * grid.areaik[i, k] / (grid.yv[j] - grid.yv[j - 1])
                de = vise * grid.areavjk[j, k] * grid.dxpwinv[i + 1]
                dw = visw * grid.areavjk[j, k] * grid.dxpwinv[i]
                dt = vist * grid.areavij[i, j] * grid.dzpbinv[k + 1]
                db = visb * grid.areavij[i, j] * grid.dzpbinv[k]

                coeffs.an[i, j, k] = _power_law_coeff(dn, fn) + ti.max(0.0, -fn)
                coeffs.as_[i, j, k] = _power_law_coeff(ds, fs) + ti.max(0.0, fs)
                coeffs.ae[i, j, k] = _power_law_coeff(de, fe) + ti.max(0.0, -fe)
                coeffs.aw[i, j, k] = _power_law_coeff(dw, fw) + ti.max(0.0, fw)
                coeffs.at[i, j, k] = _power_law_coeff(dt, ft) + ti.max(0.0, -ft)
                coeffs.ab[i, j, k] = _power_law_coeff(db, fb) + ti.max(0.0, fb)
                coeffs.apnot[i, j, k] = mat_props.den[i, j, k] * grid.volume_v[i, j, k] / simu_params.delt

                delf = fn - fs + fe - fw + ft - fb
                cp0 = ti.max(0.0, delf)
                cp1 = ti.min(0.0, delf)
                coeffs.su[i, j, k] = -cp1 * state.vVel[i, j, k]
                coeffs.su[i, j, k] += grid.areaik[i, k] * (
                    state.pressure[i, j - 1, k] - state.pressure[i, j, k]
                )
                coeffs.sp[i, j, k] = -cp0
                coeffs.su[i, j, k] += coeffs.apnot[i, j, k] * state_prev.vnot[i, j, k]

                dudyp = (state.uVel[i + 1, j, k] - state.uVel[i + 1, j - 1, k]) * grid.dypsinv[j]
                dudym = (state.uVel[i, j, k] - state.uVel[i, j - 1, k]) * grid.dypsinv[j]
                coeffs.su[i, j, k] += (vise * dudyp - visw * dudym) * grid.areavjk[j, k]

                dvdyp = (state.vVel[i, j + 1, k] - state.vVel[i, j, k]) / (grid.yv[j + 1] - grid.yv[j])
                dvdym = (state.vVel[i, j, k] - state.vVel[i, j - 1, k]) / (grid.yv[j] - grid.yv[j - 1])
                coeffs.su[i, j, k] += (visn * dvdyp - viss * dvdym) * grid.areaik[i, k]

                dwdyp = (state.wVel[i, j, k + 1] - state.wVel[i, j - 1, k + 1]) * grid.dypsinv[j]
                dwdym = (state.wVel[i, j, k] - state.wVel[i, j - 1, k]) * grid.dypsinv[j]
                coeffs.su[i, j, k] += (vist * dwdyp - visb * dwdym) * grid.areavij[i, j]


@ti.kernel
def _discretize_w(
    state: ti.template(),
    state_prev: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
):
    ni = ti.static(simu_params.ni)
    nj = ti.static(simu_params.nj)
    nk = ti.static(simu_params.nk)
    i_start = ti.static(simu_params.istatp1)
    i_end = ti.static(min(simu_params.iendm1, simu_params.ni - 2))
    j_start = ti.static(simu_params.jstatp1)
    j_end = ti.static(min(simu_params.jendm1, simu_params.nj - 2))
    k_start = ti.static(simu_params.kstatp1)
    k_end = ti.static(min(simu_params.nkm1, simu_params.nk - 2))

    for k in range(k_start, k_end + 1):
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                fracz = grid.fracz[k - 1]
                vn = state.vVel[i, j + 1, k] * (1.0 - fracz) + state.vVel[i, j + 1, k - 1] * fracz
                vs = state.vVel[i, j, k] * (1.0 - fracz) + state.vVel[i, j, k - 1] * fracz
                ue = state.uVel[i + 1, j, k] * (1.0 - fracz) + state.uVel[i + 1, j, k - 1] * fracz
                uw = state.uVel[i, j, k] * (1.0 - fracz) + state.uVel[i, j, k - 1] * fracz
                wt = 0.5 * (state.wVel[i, j, k] + state.wVel[i, j, k + 1])
                wb = 0.5 * (state.wVel[i, j, k] + state.wVel[i, j, k - 1])

                fn = vn * mat_props.den[i, j, k] * grid.areawik[i, k]
                fs = vs * mat_props.den[i, j, k] * grid.areawik[i, k]
                fe = ue * mat_props.den[i, j, k] * grid.areawjk[j, k]
                fw = uw * mat_props.den[i, j, k] * grid.areawjk[j, k]
                ft = wt * mat_props.den[i, j, k] * grid.areaij[i, j]
                fb = wb * mat_props.den[i, j, k] * grid.areaij[i, j]

                vis_w = mat_props.vis[i, j, k] * mat_props.vis[i, j, k - 1] / (
                    fracz * mat_props.vis[i, j, k - 1] + (1.0 - fracz) * mat_props.vis[i, j, k]
                )
                viswn = mat_props.vis[i, j + 1, k] * mat_props.vis[i, j + 1, k - 1] / (
                    fracz * mat_props.vis[i, j + 1, k - 1] + (1.0 - fracz) * mat_props.vis[i, j + 1, k]
                )
                visws = mat_props.vis[i, j - 1, k] * mat_props.vis[i, j - 1, k - 1] / (
                    fracz * mat_props.vis[i, j - 1, k - 1] + (1.0 - fracz) * mat_props.vis[i, j - 1, k]
                )
                visn = vis_w
                viss = vis_w
                if j != j_end:
                    visn = vis_w * viswn / ((1.0 - grid.fracy[j]) * viswn + grid.fracy[j] * vis_w)
                if j != j_start:
                    viss = vis_w * visws / (grid.fracy[j - 1] * visws + (1.0 - grid.fracy[j - 1]) * vis_w)

                viswe = mat_props.vis[i + 1, j, k] * mat_props.vis[i + 1, j, k - 1] / (
                    fracz * mat_props.vis[i + 1, j, k - 1] + (1.0 - fracz) * mat_props.vis[i + 1, j, k]
                )
                visww = mat_props.vis[i - 1, j, k] * mat_props.vis[i - 1, j, k - 1] / (
                    fracz * mat_props.vis[i - 1, j, k - 1] + (1.0 - fracz) * mat_props.vis[i - 1, j, k]
                )
                vise = vis_w
                visw = vis_w
                if i != i_end:
                    vise = vis_w * viswe / ((1.0 - grid.fracx[i]) * viswe + grid.fracx[i] * vis_w)
                if i != i_start:
                    visw = vis_w * visww / (grid.fracx[i - 1] * visww + (1.0 - grid.fracx[i - 1]) * vis_w)

                vist = mat_props.vis[i, j, k]
                visb = mat_props.vis[i, j, k - 1]

                dn = visn * grid.areawik[i, k] * grid.dypsinv[j + 1]
                ds = viss * grid.areawik[i, k] * grid.dypsinv[j]
                de = vise * grid.areawjk[j, k] * grid.dxpwinv[i + 1]
                dw = visw * grid.areawjk[j, k] * grid.dxpwinv[i]
                dt = vist * grid.areaij[i, j] / (grid.zw[k + 1] - grid.zw[k])
                db = visb * grid.areaij[i, j] / (grid.zw[k] - grid.zw[k - 1])

                coeffs.an[i, j, k] = _power_law_coeff(dn, fn) + ti.max(0.0, -fn)
                coeffs.as_[i, j, k] = _power_law_coeff(ds, fs) + ti.max(0.0, fs)
                coeffs.ae[i, j, k] = _power_law_coeff(de, fe) + ti.max(0.0, -fe)
                coeffs.aw[i, j, k] = _power_law_coeff(dw, fw) + ti.max(0.0, fw)
                coeffs.at[i, j, k] = _power_law_coeff(dt, ft) + ti.max(0.0, -ft)
                coeffs.ab[i, j, k] = _power_law_coeff(db, fb) + ti.max(0.0, fb)
                coeffs.apnot[i, j, k] = mat_props.den[i, j, k] * grid.volume_w[i, j, k] / simu_params.delt

                delf = fn - fs + fe - fw + ft - fb
                cp0 = ti.max(0.0, delf)
                cp1 = ti.min(0.0, delf)
                coeffs.su[i, j, k] = -cp1 * state.wVel[i, j, k]
                coeffs.su[i, j, k] += grid.areaij[i, j] * (
                    state.pressure[i, j, k - 1] - state.pressure[i, j, k]
                )
                coeffs.sp[i, j, k] = -cp0
                coeffs.su[i, j, k] += coeffs.apnot[i, j, k] * state_prev.wnot[i, j, k]

                dudzp = (state.uVel[i + 1, j, k] - state.uVel[i + 1, j, k - 1]) * grid.dzpbinv[k]
                dudzm = (state.uVel[i, j, k] - state.uVel[i, j, k - 1]) * grid.dzpbinv[k]
                coeffs.su[i, j, k] += (vise * dudzp - visw * dudzm) * grid.areawjk[j, k]

                dvdzp = (state.vVel[i, j + 1, k] - state.vVel[i, j + 1, k - 1]) * grid.dzpbinv[k]
                dvdzm = (state.vVel[i, j, k] - state.vVel[i, j, k - 1]) * grid.dzpbinv[k]
                coeffs.su[i, j, k] += (visn * dvdzp - viss * dvdzm) * grid.areawik[i, k]

                dwdzp = (state.wVel[i, j, k + 1] - state.wVel[i, j, k]) / (grid.zw[k + 1] - grid.zw[k])
                dwdzm = (state.wVel[i, j, k] - state.wVel[i, j, k - 1]) / (grid.zw[k] - grid.zw[k - 1])
                coeffs.su[i, j, k] += (vist * dwdzp - visb * dwdzm) * grid.areaij[i, j]


@ti.kernel
def _discretize_p(
    state: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
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
                coeffs.an[i, j, k] = grid.areaik[i, k] * coeffs.dvy[i, j + 1, k] * mat_props.den[i, j, k]
                coeffs.as_[i, j, k] = grid.areaik[i, k] * coeffs.dvy[i, j, k] * mat_props.den[i, j, k]
                coeffs.ae[i, j, k] = grid.areajk[j, k] * coeffs.dux[i + 1, j, k] * mat_props.den[i, j, k]
                coeffs.aw[i, j, k] = grid.areajk[j, k] * coeffs.dux[i, j, k] * mat_props.den[i, j, k]
                coeffs.at[i, j, k] = grid.areaij[i, j] * coeffs.dwz[i, j, k + 1] * mat_props.den[i, j, k]
                coeffs.ab[i, j, k] = grid.areaij[i, j] * coeffs.dwz[i, j, k] * mat_props.den[i, j, k]

                vn = state.vVel[i, j + 1, k]
                vs = state.vVel[i, j, k]
                ue = state.uVel[i + 1, j, k]
                uw = state.uVel[i, j, k]
                wt = state.wVel[i, j, k + 1]
                wb = state.wVel[i, j, k]

                fn = vn * grid.areaik[i, k] * mat_props.den[i, j, k]
                fs = vs * grid.areaik[i, k] * mat_props.den[i, j, k]
                fe = ue * grid.areajk[j, k] * mat_props.den[i, j, k]
                fw = uw * grid.areajk[j, k] * mat_props.den[i, j, k]
                ft = wt * grid.areaij[i, j] * mat_props.den[i, j, k]
                fb = wb * grid.areaij[i, j] * mat_props.den[i, j, k]

                delf = fn - fs + fe - fw + ft - fb
                coeffs.sp[i, j, k] = 0.0
                coeffs.su[i, j, k] = -delf


@ti.kernel
def _discretize_h(
    state: ti.template(),
    state_prev: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
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
                vn = state.vVel[i, j + 1, k]
                ue = state.uVel[i + 1, j, k]
                wt = state.wVel[i, j, k + 1]

                fn = mat_props.den[i, j, k] * vn * grid.areaik[i, k]
                fe = mat_props.den[i, j, k] * ue * grid.areajk[j, k]
                ft = mat_props.den[i, j, k] * wt * grid.areaij[i, j]

                difn = mat_props.diff[i, j, k]
                dife = mat_props.diff[i, j, k]
                dift = mat_props.diff[i, j, k]
                if j != j_end:
                    difn = mat_props.diff[i, j, k] * mat_props.diff[i, j + 1, k] / (
                        (1.0 - grid.fracy[j]) * mat_props.diff[i, j + 1, k] + grid.fracy[j] * mat_props.diff[i, j, k]
                    )
                if i != i_end:
                    dife = mat_props.diff[i, j, k] * mat_props.diff[i + 1, j, k] / (
                        (1.0 - grid.fracx[i]) * mat_props.diff[i + 1, j, k] + grid.fracx[i] * mat_props.diff[i, j, k]
                    )
                if k != k_end:
                    dift = mat_props.diff[i, j, k] * mat_props.diff[i, j, k + 1] / (
                        (1.0 - grid.fracz[k]) * mat_props.diff[i, j, k + 1] + grid.fracz[k] * mat_props.diff[i, j, k]
                    )

                dn = difn * grid.areaik[i, k] * grid.dypsinv[j + 1]
                de = dife * grid.areajk[j, k] * grid.dxpwinv[i + 1]
                dt = dift * grid.areaij[i, j] * grid.dzpbinv[k + 1]

                tmp1 = _power_law_coeff(dn, fn)
                coeffs.an[i, j, k] = tmp1 + ti.max(0.0, -fn)
                coeffs.as_[i, j + 1, k] = tmp1 + ti.max(0.0, fn)

                tmp1 = _power_law_coeff(de, fe)
                coeffs.ae[i, j, k] = tmp1 + ti.max(0.0, -fe)
                coeffs.aw[i + 1, j, k] = tmp1 + ti.max(0.0, fe)

                tmp1 = _power_law_coeff(dt, ft)
                coeffs.at[i, j, k] = tmp1 + ti.max(0.0, -ft)
                coeffs.ab[i, j, k + 1] = tmp1 + ti.max(0.0, ft)

                coeffs.apnot[i, j, k] = mat_props.den[i, j, k] * grid.vol[i, j, k] / simu_params.delt
                coeffs.sp[i, j, k] = 0.0
                coeffs.su[i, j, k] = coeffs.apnot[i, j, k] * state_prev.hnot[i, j, k]


@ti.kernel
def _discretize_h_boundaries(
    state: ti.template(),
    mat_props: ti.template(),
    grid: ti.template(),
    coeffs: ti.template(),
    simu_params: ti.template(),
):
    j0 = ti.static(simu_params.jstatp1)
    i0 = ti.static(simu_params.istatp1)
    k0 = ti.static(simu_params.kstatp1)
    i_end = ti.static(min(simu_params.iendm1, simu_params.ni - 2))
    j_end = ti.static(min(simu_params.jendm1, simu_params.nj - 2))
    k_end = ti.static(min(simu_params.nkm1, simu_params.nk - 2))

    # j = 2 plane (south)
    j_plane = j0
    for k in range(k0, k_end + 1):
        for i in range(i0, i_end + 1):
            vs = state.vVel[i, j_plane, k]
            fs = mat_props.den[i, j_plane, k] * vs * grid.areaik[i, k]
            ds = mat_props.diff[i, j_plane, k] * grid.areaik[i, k] * grid.dypsinv[j_plane]
            coeffs.as_[i, j_plane, k] = _power_law_coeff(ds, fs) + ti.max(0.0, fs)

    # i = 2 plane (west)
    i_plane = i0
    for k in range(k0, k_end + 1):
        for j in range(j0, j_end + 1):
            uw = state.uVel[i_plane, j, k]
            fw = mat_props.den[i_plane, j, k] * uw * grid.areajk[j, k]
            dw = mat_props.diff[i_plane, j, k] * grid.areajk[j, k] * grid.dxpwinv[i_plane]
            coeffs.aw[i_plane, j, k] = _power_law_coeff(dw, fw) + ti.max(0.0, fw)

    # k = 2 plane (bottom)
    k_plane = k0
    for j in range(j0, j_end + 1):
        for i in range(i0, i_end + 1):
            wb = state.wVel[i, j, k_plane]
            fb = mat_props.den[i, j, k_plane] * wb * grid.areaij[i, j]
            db = mat_props.diff[i, j, k_plane] * grid.areaij[i, j] * grid.dzpbinv[k_plane]
            coeffs.ab[i, j, k_plane] = _power_law_coeff(db, fb) + ti.max(0.0, fb)


def discretize(
    ivar: int,
    state: State,
    state_prev: StatePrev,
    grid: GridParams,
    coeffs: DiscretCoeffs,
    mat_props: MaterialProps,
    sim: SimulationParams,
    physics: PhysicsParams,
) -> DiscretCoeffs:
    """Build discretization coefficients. (From mod_discret.f90)"""
    if ivar == 1:
        _discretize_u(state, state_prev, mat_props, grid, coeffs, sim)
    elif ivar == 2:
        _discretize_v(state, state_prev, mat_props, grid, coeffs, sim)
    elif ivar == 3:
        _discretize_w(state, state_prev, mat_props, grid, coeffs, sim)
    elif ivar == 4:
        _discretize_p(state, mat_props, grid, coeffs, sim)
    elif ivar == 5:
        _discretize_h(state, state_prev, mat_props, grid, coeffs, sim)
        _discretize_h_boundaries(state, mat_props, grid, coeffs, sim)
    else:
        raise ValueError(f"Invalid ivar={ivar}. Must be 1-5.")
    return coeffs

