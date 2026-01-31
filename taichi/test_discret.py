"""
Lightweight test for discretization module.
"""

import taichi as ti
import numpy as np
from pathlib import Path

from discret import discretize
from geom import get_gridparams
from data_structures import GridParams
from data_structures import (
    State,
    StatePrev,
    MaterialProps,
    DiscretCoeffs,
    SimulationParams,
    PhysicsParams,
)


def initialize_fields(state, state_prev, mat_props, grid):
    ni, nj, nk = grid.ni, grid.nj, grid.nk
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                state.uVel[i, j, k] = 0.1
                state.vVel[i, j, k] = 0.05
                state.wVel[i, j, k] = 0.02
                state.pressure[i, j, k] = 1.0
                state.enthalpy[i, j, k] = 1.0
                state.temp[i, j, k] = 300.0
                state.fracl[i, j, k] = 0.0

                state_prev.unot[i, j, k] = state.uVel[i, j, k]
                state_prev.vnot[i, j, k] = state.vVel[i, j, k]
                state_prev.wnot[i, j, k] = state.wVel[i, j, k]
                state_prev.hnot[i, j, k] = state.enthalpy[i, j, k]

                mat_props.vis[i, j, k] = 0.006
                mat_props.diff[i, j, k] = 3.2e-6
                mat_props.den[i, j, k] = 7800.0


def create_uniform_grid(ni=10, nj=10, nk=10):
    """Fallback uniform grid when YAML is unavailable."""
    grid = GridParams(ni, nj, nk)

    xlen, ylen, zlen = 1.0e-3, 1.0e-3, 0.5e-3
    dx = xlen / (ni - 2)
    dy = ylen / (nj - 2)
    dz = zlen / (nk - 2)

    # Faces
    for i in range(ni):
        grid.xu[i] = max(0.0, (i - 1) * dx)
    for j in range(nj):
        grid.yv[j] = max(0.0, (j - 1) * dy)
    for k in range(nk):
        grid.zw[k] = max(0.0, (k - 1) * dz)

    # Centers
    for i in range(ni - 1):
        grid.x[i] = 0.5 * (grid.xu[i + 1] + grid.xu[i])
    grid.x[ni - 1] = grid.xu[ni - 1]
    for j in range(nj - 1):
        grid.y[j] = 0.5 * (grid.yv[j + 1] + grid.yv[j])
    grid.y[nj - 1] = grid.yv[nj - 1]
    for k in range(nk - 1):
        grid.z[k] = 0.5 * (grid.zw[k + 1] + grid.zw[k])
    grid.z[nk - 1] = grid.zw[nk - 1]

    # Inverses and spacing
    for i in range(1, ni):
        grid.dxpwinv[i] = 1.0 / (grid.x[i] - grid.x[i - 1])
    for i in range(0, ni - 1):
        grid.dxpeinv[i] = 1.0 / (grid.x[i + 1] - grid.x[i])
    for j in range(1, nj):
        grid.dypsinv[j] = 1.0 / (grid.y[j] - grid.y[j - 1])
    for j in range(0, nj - 1):
        grid.dypninv[j] = 1.0 / (grid.y[j + 1] - grid.y[j])
    for k in range(1, nk):
        grid.dzpbinv[k] = 1.0 / (grid.z[k] - grid.z[k - 1])
    for k in range(0, nk - 1):
        grid.dzptinv[k] = 1.0 / (grid.z[k + 1] - grid.z[k])

    # Interpolation fractions
    for i in range(0, ni - 1):
        grid.fracx[i] = 0.5
    for j in range(0, nj - 1):
        grid.fracy[j] = 0.5
    for k in range(0, nk - 1):
        grid.fracz[k] = 0.5

    # Areas and volumes
    for j in range(1, nj - 1):
        for i in range(1, ni - 1):
            grid.areaij[i, j] = dx * dy
            grid.areauij[i, j] = dx * dy
            grid.areavij[i, j] = dx * dy
    for k in range(1, nk - 1):
        for i in range(1, ni - 1):
            grid.areaik[i, k] = dx * dz
            grid.areauik[i, k] = dx * dz
            grid.areawik[i, k] = dx * dz
    for k in range(1, nk - 1):
        for j in range(1, nj - 1):
            grid.areajk[j, k] = dy * dz
            grid.areavjk[j, k] = dy * dz
            grid.areawjk[j, k] = dy * dz
    for k in range(1, nk - 1):
        for j in range(1, nj - 1):
            for i in range(1, ni - 1):
                grid.vol[i, j, k] = dx * dy * dz
                grid.volume_u[i, j, k] = dx * dy * dz
                grid.volume_v[i, j, k] = dx * dy * dz
                grid.volume_w[i, j, k] = dx * dy * dz

    return grid


def run_discret_test():
    ti.init(arch=ti.cpu)
    # Use repo-relative YAML so it works across OS paths
    repo_root = Path(__file__).resolve().parent.parent
    yaml_path = repo_root / "inputfile" / "input_param.yaml"
    if yaml_path.exists():
        grid = get_gridparams(str(yaml_path))
    else:
        print(f"Config not found: {yaml_path}. Using uniform fallback grid.")
        grid = create_uniform_grid()
    ni, nj, nk = grid.ni, grid.nj, grid.nk

    state = State(ni, nj, nk)
    state_prev = StatePrev(ni, nj, nk)
    mat_props = MaterialProps(ni, nj, nk)
    coeffs = DiscretCoeffs(ni, nj, nk)
    sim = SimulationParams(ni=ni, nj=nj, nk=nk)
    physics = PhysicsParams()

    initialize_fields(state, state_prev, mat_props, grid)

    # Run a couple of ivar cases to ensure kernels compile
    discretize(1, state, state_prev, grid, coeffs, mat_props, sim, physics)
    discretize(5, state, state_prev, grid, coeffs, mat_props, sim, physics)

    print("discretize() executed for ivar=1 and ivar=5")


if __name__ == "__main__":
    run_discret_test()

