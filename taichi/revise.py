"""
AM-CFD Taichi Implementation - Pressure-Velocity Coupling

Translated from Fortran mod_revise.f90
Version 1.0 - Python/Taichi port

This module implements the SIMPLE algorithm for pressure-velocity coupling.
It corrects velocities based on pressure gradients and updates the pressure field.
"""

import taichi as ti
from typing import Tuple

from data_structures import (
    State, GridParams, DiscretCoeffs, SimulationParams, PhysicsParams
)


def revision_p(ivar: int, state: State, coeffs: DiscretCoeffs, 
               grid: GridParams, sim: SimulationParams, physics: PhysicsParams,
               istat: int, iend: int, jstat: int, jend: int, 
               kstat: int, kend: int) -> State:
    """
    Pressure correction using SIMPLE algorithm.
    
    This function implements the pressure-velocity coupling:
    1. Corrects velocities based on pressure correction gradients (only in liquid)
    2. Updates absolute pressure field with under-relaxation
    3. Resets pressure correction to zero
    
    Corresponds to subroutine revision_p in mod_revise.f90
    
    Args:
        ivar: Variable index (only executes for ivar=4, pressure)
        state: Current state with velocity, pressure, temperature fields
        coeffs: Discretization coefficients (contains dux, dvy, dwz)
        grid: Grid parameters
        sim: Simulation parameters (contains urf_p)
        physics: Physical parameters (contains tsolid)
        istat, iend: Start/end indices in x-direction
        jstat, jend: Start/end indices in y-direction
        kstat, kend: Start/end indices in z-direction
        
    Returns:
        state: Updated state with corrected velocities and pressure
    """
    # Only execute for pressure variable (ivar=4)
    # Fortran: goto (500,500,500,400,500)ivar
    if ivar != 4:
        return state
    
    # Python indexing: Fortran's kstat:nkm1 becomes kstat:nk in Python
    nkm1 = grid.nk - 1
    
    # Step 1: Correct velocities based on pressure correction gradients
    correct_velocities_kernel(
        state.uVel, state.vVel, state.wVel,
        state.pp, state.temp,
        coeffs.dux, coeffs.dvy, coeffs.dwz,
        istat, iend, jstat, jend, kstat, nkm1,
        physics.tsolid
    )
    
    # Step 2: Update pressure and reset pressure correction
    update_pressure_kernel(
        state.pressure, state.pp, state.temp,
        istat, iend, jstat, jend, kstat, nkm1,
        sim.urf_p, physics.tsolid
    )
    
    return state


@ti.kernel
def correct_velocities_kernel(
    uVel: ti.template(), vVel: ti.template(), wVel: ti.template(),
    pp: ti.template(), temp: ti.template(),
    dux: ti.template(), dvy: ti.template(), dwz: ti.template(),
    istat: int, iend: int, jstat: int, jend: int, kstat: int, nkm1: int,
    tsolid: ti.f64
):
    """
    Taichi kernel to correct velocities based on pressure gradients.
    
    Velocity corrections (only in liquid regions):
    - u: uVel += dux * (pp[i-1] - pp[i])
    - v: vVel += dvy * (pp[j-1] - pp[j])
    - w: wVel += dwz * (pp[k-1] - pp[k])
    
    Args:
        uVel, vVel, wVel: Velocity fields
        pp: Pressure correction field
        temp: Temperature field
        dux, dvy, dwz: Velocity correction coefficients
        istat, iend, jstat, jend, kstat, nkm1: Loop bounds
        tsolid: Solidus temperature [K]
    """
    # Fortran loops: k=kstat,nkm1; j=jstat,jend; i=istatp1,iendm1
    # istatp1 = istat+1, iendm1 = iend-1
    for k in range(kstat, nkm1 + 1):
        for j in range(jstat, jend + 1):
            for i in range(istat + 1, iend):
                # u-velocity correction (uses temp at i and i-1)
                tulc = ti.min(temp[i, j, k], temp[i - 1, j, k])
                if tulc > tsolid:
                    uVel[i, j, k] += dux[i, j, k] * (pp[i - 1, j, k] - pp[i, j, k])
                
                # v-velocity correction (uses temp at j and j-1)
                tvlc = ti.min(temp[i, j, k], temp[i, j - 1, k])
                if tvlc > tsolid:
                    vVel[i, j, k] += dvy[i, j, k] * (pp[i, j - 1, k] - pp[i, j, k])
                
                # w-velocity correction (uses temp at k and k-1)
                twlc = ti.min(temp[i, j, k], temp[i, j, k - 1])
                if twlc > tsolid:
                    wVel[i, j, k] += dwz[i, j, k] * (pp[i, j, k - 1] - pp[i, j, k])


@ti.kernel
def update_pressure_kernel(
    pressure: ti.template(), pp: ti.template(), temp: ti.template(),
    istat: int, iend: int, jstat: int, jend: int, kstat: int, nkm1: int,
    urf_p: ti.f64, tsolid: ti.f64
):
    """
    Taichi kernel to update pressure and reset pressure correction.
    
    Updates pressure with under-relaxation (only in liquid):
    - pressure += urf_p * pp
    - pp = 0
    
    Args:
        pressure: Absolute pressure field
        pp: Pressure correction field
        temp: Temperature field
        istat, iend, jstat, jend, kstat, nkm1: Loop bounds
        urf_p: Pressure under-relaxation factor
        tsolid: Solidus temperature [K]
    """
    # Fortran loops: k=kstat,nkm1; j=jstat,jend; i=istatp1,iendm1
    for k in range(kstat, nkm1 + 1):
        for j in range(jstat, jend + 1):
            for i in range(istat + 1, iend):
                if temp[i, j, k] > tsolid:
                    pressure[i, j, k] += urf_p * pp[i, j, k]
                    pp[i, j, k] = 0.0


# ============================================================================
# Utility Functions for Testing
# ============================================================================

def compute_velocity_divergence(state: State, grid: GridParams,
                                istat: int, iend: int, 
                                jstat: int, jend: int,
                                kstat: int, kend: int) -> float:
    """
    Compute maximum velocity divergence in the solution domain.
    
    Used to check mass conservation after pressure correction.
    
    Args:
        state: Current state with velocity fields
        grid: Grid parameters
        istat, iend, jstat, jend, kstat, kend: Solution domain bounds
        
    Returns:
        Maximum absolute divergence [1/s]
    """
    uVel_np = state.uVel.to_numpy()
    vVel_np = state.vVel.to_numpy()
    wVel_np = state.wVel.to_numpy()
    
    dx_np = grid.dx.to_numpy()
    dy_np = grid.dy.to_numpy()
    dz_np = grid.dz.to_numpy()
    
    max_div = 0.0
    
    for k in range(kstat, min(kend, grid.nk - 1)):
        for j in range(jstat, min(jend, grid.nj - 1)):
            for i in range(istat, min(iend, grid.ni - 1)):
                # Compute divergence: ∂u/∂x + ∂v/∂y + ∂w/∂z
                du_dx = (uVel_np[i + 1, j, k] - uVel_np[i, j, k]) / dx_np[i]
                dv_dy = (vVel_np[i, j + 1, k] - vVel_np[i, j, k]) / dy_np[j]
                dw_dz = (wVel_np[i, j, k + 1] - wVel_np[i, j, k]) / dz_np[k]
                
                div = abs(du_dx + dv_dy + dw_dz)
                max_div = max(max_div, div)
    
    return max_div
