"""
AM-CFD Taichi Implementation - Melt Pool Dimensions

Translated from Fortran mod_dimen.f90
Version 1.0 - Python/Taichi port

This module computes melt pool dimensions (length, width, depth) and defines
the solution domain for momentum equations to optimize computation.
"""

import taichi as ti
import numpy as np
from typing import Tuple

from data_structures import (
    State, GridParams, PhysicsParams, PoolDimensions
)


@ti.kernel
def find_max_temp(temp: ti.template(), ni: int, nj: int, nk: int) -> ti.f64:
    """
    Find maximum temperature in the domain.
    
    Args:
        temp: Temperature field [ni, nj, nk]
        ni, nj, nk: Grid dimensions
        
    Returns:
        Maximum temperature [K]
    """
    tmax = 0.0
    for i, j, k in ti.ndrange(ni, nj, nk):
        ti.atomic_max(tmax, temp[i, j, k])
    return tmax


def pool_size(state: State, grid: GridParams, physics: PhysicsParams,
              pool: PoolDimensions) -> Tuple[int, int, int, int, int, int]:
    """
    Compute melt pool dimensions and define momentum solution domain.
    
    This function:
    1. Finds peak temperature (tpeak) in the domain
    2. Computes melt pool length, depth, and width based on solidus isotherm
    3. Defines index bounds for solving momentum equations (with buffer zone)
    
    Corresponds to subroutine pool_size in mod_dimen.f90
    
    Args:
        state: Current state with temperature field
        grid: Grid parameters
        physics: Physical parameters (tsolid)
        pool: Pool dimensions structure (updated in-place)
        
    Returns:
        Tuple of (istat, iend, jstat, jend, kstat, kend):
            - istat, iend: Start/end indices in x-direction for momentum solve
            - jstat, jend: Start/end indices in y-direction for momentum solve
            - kstat, kend: Start/end indices in z-direction for momentum solve
    """
    ni, nj, nk = grid.ni, grid.nj, grid.nk
    
    # Convert Taichi fields to numpy for CPU processing
    # (Could be optimized with Taichi kernels for GPU, but keeping structure clear)
    temp_np = state.temp.to_numpy()
    x_np = grid.x.to_numpy()
    y_np = grid.y.to_numpy()
    z_np = grid.z.to_numpy()
    
    # Find peak temperature
    tpeak = np.max(temp_np)
    pool.max_temp = tpeak
    
    # If peak temperature is below solidus, no melt pool exists
    if tpeak <= physics.tsolid:
        pool.length = 0.0
        pool.depth = 0.0
        pool.width = 0.0
        # Return minimal domain (skip momentum solve)
        return (1, 1, 1, 1, 1, 1)
    
    # Grid parameters (Fortran uses 1-based indexing)
    # In Python: istart=0 corresponds to Fortran istart
    # Surface is at k=nk-1 (top), symmetry plane at j=0 (jstart)
    istart = 0
    jstart = 0
    nim1 = ni - 1
    njm1 = nj - 1
    nkm1 = nk - 1
    
    # =========================================================================
    # MELT POOL LENGTH (x-direction)
    # =========================================================================
    # Search along surface centerline: temp[i, jstart, nk-1]
    
    imax = istart
    imin = istart
    alen = 0.0
    
    # Find maximum i where temp > tsolid (forward from center)
    for i in range(istart, nim1):
        imax = i
        if temp_np[i, jstart, nkm1] <= physics.tsolid:
            break
    
    # Interpolate to find exact position where temp = tsolid
    if imax < nim1 - 1:
        temp_diff = temp_np[imax, jstart, nkm1] - temp_np[imax + 1, jstart, nkm1]
        if abs(temp_diff) > 1.0e-10:
            dtdxxinv = (x_np[imax] - x_np[imax + 1]) / temp_diff
            xxmax = x_np[imax] + (physics.tsolid - temp_np[imax, jstart, nkm1]) * dtdxxinv
        else:
            xxmax = x_np[imax]
    else:
        xxmax = x_np[imax]
    
    # Find minimum i where temp > tsolid (backward from center)
    for i in range(istart, 1, -1):
        imin = i
        if temp_np[i, jstart, nkm1] < physics.tsolid:
            break
    
    # Interpolate to find exact position
    if imin > 1:
        temp_diff = temp_np[imin, jstart, nkm1] - temp_np[imin - 1, jstart, nkm1]
        if abs(temp_diff) > 1.0e-10:
            dtdxxinv = (x_np[imin] - x_np[imin - 1]) / temp_diff
            xxmin = x_np[imin] + (physics.tsolid - temp_np[imin, jstart, nkm1]) * dtdxxinv
        else:
            xxmin = x_np[imin]
    else:
        xxmin = x_np[imin]
    
    alen = xxmax - xxmin
    pool.length = alen
    
    # =========================================================================
    # MELT POOL DEPTH (z-direction)
    # =========================================================================
    # Search downward from surface to find deepest penetration
    
    kmin = nkm1
    depth = 0.0
    
    # Find minimum k index where liquid exists
    for i in range(1, nim1):
        for k in range(nkm1, 1, -1):
            if temp_np[i, jstart, k] < physics.tsolid:
                break
            kmin = min(kmin, k)
    
    kmin = kmin - 1
    
    # If melt reaches bottom, depth is full domain height
    if kmin == 1:
        depth = z_np[nkm1] - z_np[0]
    else:
        # Interpolate to find exact depth at each x-location
        for i in range(1, nim1):
            if temp_np[i, jstart, kmin + 1] >= physics.tsolid:
                temp_diff = temp_np[i, jstart, kmin] - temp_np[i, jstart, kmin - 1]
                if abs(temp_diff) > 1.0e-10:
                    dtdzzinv = (z_np[kmin] - z_np[kmin - 1]) / temp_diff
                    dep = z_np[nkm1] - z_np[kmin] + \
                          (temp_np[i, jstart, kmin] - physics.tsolid) * dtdzzinv
                    depth = max(dep, depth)
    
    pool.depth = depth
    kmax = nkm1
    
    # =========================================================================
    # MELT POOL WIDTH (y-direction)
    # =========================================================================
    # Search laterally from symmetry plane
    
    jmax = jstart
    jmin = jstart
    width = 0.0
    yymax = y_np[jstart]
    yymin = y_np[jstart]
    
    # Find maximum j where liquid exists
    for i in range(1, nim1):
        for j in range(jstart, njm1):
            if temp_np[i, j, nkm1] >= physics.tsolid:
                jmax = max(jmax, j)
    
    jmax = jmax + 1
    
    # Interpolate to find exact y-position at maximum extent
    if jmax > jstart:
        for i in range(1, nim1):
            if temp_np[i, jmax - 1, nkm1] >= physics.tsolid:
                if jmax < njm1:
                    temp_diff = temp_np[i, jmax, nkm1] - temp_np[i, jmax + 1, nkm1]
                    if abs(temp_diff) > 1.0e-10:
                        dtdyyinv = (y_np[jmax] - y_np[jmax + 1]) / temp_diff
                        wid = y_np[jmax] + (physics.tsolid - temp_np[i, jmax, nkm1]) * dtdyyinv
                        yymax = max(wid, yymax)
    
    # Find minimum j where liquid exists (backward from symmetry plane)
    for i in range(1, nim1):
        for j in range(jstart, 1, -1):
            if temp_np[i, j, nkm1] >= physics.tsolid:
                jmin = min(jmin, j)
    
    jmin = jmin - 1
    
    # Interpolate to find exact y-position at minimum extent
    if jmin > jstart:
        for i in range(1, nim1):
            if temp_np[i, jmin + 1, nkm1] >= physics.tsolid:
                if jmin > 1:
                    temp_diff = temp_np[i, jmin, nkm1] - temp_np[i, jmin - 1, nkm1]
                    if abs(temp_diff) > 1.0e-10:
                        dtdyyinv = (y_np[jmin] - y_np[jmin - 1]) / temp_diff
                        wid = y_np[jmin] + (physics.tsolid - temp_np[i, jmin, nkm1]) * dtdyyinv
                        yymin = min(wid, yymin)
    
    width = yymax - yymin
    pool.width = width
    
    # =========================================================================
    # DEFINE MOMENTUM SOLUTION DOMAIN
    # =========================================================================
    # Add buffer cells around liquid region to define where to solve momentum
    
    istat = max(imin - 3, 1)
    iend = min(imax + 3, nim1)
    
    jstat = max(jmin - 3, 1)
    jend = min(jmax + 2, njm1)
    
    kstat = max(kmin - 2, 2)
    kend = nkm1
    
    # Store additional indices for convenience (if needed elsewhere)
    istatp1 = istat + 1
    iendm1 = iend - 1
    
    return (istat, iend, jstat, jend, kstat, kend)


def clean_uvw(state: State, grid: GridParams, physics: PhysicsParams) -> State:
    """
    Zero velocity outside melt pool (in solid region).
    
    Sets velocity to zero for cells where temp <= tsolid.
    This ensures that solid material has no velocity.
    
    Args:
        state: Current state with temperature and velocity fields
        grid: Grid parameters
        physics: Physical parameters (tsolid)
        
    Returns:
        state: Updated state with velocities zeroed in solid region
    """
    clean_uvw_kernel(state.temp, state.uVel, state.vVel, state.wVel,
                     grid.ni, grid.nj, grid.nk, physics.tsolid)
    return state


@ti.kernel
def clean_uvw_kernel(temp: ti.template(), 
                     uVel: ti.template(), 
                     vVel: ti.template(), 
                     wVel: ti.template(),
                     ni: int, nj: int, nk: int, tsolid: ti.f64):
    """
    Taichi kernel to zero velocities in solid regions.
    
    Args:
        temp: Temperature field [ni, nj, nk]
        uVel, vVel, wVel: Velocity fields [ni, nj, nk]
        ni, nj, nk: Grid dimensions
        tsolid: Solidus temperature [K]
    """
    for i, j, k in ti.ndrange(ni, nj, nk):
        if temp[i, j, k] <= tsolid:
            uVel[i, j, k] = 0.0
            vVel[i, j, k] = 0.0
            wVel[i, j, k] = 0.0


# ============================================================================
# Utility Functions
# ============================================================================

def compute_pool_volume(state: State, grid: GridParams, 
                        physics: PhysicsParams) -> float:
    """
    Compute melt pool volume by integrating over liquid cells.
    
    Args:
        state: Current state with temperature field
        grid: Grid parameters (cell volumes)
        physics: Physical parameters (tsolid)
        
    Returns:
        Melt pool volume [m³]
    """
    temp_np = state.temp.to_numpy()
    vol_np = grid.vol.to_numpy()
    
    volume = 0.0
    for i in range(grid.ni):
        for j in range(grid.nj):
            for k in range(grid.nk):
                if temp_np[i, j, k] > physics.tsolid:
                    volume += vol_np[i, j, k]
    
    return volume


def find_max_velocity(state: State, grid: GridParams) -> Tuple[float, float, float]:
    """
    Find maximum velocity magnitudes in each direction.
    
    Args:
        state: Current state with velocity fields
        grid: Grid parameters
        
    Returns:
        Tuple of (umax, vmax, wmax): Maximum velocity magnitudes [m/s]
    """
    uVel_np = state.uVel.to_numpy()
    vVel_np = state.vVel.to_numpy()
    wVel_np = state.wVel.to_numpy()
    
    umax = np.max(np.abs(uVel_np))
    vmax = np.max(np.abs(vVel_np))
    wmax = np.max(np.abs(wVel_np))
    
    return (umax, vmax, wmax)
