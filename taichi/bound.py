"""
AM-CFD Taichi Implementation - Boundary Conditions

Converted from Fortran module: mod_bound.f90

This module applies boundary conditions for velocity and enthalpy fields.
The boundary conditions include:
- Marangoni stress on the top surface (velocity)
- No-slip conditions at solid boundaries
- Thermal boundary conditions (radiation, convection) at all surfaces
"""

import taichi as ti
from data_structures import State, GridParams, MaterialProps, PhysicsParams, SimulationParams, LaserState


@ti.kernel
def apply_boundary_conditions_u(
    state: ti.template(),
    grid: ti.template(),
    mat_props: ti.template(),
    physics: ti.template(),
    simu_params: ti.template(),
):
    """Apply boundary conditions for u-velocity (ivar=1).
    
    - Top surface: Marangoni stress (temperature gradient driven)
    - Solid boundaries: No-slip (u=0)
    """
    ni = ti.static(simu_params.ni)
    nj = ti.static(simu_params.nj)
    nk = ti.static(simu_params.nk)
    istat = ti.static(simu_params.istat)
    iend = ti.static(simu_params.iend)
    jstat = ti.static(simu_params.jstat)
    jend = ti.static(simu_params.jendm1)
    istatp1 = ti.static(simu_params.istatp1)
    iendm1 = ti.static(simu_params.iendm1)
    nkm1 = ti.static(simu_params.nkm1)

    # Top surface (k=nk): Marangoni stress
    for j in range(jstat, jend + 1):
        for i in range(istatp1, iendm1 + 1):
            # Temperature gradient in x-direction
            dtdx = (state.temp[i, j, nk - 1] - state.temp[i - 1, j, nk - 1]) * grid.dxpwinv[i]
            
            # Interpolate liquid fraction to u-face
            # fracx[i-1] is the interpolation factor
            fracx = (grid.x[i] - grid.xu[i]) / (grid.x[i] - grid.x[i - 1])
            fraclu = state.fracl[i, j, nk - 1] * (1.0 - fracx) + state.fracl[i - 1, j, nk - 1] * fracx
            
            # Harmonic mean of viscosity at u-face
            visu1 = mat_props.vis[i, j, nkm1 - 1] * mat_props.vis[i - 1, j, nkm1 - 1] / (
                mat_props.vis[i - 1, j, nkm1 - 1] * (1.0 - fracx) + mat_props.vis[i, j, nkm1 - 1] * fracx
            )
            
            # Marangoni term: fraclu * dgdt * dtdx / (visu1 * dzpbinv[nk])
            term1 = fraclu * physics.dgdt * dtdx / (visu1 * grid.dzpbinv[nk - 1])
            
            # Apply: u[nk] = u[nkm1] + term1
            state.uVel[i, j, nk - 1] = state.uVel[i, j, nkm1 - 1] + term1
    
    # Solid boundaries: u=0
    for j in range(jstat, jend + 1):
        for k in range(0, nkm1):
            state.uVel[istat, j, k] = 0.0  # Left boundary
            state.uVel[iend - 1, j, k] = 0.0  # Right boundary


@ti.kernel
def apply_boundary_conditions_v(
    state: ti.template(),
    grid: ti.template(),
    mat_props: ti.template(),
    physics: ti.template(),
    simu_params: ti.template(),
):
    """Apply boundary conditions for v-velocity (ivar=2).
    
    - Top surface: Marangoni stress (temperature gradient driven)
    - Solid boundaries: No-slip (v=0)
    """
    ni = ti.static(simu_params.ni)
    nj = ti.static(simu_params.nj)
    nk = ti.static(simu_params.nk)
    istat = ti.static(simu_params.istat)
    iend = ti.static(simu_params.iend)
    jstat = ti.static(simu_params.jstat)
    jend = ti.static(simu_params.jendm1)
    istatp1 = ti.static(simu_params.istatp1)
    iendm1 = ti.static(simu_params.iendm1)
    nkm1 = ti.static(simu_params.nkm1)

    # Top surface (k=nk): Marangoni stress
    for j in range(jstat, jend + 1):
        for i in range(istatp1, iendm1 + 1):
            # Temperature gradient in y-direction
            dtdy = (state.temp[i, j, nk - 1] - state.temp[i, j - 1, nk - 1]) * grid.dypsinv[j]
            
            # Interpolate liquid fraction to v-face
            fracy = (grid.y[j] - grid.yv[j]) / (grid.y[j] - grid.y[j - 1])
            fraclv = state.fracl[i, j, nk - 1] * (1.0 - fracy) + state.fracl[i, j - 1, nk - 1] * fracy
            
            # Harmonic mean of viscosity at v-face
            visv1 = mat_props.vis[i, j, nkm1 - 1] * mat_props.vis[i, j - 1, nkm1 - 1] / (
                mat_props.vis[i, j - 1, nkm1 - 1] * (1.0 - fracy) + mat_props.vis[i, j, nkm1 - 1] * fracy
            )
            
            # Marangoni term
            term1 = fraclv * physics.dgdt * dtdy / (visv1 * grid.dzpbinv[nk - 1])
            
            # Apply: v[nk] = v[nkm1] + term1
            state.vVel[i, j, nk - 1] = state.vVel[i, j, nkm1 - 1] + term1
    
    # Solid boundaries: v=0
    for j in range(jstat, jend + 1):
        for k in range(0, nkm1):
            state.vVel[istat, j, k] = 0.0  # Left boundary
            state.vVel[iend - 1, j, k] = 0.0  # Right boundary


@ti.kernel
def apply_boundary_conditions_w(
    state: ti.template(),
    simu_params: ti.template(),
):
    """Apply boundary conditions for w-velocity (ivar=3).
    
    - Solid boundaries: No-slip (w=0)
    """
    ni = ti.static(simu_params.ni)
    nj = ti.static(simu_params.nj)
    nk = ti.static(simu_params.nk)
    istat = ti.static(simu_params.istat)
    iend = ti.static(simu_params.iend)
    jstat = ti.static(simu_params.jstat)
    jend = ti.static(simu_params.jendm1)
    nkm1 = ti.static(simu_params.nkm1)

    # Solid boundaries: w=0
    for j in range(jstat, jend + 1):
        for k in range(0, nkm1):
            state.wVel[istat, j, k] = 0.0  # Left boundary
            state.wVel[iend - 1, j, k] = 0.0  # Right boundary


@ti.kernel
def apply_boundary_conditions_p(
    state: ti.template(),
    simu_params: ti.template(),
):
    """Apply boundary conditions for pressure correction (ivar=4).
    
    - Solid boundaries: pp=0
    """
    ni = ti.static(simu_params.ni)
    nj = ti.static(simu_params.nj)
    nk = ti.static(simu_params.nk)
    istat = ti.static(simu_params.istat)
    iend = ti.static(simu_params.iend)
    jstat = ti.static(simu_params.jstat)
    jend = ti.static(simu_params.jendm1)
    nkm1 = ti.static(simu_params.nkm1)

    # Solid boundaries: pp=0
    for j in range(jstat, jend + 1):
        for k in range(0, nkm1):
            state.pp[istat, j, k] = 0.0  # Left boundary
            state.pp[iend - 1, j, k] = 0.0  # Right boundary


@ti.kernel
def apply_boundary_conditions_enthalpy(
    state: ti.template(),
    grid: ti.template(),
    mat_props: ti.template(),
    physics: ti.template(),
    laser: ti.template(),
    simu_params: ti.template(),
) -> ti.f64:
    """Apply boundary conditions for enthalpy (ivar=5).
    
    - Top surface: Laser heating, radiation, and convection
    - Bottom surface: Convection and radiation
    - Side surfaces: Convection
    
    Returns:
        ahtoploss: Total heat loss at top surface [W]
    """
    ni = ti.static(simu_params.ni)
    nj = ti.static(simu_params.nj)
    nk = ti.static(simu_params.nk)
    nim1 = ti.static(simu_params.nim1)
    njm1 = ti.static(simu_params.njm1)
    nkm1 = ti.static(simu_params.nkm1)

    ahtoploss = 0.0
    
    # Top surface (k=nk)
    for j in range(1, njm1):
        for i in range(1, nim1):
            # Radiation loss: emiss * sigma * (T^4 - T_amb^4)
            hlossradia = physics.emiss * physics.sigma * (
                state.temp[i, j, nk - 1] ** 4 - physics.tenv ** 4
            )
            
            # Convection loss: htckn * (T - T_amb)
            hlossconvec = physics.hconv * (state.temp[i, j, nk - 1] - physics.tenv)
            
            # Thermal diffusion coefficient
            ctmp1 = mat_props.diff[i, j, nkm1 - 1] * grid.dzpbinv[nk - 1]
            
            # Enthalpy BC: H[nk] = H[nkm1] + (Q_in - Q_loss) / (k/dz)
            # Note: Original Fortran has hlossradia twice (bug?), keeping as is
            state.enthalpy[i, j, nk - 1] = (
                state.enthalpy[i, j, nkm1 - 1]
                + (laser.heatin[i, j] - hlossradia - hlossconvec) / ctmp1
            )
            
            # Accumulate total top heat loss
            ahtoploss += (hlossradia + hlossconvec) * grid.areaij[i, j]
    
    # Bottom surface (k=1)
    for j in range(1, njm1):
        for i in range(1, nim1):
            hlossconvec = physics.hconv * (state.temp[i, j, 0] - physics.tenv) + physics.emiss * physics.sigma * (
                state.temp[i, j, 0] ** 4 - physics.tenv ** 4
            )
            ctmp1 = mat_props.diff[i, j, 1] * grid.dzpbinv[1]
            state.enthalpy[i, j, 0] = state.enthalpy[i, j, 1] - hlossconvec / ctmp1
    
    # West and East boundaries
    for j in range(1, njm1):
        for k in range(1, nkm1):
            # West (i=1)
            hlossconvec = physics.hconv * (state.temp[0, j, k] - physics.tenv)
            ctmp1 = mat_props.diff[1, j, k] * grid.dxpwinv[1]
            state.enthalpy[0, j, k] = state.enthalpy[1, j, k] - hlossconvec / ctmp1
            
            # East (i=ni)
            hlossconvec = physics.hconv * (state.temp[ni - 1, j, k] - physics.tenv)
            ctmp1 = mat_props.diff[nim1 - 1, j, k] * grid.dxpwinv[nim1 - 1]
            state.enthalpy[ni - 1, j, k] = state.enthalpy[nim1 - 1, j, k] - hlossconvec / ctmp1
    
    # North and South boundaries
    for i in range(1, nim1):
        for k in range(1, nkm1):
            # South (j=1)
            hlossconvec = physics.hconv * (state.temp[i, 0, k] - physics.tenv)
            ctmp1 = mat_props.diff[i, 1, k] * grid.dypsinv[1]
            state.enthalpy[i, 0, k] = state.enthalpy[i, 1, k] - hlossconvec / ctmp1
            
            # North (j=nj)
            hlossconvec = physics.hconv * (state.temp[i, nj - 1, k] - physics.tenv)
            ctmp1 = mat_props.diff[i, njm1 - 1, k] * grid.dypsinv[njm1 - 1]
            state.enthalpy[i, nj - 1, k] = state.enthalpy[i, njm1 - 1, k] - hlossconvec / ctmp1
    
    return ahtoploss


def bound_condition(
    ivar: int,
    state: State,
    grid: GridParams,
    mat_props: MaterialProps,
    physics: PhysicsParams,
    simu_params: SimulationParams,
    laser: LaserState = None,
) -> float:
    """Apply boundary conditions based on variable type.
    
    This is the main entry point that dispatches to the appropriate
    boundary condition function based on ivar.
    
    Args:
        ivar: Variable index
              1 = u-velocity
              2 = v-velocity
              3 = w-velocity
              4 = pressure correction (pp)
              5 = enthalpy
        state: Flow field state
        grid: Grid parameters
        mat_props: Material properties
        physics: Physical parameters
        simu_params: Simulation parameters (contains index bounds)
        laser: Laser state (required for ivar=5)
    
    Returns:
        ahtoploss: Total heat loss at top surface [W] (only for ivar=5, else 0.0)
    """
    if ivar == 1:
        # u-velocity
        apply_boundary_conditions_u(
            state, grid, mat_props, physics, simu_params
        )
        return 0.0
    
    elif ivar == 2:
        # v-velocity
        apply_boundary_conditions_v(
            state, grid, mat_props, physics, simu_params
        )
        return 0.0
    
    elif ivar == 3:
        # w-velocity
        apply_boundary_conditions_w(
            state, simu_params
        )
        return 0.0
    
    elif ivar == 4:
        # pressure correction
        apply_boundary_conditions_p(
            state, simu_params
        )
        return 0.0
    
    elif ivar == 5:
        # enthalpy
        if laser is None:
            raise ValueError("LaserState required for enthalpy boundary conditions (ivar=5)")
        
        ahtoploss = apply_boundary_conditions_enthalpy(
            state, grid, mat_props, physics, laser, simu_params
        )
        return ahtoploss
    
    else:
        raise ValueError(f"Invalid ivar={ivar}. Must be 1-5.")
