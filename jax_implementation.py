"""
AM-CFD JAX Implementation Guide

This file contains:
1. State container definitions (NamedTuples)
2. File organization and module layout
3. Function patterns and JIT compilation strategy

Usage:
    from jax_implementation import FluidState, GridParams, ...
    
    state = FluidState(uVel=..., vVel=..., ...)
    new_state = state._replace(enthalpy=new_enthalpy)
"""

from typing import NamedTuple
import jax
import jax.numpy as jnp

# =============================================================================
# FILE ORGANIZATION
# =============================================================================
"""
Proposed structure:

amcfd_jax/
├── types.py          # NamedTuples (copy from this file)
├── grid.py           # Grid initialization
├── initial.py        # Initial conditions
├── boundary.py       # Boundary condition handlers
├── properties.py     # Material properties, H↔T conversion
├── discretization.py # FVM coefficients (ap, ae, aw, ...)
├── source.py         # Source terms (laser, Darcy, buoyancy)
├── solver.py         # TDMA solver
├── convergence.py    # Residuals, convergence checks
├── laser.py          # Laser heat source, toolpath
├── main.py           # Main solver orchestration
└── io.py             # Input parsing, output routines
"""

# =============================================================================
# STATE CONTAINER DEFINITIONS
# =============================================================================

class GridParams(NamedTuple):
    """Computational grid (immutable after initialization)"""
    x: jnp.ndarray       # Cell centers [ni]
    y: jnp.ndarray       # Cell centers [nj]
    z: jnp.ndarray       # Cell centers [nk]
    xu: jnp.ndarray      # u-velocity faces [ni]
    yv: jnp.ndarray      # v-velocity faces [nj]
    zw: jnp.ndarray      # w-velocity faces [nk]
    vol: jnp.ndarray     # Cell volumes [ni, nj, nk]
    areaij: jnp.ndarray  # xy-face areas [ni, nj]
    areaik: jnp.ndarray  # xz-face areas [ni, nk]
    areajk: jnp.ndarray  # yz-face areas [nj, nk]
    # Inverse distances for discretization
    dxpwinv: jnp.ndarray
    dxpeinv: jnp.ndarray
    dypsinv: jnp.ndarray
    dypninv: jnp.ndarray
    dzpbinv: jnp.ndarray
    dzptinv: jnp.ndarray
    ni: int
    nj: int
    nk: int


class FluidState(NamedTuple):
    """Primary flow field variables (updated each timestep)"""
    uVel: jnp.ndarray      # x-velocity [ni, nj, nk]
    vVel: jnp.ndarray      # y-velocity [ni, nj, nk]
    wVel: jnp.ndarray      # z-velocity [ni, nj, nk]
    pressure: jnp.ndarray  # Pressure [ni, nj, nk]
    pp: jnp.ndarray        # Pressure correction [ni, nj, nk]
    enthalpy: jnp.ndarray  # Enthalpy [ni, nj, nk]
    temp: jnp.ndarray      # Temperature [ni, nj, nk]
    fracl: jnp.ndarray     # Liquid fraction [ni, nj, nk]


class FluidStatePrev(NamedTuple):
    """Previous timestep values for transient terms"""
    unot: jnp.ndarray
    vnot: jnp.ndarray
    wnot: jnp.ndarray
    hnot: jnp.ndarray


class MaterialProps(NamedTuple):
    """Material properties (may vary spatially)"""
    vis: jnp.ndarray   # Viscosity [ni, nj, nk]
    diff: jnp.ndarray  # Diffusivity [ni, nj, nk]
    den: jnp.ndarray   # Density [ni, nj, nk]


class DiscretCoeffs(NamedTuple):
    """FVM discretization coefficients (transient, rebuilt each solve)"""
    ap: jnp.ndarray   # Center coefficient [ni, nj, nk]
    ae: jnp.ndarray   # East neighbor
    aw: jnp.ndarray   # West neighbor
    an: jnp.ndarray   # North neighbor
    as_: jnp.ndarray  # South neighbor (as_ to avoid Python keyword)
    at: jnp.ndarray   # Top neighbor
    ab: jnp.ndarray   # Bottom neighbor
    su: jnp.ndarray   # Source term
    sp: jnp.ndarray   # Linearized source coefficient


class LaserState(NamedTuple):
    """Laser and toolpath state"""
    beam_pos: float       # Current x position
    beam_posy: float      # Current y position
    heatin: jnp.ndarray   # Heat flux [ni, nj]
    laser_on: bool
    scanvelx: float
    scanvely: float


class PhysicsParams(NamedTuple):
    """Physical constants (immutable)"""
    acpa: float       # Solid Cp coefficient a
    acpb: float       # Solid Cp coefficient b
    acpl: float       # Liquid Cp
    tsolid: float     # Solidus temperature
    tliquid: float    # Liquidus temperature
    hsmelt: float     # Solidus enthalpy
    hlcal: float      # Liquidus enthalpy
    dgdt: float       # Surface tension gradient
    rho: float        # Reference density
    emiss: float      # Emissivity
    sigma: float      # Stefan-Boltzmann constant
    hconv: float      # Convection coefficient


class SimulationParams(NamedTuple):
    """Numerical parameters (immutable)"""
    delt: float       # Timestep
    timax: float      # Max simulation time
    urf_vel: float    # Under-relaxation for velocity
    urf_p: float      # Under-relaxation for pressure
    urf_h: float      # Under-relaxation for enthalpy


class TimeState(NamedTuple):
    """Time stepping information"""
    timet: float      # Current time
    iter: int         # Current iteration
    step: int         # Current timestep number


# =============================================================================
# VARIABLE LIFETIME CATEGORIES
# =============================================================================
"""
| Category              | Examples              | Stored In        | Updated           |
|-----------------------|-----------------------|------------------|-------------------|
| Persistent            | uVel, enthalpy, temp  | FluidState       | Every timestep    |
| Grid                  | x, y, z, vol          | GridParams       | Never (immutable) |
| Transient coefficients| ap, ae, su, sp        | DiscretCoeffs    | Each solver call  |
| Previous timestep     | unot, hnot            | FluidStatePrev   | Start of timestep |
| Derived on-demand     | fracl, vis            | Computed locally | As needed         |
"""

# =============================================================================
# FUNCTION PATTERNS
# =============================================================================

# Timestep function signature
@jax.jit
def time_step(state: FluidState, state_prev: FluidStatePrev,
              grid: GridParams, physics: PhysicsParams,
              laser: LaserState, sim: SimulationParams) -> FluidState:
    """
    Advance solution by one timestep.
    All intermediate variables are local and discarded after return.
    """
    # 1. Apply boundary conditions
    # 2. Compute discretization coefficients
    # 3. Solve energy equation
    # 4. Convert H → T
    # 5. Solve momentum if melted
    # 6. Pressure correction
    # Returns new FluidState
    pass  # Implementation in solver.py


@jax.jit
def solve_enthalpy(state: FluidState, coeffs: DiscretCoeffs,
                   grid: GridParams) -> jnp.ndarray:
    """TDMA solve for enthalpy field. Returns new enthalpy array."""
    pass  # Implementation in solver.py


# Updating state (immutable pattern)
def update_example(state: FluidState, new_enthalpy: jnp.ndarray, 
                   new_temp: jnp.ndarray) -> FluidState:
    """Example of immutable state update."""
    return state._replace(enthalpy=new_enthalpy, temp=new_temp)


# =============================================================================
# JIT COMPILATION STRATEGY
# =============================================================================

def run_simulation(state: FluidState, state_prev: FluidStatePrev,
                   grid: GridParams, physics: PhysicsParams,
                   laser: LaserState, sim: SimulationParams,
                   n_steps: int) -> FluidState:
    """
    Run multiple timesteps efficiently using lax.scan.
    This compiles the entire loop for GPU execution.
    """
    def body_fn(carry, _):
        state, state_prev, laser = carry
        # Save current state as previous
        new_prev = FluidStatePrev(
            unot=state.uVel,
            vnot=state.vVel,
            wnot=state.wVel,
            hnot=state.enthalpy
        )
        # Advance one step
        new_state = time_step(state, state_prev, grid, physics, laser, sim)
        # Update laser position (simplified)
        new_laser = laser._replace(
            beam_pos=laser.beam_pos + laser.scanvelx * sim.delt
        )
        return (new_state, new_prev, new_laser), None
    
    (final_state, _, _), _ = jax.lax.scan(
        body_fn, (state, state_prev, laser), None, length=n_steps
    )
    return final_state


# =============================================================================
# CONDITIONAL LOGIC IN JIT
# =============================================================================

@jax.jit
def enthalpy_to_temp(H: jnp.ndarray, physics: PhysicsParams) -> tuple:
    """
    Convert enthalpy to temperature (3 regions).
    Uses jnp.where for JIT-compatible branching.
    """
    # Solid region
    T_solid = (jnp.sqrt(physics.acpb**2 + 2*physics.acpa*H) - physics.acpb) / physics.acpa
    fracl_solid = 0.0
    
    # Mushy region
    deltemp = physics.tliquid - physics.tsolid
    fracl_mushy = (H - physics.hsmelt) / (physics.hlcal - physics.hsmelt)
    T_mushy = deltemp * fracl_mushy + physics.tsolid
    
    # Liquid region
    T_liquid = (H - physics.hlcal) / physics.acpl + physics.tliquid
    fracl_liquid = 1.0
    
    # Select based on enthalpy level
    is_solid = H < physics.hsmelt
    is_liquid = H > physics.hlcal
    
    T = jnp.where(is_solid, T_solid,
                  jnp.where(is_liquid, T_liquid, T_mushy))
    fracl = jnp.where(is_solid, fracl_solid,
                      jnp.where(is_liquid, fracl_liquid, fracl_mushy))
    
    return T, fracl


# =============================================================================
# DEBUGGING TIPS
# =============================================================================
"""
1. Disable JIT temporarily for better error messages:
   with jax.disable_jit():
       result = my_function(args)

2. Check array shapes:
   print(f"u shape: {state.uVel.shape}, expected: ({grid.ni}, {grid.nj}, {grid.nk})")

3. Use jax.debug.print inside JIT-compiled functions:
   @jax.jit
   def my_func(x):
       jax.debug.print("x = {}", x)
       return x * 2

4. For NaN debugging:
   jax.config.update("jax_debug_nans", True)
"""
