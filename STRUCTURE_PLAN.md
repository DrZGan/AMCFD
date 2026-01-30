# AM-CFD Structure Plan

## Proposed File Organization

```
jax/                          # JAX implementation
├── data_structures.py   # NamedTuples/dataclasses for state containers (from mod_const.f90)
├── param.py             # Input parsing, output routines (from mod_param.f90, mod_print.f90)
├── geom.py              # Non-uniform staggered grid generation (from mod_geom.f90)
├── init.py              # State initialization (from mod_init.f90)
├── bound.py             # Boundary condition handlers (from mod_bound.f90)
├── prop.py              # Material properties (from mod_prop.f90)
├── entot.py             # H↔T conversion (from mod_entot.f90)
├── discret.py           # FVM coefficients (from mod_discret.f90)
├── sour.py              # Source terms (from mod_sour.f90)
├── solve.py             # TDMA solver (from mod_solve.f90)
├── converge.py          # Convergence checks (from mod_converge.f90)
├── laser.py             # Laser heat source (from mod_laser.f90)
├── toolpath.py          # Toolpath loading (from mod_toolpath.f90)
├── dimen.py             # Melt pool dimensions (from mod_dimen.f90)
├── main.py              # Main solver orchestration (from main.f90)
└── inputfile/           # Input files (input_param.txt, *.crs)
```

## Data Flow Diagram

```
┌──────────────────────────┐
│  param.parse_input       │ → SimulationParams, PhysicsParams
└────────┬─────────────────┘
         ↓
┌──────────────────┐
│  geom.generate   │ → GridParams (x, y, z, xu, yv, zw, vol, areas)
└────────┬─────────┘
         ↓
┌──────────────────┐
│ init.create      │ → State (uVel, vVel, wVel, p, enthalpy, temp)
└────────┬─────────┘
         ↓
┌──────────────────────┐
│ toolpath.load        │ → ToolPath (time, x, y, z, laser_on)
└────────┬─────────────┘
         ↓
┌──────────────────────────────────────────────────────────────┐
│                      TIME LOOP (timet < timax)                │
│  ┌─────────────────────────────────────────────────────────┐ │
│  │              ITERATION LOOP (until converged)            │ │
│  │  ┌───────────┐ ┌────────────┐ ┌────────────┐ ┌─────────┐│ │
│  │  │  laser    │→│   bound    │→│  discret   │→│  solve  ││ │
│  │  └───────────┘ └────────────┘ └────────────┘ └────┬────┘│ │
│  │                                                    ↓     │ │
│  │  ┌───────────┐ ┌────────────┐ ┌────────────┐ ┌─────────┐│ │
│  │  │  dimen    │←│   prop     │←│   entot    │←│  state  ││ │
│  │  └───────────┘ └────────────┘ └────────────┘ └─────────┘│ │
│  │                       ↓                                  │ │
│  │              converge.check() → continue/break           │ │
│  └─────────────────────────────────────────────────────────┘ │
│                              ↓                                │
│                    param.write_output()                       │
└──────────────────────────────────────────────────────────────┘
```

---

## State Container Definitions (data_structures.py)

```python
from typing import NamedTuple
import jax.numpy as jnp  # or ti.types for Taichi

class State(NamedTuple):
    """Primary flow field variables (updated each timestep)"""
    uVel: jnp.ndarray      # x-velocity [ni, nj, nk]
    vVel: jnp.ndarray      # y-velocity [ni, nj, nk]
    wVel: jnp.ndarray      # z-velocity [ni, nj, nk]
    pressure: jnp.ndarray  # Pressure [ni, nj, nk]
    pp: jnp.ndarray        # Pressure correction [ni, nj, nk]
    enthalpy: jnp.ndarray  # Enthalpy [ni, nj, nk]
    temp: jnp.ndarray      # Temperature [ni, nj, nk]
    fracl: jnp.ndarray     # Liquid fraction [ni, nj, nk]

class StatePrev(NamedTuple):
    """Previous timestep values for transient terms"""
    unot: jnp.ndarray
    vnot: jnp.ndarray
    wnot: jnp.ndarray
    hnot: jnp.ndarray

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
    dgdt: float       # Surface tension gradient (Marangoni)
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
```

---

## Variable Lifetime Categories

| Category | Examples | Stored In | Updated | Scope |
|----------|----------|-----------|---------|-------|
| **Persistent** | `uVel, enthalpy, temp` | `State` | Every timestep | Global |
| **Grid** | `x, y, z, vol, areas` | `GridParams` | Never | Global (immutable) |
| **Physics** | `acpa, tsolid, dgdt` | `PhysicsParams` | Never | Global (immutable) |
| **Simulation** | `delt, urf_vel` | `SimulationParams` | Never | Global (immutable) |
| **Time** | `timet, iter, step` | `TimeState` | Every timestep | Global |
| **Properties** | `vis, diff, den` | `MaterialProps` | Each iteration | Global |
| **Coefficients** | `ap, ae, su, sp` | `DiscretCoeffs` | Each solve | Transient |
| **Previous** | `unot, hnot` | `StatePrev` | Start of timestep | Global |
| **Laser** | `beam_pos, heatin` | `LaserState` | Each timestep | Global |
| **Derived** | `fracl` | Computed locally | As needed | Local |

### Transient Variables (computed and discarded)

These exist only within function scope:

| Variable | Shape | Description | Computed In |
|----------|-------|-------------|-------------|
| `ap, ae, aw, ...` | `[ni, nj, nk]` | FVM coefficients | `discretization.compute` |
| `su, sp` | `[ni, nj, nk]` | Source terms | `source.compute` |
| `heatin` | `[ni, nj]` | Laser heat flux | `laser.compute_heat` |
| `residual` | scalar | Convergence metric | `mod_converge.compute_residual` |
| `pool_length/depth/width` | scalars | Melt pool dimensions | `mod_dimen.compute_size` |

---

## Module Responsibilities

| Module | Fortran Source | Inputs | Outputs | Notes |
|--------|----------------|--------|---------|-------|
| `data_structures.py` | `mod_const` | - | NamedTuples | Data structure definitions |
| `param.py` | `mod_param`, `mod_print` | input_param.txt | Params, output files | Parse namelist format |
| `geom.py` | `mod_geom` | Params | `GridParams` | Power-law stretching |
| `init.py` | `mod_init` | Params, Grid | `State` | Initialize to preheat T |
| `prop.py`, `entot.py` | `mod_prop`, `mod_entot` | State, Params | `MaterialProps`, T↔H | 3-region phase change |
| `laser.py`, `toolpath.py` | `mod_laser`, `mod_toolpath` | .crs file, time | `LaserState` | Gaussian heat source |
| `bound.py` | `mod_bound` | State, ivar | Updated su/sp | Marangoni, radiation |
| `discret.py` | `mod_discret` | State, Props | `DiscretCoeffs` | Power-law scheme |
| `sour.py` | `mod_sour` | State, Props | Updated su/sp | Darcy, buoyancy |
| `solve.py` | `mod_solve` | Coeffs, State | New field | TDMA line-by-line |
| `converge.py`, `resid.py` | `mod_converge`, `mod_resid` | State, Coeffs | Residuals, ratio | Heat balance check |
| `dimen.py` | `mod_dimen` | State, Grid | Pool dimensions | Solidus interpolation |
| `main.py`, `revise.py` | `main.f90`, `mod_revise` | All | Final State | Time loop, pressure correction |

---

## JIT Compilation Strategy

```python
# JAX: Compile entire iteration as single unit
@jax.jit
def iteration_step(state: State, state_prev: StatePrev,
                   grid: GridParams, physics: PhysicsParams,
                   laser: LaserState, sim: SimulationParams) -> State:
    # 1. Compute laser heat input
    heatin = laser_heat(laser, grid)
    # 2. Apply boundary conditions (updates su, sp)
    su, sp = apply_boundary(state, grid, physics, heatin)
    # 3. Compute discretization coefficients
    coeffs = compute_discretization(state, state_prev, grid, sim)
    # 4. Solve energy equation (TDMA)
    new_enthalpy = solve_tdma(coeffs, state.enthalpy)
    # 5. Convert H → T, compute fracl
    new_temp, new_fracl = enthalpy_to_temp(new_enthalpy, physics)
    # 6. Update properties
    props = compute_properties(new_temp, new_fracl, physics)
    # 7. Solve momentum if melted (optional)
    # ...
    return state._replace(enthalpy=new_enthalpy, temp=new_temp, fracl=new_fracl)

# For time loop, use lax.while_loop or lax.scan
def time_loop(state, grid, physics, laser, sim):
    def cond_fn(carry):
        state, time_state = carry
        return time_state.timet < sim.timax
    
    def body_fn(carry):
        state, time_state = carry
        # Save previous state
        state_prev = StatePrev(unot=state.uVel, ...)
        # Inner iteration loop
        state = converge_loop(state, state_prev, grid, physics, laser, sim)
        # Update time
        new_time = time_state._replace(timet=time_state.timet + sim.delt)
        return (state, new_time)
    
    return jax.lax.while_loop(cond_fn, body_fn, (state, time_state))
```

---

## Immutable Update Pattern

```python
# ❌ Wrong: in-place modification (Fortran style)
state.enthalpy[i,j,k] = new_value

# ✅ Correct: create new state with _replace()
new_state = state._replace(enthalpy=new_enthalpy)

# Multiple fields at once
new_state = state._replace(
    enthalpy=new_enthalpy,
    temp=new_temp,
    fracl=new_fracl
)
```

---

## Module Dependency Graph

```
                         ┌──────────┐
                         │ types.py │
                         └────┬─────┘
                              │ (imported by all)
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
        ▼                     ▼                     ▼
   ┌─────────┐          ┌──────────┐          ┌──────────┐
   │  io.py  │          │ grid.py  │          │ laser.py │
   └────┬────┘          └────┬─────┘          └────┬─────┘
        │                    │                     │
        └──────────┬─────────┴─────────────────────┘
                   ▼
            ┌─────────────┐
            │ initial.py  │
            └──────┬──────┘
                   │
                   ▼
          ┌──────────────┐
          │ properties.py│ (H↔T conversion)
          └──────┬───────┘
                 │
        ┌────────┴────────┐
        ▼                 ▼
 ┌─────────────┐   ┌───────────┐
 │ boundary.py │   │ source.py │
 └──────┬──────┘   └─────┬─────┘
        │                │
        └───────┬────────┘
                ▼
       ┌────────────────┐
       │discretization.py│
       └───────┬────────┘
               │
               ▼
         ┌───────────┐
         │ solver.py │ (TDMA)
         └─────┬─────┘
               │
       ┌───────┴───────┐
       ▼               ▼
 ┌───────────┐   ┌──────────────┐
 │  pool.py  │   │convergence.py│
 └─────┬─────┘   └──────┬───────┘
       │                │
       └───────┬────────┘
               ▼
         ┌───────────┐
         │  main.py  │
         └───────────┘
```

### Call Flow (single iteration)

```
main.iteration()
    │
    ├──► laser.compute_heat(laser_state, grid) → heatin
    │
    ├──► boundary.apply(state, grid, physics, heatin, ivar=5) → su, sp (enthalpy)
    │
    ├──► discretization.compute(state, state_prev, grid, sim) → coeffs
    │
    ├──► solver.tdma_enthalpy(coeffs, state.enthalpy) → new_enthalpy
    │
    ├──► properties.enthalpy_to_temp(new_enthalpy, physics) → new_temp, new_fracl
    │
    ├──► pool.compute_size(new_temp, grid, physics) → length, depth, width
    │
    ├──► [if melted] momentum solve loop:
    │         │
    │         ├──► boundary.apply(..., ivar=1,2,3) → velocity BCs
    │         ├──► discretization.compute_velocity(...)
    │         ├──► solver.tdma_uvw(...) → uVel, vVel, wVel
    │         └──► main.pressure_correction(...) → pp, corrected velocities
    │
    └──► convergence.check(state, coeffs) → converged, residual, ratio
              │
              ▼
          new_state
```

### Module Import Summary

| Module | Imports From |
|--------|--------------|
| `data_structures.py` | (none) |
| `param.py` | `data_structures` |
| `geom.py` | `data_structures` |
| `init.py` | `data_structures` |
| `laser.py` | `data_structures` |
| `toolpath.py` | `data_structures` |
| `prop.py` | `data_structures` |
| `entot.py` | `data_structures` |
| `bound.py` | `data_structures`, `prop` |
| `sour.py` | `data_structures`, `prop` |
| `discret.py` | `data_structures` |
| `solve.py` | `data_structures` |
| `dimen.py` | `data_structures` |
| `converge.py` | `data_structures` |
| `main.py` | `data_structures`, `param`, `geom`, `init`, `laser`, `toolpath`, `prop`, `entot`, `bound`, `sour`, `discret`, `solve`, `dimen`, `converge` |

---

## Fortran Module → JAX/Taichi Mapping

| Fortran Module | JAX/Taichi Module | Key Functions |
|----------------|-------------------|---------------|
| `mod_const.f90` | `data_structures.py` | Physical constants in PhysicsParams |
| `mod_param.f90` | `param.py` | `parse_input()` |
| `mod_geom.f90` | `geom.py` | `generate_grid()` |
| `mod_init.f90` | `init.py` | `initialize_state()` |
| `mod_prop.f90` | `prop.py` | `compute_properties()` |
| `mod_entot.f90` | `entot.py` | `enthalpy_to_temp()`, `temp_to_enthalpy()` |
| `mod_laser.f90` | `laser.py` | `compute_heat()` |
| `mod_toolpath.f90` | `toolpath.py` | `load_toolpath()`, `update_position()` |
| `mod_bound.f90` | `bound.py` | `apply_bc()` per variable |
| `mod_discret.f90` | `discret.py` | `compute_coefficients()` |
| `mod_sour.f90` | `sour.py` | `compute_source()` |
| `mod_solve.f90` | `solve.py` | `tdma_solve()` |
| `mod_dimen.f90` | `dimen.py` | `compute_pool_size()` |
| `mod_converge.f90` | `converge.py` | `check_convergence()` |
| `mod_resid.f90` | `resid.py` | `compute_residual()` |
| `mod_revise.f90` | `revise.py` | `pressure_correction()` |
| `mod_flux.f90` | `flux.py` | `compute_heat_balance()` |
| `mod_print.f90` | `param.py` | `write_output()`, `write_tecplot()` |
| `main.f90` | `main.py` | `run_simulation()`, time/iteration loops |
