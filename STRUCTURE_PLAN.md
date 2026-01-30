# Structure Plan: 1D Euler Solver in JAX

## Proposed File Organization

```
euler_jax/
├── types.py          # NamedTuples for state/params
├── grid.py           # Grid initialization
├── initial.py        # Initial conditions (Sod, custom)
├── boundary.py       # Boundary condition handlers
├── primitives.py     # Conservative ↔ Primitive conversion
├── flux.py           # Numerical flux schemes (Roe, HLLC)
├── muscl.py          # MUSCL reconstruction + limiters
├── integrator.py     # Time stepping (Euler, RK)
├── solver.py         # Main solver orchestration
└── io.py             # Output routines
```

## Data Flow Diagram

```
┌─────────────┐
│ initialize  │ → GridParams, PhysicsParams, FluidState
└──────┬──────┘
       ↓
┌──────────────────────────────────────────────┐
│                 TIME LOOP                     │
│  ┌─────────┐   ┌───────┐   ┌────────┐        │
│  │   BC    │ → │ Flux  │ → │ Update │ → state│
│  └─────────┘   └───────┘   └────────┘        │
└──────────────────────────────────────────────┘
       ↓
┌─────────────┐
│   output    │
└─────────────┘
```

---

## State Container Definitions (types.py)

```python
from typing import NamedTuple
import jax.numpy as jnp

class FluidState(NamedTuple):
    """Conservative variables (updated each timestep)"""
    u: jnp.ndarray          # [nx, 3] conservative variables (rho, rho*vel, E)

class Primitives(NamedTuple):
    """Primitive variables (derived on-demand, discarded after use)"""
    rho: jnp.ndarray        # Density [nx]
    vel: jnp.ndarray        # Velocity [nx]
    p: jnp.ndarray          # Pressure [nx]
    a: jnp.ndarray          # Sound speed [nx]

class GridParams(NamedTuple):
    """Computational grid (immutable after initialization)"""
    x: jnp.ndarray          # Cell center coordinates [nx]
    x_face: jnp.ndarray     # Face coordinates [nx+1]
    dx: float               # Cell width (uniform grid)
    nx: int                 # Number of cells

class PhysicsParams(NamedTuple):
    """Physical constants (immutable)"""
    gamma: float            # Ratio of specific heats
    R: float                # Specific gas constant
    CFL: float              # CFL number for stability

class TimeState(NamedTuple):
    """Time stepping information"""
    t: float                # Current simulation time
    dt: float               # Current timestep
    step: int               # Current step number
```

---

## Variable Lifetime Categories

| Category | Examples | Stored In | Updated | Scope |
|----------|----------|-----------|---------|-------|
| **Persistent** | `u` (conservatives) | `FluidState` | Every timestep | Global |
| **Grid** | `x, dx, nx` | `GridParams` | Never | Global (immutable) |
| **Physics** | `gamma, CFL` | `PhysicsParams` | Never | Global (immutable) |
| **Time** | `t, dt, step` | `TimeState` | Every timestep | Global |
| **Derived** | `rho, vel, p, a` | `Primitives` | Computed on-demand | Local to function |
| **Transient** | `flux, rho_L, rho_R` | Local arrays | Each flux computation | Local to function |

### Transient Variables (computed and discarded)

These exist only within function scope:

| Variable | Shape | Description | Computed In |
|----------|-------|-------------|-------------|
| `rho_L, rho_R` | `[nx+1]` | Reconstructed interface values | `muscl_reconstruct` |
| `flux` | `[nx+1, 3]` | Interface fluxes | `compute_fluxes` |
| `wave_speeds` | `[nx+1, 3]` | Eigenvalues for Riemann solver | `compute_fluxes` |
| `limiter_phi` | `[nx]` | Slope limiter values | `apply_limiter` |

---

## Module Responsibilities

| Module | Inputs | Outputs | Notes |
|--------|--------|---------|-------|
| `types.py` | - | NamedTuples | Data structure definitions |
| `grid.py` | nx, domain | `GridParams` | Pure, called once |
| `primitives.py` | `u`, gamma | `Primitives` | Stateless conversion |
| `flux.py` | `Primitives`, grid | `flux[nx+1, 3]` | Interface fluxes |
| `muscl.py` | field array | `L, R` arrays | Reconstruction + limiters |
| `integrator.py` | state, flux, dt | new state | Forward Euler or RK |
| `boundary.py` | state, bc_type | state w/ ghosts | Handles edges |

---

## JIT Compilation Strategy

```python
# Compile entire timestep as single unit
@jax.jit
def step(state: FluidState, grid: GridParams, 
         params: PhysicsParams) -> FluidState:
    # 1. Compute primitives (transient)
    prims = compute_primitives(state.u, params.gamma)
    # 2. Reconstruct at interfaces (transient)
    rho_L, rho_R = muscl_reconstruct(prims.rho)
    # 3. Compute fluxes (transient)
    flux = compute_fluxes(rho_L, rho_R, ...)
    # 4. Update conservatives
    u_new = state.u - dt/grid.dx * (flux[1:] - flux[:-1])
    return state._replace(u=u_new)

# For multi-step simulation, use lax.scan
def simulate(state, grid, params, n_steps):
    def body(state, _):
        return step(state, grid, params), None
    final_state, _ = jax.lax.scan(body, state, None, length=n_steps)
    return final_state
```

---

## Immutable Update Pattern

```python
# ❌ Wrong: in-place modification
state.u = new_u

# ✅ Correct: create new state with _replace()
new_state = state._replace(u=new_u)

# Multiple fields at once
new_state = state._replace(u=new_u, t=state.t + dt)
```

---

## Module Dependency Graph

```
                    ┌──────────┐
                    │ types.py │
                    └────┬─────┘
                         │ (imported by all)
        ┌────────────────┼────────────────┐
        │                │                │
        ▼                ▼                ▼
  ┌──────────┐    ┌────────────┐    ┌──────────┐
  │ grid.py  │    │ initial.py │    │  io.py   │
  └────┬─────┘    └─────┬──────┘    └────┬─────┘
       │                │                │
       └────────┬───────┘                │
                ▼                        │
         ┌─────────────┐                 │
         │ boundary.py │                 │
         └──────┬──────┘                 │
                │                        │
                ▼                        │
        ┌──────────────┐                 │
        │ primitives.py│                 │
        └──────┬───────┘                 │
               │                         │
               ▼                         │
         ┌───────────┐                   │
         │ muscl.py  │                   │
         └─────┬─────┘                   │
               │                         │
               ▼                         │
          ┌─────────┐                    │
          │ flux.py │                    │
          └────┬────┘                    │
               │                         │
               ▼                         │
       ┌──────────────┐                  │
       │ integrator.py│                  │
       └──────┬───────┘                  │
              │                          │
              ▼                          │
        ┌───────────┐                    │
        │ solver.py │◄───────────────────┘
        └───────────┘
```

### Call Flow (single timestep)

```
solver.step()
    │
    ├──► boundary.apply_bc(state)
    │         │
    │         ▼
    ├──► primitives.compute(state.u, gamma)
    │         │
    │         ▼
    ├──► muscl.reconstruct(primitives)
    │         │
    │         ▼
    ├──► flux.compute(L_states, R_states)
    │         │
    │         ▼
    └──► integrator.update(state, flux, dt)
              │
              ▼
          new_state
```

### Module Import Summary

| Module | Imports From |
|--------|--------------|
| `types.py` | (none) |
| `grid.py` | `types` |
| `initial.py` | `types` |
| `boundary.py` | `types` |
| `primitives.py` | `types` |
| `muscl.py` | `types` |
| `flux.py` | `types`, `primitives` |
| `integrator.py` | `types` |
| `solver.py` | `types`, `grid`, `boundary`, `primitives`, `muscl`, `flux`, `integrator` |
| `io.py` | `types` |
