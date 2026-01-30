# AM-CFD Fortran Migration Plan

## Code Overview
**Purpose:** Powder Bed Fusion CFD simulation (Additive Manufacturing)  
**Physics:** Heat transfer + Navier-Stokes (enthalpy, velocity, pressure)  
**Solver:** Line-by-line TDMA with OpenMP | **Grid:** Non-uniform staggered (max 1200×1200×180)

**Fortran Source Files:**
| File | Purpose |
|------|----------|
| `mod_const.f90` | Physical constants |
| `mod_param.f90` | Input parameters, namelist parsing |
| `mod_geom.f90` | Grid generation, geometry |
| `mod_init.f90` | State initialization |
| `mod_prop.f90` | Material properties |
| `mod_entot.f90` | Enthalpy↔Temperature conversion |
| `mod_laser.f90` | Laser heat source |
| `mod_toolpath.f90` | Toolpath loading (.crs files) |
| `mod_bound.f90` | Boundary conditions |
| `mod_discret.f90` | FVM discretization |
| `mod_sour.f90` | Source terms |
| `mod_solve.f90` | TDMA solver |
| `mod_dimen.f90` | Melt pool dimensions |
| `mod_converge.f90` | Convergence checks |
| `mod_resid.f90` | Residual calculation |
| `mod_revise.f90` | Pressure/velocity correction |
| `mod_flux.f90` | Heat flux balance |
| `mod_print.f90` | Output routines |
| `main.f90` | Main program, time loop |

**Input:** `inputfile/input_param.txt` (geometry, material, numerics, BCs) | `ToolFiles/*.crs` (toolpath)

---

## Code Architecture Principles

1. **Pure Functions**: All functions should be side-effect free. Inputs are read-only; return new objects instead of modifying in-place.

2. **State Bundling**: Group related variables into containers. This replaces Fortran's `common` blocks and global arrays.

3. **Function Signature Pattern**: All computational functions follow `(immutable_state, params) → new_state`

4. **Explicit Updates**: Create new state objects with modified fields rather than in-place mutation.

### Language-Specific Implementations

See the implementation file for your target language:
- **JAX**: [`jax_implementation.py`](./jax_implementation.py)
- **Taichi**: [`taichi_implementation.py`](./taichi_implementation.py) *(to be created)*

These files contain:
- State container definitions (data structures)
- File organization and module layout
- Variable lifetime categories
- Function patterns and JIT/compilation strategy

### Variable Lifetime Categories

| Category | Examples | Description |
|----------|----------|-------------|
| **Persistent** | `uVel, enthalpy, temp` | Updated every timestep |
| **Grid** | `x, y, z, vol` | Immutable after initialization |
| **Transient coefficients** | `ap, ae, su, sp` | Rebuilt each solver call |
| **Previous timestep** | `unot, hnot` | Saved at start of timestep |
| **Derived on-demand** | `fracl, vis` | Computed as needed |

---

## Migration Steps

> See [`STRUCTURE_PLAN.md`](./STRUCTURE_PLAN.md) for target module organization and data flow.
>
> **Dependencies:** Steps 1-2 must complete first. Steps 3-5 can run in parallel. Steps 6-8 depend on 3-5. Steps 9-10 depend on all previous.

---

### Step 1: Types & Config ✅ COMPLETE

**Create:** `data_structures.py`, `param.py`  
**Convert:** `mod_const.f90`, `mod_param.f90`

**Prompt:**
```
Create state containers and parse input configuration.
- Define NamedTuples: State, StatePrev, GridParams, MaterialProps, 
  DiscretCoeffs, LaserState, PhysicsParams, SimulationParams, TimeState
- Parse input_param.txt (namelist format) into PhysicsParams + SimulationParams
- Physical constants: acpa, acpb, acpl, tsolid, tliquid, dgdt, emiss, sigma, hconv
- Numerical params: delt, timax, urf_vel, urf_p, urf_h
```
**Test:** Verify all parameters parsed correctly. ✅ 8 tests pass

---

### Step 2: Grid Generation

**Create:** `geom.py`  
**Convert:** `mod_geom.f90`

**Prompt:**
```
Generate non-uniform staggered grid with power-law stretching.
- Arrays: x, y, z (cell centers), xu, yv, zw (velocity faces)
- Inverse distances: dxpwinv, dxpeinv, dypsinv, dypninv, dzpbinv, dzptinv
- Cell volumes: vol[ni,nj,nk], face areas: areaij, areaik, areajk
- Use runtime allocation instead of hardcoded nx=1200, ny=1200, nz=180
- Return GridParams NamedTuple
```
**Test:** Compare grid arrays. Exact match expected.

---

### Step 3: State Initialization

**Create:** `init.py`  
**Convert:** `mod_init.f90`

**Prompt:**
```
Initialize flow field to preheat conditions.
- Fields: uVel, vVel, wVel (zeros), pressure, pp (zeros)
- enthalpy: computed from preheat temperature
- temp: preheat temperature
- fracl: 0 (solid)
- Return State NamedTuple
```
**Test:** Compare initial state. Tolerance ~1e-10.

---

### Step 4: Properties & H↔T Conversion

**Create:** `prop.py`, `entot.py`  
**Convert:** `mod_prop.f90`, `mod_entot.f90`

**Prompt:**
```
Convert material properties and enthalpy↔temperature relations.
- Properties: vis(T), diff(T), den (constant or T-dependent)
- Enthalpy→Temperature (3 regions):
  - Solid (H < hsmelt): T = (sqrt(acpb² + 2*acpa*H) - acpb) / acpa
  - Mushy (hsmelt ≤ H ≤ hlcal): fracl = (H-hsmelt)/(hlcal-hsmelt), T = deltemp*fracl + tsolid
  - Liquid (H > hlcal): T = (H - hlcal)/acpl + tliquid
- Temperature→Enthalpy (inverse)
- Use jnp.where for JIT-compatible branching
```
**Test:** Compare T↔H round-trip. Tolerance ~1e-10.

---

### Step 5: Laser & Toolpath

**Create:** `laser.py`, `toolpath.py`  
**Convert:** `mod_laser.f90`, `mod_toolpath.f90`

**Prompt:**
```
Convert laser heat source and toolpath loading.
- Load .crs file: 5 columns (time, x, y, z, laser_on)
- Track beam_pos, beam_posy; calculate scanvelx, scanvely from toolpath
- Gaussian heat: heatin[i,j] = peakhin * exp(-alasfact * dist² / rb²)
- Total input: heatinLaser = sum(areaij * heatin)
- Return LaserState NamedTuple
```
**Test:** Compare heatin array for first timesteps. Tolerance ~1e-8.

---

### Step 6: Boundary Conditions

**Create:** `bound.py`  
**Convert:** `mod_bound.f90`

**Prompt:**
```
Convert boundary conditions (selected by ivar=1,2,3,4,5).
- Velocity (ivar=1,2,3): Marangoni at top: uVel(i,j,nk) = uVel(i,j,nkm1) + fracl*dgdt*dT/dx/(vis*dz)
- Pressure (ivar=4): Zero gradient
- Enthalpy (ivar=5): Top=laser+radiation+convection; sides/bottom=convection; j=1 symmetry
- Return updated su, sp arrays
```
**Test:** Compare su, sp arrays. Tolerance ~1e-8.

---

### Step 7: Discretization & Source

**Create:** `discret.py`, `sour.py`  
**Convert:** `mod_discret.f90`, `mod_sour.f90`

**Prompt:**
```
Convert FVM discretization (performance critical).
- Power-law scheme for convection/diffusion
- 7-point stencil coefficients: ap, ae, aw, an, as_, at, ab (staggered grid)
- Source terms: mushy zone Darcy damping, buoyancy, latent heat
- Under-relaxation: ap = ap/urf, su += (1-urf)*ap*phi
- Vectorize for GPU (no explicit loops)
- Return DiscretCoeffs NamedTuple
```
**Test:** Compare coefficient arrays. Tolerance ~1e-8.

---

### Step 8: TDMA Solver

**Create:** `solve.py`  
**Convert:** `mod_solve.f90`

**Prompt:**
```
Convert TDMA solver.
- Line-by-line Thomas algorithm: solve_enthalpy, solve_uvw
- OpenMP parallel over j-lines → batch/vectorize for GPU
- Consider batched TDMA or iterative solver for better GPU utilization
- Return new field array (do not modify in place)
```
**Test:** Compare solved fields. Tolerance ~1e-6.

---

### Step 9: Convergence & Pool Size

**Create:** `converge.py`, `dimen.py`, `resid.py`, `flux.py`  
**Convert:** `mod_converge.f90`, `mod_resid.f90`, `mod_dimen.f90`, `mod_flux.f90`

**Prompt:**
```
Convert convergence checks and melt pool calculations.
- Residual: sum|ap*phi - sum(anb*phi_nb) - su|
- Heat balance ratio: (hin + heatvol) / (hout + accumulation)
- Convergence: amaxres < 5e-4 and 0.99 < ratio < 1.01 (heating); resorh < 5e-7 (cooling)
- Pool dimensions: length/depth/width by interpolating solidus isotherm
```
**Test:** Compare pool dimensions and residuals. Tolerance ~1e-6.

---

### Step 10: Main Loop & Output

**Create:** `main.py`, `revise.py`  
**Convert:** `main.f90`, `mod_revise.f90`, `mod_print.f90`

**Prompt:**
```
Convert main time loop and output.
- Outer loop: timet += delt until timax; update laser position each step
- Inner loop: iterate until convergence (solve energy → H↔T → pool_size → momentum if melted)
- Pressure correction: revision_p → velocity correction u += du*(p_west - p_east), pressure update
- Replace GOTO 10/30/41/50 with proper while loops and conditionals
- Output: time, iterations, residuals, pool dimensions, heat balance to output.txt
- Tecplot output: tecmov*.tec
- Use jax.lax.while_loop for outer time loop, lax.scan for inner iteration
```
**Test:** Compare output.txt. Tolerance ~1e-4.

---

## Key Data Structures

**3D Fields:** `uVel, vVel, wVel, pressure, pp, enthalpy, temp, fracl, vis, diff, den, ap/ae/aw/an/as/at/ab, su, sp`  
**2D:** `heatin(nx,ny)` - laser flux | `areaij, areaik, areajk` - face areas  
**Toolpath:** `toolmatrix(TOOLLINES,5)` - time,x,y,z,laser_on | `coordhistory(COORDLINES,8)`

---

## Testing

1. Add dump statements after: grid generation, initialization, each solver call, each timestep
2. Write arrays to binary files
3. **Tolerances:** Grid/init: exact | Per-iteration: ~1e-8 | Multi-step: ~1e-4

---

## Technical Notes

### Fortran → JAX/Taichi Patterns

| Fortran Pattern | JAX Solution | Taichi Solution |
|-----------------|--------------|------------------|
| OpenMP parallel j-loops | `jax.vmap` over j | `ti.kernel` with parallel for |
| In-place array update | Return new array with `_replace()` | Use `ti.field` with explicit copy |
| GOTO statements | `jax.lax.while_loop` / `cond` | Python while/if |
| EQUIVALENCE aliases | Explicit field access | Explicit field access |
| Allocatable arrays | Dynamic shapes via params | Static shapes with `ti.field` |

### TDMA Solver Strategy

- **Fortran:** Line-by-line serial TDMA, OpenMP parallel over j-lines
- **JAX Option 1:** Batched TDMA using `jax.lax.scan` over lines
- **JAX Option 2:** Iterative solver (Jacobi/Gauss-Seidel) for better GPU parallelism
- **Taichi:** Parallel TDMA with `ti.kernel`, or use sparse solver

### Control Flow Conversion

```fortran
! Fortran GOTO pattern (main.f90)
10 CONTINUE
   ... compute ...
   IF (condition) GOTO 30
   GOTO 10
30 CONTINUE
```

```python
# JAX equivalent
def iteration_loop(carry):
    state, converged = carry
    new_state = compute_step(state)
    return (new_state, check_convergence(new_state))

final_state, _ = jax.lax.while_loop(
    lambda carry: ~carry[1],  # continue while not converged
    iteration_loop,
    (initial_state, False)
)
```

### Phase Change Handling

- **Solid:** H < hsmelt → T from quadratic Cp relation
- **Mushy:** hsmelt ≤ H ≤ hlcal → Linear interpolation, fracl ∈ [0,1]
- **Liquid:** H > hlcal → T from constant Cp
- Use `jnp.where` (JAX) or `ti.select` (Taichi) for branchless computation

### Performance Considerations

1. **Memory:** 3D fields at 1200×1200×180 = ~260M cells × 8 bytes × ~20 fields ≈ 40 GB
2. **Batch size:** May need to tile domain for GPU memory limits
3. **JIT compilation:** First call slow, subsequent calls fast
4. **Precision:** Use float64 for accuracy matching Fortran
