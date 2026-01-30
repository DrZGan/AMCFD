# AM-CFD Fortran Migration Plan

## Code Overview
**Purpose:** Powder Bed Fusion CFD simulation (Additive Manufacturing)  
**Physics:** Heat transfer + Navier-Stokes (enthalpy, velocity, pressure)  
**Solver:** Line-by-line TDMA with OpenMP | **Grid:** Non-uniform staggered (max 1200×1200×180)

**Modules:** `constant`, `parameters`, `geometry`, `initialization`, `property`, `entotemp`, `laserinput`, `toolpath`, `boundary`, `discretization`, `source`, `solver`, `dimensions`, `convergence`, `residue`, `revision`, `fluxes`, `printing`, `main`

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

### Fortran Pattern Mapping

| Fortran Pattern | Target Equivalent |
|-----------------|-------------------|
| `common /block/ u,v,w` | State container/struct |
| `u(i,j,k) = ...` (in-place) | Return new array |
| Subroutine with `intent(inout)` | Pure function returning new state |
| `EQUIVALENCE(phi, uVel)` | Explicit field access |
| `allocate(u(ni,nj,nk))` | Array initialization in container |

---

## Migration Steps

---

### Step 1: Config & Grid

**Convert:** `mod_const.f90`, `mod_param.f90`, `mod_geom.f90`

**Prompt:**
```
Convert constants, parameters, and grid generation.
- Parse input_param.txt (namelist format) into config object
- Generate non-uniform staggered grid with power-law stretching
- Arrays: x, y, z, xu, yv, zw, dxpwinv, dxpeinv, dypsinv, dypninv, dzpbinv, dzptinv
- Cell volumes: vol(ni,nj,nk), face areas: areaij, areaik, areajk
- Use runtime allocation instead of hardcoded nx=1200, ny=1200, nz=180
```
**Test:** Compare grid arrays. Exact match expected.

---

### Step 2: State & Properties

**Convert:** `mod_init.f90`, `mod_prop.f90`, `mod_entot.f90`

**Prompt:**
```
Convert state initialization and material properties.
- Fields: uVel, vVel, wVel, pressure, pp, enthalpy, temp + previous timestep (unot, vnot, etc.)
- Properties: vis, diff, den | Coefficients: ap, an, as, ae, aw, at, ab | Source: su, sp
- Initialize to preheat temperature/enthalpy
- Enthalpy↔Temperature (3 regions):
  - Solid: T = (sqrt(acpb² + 2*acpa*H) - acpb) / acpa
  - Mushy: fracl = (H-hsmelt)/(hlcal-hsmelt), T = deltemp*fracl + tsolid
  - Liquid: T = (H - hlcal)/acpl + tliquid
```
**Test:** Compare initial state. Tolerance ~1e-10.

---

### Step 3: Laser & Toolpath

**Convert:** `mod_laser.f90`, `mod_toolpath.f90`

**Prompt:**
```
Convert laser heat source and toolpath.
- Load .crs file: 5 columns (time, x, y, z, laser_on)
- Track beam_pos, beam_posy; calculate scanvelx, scanvely
- Gaussian heat: heatin(i,j) = peakhin * exp(-alasfact * dist² / rb²)
- Total input: heatinLaser = sum(areaij * heatin)
```
**Test:** Compare heatin array for first timesteps. Tolerance ~1e-8.

---

### Step 4: Boundary Conditions

**Convert:** `mod_bound.f90`

**Prompt:**
```
Convert boundary conditions (selected by ivar=1,2,3,4,5).
- Velocity (ivar=1,2,3): Marangoni at top: uVel(i,j,nk) = uVel(i,j,nkm1) + fracl*dgdt*dT/dx/(vis*dz)
- Pressure (ivar=4): Zero gradient
- Enthalpy (ivar=5): Top=laser+radiation+convection; sides/bottom=convection; j=1 symmetry
```
**Test:** Compare su, sp arrays. Tolerance ~1e-8.

---

### Step 5: Discretization & Source

**Convert:** `mod_discret.f90`, `mod_sour.f90`

**Prompt:**
```
Convert FVM discretization (performance critical).
- Power-law scheme for convection/diffusion
- 7-point stencil: ap, ae, aw, an, as, at, ab (staggered grid)
- Source: mushy zone Darcy damping, buoyancy, latent heat
- Under-relaxation: ap = ap/urf, su += (1-urf)*ap*phi
- Vectorize for GPU
```
**Test:** Compare coefficient arrays. Tolerance ~1e-8.

---

### Step 6: Solver

**Convert:** `mod_solve.f90`

**Prompt:**
```
Convert TDMA solver.
- Line-by-line Thomas algorithm: solution_enthalpy, solution_uvw
- OpenMP parallel over j-lines → batch for GPU
```
**Test:** Compare solved fields. Tolerance ~1e-6.

---

### Step 7: Supporting Modules

**Convert:** `mod_dimen.f90`, `mod_converge.f90`, `mod_resid.f90`, `mod_revise.f90`

**Prompt:**
```
Convert supporting calculations.
- pool_size: melt pool length/depth/width by interpolating solidus isotherm
- residual: sum|ap*phi - sum(anb*phi_nb) - su|
- revision_p: velocity correction u += du*(p_west - p_east), pressure update
```
**Test:** Compare pool dimensions and residuals. Tolerance ~1e-6.

---

### Step 8: Main Loop & Output

**Convert:** `main.f90`, `mod_flux.f90`, `mod_print.f90`

**Prompt:**
```
Convert main time loop and output.
- Outer: timet += delt until timax; update laser position
- Inner: iterate until convergence (solve energy → H↔T → pool_size → momentum if melted)
- Convergence: amaxres < 5e-4 and 0.99 < ratio < 1.01 (heating); resorh < 5e-7 (cooling)
- Heat balance: ratio = (hin + heatvol) / (hout + accumulation)
- Output: time, iterations, residuals, pool dimensions, heat balance
- Replace GOTO 10/30/41/50 with proper loops
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

- **OpenMP → GPU:** Parallel j-loops → vectorized ops
- **TDMA:** Serial per line; consider batched TDMA or iterative for GPU
- **GOTO → loops:** main.f90 uses GOTO 10/30/41/50
- **EQUIVALENCE:** `phi` aliases u,v,w,p → make explicit
- **Phase change:** Solid/mushy/liquid with different H↔T relations
