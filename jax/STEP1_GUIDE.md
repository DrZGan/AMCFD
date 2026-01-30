# Step 1: Types & Config ✅ COMPLETE

## Files Created

| File | Lines | Purpose | Fortran Source |
|------|-------|---------|----------------|
| `data_structures.py` | 370 | All NamedTuple definitions (State, GridParams, PhysicsParams, etc.) | `mod_const.f90` |
| `param.py` | 530 | Parse Fortran input files (geometry + namelists), load toolpath, write output | `mod_param.f90`, `mod_print.f90` |
| `test_step1.py` | 280 | Unit tests validating parsing and data structures | - |
| `inputfile/` | - | Input files copied from Fortran: `input_param.txt`, `B26.crs`, `B26-1.crs` | - |

**Naming Convention**: Python files match Fortran module names without the `mod_` prefix (e.g., `mod_param.f90` → `param.py`).

## Input Files

| File | Format | Purpose |
|------|--------|---------|
| `input_param.txt` | Fortran namelist + line-by-line geometry | Simulation parameters |
| `B26.crs` | 5-column text (time, x, y, z, laser_on) | Full toolpath |
| `B26-1.crs` | 5-column text | Simplified toolpath (6 segments) |

## Running Tests

```bash
cd /home/mzomoro1/bin/Hackaton/jax
python test_step1.py
```

Expected: 8 tests pass.

## What Each Test Does

| Test | Validates |
|------|-----------|
| `test_parse_value()` | Fortran notation (`1.0d-6`), booleans (`.true.`), arrays |
| `test_enthalpy_conversion()` | `H = (acpa/2)*T² + acpb*T` calculation |
| `test_physics_params()` | Material properties, computed enthalpies (`hsmelt`, `hlcal`) |
| `test_simulation_params()` | Grid dims, timestep, under-relaxation factors |
| `test_laser_params()` | Peak flux: `q = 2*P*α/(π*rb²)` |
| `test_geometry_parsing()` | Zone-based grid (`nzx`, `ncvx`, `powrx`, etc.) |
| `test_sample_input_file()` | End-to-end parsing of actual `input_param.txt` |
| `test_toolpath_loading()` | Load `.crs` toolpath files |

## Input File Format

The actual input format from Fortran (`input_param.txt`):

```
!-----geometrical parameters---------
1                    ! nzx (number of x-zones)
4.0e-3               ! xzone lengths
200                  ! ncvx (control volumes per zone)
1                    ! powrx (stretching exponent)
...

&material_properties  dens=8440, tsolid=1563, tliquid=1623, hsmelt=861e3, .../
&volumetric_parameters alaspowvol=300, alasetavol=0.3512, sourcerad=31.71e-6/
&numerical_relax maxit=50, delt=2e-5, urfu=0.7, urfp=0.7, urfh=0.7/
```

## Key Design Choices

1. **NamedTuples**: Immutable for JAX JIT compilation
2. **Update pattern**: `state._replace(temp=new_temp)` not `state.temp[i] = x`
3. **Computed params**: `hsmelt`, `hlcal`, `peak_flux` derived from inputs
4. **Zone-based grid**: Supports multi-zone with power-law stretching

## Fortran Variable Mapping

| Fortran | Python | Description |
|---------|--------|-------------|
| `dens`, `denl` | `rho`, `rholiq` | Density (solid/liquid) |
| `thconsa`, `thconsb` | `tcond` | Thermal conductivity |
| `alaspowvol`, `alasetavol` | `power`, `absorptivity` | Laser params |
| `dgdtp` | `dgdt` | Surface tension gradient |
| `nzx`, `ncvx`, `powrx` | Zone geometry | Multi-zone grid definition |

## Optional Modules (Not in Step 1)

These Fortran modules in `other files/` are **optional physics** for future steps:

| Module | Purpose | When Needed |
|--------|---------|-------------|
| `mod_freesur.f90` | Free surface deformation (melt pool shape) | Step 5+ |
| `mod_solidification.f90` | Solidification tracking (G, R) | Post-processing |
| `mod_species.f90` | Species transport (multi-component) | Advanced features |

## Migration Plan Checklist

- [x] Define NamedTuples: State, GridParams, PhysicsParams, SimulationParams, etc.
- [x] Parse input_param.txt (geometry + namelist format)
- [x] Physical constants: acpa, acpb, tsolid, tliquid, dgdt, emiss, sigma, hconv
- [x] Numerical params: delt, timax, urf_vel, urf_p, urf_h
- [x] Load toolpath from .crs files
- [x] Copy input files to jax/inputfile/
- [x] Unit tests: All parameters parsed correctly (8 tests)

## Next: Step 2 (Grid Generation)

Convert `mod_geom.f90` → `grid.py`: Non-uniform staggered grid with power-law stretching.
