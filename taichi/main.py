"""
AM-CFD Taichi Implementation - Main Program

Translated from Fortran main.f90
Version 1.0 - Python/Taichi port

This is the main simulation driver that orchestrates all CFD computations.
"""

import os
import time as pytime
from pathlib import Path

import taichi as ti
import numpy as np

# Initialize Taichi (must be done before importing Taichi-decorated modules)
ti.init(arch=ti.cpu)

# Data structures
from data_structures import (
    PhysicsParams, SimulationParams, LaserParams, OutputConfig, ToolPath,
    GridParams, State, StatePrev, MaterialProps, DiscretCoeffs, LaserState,
    TimeState, PoolDimensions, ConvergenceState
)

# I/O and geometry
from param import parse_input, load_toolpath
from geom import get_gridparams
from bound import bound_condition
<<<<<<< HEAD
from entot import enthalpy_to_temp as entot_enthalpy_to_temp
from dimen import pool_size, clean_uvw
from revise import revision_p as revise_revision_p
=======
from discret import discretize
>>>>>>> 2390666 (adding discret.py and test_discret.py)

# Import the new flux module
from flux import HeatFluxCalculator, heat_fluxes as compute_heat_fluxes

# TODO: These modules will be created from corresponding .f90 files
# from initialization import initialize
# from source import source_term
# from residue import residual
# from solver import solution_enthalpy, solution_uvw
# from fluxes import heat_fluxes
# from properties import properties
# from entotemp import enthalpy_to_temp
# from convergence import enhance_converge_speed
# from revision import revision_p
# from laserinput import laser_beam, calc_rhf
# from toolpath import read_coordinates
# from printing import output_results, custom_output


# ============================================================================
# Placeholder functions (to be replaced by actual module imports)
# ============================================================================


def initialize(state: State, state_prev: StatePrev, mat_props: MaterialProps,
               physics: PhysicsParams, grid: GridParams) -> None:
    """Initialize all field variables. (From initialization.f90)
    
    Modifies state, state_prev, and mat_props in-place.
    """
    # Set initial temperature and enthalpy to preheat values
    for i in range(grid.ni):
        for j in range(grid.nj):
            for k in range(grid.nk):
                state.temp[i, j, k] = physics.tpreheat
                state.enthalpy[i, j, k] = physics.hpreheat
                state.fracl[i, j, k] = 0.0
                state.uVel[i, j, k] = 0.0
                state.vVel[i, j, k] = 0.0
                state.wVel[i, j, k] = 0.0
                state.pressure[i, j, k] = 0.0


def properties(state: State, mat_props: MaterialProps, physics: PhysicsParams) -> MaterialProps:
    """Update temperature-dependent properties. (From property.f90)
    
    Returns:
        mat_props: Updated material properties (also modified in-place)
    """
    # TODO: Implement temperature-dependent viscosity, diffusion, etc.
    return mat_props



def source_term(ivar: int, state: State, state_prev: StatePrev,
                grid: GridParams, coeffs: DiscretCoeffs,
                laser_state: LaserState, physics: PhysicsParams,
                sim: SimulationParams) -> DiscretCoeffs:
    """Add source terms to discretized equations. (From source.f90)
    
    Args:
        ivar: Variable index (1=u, 2=v, 3=w, 4=p, 5=enthalpy)
        
    Returns:
        coeffs: Updated coefficients with source terms (su, sp) added
    """
    # TODO: Implement source terms (laser heat, buoyancy, Darcy damping)
    return coeffs


def residual(ivar: int, state: State, coeffs: DiscretCoeffs, 
             conv: ConvergenceState) -> ConvergenceState:
    """Compute residual errors. (From residue.f90)
    
    Args:
        ivar: Variable index (1=u, 2=v, 3=w, 4=p, 5=enthalpy)
        
    Returns:
        conv: Updated convergence state with residuals (resorm/resoru/resorv/resorw/resorh)
    """
    # TODO: Implement residual computation
    return conv


def enhance_converge_speed(coeffs: DiscretCoeffs, grid: GridParams) -> DiscretCoeffs:
    """Enhance convergence using slice-wise TDMA. (From convergence.f90)
    
    Computes residual error in each x-direction slice to enhance convergence.
    
    Returns:
        coeffs: Updated discretization coefficients
    """
    # TODO: Implement slice-wise residual computation
    return coeffs


def solution_enthalpy(state: State, coeffs: DiscretCoeffs, grid: GridParams,
                      sim: SimulationParams) -> State:
    """Solve enthalpy equation using line-by-line TDMA. (From solver.f90)
    
    Returns:
        state: Updated state with new enthalpy field
    """
    # TODO: Implement line-by-line TDMA solver
    return state


def solution_uvw(ivar: int, state: State, coeffs: DiscretCoeffs, 
                 grid: GridParams, sim: SimulationParams) -> State:
    """Solve momentum equations. (From solver.f90)
    
    Args:
        ivar: Variable index (1=u, 2=v, 3=w)
        
    Returns:
        state: Updated state with new velocity components
    """
    # TODO: Implement momentum solver
    return state


def enthalpy_to_temp(state: State, physics: PhysicsParams) -> State:
    """Convert enthalpy field to temperature. (From entotemp.f90)
    
    Translates enthalpy to temperature, handling phase change.
    
    Returns:
        state: Updated state with temperature and liquid fraction fields
    """
    # TODO: Implement enthalpy-to-temperature conversion with mushy zone
    return state


def pool_size(state: State, grid: GridParams, physics: PhysicsParams,
              pool: PoolDimensions) -> tuple[int, int, int, int, int, int]:
    """Determine melt pool dimensions and bounds. (From dimensions.f90)
    
    Gets melt pool dimension, start and end index of i,j,k to determine fluid region.
    Also updates pool.max_temp (tpeak in Fortran).
    
    Returns:
        tuple: (ist, ien, jst, jen, kst, ken) - start/end indices for melt pool region
    """
    # TODO: Find liquid region bounds and compute pool dimensions
    pool.max_temp = 300.0  # Placeholder - should compute actual peak temp
    return (1, 1, 1, 1, 1, 1)  # ist, ien, jst, jen, kst, ken


def clean_uvw(state: State, grid: GridParams, physics: PhysicsParams) -> State:
    """Zero velocity outside melt pool. (From dimensions.f90)
    
    Sets velocity to zero for cells where temp <= tsolid (outside liquid region).
    
    Returns:
        state: Updated state with velocities zeroed in solid region
    """
    # TODO: Zero velocities in solid cells
    return state


def revision_p(state: State, coeffs: DiscretCoeffs, grid: GridParams,
               sim: SimulationParams) -> State:
    """Pressure correction (SIMPLE algorithm). (From revision.f90)
    
    Corrects pressure and velocities using SIMPLE algorithm.
    
    NOTE: This is now imported from revise.py module.
    The actual function signature requires additional parameters (ivar, physics, domain bounds).
    This wrapper is kept for compatibility but needs updating in the main loop.
    
    Returns:
        state: Updated state with corrected pressure and velocities
    """
    # This needs to be called with proper parameters in the main loop
    # For now, return unchanged state
    return state


def read_coordinates(time_state: TimeState, toolpath: ToolPath, 
                     laser_state: LaserState) -> LaserState:
    """Read/interpolate coordinates from toolpath. (From toolpath.f90)
    
    Reads current position along toolpath based on simulation time.
    
    Returns:
        laser_state: Updated with current beam coordinates
    """
    # TODO: Interpolate toolpath to get current coordinates
    return laser_state


def output_results(time_state: TimeState, state: State, conv: ConvergenceState,
                   pool: PoolDimensions, output: OutputConfig) -> None:
    """Output results to files and screen. (From printing.f90)
    
    Outputs parameters to file and screen (corresponds to outputres in Fortran).
    """
    print(f"  Time: {time_state.timet:.6e}s, Iter: {time_state.iteration}, "
          f"MaxRes: {conv.max_residual:.4e}, Ratio: {conv.heat_ratio:.4f}")


def custom_output(time_state: TimeState, state: State, grid: GridParams,
                  output: OutputConfig) -> None:
    """Output thermal cycles and surface data. (From printing.f90)
    
    Outputs thermal cycle, velocity and temperature at surface and symmetry plane.
    Corresponds to Cust_Out in Fortran.
    """
    # TODO: Implement custom output for thermal cycles
    pass


# ============================================================================
# Main Program
# ============================================================================

def main():
    """Main simulation driver - mirrors main.f90 structure."""
    
    # ------------------------------------------------------------------
    # Setup: Read input and initialize
    # ------------------------------------------------------------------
    
    # Paths
    base_dir = Path(__file__).parent.parent
    input_yaml = base_dir / "inputfile" / "input_param.yaml"
    toolpath_file = base_dir / "inputfile" / "toolpath.crs"
    
    # Check input file exists
    if not input_yaml.exists():
        raise FileNotFoundError(f"Input file not found: {input_yaml}")
    
    print("=" * 60)
    print("AM-CFD Taichi Implementation")
    print("=" * 60)
    
    # Read input parameters (corresponds to read_data in Fortran)
    print("Reading input parameters...")
    physics, sim, laser, output_cfg = parse_input(str(input_yaml))
  
    
    # Read toolpath (corresponds to read_toolpath in Fortran)
    if toolpath_file.exists():
        print("Reading toolpath...")
        toolpath = load_toolpath(str(toolpath_file))
    else:
        print("No toolpath file found, using default.")
        toolpath = ToolPath()
    
    # Generate grid (corresponds to generate_grid in Fortran)
    print("Generating computational grid...")
    grid = get_gridparams(str(input_yaml))
    ni, nj, nk = grid.ni, grid.nj, grid.nk
    print(f"  Grid dimensions: ni={ni}, nj={nj}, nk={nk}")
    
    # Allocate state arrays
    print("Allocating state arrays...")
    state = State(ni, nj, nk)
    state_prev = StatePrev(ni, nj, nk)
    material_props = MaterialProps(ni, nj, nk)
    discret_coeffs = DiscretCoeffs(ni, nj, nk)
    laser_state = LaserState(ni, nj)
    
    # Allocate flux calculator
    flux_calc = HeatFluxCalculator(ni, nj, nk)
    
    # Initialize variables (corresponds to initialize in Fortran)
    print("Initializing fields...")
    initialize(state, state_prev, mat_props, physics, grid)
    simulation_params = SimulationParams(ni, nj, nk)
    # Initialize tracking structures
    time_state = TimeState()
    pool_dimensions = PoolDimensions()
    convergence_state = ConvergenceState()
    
    # Record start time
    wall_start = pytime.time()
    print(f"\nSimulation started at: {pytime.strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 60)
    
    # ------------------------------------------------------------------
    # Time integration loop (GOTO 10 in Fortran)
    # ------------------------------------------------------------------
    
    small = 1.0e-7
    time_state.timet = small
    itertot = 0
    
    while time_state.timet < sim.timax:
        # Advance time
        time_state.timet += sim.delt
        time_state.timestep += 1
        niter = 0
        
        # Update laser position and compute heat flux
        # Corresponds to: call laser_beam, call read_coordinates, call calcRHF
        laser_beam(time_state, toolpath, laser_state, laser)
        read_coordinates(time_state, toolpath, laser_state)
        calc_rhf(laser_state, grid, laser)
        
        # ------------------------------------------------------------------
        # Iteration loop within timestep (GOTO 30 in Fortran)
        # ------------------------------------------------------------------
        
        converged = False
        while not converged and niter < sim.max_iter:
            niter += 1
            itertot += 1
            time_state.iteration = niter
            
            # ============================================================
            # Solve energy equation (ivar=5)
            # ============================================================
            ivar = 5
            
            properties(state, mat_props, physics)
            ahtoploss = bound_condition(ivar, state, grid, mat_props, physics, simu_params, laser_state)
            discretize(ivar, state, state_prev, grid, coeffs, mat_props, sim, physics)
            source_term(ivar, state, state_prev, grid, coeffs, 
                       laser_state, physics, sim)
            residual(ivar, state, coeffs, conv)
            
            enhance_converge_speed(coeffs, grid)
            solution_enthalpy(state, coeffs, grid, sim)
            
            enthalpy_to_temp(state, physics)
            
            # ============================================================
            # Determine melt pool and solve momentum if needed
            # ============================================================
            ist, ien, jst, jen, kst, ken = pool_size(state, grid, physics, pool)
            
            # Skip momentum if peak temp <= solidus (no liquid)
            if pool.max_temp <= physics.tsolid:
                pass  # Skip to convergence check
            else:
                clean_uvw(state, grid, physics)
                
                # Solve momentum equations (ivar=1,2,3) and pressure (ivar=4)
                for ivar in range(1, 5):
                    bound_condition(ivar, state, grid, mat_props, physics, simu_params)
                    discretize(ivar, state, state_prev, grid, coeffs, mat_props, sim, physics)
                    source_term(ivar, state, state_prev, grid, coeffs,
                               laser_state, physics, sim)
                    residual(ivar, state, coeffs, conv)
                    solution_uvw(ivar, state, coeffs, grid, sim)
                    revision_p(state, coeffs, grid, sim)
            
            # ============================================================
            # Check convergence
            # ============================================================
            # Use the new heat_fluxes function
            compute_heat_fluxes(state, state_prev, mat_props, grid,
                               physics, sim, laser_state,
                               mat_props.sourceinput, ahtoploss,
                               conv, flux_calc)
            
            conv.max_residual = max(conv.residual_h, conv.residual_u,
                                   conv.residual_v, conv.residual_w)
            
            # Convergence criteria (heating vs cooling)
            if laser_state.laser_on:
                # Heating: check residual and heat balance
                if (conv.max_residual < sim.conv_tol and 
                    0.99 <= conv.heat_ratio <= 1.01):
                    converged = True
            else:
                # Cooling: only check enthalpy residual
                if conv.residual_h < 5.0e-7:
                    converged = True
        
        # End of iteration loop
        
        # ------------------------------------------------------------------
        # Post-timestep processing
        # ------------------------------------------------------------------
        
        # Output results
        output_results(time_state, state, conv, pool, output_cfg)
        
        # Zero velocity in solid regions and store previous values
        # Corresponds to Fortran loop before Cust_Out
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    if state.temp[i, j, k] <= physics.tsolid:
                        state.uVel[i, j, k] = 0.0
                        state.vVel[i, j, k] = 0.0
                        state.wVel[i, j, k] = 0.0
                    
                    # Store for next timestep (unot, vnot, wnot, tnot, hnot, fraclnot)
                    state_prev.unot[i, j, k] = state.uVel[i, j, k]
                    state_prev.vnot[i, j, k] = state.vVel[i, j, k]
                    state_prev.wnot[i, j, k] = state.wVel[i, j, k]
                    state_prev.tnot[i, j, k] = state.temp[i, j, k]
                    state_prev.hnot[i, j, k] = state.enthalpy[i, j, k]
                    state_prev.fraclnot[i, j, k] = state.fracl[i, j, k]
        
        # Custom output (thermal cycles, surface data)
        custom_output(time_state, state, grid, output_cfg)
    
    # End of time loop
    
    # ------------------------------------------------------------------
    # Finalize
    # ------------------------------------------------------------------
    
    wall_end = pytime.time()
    elapsed = wall_end - wall_start
    
    print("-" * 60)
    print(f"Simulation completed at: {pytime.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total wall time: {elapsed:.2f} seconds")
    print(f"Total iterations: {itertot}")
    print("=" * 60)


if __name__ == "__main__":
    main()
