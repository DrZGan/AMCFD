"""
Test script for flux.py module.

Verifies that the HeatFluxCalculator and heat_fluxes function work correctly.
"""

import taichi as ti
import numpy as np

# Initialize Taichi before importing Taichi-decorated modules
ti.init(arch=ti.cpu)

from data_structures import (
    State, StatePrev, GridParams, MaterialProps, LaserState,
    PhysicsParams, SimulationParams, ConvergenceState
)
from flux import HeatFluxCalculator, heat_fluxes


def setup_test_grid(ni: int, nj: int, nk: int) -> GridParams:
    """Create a simple uniform grid for testing."""
    grid = GridParams(ni, nj, nk)
    
    # Set uniform grid spacing
    dx_val = 1.0e-5  # 10 microns
    dy_val = 1.0e-5
    dz_val = 1.0e-5
    
    for i in range(ni):
        grid.x[i] = (i + 0.5) * dx_val
        grid.xu[i] = i * dx_val
        grid.dx[i] = dx_val
        grid.dxpwinv[i] = 1.0 / dx_val if i > 0 else 0.0
        grid.dxpeinv[i] = 1.0 / dx_val if i < ni - 1 else 0.0
    
    for j in range(nj):
        grid.y[j] = (j + 0.5) * dy_val
        grid.yv[j] = j * dy_val
        grid.dy[j] = dy_val
        grid.dypsinv[j] = 1.0 / dy_val if j > 0 else 0.0
        grid.dypninv[j] = 1.0 / dy_val if j < nj - 1 else 0.0
    
    for k in range(nk):
        grid.z[k] = (k + 0.5) * dz_val
        grid.zw[k] = k * dz_val
        grid.dz[k] = dz_val
        grid.dzpbinv[k] = 1.0 / dz_val if k > 0 else 0.0
        grid.dzptinv[k] = 1.0 / dz_val if k < nk - 1 else 0.0
    
    # Set volumes and areas
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                grid.vol[i, j, k] = dx_val * dy_val * dz_val
    
    for i in range(ni):
        for j in range(nj):
            grid.areaij[i, j] = dx_val * dy_val
    
    for i in range(ni):
        for k in range(nk):
            grid.areaik[i, k] = dx_val * dz_val
    
    for j in range(nj):
        for k in range(nk):
            grid.areajk[j, k] = dy_val * dz_val
    
    return grid


def test_flux_calculator():
    """Test the HeatFluxCalculator with synthetic data."""
    print("=" * 60)
    print("Testing HeatFluxCalculator")
    print("=" * 60)
    
    # Small test grid
    ni, nj, nk = 10, 10, 10
    
    # Create data structures
    grid = setup_test_grid(ni, nj, nk)
    state = State(ni, nj, nk)
    state_prev = StatePrev(ni, nj, nk)
    mat_props = MaterialProps(ni, nj, nk)
    laser_state = LaserState(ni, nj)
    conv = ConvergenceState()
    flux_calc = HeatFluxCalculator(ni, nj, nk)
    
    # Physics and simulation params
    physics = PhysicsParams(
        tpreheat=300.0,
        hpreheat=150000.0,
        hlatent=2.7e5,
        tsolid=1563.0,
        tliquid=1623.0
    )
    sim = SimulationParams(
        ni=ni, nj=nj, nk=nk,
        delt=1.0e-6
    )
    
    # Initialize fields with a temperature gradient (hot center, cool edges)
    print("\nInitializing fields with temperature gradient...")
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                # Distance from center
                cx, cy, cz = ni // 2, nj // 2, nk - 1  # Hot spot at top center
                dist = np.sqrt((i - cx)**2 + (j - cy)**2 + (k - cz)**2)
                
                # Enthalpy decreases with distance from hot spot
                h_max = 1.0e6  # Hot center
                h_min = 1.5e5  # Cool edges
                max_dist = np.sqrt(cx**2 + cy**2 + cz**2)
                
                h_val = h_max - (h_max - h_min) * (dist / max_dist)
                state.enthalpy[i, j, k] = h_val
                
                # Previous enthalpy slightly lower (heating)
                state_prev.hnot[i, j, k] = h_val * 0.99
                
                # Liquid fraction (1 where hot, 0 where cold)
                state.fracl[i, j, k] = min(1.0, max(0.0, (h_val - 5e5) / 5e5))
                state_prev.fraclnot[i, j, k] = state.fracl[i, j, k] * 0.99
                
                # Material properties
                mat_props.diff[i, j, k] = 25.0 / 800.0  # k/Cp approximation
                mat_props.den[i, j, k] = 7800.0
                mat_props.sourceinput[i, j, k] = 0.0
    
    # Add laser heat source at top center
    cx, cy = ni // 2, nj // 2
    for i in range(ni):
        for j in range(nj):
            dist = np.sqrt((i - cx)**2 + (j - cy)**2)
            if dist < 3:
                # Volumetric source in top layer
                mat_props.sourceinput[i, j, nk - 2] = 1.0e10  # W/m³
    
    laser_state.heat_total = 100.0  # 100 W laser
    ahtoploss = 50.0  # 50 W surface loss
    
    print("\nComputing boundary fluxes...")
    
    # Test individual flux computations
    flux_calc.compute_west_east_flux(
        state.enthalpy, mat_props.diff, grid.dxpwinv, grid.areajk)
    print(f"  West flux: {flux_calc.flux_west[None]:.6e} W")
    print(f"  East flux: {flux_calc.flux_east[None]:.6e} W")
    
    flux_calc.compute_bottom_top_flux(
        state.enthalpy, mat_props.diff, grid.dzpbinv, grid.areaij)
    print(f"  Bottom flux: {flux_calc.flux_bottom[None]:.6e} W")
    print(f"  Top flux: {flux_calc.flux_top[None]:.6e} W")
    
    flux_calc.compute_south_north_flux(
        state.enthalpy, mat_props.diff, grid.dypsinv, grid.areaik)
    print(f"  South flux: {flux_calc.flux_south[None]:.6e} W")
    print(f"  North flux: {flux_calc.flux_north[None]:.6e} W")
    
    flux_calc.compute_accumulation(
        state.enthalpy, state_prev.hnot,
        state.fracl, state_prev.fraclnot,
        mat_props.den, grid.vol,
        mat_props.sourceinput,
        physics.hlatent, sim.delt)
    print(f"\n  Heat accumulation: {flux_calc.accul[None]:.6e} W")
    print(f"  Volumetric heat source: {flux_calc.heatvol[None]:.6e} W")
    
    # Compute heat balance
    ratio = flux_calc.compute_heat_balance(ahtoploss, laser_state.heat_total)
    print(f"\n  Heat balance ratio: {ratio:.6f}")
    print(f"  Total heat out: {flux_calc.heatout[None]:.6e} W")
    
    print("\n" + "-" * 60)
    print("Testing full heat_fluxes function...")
    
    # Reset convergence state
    conv = ConvergenceState()
    
    # Call the full function
    conv = heat_fluxes(
        state, state_prev, mat_props, grid,
        physics, sim, laser_state,
        mat_props.sourceinput, ahtoploss,
        conv, flux_calc
    )
    
    print(f"  Convergence heat_ratio: {conv.heat_ratio:.6f}")
    
    # Validation checks
    print("\n" + "=" * 60)
    print("Validation Results")
    print("=" * 60)
    
    passed = True
    
    # Check that fluxes are computed (non-zero for this test case)
    if flux_calc.flux_west[None] == 0 and flux_calc.flux_east[None] == 0:
        print("WARNING: Both west and east fluxes are zero")
    
    # Check that accumulation is positive (we're heating)
    if flux_calc.accul[None] <= 0:
        print("WARNING: Heat accumulation should be positive for heating case")
        passed = False
    else:
        print("PASS: Heat accumulation is positive")
    
    # Check that ratio is reasonable (not NaN or Inf)
    if np.isnan(conv.heat_ratio) or np.isinf(conv.heat_ratio):
        print("FAIL: Heat ratio is NaN or Inf")
        passed = False
    else:
        print(f"PASS: Heat ratio is finite ({conv.heat_ratio:.6f})")
    
    # Check that volumetric source matches what we set
    if flux_calc.heatvol[None] > 0:
        print("PASS: Volumetric heat source is positive")
    else:
        print("WARNING: Volumetric heat source is zero or negative")
    
    if passed:
        print("\n*** All tests PASSED ***")
    else:
        print("\n*** Some tests FAILED ***")
    
    return passed


def test_zero_flux():
    """Test with uniform field (should have zero fluxes)."""
    print("\n" + "=" * 60)
    print("Testing with uniform field (zero flux case)")
    print("=" * 60)
    
    ni, nj, nk = 5, 5, 5
    
    grid = setup_test_grid(ni, nj, nk)
    state = State(ni, nj, nk)
    state_prev = StatePrev(ni, nj, nk)
    mat_props = MaterialProps(ni, nj, nk)
    flux_calc = HeatFluxCalculator(ni, nj, nk)
    
    # Uniform enthalpy everywhere
    uniform_h = 2.0e5
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                state.enthalpy[i, j, k] = uniform_h
                state_prev.hnot[i, j, k] = uniform_h
                state.fracl[i, j, k] = 0.0
                state_prev.fraclnot[i, j, k] = 0.0
                mat_props.diff[i, j, k] = 0.03
                mat_props.den[i, j, k] = 7800.0
                mat_props.sourceinput[i, j, k] = 0.0
    
    # Compute fluxes
    flux_calc.compute_west_east_flux(
        state.enthalpy, mat_props.diff, grid.dxpwinv, grid.areajk)
    flux_calc.compute_bottom_top_flux(
        state.enthalpy, mat_props.diff, grid.dzpbinv, grid.areaij)
    flux_calc.compute_south_north_flux(
        state.enthalpy, mat_props.diff, grid.dypsinv, grid.areaik)
    
    print(f"  West flux: {flux_calc.flux_west[None]:.6e}")
    print(f"  East flux: {flux_calc.flux_east[None]:.6e}")
    print(f"  Bottom flux: {flux_calc.flux_bottom[None]:.6e}")
    print(f"  Top flux: {flux_calc.flux_top[None]:.6e}")
    print(f"  South flux: {flux_calc.flux_south[None]:.6e}")
    print(f"  North flux: {flux_calc.flux_north[None]:.6e}")
    
    # All fluxes should be zero for uniform field
    tol = 1.0e-15
    all_zero = (
        abs(flux_calc.flux_west[None]) < tol and
        abs(flux_calc.flux_east[None]) < tol and
        abs(flux_calc.flux_bottom[None]) < tol and
        abs(flux_calc.flux_top[None]) < tol and
        abs(flux_calc.flux_south[None]) < tol and
        abs(flux_calc.flux_north[None]) < tol
    )
    
    if all_zero:
        print("\nPASS: All fluxes are zero for uniform field")
        return True
    else:
        print("\nFAIL: Fluxes should be zero for uniform field")
        return False


def test_fortran_index_mapping():
    """Test that index mapping matches Fortran convention.
    
    Fortran uses 1-based indexing:
    - Loops: do k=2,nkm1 means k from 2 to nk-1 (inclusive)
    - Arrays: enthalpy(1,j,k) is the west boundary
    
    Python uses 0-based indexing:
    - Loops: range(1, nk-1) means k from 1 to nk-2 (inclusive)
    - Arrays: enthalpy[0,j,k] is the west boundary
    """
    print("\n" + "=" * 60)
    print("Testing Fortran index mapping")
    print("=" * 60)
    
    ni, nj, nk = 5, 5, 5
    grid = setup_test_grid(ni, nj, nk)
    state = State(ni, nj, nk)
    mat_props = MaterialProps(ni, nj, nk)
    flux_calc = HeatFluxCalculator(ni, nj, nk)
    
    # Set specific values to verify indexing
    # Set a gradient only in i-direction
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                state.enthalpy[i, j, k] = float(i * 100)  # 0, 100, 200, 300, 400
                mat_props.diff[i, j, k] = 0.01  # Constant diffusion
    
    # Compute west/east flux
    flux_calc.compute_west_east_flux(
        state.enthalpy, mat_props.diff, grid.dxpwinv, grid.areajk)
    
    # Expected flux calculation:
    # West: diff[1,j,k] * (h[0] - h[1]) * dxpwinv[1] = 0.01 * (0 - 100) * 1e5 = -100
    # Sum over (nj-2) * (nk-2) = 3 * 3 = 9 cells
    # Each cell has areajk = 1e-5 * 1e-5 = 1e-10
    # flux_west = 9 * 1e-10 * (-100) = -9e-8
    
    # East: diff[3,j,k] * (h[4] - h[3]) * dxpwinv[4] = 0.01 * (400 - 300) * 1e5 = 100
    # flux_east = 9 * 1e-10 * 100 = 9e-8
    
    print(f"  West flux: {flux_calc.flux_west[None]:.6e}")
    print(f"  East flux: {flux_calc.flux_east[None]:.6e}")
    
    # Verify sign convention: flux is positive when heat flows OUT
    # West boundary: h[0] < h[1], so heat flows OUT (negative flux out = into domain)
    # East boundary: h[4] > h[3], so heat flows OUT (positive flux)
    
    passed = True
    
    # West flux should be negative (heat flowing INTO domain at west boundary)
    if flux_calc.flux_west[None] < 0:
        print("PASS: West flux is negative (heat flows into domain)")
    else:
        print("FAIL: West flux should be negative")
        passed = False
    
    # East flux should be positive (heat flowing OUT of domain at east boundary)
    if flux_calc.flux_east[None] > 0:
        print("PASS: East flux is positive (heat flows out of domain)")
    else:
        print("FAIL: East flux should be positive")
        passed = False
    
    # West and east should have equal magnitude (symmetric setup)
    if abs(abs(flux_calc.flux_west[None]) - abs(flux_calc.flux_east[None])) < 1e-15:
        print("PASS: West and east flux magnitudes are equal")
    else:
        print("FAIL: West and east flux magnitudes should be equal")
        passed = False
    
    return passed


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("FLUX MODULE TEST SUITE")
    print("=" * 60 + "\n")
    
    test1_passed = test_flux_calculator()
    test2_passed = test_zero_flux()
    test3_passed = test_fortran_index_mapping()
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Test 1 (gradient field):    {'PASSED' if test1_passed else 'FAILED'}")
    print(f"  Test 2 (uniform field):     {'PASSED' if test2_passed else 'FAILED'}")
    print(f"  Test 3 (index mapping):     {'PASSED' if test3_passed else 'FAILED'}")
    
    if test1_passed and test2_passed and test3_passed:
        print("\n*** ALL TESTS PASSED ***\n")
    else:
        print("\n*** SOME TESTS FAILED ***\n")
