"""
Test script for boundary conditions module

This tests the bound_condition function with a simple setup.

Test Summary:
-----------------------------------------------------------------------------
| Test                            | Purpose                                  |
-----------------------------------------------------------------------------
| test_velocity_boundary_conditions| Tests velocity BCs at domain boundaries |
| test_enthalpy_boundary_conditions| Tests enthalpy BCs at domain boundaries |
| test_marangoni_effect           | Tests Marangoni (surface tension) effect|
-----------------------------------------------------------------------------
"""

import taichi as ti
import numpy as np
from bound import bound_condition
from data_structures import (
    State, GridParams, MaterialProps, PhysicsParams, SimulationParams, LaserState
)


def create_test_grid(ni=10, nj=10, nk=10):
    """Create a simple uniform test grid."""
    grid = GridParams(ni, nj, nk)
    
    # Simple uniform grid
    xlen, ylen, zlen = 1e-3, 1e-3, 0.5e-3
    dx = xlen / ni
    dy = ylen / nj
    dz = zlen / nk
    
    # Fill grid arrays
    for i in range(ni):
        grid.x[i] = (i + 0.5) * dx
        grid.xu[i] = i * dx
        grid.dx[i] = dx
        grid.dxpwinv[i] = 1.0 / dx
        grid.dxpeinv[i] = 1.0 / dx
    
    for j in range(nj):
        grid.y[j] = (j + 0.5) * dy
        grid.yv[j] = j * dy
        grid.dy[j] = dy
        grid.dypsinv[j] = 1.0 / dy
        grid.dypninv[j] = 1.0 / dy
    
    for k in range(nk):
        grid.z[k] = (k + 0.5) * dz
        grid.zw[k] = k * dz
        grid.dz[k] = dz
        grid.dzpbinv[k] = 1.0 / dz
        grid.dzptinv[k] = 1.0 / dz
    
    # Fill areas
    for i in range(ni):
        for j in range(nj):
            grid.areaij[i, j] = dx * dy
    
    return grid


def create_test_state(ni=10, nj=10, nk=10):
    """Create test state with some initial values."""
    state = State(ni, nj, nk)
    
    # Initialize with some test values
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                # Temperature: hotter at top
                state.temp[i, j, k] = 300.0 + 1000.0 * (k / nk)
                
                # Enthalpy
                state.enthalpy[i, j, k] = 150000.0 + 400000.0 * (k / nk)
                
                # Liquid fraction: fully solid except near top
                if k >= nk - 2:
                    state.fracl[i, j, k] = 0.5 + 0.5 * (k / nk)
                else:
                    state.fracl[i, j, k] = 0.0
                
                # Initial velocities
                state.uVel[i, j, k] = 0.001 * (i / ni)
                state.vVel[i, j, k] = 0.001 * (j / nj)
                state.wVel[i, j, k] = 0.0
                state.pp[i, j, k] = 0.0
    
    return state


def create_test_material_props(ni=10, nj=10, nk=10):
    """Create test material properties."""
    mat_props = MaterialProps(ni, nj, nk)
    
    # Initialize with constant values
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                mat_props.vis[i, j, k] = 0.006  # 6 mPa·s
                mat_props.diff[i, j, k] = 3.2e-6  # k/(rho*cp)
                mat_props.den[i, j, k] = 7800.0
                mat_props.tcond[i, j, k] = 25.0
    
    return mat_props


def create_test_laser(ni=10, nj=10):
    """Create test laser state."""
    laser = LaserState(ni, nj)
    
    # Simple Gaussian heat flux at center
    xc, yc = ni // 2, nj // 2
    radius = 2.0
    peak_flux = 1e7  # 10 MW/m²
    
    for i in range(ni):
        for j in range(nj):
            r2 = (i - xc) ** 2 + (j - yc) ** 2
            laser.heatin[i, j] = peak_flux * np.exp(-2.0 * r2 / radius ** 2)
    
    laser.laser_on = True
    return laser


def test_velocity_boundary_conditions():
    """Test velocity boundary conditions."""
    print("\n=== Testing Velocity Boundary Conditions ===\n")
    
    ti.init(arch=ti.cpu)
    
    ni, nj, nk = 10, 10, 10
    grid = create_test_grid(ni, nj, nk)
    state = create_test_state(ni, nj, nk)
    mat_props = create_test_material_props(ni, nj, nk)
    physics = PhysicsParams()
    simu_params = SimulationParams(ni=ni, nj=nj, nk=nk)
    
    # Test u-velocity BC
    print("Testing u-velocity boundary conditions (ivar=1)...")
    ahtoploss = bound_condition(1, state, grid, mat_props, physics, simu_params)
    print(f"  - Heat loss returned: {ahtoploss:.2e} W (should be 0 for velocity BC)")
    
    # Check boundary values
    u_left = state.uVel.to_numpy()[0, :, :]
    u_right = state.uVel.to_numpy()[-1, :, :]
    print(f"  - Left boundary u max: {np.max(np.abs(u_left)):.2e} (should be ~0)")
    print(f"  - Right boundary u max: {np.max(np.abs(u_right)):.2e} (should be ~0)")
    
    # Test v-velocity BC
    print("\nTesting v-velocity boundary conditions (ivar=2)...")
    ahtoploss = bound_condition(2, state, grid, mat_props, physics, simu_params)
    
    v_left = state.vVel.to_numpy()[0, :, :]
    v_right = state.vVel.to_numpy()[-1, :, :]
    print(f"  - Left boundary v max: {np.max(np.abs(v_left)):.2e} (should be ~0)")
    print(f"  - Right boundary v max: {np.max(np.abs(v_right)):.2e} (should be ~0)")
    
    # Test w-velocity BC
    print("\nTesting w-velocity boundary conditions (ivar=3)...")
    ahtoploss = bound_condition(3, state, grid, mat_props, physics, simu_params)
    
    w_left = state.wVel.to_numpy()[0, :, :]
    w_right = state.wVel.to_numpy()[-1, :, :]
    print(f"  - Left boundary w max: {np.max(np.abs(w_left)):.2e} (should be ~0)")
    print(f"  - Right boundary w max: {np.max(np.abs(w_right)):.2e} (should be ~0)")
    
    # Test pressure BC
    print("\nTesting pressure boundary conditions (ivar=4)...")
    ahtoploss = bound_condition(4, state, grid, mat_props, physics, simu_params)
    
    pp_left = state.pp.to_numpy()[0, :, :]
    pp_right = state.pp.to_numpy()[-1, :, :]
    print(f"  - Left boundary pp max: {np.max(np.abs(pp_left)):.2e} (should be ~0)")
    print(f"  - Right boundary pp max: {np.max(np.abs(pp_right)):.2e} (should be ~0)")
    
    print("\n✓ Velocity boundary conditions test completed")


def test_enthalpy_boundary_conditions():
    """Test enthalpy boundary conditions."""
    print("\n=== Testing Enthalpy Boundary Conditions ===\n")
    
    ti.init(arch=ti.cpu)
    
    ni, nj, nk = 10, 10, 10
    grid = create_test_grid(ni, nj, nk)
    state = create_test_state(ni, nj, nk)
    mat_props = create_test_material_props(ni, nj, nk)
    physics = PhysicsParams()
    simu_params = SimulationParams(ni=ni, nj=nj, nk=nk)
    laser = create_test_laser(ni, nj)
    
    print("Testing enthalpy boundary conditions (ivar=5)...")
    
    # Store initial enthalpy at boundaries
    h_top_before = state.enthalpy.to_numpy()[:, :, -1].copy()
    h_bottom_before = state.enthalpy.to_numpy()[:, :, 0].copy()
    
    # Apply BC
    ahtoploss = bound_condition(5, state, grid, mat_props, physics, simu_params, laser)
    
    # Check results
    h_top_after = state.enthalpy.to_numpy()[:, :, -1]
    h_bottom_after = state.enthalpy.to_numpy()[:, :, 0]
    
    print(f"  - Total heat loss at top surface: {ahtoploss:.2e} W")
    print(f"  - Top surface enthalpy changed: {not np.allclose(h_top_before, h_top_after)}")
    print(f"  - Bottom surface enthalpy changed: {not np.allclose(h_bottom_before, h_bottom_after)}")
    print(f"  - Top enthalpy range: [{np.min(h_top_after):.2e}, {np.max(h_top_after):.2e}]")
    print(f"  - Bottom enthalpy range: [{np.min(h_bottom_after):.2e}, {np.max(h_bottom_after):.2e}]")
    
    # Check that heat loss is positive (reasonable)
    if ahtoploss > 0:
        print(f"  ✓ Heat loss is positive (physically reasonable)")
    else:
        print(f"  ⚠ Heat loss is non-positive (check BC implementation)")
    
    print("\n✓ Enthalpy boundary conditions test completed")


def test_marangoni_effect():
    """Test Marangoni stress at top surface."""
    print("\n=== Testing Marangoni Effect ===\n")
    
    ti.init(arch=ti.cpu)
    
    ni, nj, nk = 20, 20, 10
    grid = create_test_grid(ni, nj, nk)
    state = create_test_state(ni, nj, nk)
    mat_props = create_test_material_props(ni, nj, nk)
    physics = PhysicsParams()
    simu_params = SimulationParams(ni=ni, nj=nj, nk=nk)
    
    # Create temperature gradient at top surface
    for i in range(ni):
        for j in range(nj):
            # Linear gradient: hotter in center, cooler at edges
            dx = (i - ni / 2) / (ni / 2)
            dy = (j - nj / 2) / (nj / 2)
            r = np.sqrt(dx ** 2 + dy ** 2)
            state.temp[i, j, nk - 1] = 1800.0 - 200.0 * r
            state.fracl[i, j, nk - 1] = 1.0  # Fully liquid at top
    
    # Store initial u, v at top
    u_top_before = state.uVel.to_numpy()[:, :, -1].copy()
    v_top_before = state.vVel.to_numpy()[:, :, -1].copy()
    
    # Apply boundary conditions
    bound_condition(1, state, grid, mat_props, physics, simu_params)  # u-velocity
    bound_condition(2, state, grid, mat_props, physics, simu_params)  # v-velocity
    
    # Check results
    u_top_after = state.uVel.to_numpy()[:, :, -1]
    v_top_after = state.vVel.to_numpy()[:, :, -1]
    
    u_change = np.max(np.abs(u_top_after - u_top_before))
    v_change = np.max(np.abs(v_top_after - v_top_before))
    
    print(f"  - Maximum u-velocity change at top: {u_change:.2e} m/s")
    print(f"  - Maximum v-velocity change at top: {v_change:.2e} m/s")
    print(f"  - dgdt parameter: {physics.dgdt:.2e} N/(m·K)")
    
    if u_change > 1e-10 or v_change > 1e-10:
        print(f"  ✓ Marangoni effect is active (velocities changed)")
    else:
        print(f"  ⚠ Marangoni effect seems inactive (no velocity change)")
    
    print("\n✓ Marangoni effect test completed")


if __name__ == "__main__":
    print("=" * 60)
    print("Boundary Conditions Module Test Suite")
    print("=" * 60)
    
    try:
        test_velocity_boundary_conditions()
        test_enthalpy_boundary_conditions()
        test_marangoni_effect()
        
        print("\n" + "=" * 60)
        print("All tests completed successfully!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
