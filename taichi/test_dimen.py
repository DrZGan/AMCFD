"""
Test module for dimen.py - Melt Pool Dimensions

Tests the pool_size and clean_uvw functions with synthetic data.
"""

import numpy as np
import taichi as ti

# Initialize Taichi
ti.init(arch=ti.cpu)

from data_structures import (
    State, GridParams, PhysicsParams, PoolDimensions
)
from dimen import pool_size, clean_uvw, compute_pool_volume, find_max_velocity


def test_pool_size_no_melt():
    """Test pool_size when temperature is below solidus (no melt pool)."""
    print("=" * 60)
    print("Test 1: No melt pool (all solid)")
    print("=" * 60)
    
    # Setup small grid
    ni, nj, nk = 10, 10, 10
    grid = GridParams(ni, nj, nk)
    
    # Initialize grid coordinates (simple uniform)
    for i in range(ni):
        grid.x[i] = i * 0.1e-3
    for j in range(nj):
        grid.y[j] = j * 0.1e-3
    for k in range(nk):
        grid.z[k] = k * 0.05e-3
    
    # Initialize state with low temperature (all solid)
    state = State(ni, nj, nk)
    state.temp.fill(300.0)  # Below solidus
    
    # Physics parameters
    physics = PhysicsParams()
    physics.tsolid = 1563.0
    physics.tliquid = 1623.0
    
    # Pool dimensions
    pool = PoolDimensions()
    
    # Call pool_size
    istat, iend, jstat, jend, kstat, kend = pool_size(state, grid, physics, pool)
    
    # Verify results
    print(f"  Peak temperature: {pool.max_temp:.2f} K")
    print(f"  Pool length: {pool.length*1e3:.6f} mm")
    print(f"  Pool width: {pool.width*1e3:.6f} mm")
    print(f"  Pool depth: {pool.depth*1e3:.6f} mm")
    print(f"  Solution domain: i=[{istat},{iend}], j=[{jstat},{jend}], k=[{kstat},{kend}]")
    
    assert pool.length == 0.0, "Pool length should be zero with no melt"
    assert pool.width == 0.0, "Pool width should be zero with no melt"
    assert pool.depth == 0.0, "Pool depth should be zero with no melt"
    print("  ✓ Test passed: No melt pool detected correctly\n")


def test_pool_size_with_melt():
    """Test pool_size with a synthetic melt pool."""
    print("=" * 60)
    print("Test 2: With melt pool")
    print("=" * 60)
    
    # Setup grid
    ni, nj, nk = 20, 15, 15
    grid = GridParams(ni, nj, nk)
    
    # Initialize grid coordinates
    for i in range(ni):
        grid.x[i] = i * 0.05e-3
    for j in range(nj):
        grid.y[j] = j * 0.05e-3
    for k in range(nk):
        grid.z[k] = k * 0.05e-3
    
    # Initialize state
    state = State(ni, nj, nk)
    
    # Create synthetic Gaussian melt pool at surface center
    temp_np = np.zeros((ni, nj, nk))
    x_center = 0.5e-3  # Center position
    y_center = 0.0     # Symmetry plane
    z_surface = (nk - 1) * 0.05e-3
    
    sigma_x = 0.3e-3   # Pool extent in x
    sigma_y = 0.15e-3  # Pool extent in y
    sigma_z = 0.2e-3   # Pool depth
    
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                x = i * 0.05e-3
                y = j * 0.05e-3
                z = k * 0.05e-3
                
                # Gaussian temperature distribution
                r2 = ((x - x_center) / sigma_x)**2 + \
                     ((y - y_center) / sigma_y)**2 + \
                     ((z - z_surface) / sigma_z)**2
                
                temp_np[i, j, k] = 300.0 + 2000.0 * np.exp(-r2)
    
    # Copy to Taichi field
    state.temp.from_numpy(temp_np)
    
    # Physics parameters
    physics = PhysicsParams()
    physics.tsolid = 1563.0
    physics.tliquid = 1623.0
    
    # Pool dimensions
    pool = PoolDimensions()
    
    # Call pool_size
    istat, iend, jstat, jend, kstat, kend = pool_size(state, grid, physics, pool)
    
    # Print results
    print(f"  Peak temperature: {pool.max_temp:.2f} K")
    print(f"  Pool length: {pool.length*1e6:.2f} µm")
    print(f"  Pool width: {pool.width*1e6:.2f} µm")
    print(f"  Pool depth: {pool.depth*1e6:.2f} µm")
    print(f"  Solution domain: i=[{istat},{iend}], j=[{jstat},{jend}], k=[{kstat},{kend}]")
    
    # Verify that pool dimensions are reasonable
    assert pool.length > 0.0, "Pool length should be positive"
    assert pool.width >= 0.0, "Pool width should be non-negative"
    assert pool.depth >= 0.0, "Pool depth should be non-negative"
    assert pool.max_temp > physics.tsolid, "Peak temp should exceed solidus"
    print("  ✓ Test passed: Melt pool detected correctly\n")


def test_clean_uvw():
    """Test clean_uvw function to zero velocities in solid regions."""
    print("=" * 60)
    print("Test 3: Clean velocities in solid region")
    print("=" * 60)
    
    # Setup grid
    ni, nj, nk = 10, 10, 10
    grid = GridParams(ni, nj, nk)
    
    # Initialize state
    state = State(ni, nj, nk)
    
    # Set temperature field: half liquid, half solid
    temp_np = np.zeros((ni, nj, nk))
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                if i < 5:
                    temp_np[i, j, k] = 1700.0  # Liquid
                else:
                    temp_np[i, j, k] = 1000.0  # Solid
    
    state.temp.from_numpy(temp_np)
    
    # Set velocities everywhere (non-zero)
    state.uVel.fill(1.0)
    state.vVel.fill(0.5)
    state.wVel.fill(-0.3)
    
    # Physics parameters
    physics = PhysicsParams()
    physics.tsolid = 1563.0
    
    # Clean velocities
    state = clean_uvw(state, grid, physics)
    
    # Check that velocities are zeroed in solid region
    uVel_np = state.uVel.to_numpy()
    vVel_np = state.vVel.to_numpy()
    wVel_np = state.wVel.to_numpy()
    
    # Verify liquid region still has velocity
    assert np.any(uVel_np[:5, :, :] != 0.0), "Liquid region should have non-zero velocity"
    
    # Verify solid region has zero velocity
    u_solid_max = np.max(np.abs(uVel_np[5:, :, :]))
    v_solid_max = np.max(np.abs(vVel_np[5:, :, :]))
    w_solid_max = np.max(np.abs(wVel_np[5:, :, :]))
    
    print(f"  Max |u| in solid region: {u_solid_max:.6e}")
    print(f"  Max |v| in solid region: {v_solid_max:.6e}")
    print(f"  Max |w| in solid region: {w_solid_max:.6e}")
    
    assert u_solid_max == 0.0, "Solid region should have zero u-velocity"
    assert v_solid_max == 0.0, "Solid region should have zero v-velocity"
    assert w_solid_max == 0.0, "Solid region should have zero w-velocity"
    print("  ✓ Test passed: Velocities zeroed in solid region correctly\n")


def test_utility_functions():
    """Test utility functions: compute_pool_volume and find_max_velocity."""
    print("=" * 60)
    print("Test 4: Utility functions")
    print("=" * 60)
    
    # Setup grid
    ni, nj, nk = 10, 10, 10
    grid = GridParams(ni, nj, nk)
    
    # Initialize volumes (uniform cells)
    vol_np = np.ones((ni, nj, nk)) * 1.0e-12  # 1 µm³ per cell
    grid.vol.from_numpy(vol_np)
    
    # Initialize state
    state = State(ni, nj, nk)
    
    # Set temperature: 200 cells liquid, 800 cells solid
    temp_np = np.zeros((ni, nj, nk))
    temp_np[:5, :, :4] = 1700.0  # 5*10*4 = 200 cells liquid
    temp_np[5:, :, :] = 1000.0   # Solid
    temp_np[:5, :, 4:] = 1000.0  # Solid
    state.temp.from_numpy(temp_np)
    
    # Physics
    physics = PhysicsParams()
    physics.tsolid = 1563.0
    
    # Compute pool volume
    volume = compute_pool_volume(state, grid, physics)
    expected_volume = 200 * 1.0e-12  # 200 liquid cells
    
    print(f"  Computed pool volume: {volume:.6e} m³")
    print(f"  Expected volume: {expected_volume:.6e} m³")
    assert abs(volume - expected_volume) < 1e-15, "Pool volume calculation error"
    
    # Set velocities
    state.uVel.fill(1.5)
    state.vVel.fill(-0.8)
    state.wVel.fill(2.3)
    
    # Find max velocities
    umax, vmax, wmax = find_max_velocity(state, grid)
    
    print(f"  Max velocities: u={umax:.6f}, v={vmax:.6f}, w={wmax:.6f} m/s")
    print(f"  Expected: u=1.5, v=0.8, w=2.3 m/s")
    assert abs(umax - 1.5) < 1e-6, f"Max u-velocity incorrect: {umax}"
    assert abs(vmax - 0.8) < 1e-6, f"Max v-velocity incorrect: {vmax}"
    assert abs(wmax - 2.3) < 1e-6, f"Max w-velocity incorrect: {wmax}"
    print("  ✓ Test passed: Utility functions working correctly\n")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Testing dimen.py - Melt Pool Dimensions Module")
    print("=" * 60 + "\n")
    
    test_pool_size_no_melt()
    test_pool_size_with_melt()
    test_clean_uvw()
    test_utility_functions()
    
    print("=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)
