"""
Tests for data_structures.py - Taichi field creation and data classes.
"""

import sys
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

# Initialize Taichi before importing data_structures
import taichi as ti
ti.init(arch=ti.cpu)  # Use CPU for testing

from data_structures import (
    PhysicsParams, SimulationParams, LaserParams, OutputConfig, ToolPath,
    GridParams, State, StatePrev, MaterialProps, DiscretCoeffs, LaserState,
    TimeState, PoolDimensions, ConvergenceState
)


def test_physics_params_defaults():
    """Test PhysicsParams dataclass with defaults."""
    physics = PhysicsParams()
    
    assert physics.tsolid == 1563.0
    assert physics.tliquid == 1623.0
    assert physics.tpreheat == 300.0
    assert physics.sigma == 5.67e-8  # Stefan-Boltzmann constant
    assert physics.hlatent == 2.7e5
    
    print("✓ PhysicsParams defaults tests passed")


def test_simulation_params_defaults():
    """Test SimulationParams dataclass with defaults."""
    sim = SimulationParams()
    
    assert sim.delt == 1.0e-6
    assert sim.timax == 1.0e-3
    assert sim.ni == 100
    assert sim.nj == 100
    assert sim.nk == 50
    assert sim.max_iter == 1000
    
    print("✓ SimulationParams defaults tests passed")


def test_laser_params_peak_flux():
    """Test LaserParams peak flux calculation."""
    import math
    
    laser = LaserParams(power=200.0, radius=50.0e-6, absorptivity=0.35)
    
    expected_peak = 2.0 * 200.0 * 0.35 / (math.pi * (50.0e-6)**2)
    assert abs(laser.peak_flux - expected_peak) / expected_peak < 1e-10
    
    print("✓ LaserParams peak flux tests passed")


def test_toolpath_defaults():
    """Test ToolPath dataclass initialization."""
    toolpath = ToolPath()
    
    assert toolpath.n_segments == 1
    assert len(toolpath.time) == 1
    assert toolpath.time[0] == 0.0
    assert toolpath.laser_on[0] == 0
    
    print("✓ ToolPath defaults tests passed")


def test_grid_params_creation():
    """Test GridParams Taichi field creation."""
    ni, nj, nk = 10, 10, 5
    grid = GridParams(ni, nj, nk)
    
    assert grid.ni == ni
    assert grid.nj == nj
    assert grid.nk == nk
    
    # Verify 1D field shapes
    assert grid.x.shape == (ni,)
    assert grid.y.shape == (nj,)
    assert grid.z.shape == (nk,)
    assert grid.dx.shape == (ni,)
    assert grid.dy.shape == (nj,)
    assert grid.dz.shape == (nk,)
    
    # Verify 3D field shapes
    assert grid.vol.shape == (ni, nj, nk)
    
    # Verify 2D field shapes
    assert grid.areaij.shape == (ni, nj)
    assert grid.areaik.shape == (ni, nk)
    assert grid.areajk.shape == (nj, nk)
    
    print("✓ GridParams creation tests passed")


def test_state_creation():
    """Test State Taichi field creation."""
    ni, nj, nk = 10, 10, 5
    state = State(ni, nj, nk)
    
    assert state.ni == ni
    assert state.nj == nj
    assert state.nk == nk
    
    # Verify velocity field shapes
    assert state.uVel.shape == (ni, nj, nk)
    assert state.vVel.shape == (ni, nj, nk)
    assert state.wVel.shape == (ni, nj, nk)
    
    # Verify thermal field shapes
    assert state.temp.shape == (ni, nj, nk)
    assert state.enthalpy.shape == (ni, nj, nk)
    assert state.fracl.shape == (ni, nj, nk)
    
    # Verify pressure field shapes
    assert state.pressure.shape == (ni, nj, nk)
    assert state.pp.shape == (ni, nj, nk)
    
    print("✓ State creation tests passed")


def test_state_prev_creation():
    """Test StatePrev Taichi field creation."""
    ni, nj, nk = 10, 10, 5
    state_prev = StatePrev(ni, nj, nk)
    
    # Verify all previous-timestep field shapes
    assert state_prev.unot.shape == (ni, nj, nk)
    assert state_prev.vnot.shape == (ni, nj, nk)
    assert state_prev.wnot.shape == (ni, nj, nk)
    assert state_prev.tnot.shape == (ni, nj, nk)
    assert state_prev.hnot.shape == (ni, nj, nk)
    assert state_prev.fraclnot.shape == (ni, nj, nk)
    
    print("✓ StatePrev creation tests passed")


def test_material_props_creation():
    """Test MaterialProps Taichi field creation."""
    ni, nj, nk = 10, 10, 5
    mat_props = MaterialProps(ni, nj, nk)
    
    assert mat_props.vis.shape == (ni, nj, nk)
    assert mat_props.diff.shape == (ni, nj, nk)
    assert mat_props.den.shape == (ni, nj, nk)
    assert mat_props.tcond.shape == (ni, nj, nk)
    
    print("✓ MaterialProps creation tests passed")


def test_discret_coeffs_creation():
    """Test DiscretCoeffs Taichi field creation."""
    ni, nj, nk = 10, 10, 5
    coeffs = DiscretCoeffs(ni, nj, nk)
    
    # Check all coefficient fields
    assert coeffs.ap.shape == (ni, nj, nk)
    assert coeffs.ae.shape == (ni, nj, nk)
    assert coeffs.aw.shape == (ni, nj, nk)
    assert coeffs.an.shape == (ni, nj, nk)
    assert coeffs.as_.shape == (ni, nj, nk)  # as_ to avoid keyword
    assert coeffs.at.shape == (ni, nj, nk)
    assert coeffs.ab.shape == (ni, nj, nk)
    assert coeffs.su.shape == (ni, nj, nk)
    assert coeffs.sp.shape == (ni, nj, nk)
    
    print("✓ DiscretCoeffs creation tests passed")


def test_laser_state_creation():
    """Test LaserState creation."""
    ni, nj = 10, 10
    laser_state = LaserState(ni, nj)
    
    assert laser_state.beam_x == 0.0
    assert laser_state.beam_y == 0.0
    assert laser_state.laser_on == False
    assert laser_state.heatin.shape == (ni, nj)
    
    print("✓ LaserState creation tests passed")


def test_time_state_defaults():
    """Test TimeState dataclass defaults."""
    time_state = TimeState()
    
    assert time_state.timet == 0.0
    assert time_state.iteration == 0
    assert time_state.timestep == 0
    assert time_state.converged == False
    
    print("✓ TimeState defaults tests passed")


def test_pool_dimensions_defaults():
    """Test PoolDimensions dataclass defaults."""
    pool = PoolDimensions()
    
    assert pool.length == 0.0
    assert pool.width == 0.0
    assert pool.depth == 0.0
    assert pool.volume == 0.0
    assert pool.max_temp == 0.0
    
    print("✓ PoolDimensions defaults tests passed")


def test_convergence_state_defaults():
    """Test ConvergenceState dataclass defaults."""
    conv = ConvergenceState()
    
    assert conv.residual_h == 0.0
    assert conv.residual_u == 0.0
    assert conv.residual_v == 0.0
    assert conv.residual_w == 0.0
    assert conv.heat_ratio == 0.0
    assert conv.max_residual == 0.0
    
    print("✓ ConvergenceState defaults tests passed")


def test_taichi_field_operations():
    """Test basic Taichi field operations (fill and read)."""
    ni, nj, nk = 5, 5, 3
    
    # Create state and fill with preheat temperature
    state = State(ni, nj, nk)
    
    @ti.kernel
    def fill_temp(temp: ti.template(), value: ti.f64):
        for i, j, k in temp:
            temp[i, j, k] = value
    
    fill_temp(state.temp, 300.0)
    
    # Verify values using numpy conversion
    temp_np = state.temp.to_numpy()
    assert temp_np.shape == (ni, nj, nk)
    assert abs(temp_np[0, 0, 0] - 300.0) < 1e-10
    assert abs(temp_np[ni-1, nj-1, nk-1] - 300.0) < 1e-10
    
    print("✓ Taichi field operations tests passed")


def run_all_tests():
    """Run all data_structures tests."""
    print("\n" + "=" * 60)
    print("Testing data_structures.py - Type Definitions & Taichi Fields")
    print("=" * 60 + "\n")
    
    test_physics_params_defaults()
    test_simulation_params_defaults()
    test_laser_params_peak_flux()
    test_toolpath_defaults()
    test_grid_params_creation()
    test_state_creation()
    test_state_prev_creation()
    test_material_props_creation()
    test_discret_coeffs_creation()
    test_laser_state_creation()
    test_time_state_defaults()
    test_pool_dimensions_defaults()
    test_convergence_state_defaults()
    test_taichi_field_operations()
    
    print("\n" + "=" * 60)
    print("All data_structures.py tests passed! ✓")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    run_all_tests()
