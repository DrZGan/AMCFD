"""
Tests for toolpath.py - Toolpath loading and coordinate tracking.
"""

import sys
import numpy as np
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

# Initialize Taichi before importing data_structures
import taichi as ti
ti.init(arch=ti.cpu)

from data_structures import ToolPath, LaserState, TimeState
from toolpath import (
    load_toolpath, get_current_segment, interpolate_position,
    read_coordinates, calc_rhf, CoordinateHistory
)

# Path to toolpath files
TOOLPATH_FILE = Path(__file__).parent / "inputfile" / "B26-1.crs"


def test_load_toolpath():
    """Test loading toolpath from .crs file."""
    if not TOOLPATH_FILE.exists():
        print("⚠ Skipping toolpath loading test - B26-1.crs not found")
        return
    
    toolpath = load_toolpath(str(TOOLPATH_FILE))
    
    # Check that data was loaded
    assert toolpath.n_segments == 6, f"Expected 6 segments, got {toolpath.n_segments}"
    
    # Check first segment (all zeros)
    assert float(toolpath.time[0]) == 0.0, f"Expected time[0]=0, got {toolpath.time[0]}"
    assert float(toolpath.x[0]) == 0.0, f"Expected x[0]=0, got {toolpath.x[0]}"
    assert int(toolpath.laser_on[0]) == 0, f"Expected laser_on[0]=0, got {toolpath.laser_on[0]}"
    
    # Check segment with laser on (index 2)
    assert int(toolpath.laser_on[2]) == 1, f"Expected laser_on[2]=1, got {toolpath.laser_on[2]}"
    
    # Check that x, y, z coordinates are reasonable (in mm range)
    assert float(toolpath.x[2]) > 0, f"Expected positive x[2], got {toolpath.x[2]}"
    
    print("✓ Toolpath loading tests passed")


def test_empty_toolpath():
    """Test handling of empty/missing toolpath file."""
    # Create a ToolPath with defaults
    toolpath = ToolPath()
    
    assert toolpath.n_segments == 1
    assert len(toolpath.time) == 1
    assert toolpath.time[0] == 0.0
    
    print("✓ Empty toolpath tests passed")


def test_get_current_segment():
    """Test finding current segment from time."""
    # Create a simple toolpath
    toolpath = ToolPath(
        time=np.array([0.0, 1.0, 2.0, 3.0]),
        x=np.array([0.0, 1.0, 2.0, 3.0]),
        y=np.array([0.0, 0.0, 0.0, 0.0]),
        z=np.array([0.0, 0.0, 0.0, 0.0]),
        laser_on=np.array([0, 1, 1, 0], dtype=np.int32),
        n_segments=4
    )
    
    # Test various times
    assert get_current_segment(0.0, toolpath) == 0
    assert get_current_segment(0.5, toolpath) == 0
    assert get_current_segment(1.0, toolpath) == 1
    assert get_current_segment(1.5, toolpath) == 1
    assert get_current_segment(2.5, toolpath) == 2
    assert get_current_segment(5.0, toolpath) == 3  # Past end
    
    print("✓ Get current segment tests passed")


def test_interpolate_position():
    """Test position interpolation along toolpath."""
    toolpath = ToolPath(
        time=np.array([0.0, 1.0, 2.0]),
        x=np.array([0.0, 1.0, 3.0]),
        y=np.array([0.0, 2.0, 2.0]),
        z=np.array([0.0, 0.0, 0.0]),
        laser_on=np.array([0, 1, 1], dtype=np.int32),
        n_segments=3
    )
    
    # At start
    x, y, z, vx, vy = interpolate_position(0.0, toolpath)
    assert x == 0.0
    assert y == 0.0
    
    # Midway through first segment
    x, y, z, vx, vy = interpolate_position(0.5, toolpath)
    assert abs(x - 0.5) < 1e-10
    assert abs(y - 1.0) < 1e-10
    assert abs(vx - 1.0) < 1e-10  # dx/dt = 1.0/1.0
    assert abs(vy - 2.0) < 1e-10  # dy/dt = 2.0/1.0
    
    # At segment boundary
    x, y, z, vx, vy = interpolate_position(1.0, toolpath)
    assert abs(x - 1.0) < 1e-10
    assert abs(y - 2.0) < 1e-10
    
    print("✓ Position interpolation tests passed")


def test_read_coordinates():
    """Test read_coordinates function updates LaserState."""
    toolpath = ToolPath(
        time=np.array([0.0, 1.0, 2.0]),
        x=np.array([0.0, 1.0e-3, 2.0e-3]),
        y=np.array([0.0, 0.0, 0.0]),
        z=np.array([0.0, 0.0, 0.0]),
        laser_on=np.array([0, 1, 0], dtype=np.int32),
        n_segments=3
    )
    
    time_state = TimeState(timet=0.5)
    laser_state = LaserState(ni=10, nj=10)
    
    result = read_coordinates(time_state, toolpath, laser_state)
    
    assert result.current_segment == 0
    assert abs(result.beam_x - 0.5e-3) < 1e-10
    assert result.laser_on == False  # laser_on[0] = 0
    
    # Move to segment with laser on
    time_state.timet = 1.5
    result = read_coordinates(time_state, toolpath, laser_state)
    
    assert result.current_segment == 1
    assert result.laser_on == True  # laser_on[1] = 1
    
    print("✓ Read coordinates tests passed")


def test_coordinate_history():
    """Test CoordinateHistory class for RHF calculation."""
    history = CoordinateHistory(max_entries=5)
    
    # Add some entries
    history.add_entry(0.0, 0.0, 0.0, 0.0, 200.0, 1.0, 1.0, 0.0)
    history.add_entry(0.1, 0.1e-3, 0.0, 0.0, 200.0, 1.0, 1.0, 0.0)
    history.add_entry(0.2, 0.2e-3, 0.0, 0.0, 200.0, 1.0, 1.0, 0.0)
    
    assert history.count == 3
    
    # Most recent should be at index 0
    assert history.history[0, 0] == 0.2  # time
    assert history.history[0, 1] == 0.2e-3  # x
    
    # Older entries shifted down
    assert history.history[1, 0] == 0.1
    assert history.history[2, 0] == 0.0
    
    print("✓ Coordinate history tests passed")


def test_calc_rhf():
    """Test RHF (Reheating Factor) calculation."""
    # Create a simple history
    history = np.array([
        [0.1, 0.1e-3, 0.0, 0.0, 300.0, 1.0, 1.0, 0.0],  # Most recent
        [0.0, 0.0, 0.0, 0.0, 300.0, 1.0, 1.0, 0.0],     # Previous
        [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],  # Empty
    ])
    
    # Current position at (0.1e-3, 0.0), current time 0.1
    rhf = calc_rhf(
        history,
        current_pos=(0.1e-3, 0.0),
        current_time=0.1,
        power=300.0,
        R=0.2e-3,  # 200 um radius
        T=2e-3,    # 2 ms window
        P0=300.0   # Reference power
    )
    
    # RHF should be non-negative
    assert rhf >= 0.0
    
    print("✓ RHF calculation tests passed")


def run_all_tests():
    """Run all toolpath tests."""
    print("\n" + "=" * 60)
    print("Testing toolpath.py - Toolpath Loading & Tracking")
    print("=" * 60 + "\n")
    
    test_empty_toolpath()
    test_get_current_segment()
    test_interpolate_position()
    test_load_toolpath()
    test_read_coordinates()
    test_coordinate_history()
    test_calc_rhf()
    
    print("\n" + "=" * 60)
    print("All toolpath.py tests passed! ✓")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    run_all_tests()
