"""
Test Step 1: Types & Config

Verifies that input parsing works correctly and all parameters are extracted.
Uses the actual input_param.txt from the Fortran code.
"""

import sys
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

from data_structures import PhysicsParams, SimulationParams, LaserParams
from param import parse_input, load_toolpath, _read_namelist, _temp_to_enthalpy_solid

# Path to actual input file (copied from fortran/inputfile/)
INPUT_FILE = Path(__file__).parent / "inputfile" / "input_param.txt"


def test_parse_value():
    """Test parsing individual values."""
    from param import _parse_value, _parse_single_value
    
    # Test integers
    assert _parse_single_value("42") == 42
    assert _parse_single_value("  100  ") == 100
    
    # Test floats
    assert abs(_parse_single_value("3.14") - 3.14) < 1e-10
    assert abs(_parse_single_value("1.0e-6") - 1.0e-6) < 1e-16
    
    # Test Fortran double precision
    assert abs(_parse_single_value("1.0d-6") - 1.0e-6) < 1e-16
    assert abs(_parse_single_value("2.5D+3") - 2500.0) < 1e-10
    
    # Test booleans
    assert _parse_value(".true.") == True
    assert _parse_value(".false.") == False
    assert _parse_value(".T.") == True
    assert _parse_value(".F.") == False
    
    # Test strings
    assert _parse_value("'hello'") == "hello"
    assert _parse_value('"world"') == "world"
    
    # Test arrays
    result = _parse_value("1.0, 2.0, 3.0")
    assert len(result) == 3
    assert abs(result[0] - 1.0) < 1e-10
    
    print("✓ Value parsing tests passed")


def test_enthalpy_conversion():
    """Test enthalpy calculation."""
    # Test with typical steel values
    acpa = 0.0
    acpb = 500.0
    T = 1000.0
    
    H = _temp_to_enthalpy_solid(T, acpa, acpb)
    expected = 0.5 * acpa * T * T + acpb * T
    
    assert abs(H - expected) < 1e-10
    assert abs(H - 500000.0) < 1e-10  # 500 * 1000
    
    # Test with quadratic term
    acpa = 0.5
    H = _temp_to_enthalpy_solid(T, acpa, acpb)
    expected = 0.5 * 0.5 * 1000.0 * 1000.0 + 500.0 * 1000.0
    assert abs(H - expected) < 1e-10
    
    print("✓ Enthalpy conversion tests passed")


def test_physics_params():
    """Test PhysicsParams creation with defaults."""
    # Create with minimal parameters
    params = {
        'tsolid': 1563.0,
        'tliquid': 1623.0,
        'tpreheat': 300.0,
        'rho': 7800.0,
        'hlatent': 2.7e5,
    }
    
    from param import _create_physics_params
    physics = _create_physics_params(params)
    
    assert physics.tsolid == 1563.0
    assert physics.tliquid == 1623.0
    assert physics.rho == 7800.0
    assert physics.hlatent == 2.7e5
    assert physics.sigma == 5.67e-8  # Stefan-Boltzmann constant
    
    # Check computed enthalpy values
    assert physics.hsmelt > 0
    assert physics.hlcal > physics.hsmelt
    assert physics.hlcal == physics.hsmelt + physics.hlatent
    
    print("✓ PhysicsParams tests passed")


def test_simulation_params():
    """Test SimulationParams creation."""
    params = {
        'delt': 1.0e-6,
        'timax': 1.0e-3,
        'nx': 100,
        'ny': 100,
        'nz': 50,
    }
    
    from param import _create_simulation_params
    sim = _create_simulation_params(params)
    
    assert sim.delt == 1.0e-6
    assert sim.timax == 1.0e-3
    assert sim.ni == 100
    assert sim.nj == 100
    assert sim.nk == 50
    
    print("✓ SimulationParams tests passed")


def test_laser_params():
    """Test LaserParams creation and peak flux calculation."""
    import math
    
    params = {
        'power': 200.0,
        'rb': 50.0e-6,
        'absorp': 0.35,
    }
    
    from param import _create_laser_params
    laser = _create_laser_params(params)
    
    assert laser.power == 200.0
    assert laser.radius == 50.0e-6
    assert laser.absorptivity == 0.35
    
    # Verify peak flux calculation
    expected_peak = 2.0 * 200.0 * 0.35 / (math.pi * (50.0e-6)**2)
    assert abs(laser.peak_flux - expected_peak) / expected_peak < 1e-10
    
    print("✓ LaserParams tests passed")


def test_sample_input_file():
    """Test parsing the actual input_param.txt from Fortran code."""
    
    # Verify input file exists
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}\n"
                                "Please copy from fortran/inputfile/input_param.txt")
    
    # Parse the actual file
    physics, simulation, laser, output = parse_input(str(INPUT_FILE))
    
    # Verify grid dimensions from geometry section
    # From input file: ncvx=[200], ncvy=[200], ncvz=[10,20] -> total 200, 200, 30
    assert simulation.ni == 200, f"Expected ni=200, got {simulation.ni}"
    assert simulation.nj == 200, f"Expected nj=200, got {simulation.nj}"
    assert simulation.nk == 30, f"Expected nk=30, got {simulation.nk}"
    
    # Verify domain lengths
    # xzone=[4.0e-3], yzone=[4.0e-3], zzone=[0.5e-3, 0.2e-3] -> 4mm, 4mm, 0.7mm
    assert abs(simulation.xlen - 4.0e-3) < 1e-10, f"Expected xlen=4mm, got {simulation.xlen}"
    assert abs(simulation.ylen - 4.0e-3) < 1e-10, f"Expected ylen=4mm, got {simulation.ylen}"
    assert abs(simulation.zlen - 0.7e-3) < 1e-10, f"Expected zlen=0.7mm, got {simulation.zlen}"
    
    # Verify material properties from &material_properties namelist
    assert physics.rho == 8440.0, f"Expected rho=8440, got {physics.rho}"
    assert physics.tsolid == 1563.0, f"Expected tsolid=1563, got {physics.tsolid}"
    assert physics.tliquid == 1623.0, f"Expected tliquid=1623, got {physics.tliquid}"
    assert physics.hsmelt == 861e3, f"Expected hsmelt=861e3, got {physics.hsmelt}"
    
    # Verify laser parameters from &volumetric_parameters namelist
    assert laser.power == 300.0, f"Expected power=300, got {laser.power}"
    assert abs(laser.absorptivity - 0.3512) < 1e-10, f"Expected absorptivity=0.3512, got {laser.absorptivity}"
    
    # Verify numerical parameters from &numerical_relax namelist
    assert abs(simulation.delt - 2e-5) < 1e-10, f"Expected delt=2e-5, got {simulation.delt}"
    assert simulation.urf_vel == 0.7, f"Expected urfu=0.7, got {simulation.urf_vel}"
    assert simulation.max_iter == 50, f"Expected maxit=50, got {simulation.max_iter}"
    
    print("✓ Actual input file (input_param.txt) parsing tests passed")


def test_geometry_parsing():
    """Test parsing of geometry section (zone-based grid)."""
    from param import _read_namelist
    
    if not INPUT_FILE.exists():
        print("⚠ Skipping geometry test - input file not found")
        return
    
    params = _read_namelist(str(INPUT_FILE))
    
    # Debug: print available keys
    geom_keys = ['nzx', 'nzy', 'nzz', 'xzone', 'yzone', 'zzone', 'ncvx', 'ncvy', 'ncvz', 'nx', 'ny', 'nz']
    missing = [k for k in geom_keys if k not in params]
    if missing:
        print(f"  Available keys: {sorted(params.keys())}")
        raise KeyError(f"Missing geometry keys: {missing}")
    
    # Check zone counts
    assert params['nzx'] == 1, f"Expected nzx=1, got {params['nzx']}"
    assert params['nzy'] == 1, f"Expected nzy=1, got {params['nzy']}"
    assert params['nzz'] == 2, f"Expected nzz=2, got {params['nzz']}"
    
    # Check zone lengths
    assert params['xzone'] == [4.0e-3], f"Expected xzone=[4e-3], got {params['xzone']}"
    assert params['yzone'] == [4.0e-3], f"Expected yzone=[4e-3], got {params['yzone']}"
    assert params['zzone'] == [0.5e-3, 0.2e-3], f"Expected zzone=[0.5e-3, 0.2e-3], got {params['zzone']}"
    
    # Check control volumes per zone
    assert params['ncvx'] == [200], f"Expected ncvx=[200], got {params['ncvx']}"
    assert params['ncvy'] == [200], f"Expected ncvy=[200], got {params['ncvy']}"
    assert params['ncvz'] == [10, 20], f"Expected ncvz=[10,20], got {params['ncvz']}"
    
    # Check computed totals
    assert params['nx'] == 200, f"Expected nx=200, got {params['nx']}"
    assert params['ny'] == 200, f"Expected ny=200, got {params['ny']}"
    assert params['nz'] == 30, f"Expected nz=30, got {params['nz']}"
    
    print("✓ Geometry section parsing tests passed")


def test_toolpath_loading():
    """Test loading toolpath from .crs file."""
    toolpath_file = Path(__file__).parent / "inputfile" / "B26-1.crs"
    
    if not toolpath_file.exists():
        print("⚠ Skipping toolpath test - B26-1.crs not found")
        return
    
    toolpath = load_toolpath(str(toolpath_file))
    
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


def run_all_tests():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("Testing Step 1: Types & Config")
    print("=" * 60 + "\n")
    
    test_parse_value()
    test_enthalpy_conversion()
    test_physics_params()
    test_simulation_params()
    test_laser_params()
    test_geometry_parsing()
    test_sample_input_file()
    test_toolpath_loading()
    
    print("\n" + "=" * 60)
    print("All Step 1 tests passed! ✓")
    print(f"Input file used: {INPUT_FILE}")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    run_all_tests()
