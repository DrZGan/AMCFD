"""
Tests for geom.py - Grid generation and geometry.
"""

import sys
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

# Initialize Taichi before importing
import taichi as ti
ti.init(arch=ti.cpu)

from data_structures import GridParams
from geom import get_gridparams

# Path to YAML input file
INPUT_FILE = Path(__file__).parent / "inputfile" / "input_param.yaml"


def test_get_gridparams():
    """Test grid parameter extraction from YAML."""
    if not INPUT_FILE.exists():
        print("⚠ Skipping geom test - input file not found")
        return
    
    grid = get_gridparams(str(INPUT_FILE))
    
    # Check dimensions from YAML (includes +2 ghost cells for staggered grid)
    # cv_per_zone: x=200, y=200, z=[10,20] -> 202, 202, 32 (with ghost cells)
    assert grid.ni == 202, f"Expected ni=202 (200+2 ghost), got {grid.ni}"
    assert grid.nj == 202, f"Expected nj=202 (200+2 ghost), got {grid.nj}"
    assert grid.nk == 32, f"Expected nk=32 (30+2 ghost), got {grid.nk}"
    
    # Check field shapes
    assert grid.x.shape == (202,)
    assert grid.y.shape == (202,)
    assert grid.z.shape == (32,)
    assert grid.vol.shape == (202, 202, 32)
    
    print("✓ Grid parameters extraction tests passed")


def test_grid_field_shapes():
    """Test that grid fields have correct shapes."""
    ni, nj, nk = 20, 20, 10
    grid = GridParams(ni, nj, nk)
    
    # 1D coordinate arrays
    assert grid.x.shape == (ni,)
    assert grid.y.shape == (nj,)
    assert grid.z.shape == (nk,)
    
    # Staggered grid locations
    assert grid.xu.shape == (ni,)
    assert grid.yv.shape == (nj,)
    assert grid.zw.shape == (nk,)
    
    # Cell dimensions
    assert grid.dx.shape == (ni,)
    assert grid.dy.shape == (nj,)
    assert grid.dz.shape == (nk,)
    
    # Inverse distances
    assert grid.dxpwinv.shape == (ni,)
    assert grid.dxpeinv.shape == (ni,)
    assert grid.dypsinv.shape == (nj,)
    assert grid.dypninv.shape == (nj,)
    assert grid.dzpbinv.shape == (nk,)
    assert grid.dzptinv.shape == (nk,)
    
    # 3D volume array
    assert grid.vol.shape == (ni, nj, nk)
    
    # 2D area arrays
    assert grid.areaij.shape == (ni, nj)
    assert grid.areaik.shape == (ni, nk)
    assert grid.areajk.shape == (nj, nk)
    
    print("✓ Grid field shapes tests passed")


def run_all_tests():
    """Run all geom tests."""
    print("\n" + "=" * 60)
    print("Testing geom.py - Grid Generation")
    print("=" * 60 + "\n")
    
    test_grid_field_shapes()
    test_get_gridparams()
    
    print("\n" + "=" * 60)
    print("All geom.py tests passed! ✓")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    run_all_tests()
