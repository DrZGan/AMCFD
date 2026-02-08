"""
Tests for entot.py - Enthalpy to Temperature conversion.

Test Summary:
-----------------------------------------------------------------------------
| Test                        | Purpose                                      |
-----------------------------------------------------------------------------
| test_temp_to_enthalpy_solid | Tests solid phase enthalpy calculation       |
| test_temp_to_enthalpy_all_  | Tests enthalpy calculation for all phases    |
|   phases                    | (solid, mushy, liquid)                       |
| test_enthalpy_to_temp_kernel| Tests Taichi kernel for H->T conversion      |
| test_roundtrip_conversion   | Tests T->H->T roundtrip consistency          |
-----------------------------------------------------------------------------
"""

import sys
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

# Initialize Taichi before importing
import taichi as ti
ti.init(arch=ti.cpu)

from data_structures import State, PhysicsParams
from entot import (
    enthalpy_to_temp, temp_to_enthalpy_solid, temp_to_enthalpy,
    enthalpy_to_temp_kernel
)


def test_temp_to_enthalpy_solid():
    """Test solid phase enthalpy calculation."""
    # Linear Cp case (acpa = 0)
    acpa = 0.0
    acpb = 500.0
    T = 1000.0
    
    H = temp_to_enthalpy_solid(T, acpa, acpb)
    expected = acpb * T  # 500 * 1000 = 500000
    assert abs(H - expected) < 1e-10, f"Expected {expected}, got {H}"
    
    # Quadratic Cp case
    acpa = 0.5
    acpb = 500.0
    T = 1000.0
    
    H = temp_to_enthalpy_solid(T, acpa, acpb)
    expected = 0.5 * acpa * T * T + acpb * T  # 0.25*1e6 + 500*1000 = 750000
    assert abs(H - expected) < 1e-10, f"Expected {expected}, got {H}"
    
    print("✓ temp_to_enthalpy_solid tests passed")


def test_temp_to_enthalpy_all_phases():
    """Test enthalpy calculation for all phases."""
    physics = PhysicsParams(
        acpa=0.0,
        acpb=500.0,
        acpl=800.0,
        tsolid=1563.0,
        tliquid=1623.0,
        hsmelt=781500.0,  # 500 * 1563
        hlcal=1051500.0,  # hsmelt + latent heat (270000)
        hlatent=270000.0
    )
    
    # Solid region
    T_solid = 1000.0
    H = temp_to_enthalpy(T_solid, physics)
    expected = 500.0 * 1000.0
    assert abs(H - expected) < 1e-6, f"Solid: expected {expected}, got {H}"
    
    # Liquid region
    T_liquid = 1700.0
    H = temp_to_enthalpy(T_liquid, physics)
    expected = physics.hlcal + physics.acpl * (T_liquid - physics.tliquid)
    assert abs(H - expected) < 1e-6, f"Liquid: expected {expected}, got {H}"
    
    # Mushy zone (midpoint)
    T_mushy = (physics.tsolid + physics.tliquid) / 2.0
    H = temp_to_enthalpy(T_mushy, physics)
    expected = (physics.hsmelt + physics.hlcal) / 2.0  # Midpoint enthalpy
    assert abs(H - expected) < 1e-6, f"Mushy: expected {expected}, got {H}"
    
    print("✓ temp_to_enthalpy all phases tests passed")


def test_enthalpy_to_temp_kernel():
    """Test the Taichi kernel for enthalpy to temperature conversion."""
    ni, nj, nk = 5, 5, 3
    
    physics = PhysicsParams(
        acpa=0.0,
        acpb=500.0,
        acpl=800.0,
        tsolid=1563.0,
        tliquid=1623.0,
        hsmelt=781500.0,
        hlcal=1051500.0,
        hlatent=270000.0
    )
    
    state = State(ni, nj, nk)
    
    # Fill with known enthalpy values
    @ti.kernel
    def fill_enthalpy(enthalpy: ti.template(), value: ti.f64):
        for i, j, k in enthalpy:
            enthalpy[i, j, k] = value
    
    # Test solid region
    H_solid = 500000.0  # Should give T = 1000 K
    fill_enthalpy(state.enthalpy, H_solid)
    enthalpy_to_temp(state, physics)
    
    temp_np = state.temp.to_numpy()
    fracl_np = state.fracl.to_numpy()
    
    expected_T = 1000.0  # H / acpb = 500000 / 500
    assert abs(temp_np[0, 0, 0] - expected_T) < 1.0, f"Solid T: expected ~{expected_T}, got {temp_np[0,0,0]}"
    assert fracl_np[0, 0, 0] == 0.0, f"Solid fracl: expected 0, got {fracl_np[0,0,0]}"
    
    # Test liquid region
    H_liquid = physics.hlcal + 800.0 * 100.0  # 100K above liquidus
    fill_enthalpy(state.enthalpy, H_liquid)
    enthalpy_to_temp(state, physics)
    
    temp_np = state.temp.to_numpy()
    fracl_np = state.fracl.to_numpy()
    
    expected_T = physics.tliquid + 100.0
    assert abs(temp_np[0, 0, 0] - expected_T) < 1.0, f"Liquid T: expected ~{expected_T}, got {temp_np[0,0,0]}"
    assert fracl_np[0, 0, 0] == 1.0, f"Liquid fracl: expected 1, got {fracl_np[0,0,0]}"
    
    # Test mushy zone (midpoint)
    H_mushy = (physics.hsmelt + physics.hlcal) / 2.0
    fill_enthalpy(state.enthalpy, H_mushy)
    enthalpy_to_temp(state, physics)
    
    temp_np = state.temp.to_numpy()
    fracl_np = state.fracl.to_numpy()
    
    expected_fracl = 0.5
    expected_T = physics.tsolid + 0.5 * (physics.tliquid - physics.tsolid)
    assert abs(fracl_np[0, 0, 0] - expected_fracl) < 0.01, f"Mushy fracl: expected ~{expected_fracl}, got {fracl_np[0,0,0]}"
    assert abs(temp_np[0, 0, 0] - expected_T) < 1.0, f"Mushy T: expected ~{expected_T}, got {temp_np[0,0,0]}"
    
    print("✓ enthalpy_to_temp kernel tests passed")


def test_roundtrip_conversion():
    """Test that temp->enthalpy->temp gives back original temperature."""
    physics = PhysicsParams(
        acpa=0.0,
        acpb=500.0,
        acpl=800.0,
        tsolid=1563.0,
        tliquid=1623.0,
        hsmelt=781500.0,
        hlcal=1051500.0,
        hlatent=270000.0
    )
    
    ni, nj, nk = 3, 3, 3
    state = State(ni, nj, nk)
    
    # Test temperatures in each region
    test_temps = [300.0, 1000.0, 1563.0, 1593.0, 1623.0, 1800.0]
    
    for T_original in test_temps:
        # Convert T -> H
        H = temp_to_enthalpy(T_original, physics)
        
        # Set enthalpy field
        @ti.kernel
        def fill_h(enthalpy: ti.template(), val: ti.f64):
            for i, j, k in enthalpy:
                enthalpy[i, j, k] = val
        
        fill_h(state.enthalpy, H)
        
        # Convert H -> T
        enthalpy_to_temp(state, physics)
        
        T_result = state.temp.to_numpy()[0, 0, 0]
        
        # Check roundtrip accuracy
        error = abs(T_result - T_original)
        assert error < 1.0, f"Roundtrip failed for T={T_original}: got {T_result}, error={error}"
    
    print("✓ Roundtrip conversion tests passed")


def run_all_tests():
    """Run all entot tests."""
    print("\n" + "=" * 60)
    print("Testing entot.py - Enthalpy/Temperature Conversion")
    print("=" * 60 + "\n")
    
    test_temp_to_enthalpy_solid()
    test_temp_to_enthalpy_all_phases()
    test_enthalpy_to_temp_kernel()
    test_roundtrip_conversion()
    
    print("\n" + "=" * 60)
    print("All entot.py tests passed! ✓")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    run_all_tests()
