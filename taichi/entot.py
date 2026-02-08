"""
AM-CFD Taichi Implementation - Enthalpy to Temperature Conversion

Translated from Fortran mod_entot.f90
Converts enthalpy field to temperature and computes liquid fraction.
"""

import taichi as ti
import math

from data_structures import State, PhysicsParams, GridParams


@ti.kernel
def enthalpy_to_temp_kernel(
    enthalpy: ti.template(),
    temp: ti.template(),
    fracl: ti.template(),
    ni: ti.i32, nj: ti.i32, nk: ti.i32,
    hlcal: ti.f64, hsmelt: ti.f64,
    tliquid: ti.f64, tsolid: ti.f64,
    acpl: ti.f64, acpa: ti.f64, acpb: ti.f64
):
    """
    Taichi kernel to convert enthalpy to temperature.
    
    Three regions:
    1. Liquid (H >= hlcal): T = (H - hlcal)/acpl + tliquid, fracl = 1
    2. Solid (H <= hsmelt): T from inverse Cp integral, fracl = 0
    3. Mushy (between): Linear interpolation of fracl and T
    """
    deltemp = tliquid - tsolid
    
    for i, j, k in ti.ndrange(ni, nj, nk):
        H = enthalpy[i, j, k]
        
        if H >= hlcal:
            # Liquid region: fully melted
            fracl[i, j, k] = 1.0
            temp[i, j, k] = (H - hlcal) / acpl + tliquid
            
        elif H <= hsmelt:
            # Solid region: fully solid
            fracl[i, j, k] = 0.0
            
            # Inverse of enthalpy integral: H = acpa/2 * T^2 + acpb * T
            # Solving quadratic: acpa/2 * T^2 + acpb * T - H = 0
            # T = (-acpb + sqrt(acpb^2 + 2*acpa*H)) / acpa
            if ti.abs(acpa) > 1.0e-12:
                discriminant = acpb * acpb + 2.0 * acpa * H
                if discriminant > 0.0:
                    temp[i, j, k] = (ti.sqrt(discriminant) - acpb) / acpa
                else:
                    temp[i, j, k] = H / acpb  # Fallback for edge case
            else:
                # Linear Cp (acpa = 0): H = acpb * T => T = H / acpb
                temp[i, j, k] = H / acpb
                
        else:
            # Mushy zone: phase change region
            # Linear interpolation of liquid fraction
            fracl[i, j, k] = (H - hsmelt) / (hlcal - hsmelt)
            temp[i, j, k] = deltemp * fracl[i, j, k] + tsolid


def enthalpy_to_temp(state: State, physics: PhysicsParams) -> State:
    """
    Convert enthalpy field to temperature. (From mod_entot.f90)
    
    Translates enthalpy to temperature, handling phase change in the mushy zone.
    Updates both state.temp and state.fracl fields.
    
    Args:
        state: Current state containing enthalpy field
        physics: Physics parameters with phase change properties
        
    Returns:
        state: Updated state with temperature and liquid fraction fields
    """
    ni, nj, nk = state.ni, state.nj, state.nk
    
    enthalpy_to_temp_kernel(
        state.enthalpy,
        state.temp,
        state.fracl,
        ni, nj, nk,
        physics.hlcal,
        physics.hsmelt,
        physics.tliquid,
        physics.tsolid,
        physics.acpl,
        physics.acpa,
        physics.acpb
    )
    
    return state


def temp_to_enthalpy_solid(T: float, acpa: float, acpb: float) -> float:
    """
    Convert temperature to enthalpy for solid phase.
    
    For solid: Cp = acpa*T + acpb
    H = integral(Cp dT) = acpa/2 * T^2 + acpb * T
    
    Args:
        T: Temperature [K]
        acpa: Cp coefficient a [J/(kg·K²)]
        acpb: Cp coefficient b [J/(kg·K)]
        
    Returns:
        Enthalpy [J/kg]
    """
    return 0.5 * acpa * T * T + acpb * T


def temp_to_enthalpy(T: float, physics: PhysicsParams) -> float:
    """
    Convert temperature to enthalpy, handling all phases.
    
    Args:
        T: Temperature [K]
        physics: Physics parameters
        
    Returns:
        Enthalpy [J/kg]
    """
    if T <= physics.tsolid:
        # Solid region
        return temp_to_enthalpy_solid(T, physics.acpa, physics.acpb)
    elif T >= physics.tliquid:
        # Liquid region
        return physics.hlcal + physics.acpl * (T - physics.tliquid)
    else:
        # Mushy zone - linear interpolation
        fracl = (T - physics.tsolid) / (physics.tliquid - physics.tsolid)
        return physics.hsmelt + fracl * (physics.hlcal - physics.hsmelt)


@ti.kernel
def temp_to_enthalpy_kernel(
    temp: ti.template(),
    enthalpy: ti.template(),
    fracl: ti.template(),
    ni: ti.i32, nj: ti.i32, nk: ti.i32,
    hlcal: ti.f64, hsmelt: ti.f64,
    tliquid: ti.f64, tsolid: ti.f64,
    acpl: ti.f64, acpa: ti.f64, acpb: ti.f64
):
    """
    Taichi kernel to convert temperature to enthalpy (inverse operation).
    """
    deltemp = tliquid - tsolid
    
    for i, j, k in ti.ndrange(ni, nj, nk):
        T = temp[i, j, k]
        
        if T >= tliquid:
            # Liquid region
            fracl[i, j, k] = 1.0
            enthalpy[i, j, k] = hlcal + acpl * (T - tliquid)
            
        elif T <= tsolid:
            # Solid region
            fracl[i, j, k] = 0.0
            enthalpy[i, j, k] = 0.5 * acpa * T * T + acpb * T
            
        else:
            # Mushy zone
            fracl[i, j, k] = (T - tsolid) / deltemp
            enthalpy[i, j, k] = hsmelt + fracl[i, j, k] * (hlcal - hsmelt)
