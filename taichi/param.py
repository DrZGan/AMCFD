"""
AM-CFD Taichi Implementation - Input/Output

Parses YAML input files and handles output.
"""

import yaml
from pathlib import Path
from typing import Tuple, Dict, Any, List, Union
import numpy as np

from data_structures import (
    PhysicsParams, SimulationParams, LaserParams, OutputConfig, ToolPath
)


def parse_input(filepath: str) -> Tuple[PhysicsParams, SimulationParams, LaserParams, OutputConfig]:
    """
    Parse YAML input parameter file.
    
    Args:
        filepath: Path to input_param.yaml
        
    Returns:
        Tuple of (PhysicsParams, SimulationParams, LaserParams, OutputConfig)
    """
    params = _read_yaml(filepath)
    
    # Extract physics parameters
    physics = _create_physics_params(params)
    
    # Extract simulation parameters
    simulation = _create_simulation_params(params)
    
    # Extract laser parameters
    laser = _create_laser_params(params)
    
    # Extract output configuration
    output = _create_output_config(params)
    
    return physics, simulation, laser, output


def _convert_value(val: Any) -> Any:
    """Convert YAML value to proper numeric type if possible."""
    if isinstance(val, str):
        # Try to convert scientific notation strings to float
        try:
            return float(val)
        except ValueError:
            return val
    elif isinstance(val, list):
        return [_convert_value(v) for v in val]
    elif isinstance(val, dict):
        return {k: _convert_value(v) for k, v in val.items()}
    return val


def _read_yaml(filepath: str) -> Dict[str, Any]:
    """
    Read YAML input file and flatten to parameter dictionary.
    
    Args:
        filepath: Path to YAML file
        
    Returns:
        Flattened dictionary of parameters
    """
    with open(filepath, 'r') as f:
        config = yaml.safe_load(f)
    
    # Convert any string numbers to proper numeric types
    config = _convert_value(config)
    
    # Flatten nested structure for easier access
    params = {}
    
    # Parse geometry section
    if 'geometry' in config:
        geom = config['geometry']
        params.update(_parse_geometry_yaml(geom))
    
    # Parse process parameters
    if 'process_parameters' in config:
        params.update(config['process_parameters'])
    
    # Parse volumetric parameters
    if 'volumetric_parameters' in config:
        params.update(config['volumetric_parameters'])
    
    # Parse material properties
    if 'material_properties' in config:
        params.update(config['material_properties'])
    
    # Parse powder properties
    if 'powder_properties' in config:
        params.update(config['powder_properties'])
    
    # Parse numerical relaxation
    if 'numerical_relax' in config:
        params.update(config['numerical_relax'])
    
    # Parse boundary conditions
    if 'boundary_conditions' in config:
        bc = config['boundary_conditions']
        params.update(bc)
        # Map boundary condition names
        if 'tempPreheat' in bc:
            params['temppreheat'] = bc['tempPreheat']
        if 'tempAmb' in bc:
            params['tempamb'] = bc['tempAmb']
        if 'htckn' in bc:
            params['htckn'] = bc['htckn']
    
    # Parse output config if present
    if 'output' in config:
        params.update(config['output'])
    
    return params


def _parse_geometry_yaml(geom: Dict[str, Any]) -> Dict[str, Any]:
    """
    Parse geometry section from YAML structure.
    
    Args:
        geom: Geometry dictionary from YAML
        
    Returns:
        Flattened geometry parameters
    """
    params = {}
    
    # Helper to ensure list format
    def to_list(val):
        if isinstance(val, list):
            return val
        return [val]
    
    # X-direction
    if 'x' in geom:
        x = geom['x']
        params['nzx'] = x.get('zones', 1)
        params['xzone'] = to_list(x.get('zone_length_m', 1.0e-3))
        params['ncvx'] = to_list(x.get('cv_per_zone', 100))
        params['powrx'] = to_list(x.get('cv_boundary_exponent', 1.0))
    
    # Y-direction
    if 'y' in geom:
        y = geom['y']
        params['nzy'] = y.get('zones', 1)
        params['yzone'] = to_list(y.get('zone_length_m', 1.0e-3))
        params['ncvy'] = to_list(y.get('cv_per_zone', 100))
        params['powry'] = to_list(y.get('cv_boundary_exponent', 1.0))
    
    # Z-direction
    if 'z' in geom:
        z = geom['z']
        params['nzz'] = z.get('zones', 1)
        params['zzone'] = to_list(z.get('zone_length_m', 0.5e-3))
        params['ncvz'] = to_list(z.get('cv_per_zone', 50))
        params['powrz'] = to_list(z.get('cv_boundary_exponent', 1.0))
    
    # Compute total grid dimensions from zones
    params['nx'] = sum(params.get('ncvx', [100]))
    params['ny'] = sum(params.get('ncvy', [100]))
    params['nz'] = sum(params.get('ncvz', [50]))
    
    # Compute total domain lengths
    params['xl'] = sum(params.get('xzone', [1.0e-3]))
    params['yl'] = sum(params.get('yzone', [1.0e-3]))
    params['zl'] = sum(params.get('zzone', [0.5e-3]))
    
    return params


def _create_physics_params(params: Dict[str, Any]) -> PhysicsParams:
    """Create PhysicsParams from parsed dictionary.
    
    Maps YAML variable names to PhysicsParams fields.
    """
    
    # Specific heat coefficients (Cp = acpa*T + acpb for solid)
    acpa = _get_scalar(params, 'acpa', default=0.0)
    acpb = _get_scalar(params, 'acpb', default=500.0)
    acpl = _get_scalar(params, 'acpl', default=800.0)  # liquid specific heat
    
    # Temperatures
    tsolid = _get_scalar(params, 'tsolid', default=1563.0)
    tliquid = _get_scalar(params, 'tliquid', default=1623.0)
    tpreheat = _get_scalar(params, 'temppreheat', 'tpreheat', 'tempPreheat', default=300.0)
    tvapor = _get_scalar(params, 'tboiling', 'tvapor', default=3000.0)
    
    # Enthalpy at phase boundaries
    if 'hsmelt' in params and 'hlfriz' in params:
        hsmelt = _get_scalar(params, 'hsmelt', default=0.0)
        hlcal = _get_scalar(params, 'hlfriz', default=0.0)
        hlatent = hlcal - hsmelt
    else:
        # Compute from temperature if not provided
        hlatent = _get_scalar(params, 'hlatent', default=2.7e5)
        hsmelt = _temp_to_enthalpy_solid(tsolid, acpa, acpb)
        hlcal = hsmelt + hlatent
    
    hpreheat = _temp_to_enthalpy_solid(tpreheat, acpa, acpb)
    
    # Density
    rho = _get_scalar(params, 'dens', 'rho', default=7800.0)
    rholiq = _get_scalar(params, 'denl', 'rholiq', default=rho)
    
    # Thermal conductivity
    thconsa = _get_scalar(params, 'thconsa', default=0.0)
    thconsb = _get_scalar(params, 'thconsb', default=25.0)
    tcond = thconsb  # Use constant term as representative solid conductivity
    tcondl = _get_scalar(params, 'thconl', 'tcondl', default=tcond)
    
    # Surface tension
    dgdt = _get_scalar(params, 'dgdtp', 'dgdt', default=-4.3e-4)
    sigma_surf = _get_scalar(params, 'sigma_surf', default=1.8)
    
    # Radiation/convection
    emiss = _get_scalar(params, 'emiss', default=0.4)
    sigma = 5.67e-8  # Stefan-Boltzmann constant
    hconv = _get_scalar(params, 'htckn', 'hconv', default=10.0)
    tenv = _get_scalar(params, 'tempamb', 'tempAmb', 'tenv', default=300.0)
    
    # Viscosity
    vis0 = _get_scalar(params, 'viscos', 'vis0', default=6.0e-3)
    
    # Gravity/buoyancy
    grav = params.get('grav', 9.81)
    beta = params.get('beta', 1.0e-4)
    
    # Mushy zone
    darcy_coeff = params.get('darcy_coeff', 1.0e6)
    
    return PhysicsParams(
        acpa=acpa,
        acpb=acpb,
        acpl=acpl,
        tsolid=tsolid,
        tliquid=tliquid,
        tpreheat=tpreheat,
        tvapor=tvapor,
        hsmelt=hsmelt,
        hlcal=hlcal,
        hpreheat=hpreheat,
        hlatent=hlatent,
        rho=rho,
        rholiq=rholiq,
        tcond=tcond,
        tcondl=tcondl,
        dgdt=dgdt,
        sigma_surf=sigma_surf,
        emiss=emiss,
        sigma=sigma,
        hconv=hconv,
        tenv=tenv,
        vis0=vis0,
        grav=grav,
        beta=beta,
        darcy_coeff=darcy_coeff,
    )


def _temp_to_enthalpy_solid(T: float, acpa: float, acpb: float) -> float:
    """Convert temperature to enthalpy for solid phase (quadratic Cp)."""
    # Cp = acpa*T + acpb
    # H = integral(Cp dT) = acpa/2 * T^2 + acpb * T
    return 0.5 * acpa * T * T + acpb * T


def _create_simulation_params(params: Dict[str, Any]) -> SimulationParams:
    """Create SimulationParams from parsed dictionary."""
    return SimulationParams(
        delt=params.get('delt', 1.0e-6),
        timax=params.get('timax', 1.0e-3),
        urf_vel=params.get('urfu', params.get('urf_vel', 0.5)),
        urf_p=params.get('urfp', params.get('urf_p', 0.3)),
        urf_h=params.get('urfh', params.get('urf_h', 0.8)),
        max_iter=params.get('maxit', params.get('max_iter', params.get('niter', 1000))),
        conv_tol=params.get('conv_tol', 5.0e-4),
        ni=params.get('nx', params.get('ni', 100)),
        nj=params.get('ny', params.get('nj', 100)),
        nk=params.get('nz', params.get('nk', 50)),
        xlen=params.get('xl', params.get('xlen', 1.0e-3)),
        ylen=params.get('yl', params.get('ylen', 1.0e-3)),
        zlen=params.get('zl', params.get('zlen', 0.5e-3)),
        stretch_x=_get_scalar(params, 'powrx', 'stretch_x', 'stretchx', default=1.0),
        stretch_y=_get_scalar(params, 'powry', 'stretch_y', 'stretchy', default=1.0),
        stretch_z=_get_scalar(params, 'powrz', 'stretch_z', 'stretchz', default=1.0),
    )


def _get_scalar(params: Dict[str, Any], *keys, default=None):
    """Get a scalar value from params, handling lists by taking first element."""
    for key in keys:
        if key in params:
            val = params[key]
            if isinstance(val, list):
                return val[0] if val else default
            return val
    return default


def _create_laser_params(params: Dict[str, Any]) -> LaserParams:
    """Create LaserParams from parsed dictionary."""
    import math
    
    # Try volumetric parameters first (more commonly used)
    power = _get_scalar(params, 'alaspowvol', 'alaspow', 'power', 'plaser', default=200.0)
    absorptivity = _get_scalar(params, 'alasetavol', 'alaseta', 'absorptivity', 'absorp', default=0.35)
    radius = _get_scalar(params, 'sourcerad', 'alasrb', 'radius', 'rb', default=50.0e-6)
    efficiency = _get_scalar(params, 'alasfact', 'efficiency', 'eff', default=1.0)
    
    # Compute peak heat flux for Gaussian distribution
    peak_flux = 2.0 * power * absorptivity * efficiency / (math.pi * radius * radius)
    
    return LaserParams(
        power=power,
        radius=radius,
        absorptivity=absorptivity,
        efficiency=efficiency,
        peak_flux=peak_flux,
    )


def _create_output_config(params: Dict[str, Any]) -> OutputConfig:
    """Create OutputConfig from parsed dictionary."""
    return OutputConfig(
        output_interval=params.get('output_interval', params.get('nprint', 100)),
        tecplot_interval=params.get('tecplot_interval', params.get('ntec', 1000)),
        output_dir=params.get('output_dir', './output'),
        case_name=params.get('case_name', 'amcfd'),
    )


# Note: load_toolpath has been moved to toolpath.py
# For backwards compatibility, re-export it here
from toolpath import load_toolpath
