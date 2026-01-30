"""
AM-CFD JAX Implementation - Input/Output

Converted from Fortran modules: mod_param.f90, mod_print.f90
Parses namelist-style input files and handles output.
"""

import re
from pathlib import Path
from typing import Tuple, Dict, Any, Optional
import jax.numpy as jnp

from data_structures import (
    PhysicsParams, SimulationParams, LaserParams, OutputConfig, ToolPath
)


def parse_input(filepath: str) -> Tuple[PhysicsParams, SimulationParams, LaserParams, OutputConfig]:
    """
    Parse input parameter file (Fortran namelist format).
    
    Args:
        filepath: Path to input_param.txt
        
    Returns:
        Tuple of (PhysicsParams, SimulationParams, LaserParams, OutputConfig)
    """
    params = _read_namelist(filepath)
    
    # Extract physics parameters
    physics = _create_physics_params(params)
    
    # Extract simulation parameters
    simulation = _create_simulation_params(params)
    
    # Extract laser parameters
    laser = _create_laser_params(params)
    
    # Extract output configuration
    output = _create_output_config(params)
    
    return physics, simulation, laser, output


def _read_namelist(filepath: str) -> Dict[str, Any]:
    """
    Read Fortran input file (AM-CFD format).
    
    Handles the actual format from mod_param.f90:
        - First section: Line-by-line geometry data
        - Then namelists: &process_parameters, &material_properties, etc.
    """
    params = {}
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Parse geometry section (lines 1-13, before namelists)
    # Format from mod_param.f90 read_data subroutine
    line_idx = 0
    geometry_lines = []
    namelist_content = []
    
    for line in lines:
        stripped = line.strip()
        # Check if we've reached a namelist block
        if stripped.startswith('&'):
            namelist_content.append(line)
        elif namelist_content:
            namelist_content.append(line)
        else:
            geometry_lines.append(stripped)
    
    # Parse geometry data (line-by-line format)
    try:
        params.update(_parse_geometry_section(geometry_lines))
    except Exception as e:
        print(f"Warning: Could not parse geometry section: {e}")
    
    # Parse namelist blocks
    content = ''.join(namelist_content)
    
    # Remove inline comments
    content_lines = []
    for line in content.split('\n'):
        if '!' in line:
            line = line[:line.index('!')]
        content_lines.append(line)
    content = '\n'.join(content_lines)
    
    # Parse variable assignments from namelists
    # Fortran namelist format: &name var1=val1, var2=val2 /
    # Split by comma or whitespace, then parse each assignment
    # First remove namelist markers
    content = re.sub(r'&\w+', ' ', content)  # Remove &namelist_name
    content = re.sub(r'/', ' ', content)      # Remove /
    
    # Split into tokens by comma and whitespace, but keep = together with name/value
    # Match pattern: name=value where value can have scientific notation
    pattern = r'(\w+)\s*=\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?)'
    
    for match in re.finditer(pattern, content):
        name = match.group(1).lower()
        value_str = match.group(2).strip()
        
        # Try to parse the value
        params[name] = _parse_value(value_str)
    
    return params


def _parse_geometry_section(lines: list) -> Dict[str, Any]:
    """
    Parse the geometry section of input_param.txt.
    
    Based on mod_param.f90 read_data subroutine:
        Line 0: comment (skipped by read(10,*))
        Line 1: nzx (number of x-zones)
        Line 2: xzone values
        Line 3: ncvx values (control volumes per zone)
        Line 4: powrx values (exponents)
        Line 5: nzy
        Line 6: yzone values
        Line 7: ncvy values
        Line 8: powry values
        Line 9: nzz
        Line 10: zzone values
        Line 11: ncvz values
        Line 12: powrz values
    """
    params = {}
    
    # Filter out empty lines but keep comments for proper indexing
    # The Fortran code skips comment line with read(10,*)
    data_lines = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # Remove inline comments but keep the data part
        if '!' in line:
            # Check if there's data before the comment
            data_part = line[:line.index('!')].strip()
            if data_part:
                line = data_part
            else:
                continue  # Pure comment line, skip
        if line:
            data_lines.append(line)
    
    if len(data_lines) < 12:
        return params  # Not enough lines for geometry (need 12 data lines)
    
    idx = 0
    
    # X-direction
    params['nzx'] = int(data_lines[idx]); idx += 1
    params['xzone'] = [_parse_single_value(v) for v in data_lines[idx].split()]; idx += 1
    params['ncvx'] = [int(_parse_single_value(v)) for v in data_lines[idx].split()]; idx += 1
    params['powrx'] = [_parse_single_value(v) for v in data_lines[idx].split()]; idx += 1
    
    # Y-direction
    params['nzy'] = int(data_lines[idx]); idx += 1
    params['yzone'] = [_parse_single_value(v) for v in data_lines[idx].split()]; idx += 1
    params['ncvy'] = [int(_parse_single_value(v)) for v in data_lines[idx].split()]; idx += 1
    params['powry'] = [_parse_single_value(v) for v in data_lines[idx].split()]; idx += 1
    
    # Z-direction
    params['nzz'] = int(data_lines[idx]); idx += 1
    params['zzone'] = [_parse_single_value(v) for v in data_lines[idx].split()]; idx += 1
    params['ncvz'] = [int(_parse_single_value(v)) for v in data_lines[idx].split()]; idx += 1
    params['powrz'] = [_parse_single_value(v) for v in data_lines[idx].split()]; idx += 1
    
    # Compute total grid dimensions from zones
    params['nx'] = sum(params['ncvx'])
    params['ny'] = sum(params['ncvy'])
    params['nz'] = sum(params['ncvz'])
    
    # Compute total domain lengths
    params['xl'] = sum(params['xzone'])
    params['yl'] = sum(params['yzone'])
    params['zl'] = sum(params['zzone'])
    
    return params


def _parse_value(value_str: str) -> Any:
    """Parse a value string into appropriate Python type."""
    value_str = value_str.strip().rstrip(',')
    
    # Check for boolean
    if value_str.lower() in ['.true.', 'true', '.t.', 't']:
        return True
    if value_str.lower() in ['.false.', 'false', '.f.', 'f']:
        return False
    
    # Check for string (quoted)
    if (value_str.startswith("'") and value_str.endswith("'")) or \
       (value_str.startswith('"') and value_str.endswith('"')):
        return value_str[1:-1]
    
    # Check for array (comma-separated)
    if ',' in value_str:
        values = [_parse_single_value(v.strip()) for v in value_str.split(',')]
        return values
    
    return _parse_single_value(value_str)


def _parse_single_value(value_str: str) -> Any:
    """Parse a single value (int, float, or string)."""
    value_str = value_str.strip()
    
    # Handle Fortran double precision notation (1.0d0 -> 1.0e0)
    value_str = re.sub(r'([0-9.]+)[dD]([+-]?[0-9]+)', r'\1e\2', value_str)
    
    try:
        # Try integer first
        if '.' not in value_str and 'e' not in value_str.lower():
            return int(value_str)
        # Then float
        return float(value_str)
    except ValueError:
        return value_str


def _create_physics_params(params: Dict[str, Any]) -> PhysicsParams:
    """Create PhysicsParams from parsed dictionary.
    
    Maps Fortran variable names from mod_param.f90:
        dens -> rho, denl -> rholiq
        thconsa, thconsb -> thermal conductivity coefficients
        thconl -> liquid thermal conductivity
        hsmelt, hlfriz -> enthalpy at solidus/liquidus (precomputed in Fortran)
        dgdtp -> dgdt (surface tension temperature gradient)
    """
    
    # Specific heat coefficients (Cp = acpa*T + acpb for solid)
    acpa = _get_scalar(params, 'acpa', default=0.0)
    acpb = _get_scalar(params, 'acpb', default=500.0)
    acpl = _get_scalar(params, 'acpl', default=800.0)  # liquid specific heat
    
    # Temperatures
    tsolid = _get_scalar(params, 'tsolid', default=1563.0)
    tliquid = _get_scalar(params, 'tliquid', default=1623.0)
    tpreheat = _get_scalar(params, 'temppreheat', 'tpreheat', default=300.0)
    tvapor = _get_scalar(params, 'tboiling', 'tvapor', default=3000.0)
    
    # Enthalpy at phase boundaries
    # hsmelt and hlfriz may be provided directly from Fortran input
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
    
    # Density (Fortran uses dens/denl)
    rho = _get_scalar(params, 'dens', 'rho', default=7800.0)
    rholiq = _get_scalar(params, 'denl', 'rholiq', default=rho)
    
    # Thermal conductivity (Fortran uses thconsa, thconsb for solid, thconl for liquid)
    # k_solid = thconsa*T + thconsb
    thconsa = _get_scalar(params, 'thconsa', default=0.0)
    thconsb = _get_scalar(params, 'thconsb', default=25.0)
    tcond = thconsb  # Use constant term as representative solid conductivity
    tcondl = _get_scalar(params, 'thconl', 'tcondl', default=tcond)
    
    # Surface tension (Fortran uses dgdtp)
    dgdt = _get_scalar(params, 'dgdtp', 'dgdt', default=-4.3e-4)
    sigma_surf = _get_scalar(params, 'sigma_surf', default=1.8)
    
    # Radiation/convection
    emiss = _get_scalar(params, 'emiss', default=0.4)
    sigma = 5.67e-8  # Stefan-Boltzmann constant
    # Use htck1 (top surface) as representative convection coefficient
    hconv = _get_scalar(params, 'htckn', 'hconv', default=10.0)
    tenv = _get_scalar(params, 'tempamb', 'tenv', default=300.0)
    
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
    """Create SimulationParams from parsed dictionary.
    
    Uses geometry computed from zone data and Fortran namelist values.
    """
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
        stretch_x=params.get('powrx', params.get('stretch_x', params.get('stretchx', 1.0))),
        stretch_y=params.get('powry', params.get('stretch_y', params.get('stretchy', 1.0))),
        stretch_z=params.get('powrz', params.get('stretch_z', params.get('stretchz', 1.0))),
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
    """Create LaserParams from parsed dictionary.
    
    Supports both surface and volumetric laser sources from Fortran:
        Surface: alaspow, alaseta, alasrb
        Volumetric: alaspowvol, alasetavol, sourcerad, sourcedepth
    """
    # Try volumetric parameters first (more commonly used in input file)
    power = _get_scalar(params, 'alaspowvol', 'alaspow', 'power', 'plaser', default=200.0)
    absorptivity = _get_scalar(params, 'alasetavol', 'alaseta', 'absorptivity', 'absorp', default=0.35)
    radius = _get_scalar(params, 'sourcerad', 'alasrb', 'radius', 'rb', default=50.0e-6)
    efficiency = _get_scalar(params, 'alasfact', 'efficiency', 'eff', default=1.0)
    
    # Compute peak heat flux for Gaussian distribution
    # q(r) = q_peak * exp(-2*r^2/rb^2)
    # Total power = integral(q * 2*pi*r dr) = q_peak * pi * rb^2 / 2
    # q_peak = 2 * P * absorptivity / (pi * rb^2)
    import math
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


def load_toolpath(filepath: str) -> ToolPath:
    """
    Load toolpath from .crs file.
    
    Format: 5 columns - time, x, y, z, laser_on
    
    Args:
        filepath: Path to .crs toolpath file
        
    Returns:
        ToolPath NamedTuple
        
    Note: Uses numpy arrays for fast loading. JAX can use these directly.
    """
    import numpy as np
    
    data = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            
            parts = line.split()
            if len(parts) >= 5:
                try:
                    row = [float(p) for p in parts[:5]]
                    data.append(row)
                except ValueError:
                    continue
    
    if not data:
        # Return empty toolpath
        return ToolPath(
            time=np.array([0.0]),
            x=np.array([0.0]),
            y=np.array([0.0]),
            z=np.array([0.0]),
            laser_on=np.array([0], dtype=np.int32),
            n_segments=1,
        )
    
    data = np.array(data)
    
    return ToolPath(
        time=data[:, 0],
        x=data[:, 1],
        y=data[:, 2],
        z=data[:, 3],
        laser_on=data[:, 4].astype(np.int32),
        n_segments=len(data),
    )


def write_output_header(filepath: str, physics: PhysicsParams, 
                        simulation: SimulationParams, laser: LaserParams) -> None:
    """Write output file header with simulation parameters."""
    with open(filepath, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("AM-CFD Simulation Output (JAX Implementation)\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("SIMULATION PARAMETERS:\n")
        f.write(f"  Grid: {simulation.ni} x {simulation.nj} x {simulation.nk}\n")
        f.write(f"  Domain: {simulation.xlen*1e3:.3f} x {simulation.ylen*1e3:.3f} x {simulation.zlen*1e3:.3f} mm\n")
        f.write(f"  Timestep: {simulation.delt*1e6:.3f} µs\n")
        f.write(f"  Max time: {simulation.timax*1e3:.3f} ms\n\n")
        
        f.write("MATERIAL PROPERTIES:\n")
        f.write(f"  Density: {physics.rho:.1f} kg/m³\n")
        f.write(f"  Solidus: {physics.tsolid:.1f} K\n")
        f.write(f"  Liquidus: {physics.tliquid:.1f} K\n")
        f.write(f"  Latent heat: {physics.hlatent/1e3:.1f} kJ/kg\n\n")
        
        f.write("LASER PARAMETERS:\n")
        f.write(f"  Power: {laser.power:.1f} W\n")
        f.write(f"  Radius: {laser.radius*1e6:.1f} µm\n")
        f.write(f"  Absorptivity: {laser.absorptivity:.2f}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write(f"{'Time':>12} {'Iter':>6} {'Residual':>12} {'MaxT':>10} "
                f"{'Length':>10} {'Width':>10} {'Depth':>10} {'Ratio':>8}\n")
        f.write(f"{'[ms]':>12} {'':>6} {'':>12} {'[K]':>10} "
                f"{'[µm]':>10} {'[µm]':>10} {'[µm]':>10} {'':>8}\n")
        f.write("-" * 80 + "\n")


def write_output_line(filepath: str, time: float, iteration: int,
                      residual: float, max_temp: float,
                      pool_length: float, pool_width: float, pool_depth: float,
                      heat_ratio: float) -> None:
    """Append a line to the output file."""
    with open(filepath, 'a') as f:
        f.write(f"{time*1e3:12.6f} {iteration:6d} {residual:12.4e} {max_temp:10.1f} "
                f"{pool_length*1e6:10.2f} {pool_width*1e6:10.2f} {pool_depth*1e6:10.2f} "
                f"{heat_ratio:8.4f}\n")


def write_tecplot(filepath: str, state, grid, time: float) -> None:
    """
    Write Tecplot-format output file.
    
    Args:
        filepath: Output file path
        state: State NamedTuple
        grid: GridParams NamedTuple  
        time: Current simulation time
    """
    ni, nj, nk = grid.ni, grid.nj, grid.nk
    
    with open(filepath, 'w') as f:
        f.write(f'TITLE = "AM-CFD t={time*1e6:.2f}us"\n')
        f.write('VARIABLES = "X" "Y" "Z" "T" "H" "fracl" "U" "V" "W" "P"\n')
        f.write(f'ZONE T="Zone1", I={ni}, J={nj}, K={nk}, F=POINT\n')
        
        # Convert JAX arrays to numpy for iteration
        x = jnp.array(grid.x)
        y = jnp.array(grid.y)
        z = jnp.array(grid.z)
        temp = jnp.array(state.temp)
        enthalpy = jnp.array(state.enthalpy)
        fracl = jnp.array(state.fracl)
        uVel = jnp.array(state.uVel)
        vVel = jnp.array(state.vVel)
        wVel = jnp.array(state.wVel)
        pressure = jnp.array(state.pressure)
        
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    f.write(f"{float(x[i]):14.6e} {float(y[j]):14.6e} {float(z[k]):14.6e} "
                            f"{float(temp[i,j,k]):14.6e} {float(enthalpy[i,j,k]):14.6e} "
                            f"{float(fracl[i,j,k]):14.6e} {float(uVel[i,j,k]):14.6e} "
                            f"{float(vVel[i,j,k]):14.6e} {float(wVel[i,j,k]):14.6e} "
                            f"{float(pressure[i,j,k]):14.6e}\n")
