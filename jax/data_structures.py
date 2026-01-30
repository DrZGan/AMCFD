"""
AM-CFD JAX Implementation - Type Definitions

Converted from Fortran modules: mod_const.f90, mod_param.f90, and variable declarations
across all modules.
"""

from typing import NamedTuple, Optional
import jax.numpy as jnp
from jax import Array


class PhysicsParams(NamedTuple):
    """Physical constants and material properties (from mod_const.f90, mod_param.f90)"""
    # Specific heat coefficients (solid: Cp = acpa*T + acpb)
    acpa: float          # Solid Cp coefficient a [J/(kg·K²)]
    acpb: float          # Solid Cp coefficient b [J/(kg·K)]
    acpl: float          # Liquid specific heat [J/(kg·K)]
    
    # Phase change temperatures
    tsolid: float        # Solidus temperature [K]
    tliquid: float       # Liquidus temperature [K]
    tpreheat: float      # Preheat temperature [K]
    tvapor: float        # Vaporization temperature [K]
    
    # Enthalpy at phase boundaries (computed from T)
    hsmelt: float        # Enthalpy at solidus [J/kg]
    hlcal: float         # Enthalpy at liquidus [J/kg]
    hpreheat: float      # Enthalpy at preheat [J/kg]
    
    # Latent heat
    hlatent: float       # Latent heat of fusion [J/kg]
    
    # Density
    rho: float           # Reference density [kg/m³]
    rholiq: float        # Liquid density [kg/m³]
    
    # Thermal properties
    tcond: float         # Thermal conductivity [W/(m·K)]
    tcondl: float        # Liquid thermal conductivity [W/(m·K)]
    
    # Surface tension
    dgdt: float          # Surface tension gradient dγ/dT [N/(m·K)] (Marangoni)
    sigma_surf: float    # Surface tension [N/m]
    
    # Radiation/convection
    emiss: float         # Emissivity [-]
    sigma: float         # Stefan-Boltzmann constant [W/(m²·K⁴)]
    hconv: float         # Convection coefficient [W/(m²·K)]
    tenv: float          # Environment temperature [K]
    
    # Viscosity
    vis0: float          # Reference viscosity [Pa·s]
    
    # Gravity
    grav: float          # Gravitational acceleration [m/s²]
    beta: float          # Thermal expansion coefficient [1/K]
    
    # Mushy zone
    darcy_coeff: float   # Darcy coefficient for mushy zone damping


class SimulationParams(NamedTuple):
    """Numerical/simulation parameters (from mod_param.f90)"""
    # Time stepping
    delt: float          # Timestep [s]
    timax: float         # Maximum simulation time [s]
    
    # Under-relaxation factors
    urf_vel: float       # Velocity under-relaxation
    urf_p: float         # Pressure under-relaxation
    urf_h: float         # Enthalpy under-relaxation
    
    # Convergence criteria
    max_iter: int        # Maximum iterations per timestep
    conv_tol: float      # Convergence tolerance
    
    # Grid dimensions
    ni: int              # Number of cells in x
    nj: int              # Number of cells in y
    nk: int              # Number of cells in z
    
    # Domain size
    xlen: float          # Domain length in x [m]
    ylen: float          # Domain length in y [m]
    zlen: float          # Domain length in z [m]
    
    # Grid stretching parameters
    stretch_x: float     # Power-law exponent for x
    stretch_y: float     # Power-law exponent for y
    stretch_z: float     # Power-law exponent for z


class LaserParams(NamedTuple):
    """Laser parameters (from mod_param.f90, mod_laser.f90)"""
    power: float         # Laser power [W]
    radius: float        # Beam radius (1/e²) [m]
    absorptivity: float  # Absorptivity [-]
    efficiency: float    # Laser efficiency [-]
    
    # Derived: peak heat flux
    peak_flux: float     # Peak heat flux [W/m²]


class GridParams(NamedTuple):
    """Computational grid (immutable after initialization)"""
    # Cell centers
    x: Array             # x-coordinates of cell centers [ni]
    y: Array             # y-coordinates of cell centers [nj]
    z: Array             # z-coordinates of cell centers [nk]
    
    # Velocity face locations (staggered grid)
    xu: Array            # x-locations of u-velocity faces [ni]
    yv: Array            # y-locations of v-velocity faces [nj]
    zw: Array            # z-locations of w-velocity faces [nk]
    
    # Cell dimensions
    dx: Array            # Cell widths in x [ni]
    dy: Array            # Cell widths in y [nj]
    dz: Array            # Cell widths in z [nk]
    
    # Inverse distances for discretization (efficiency)
    dxpwinv: Array       # 1/(x[i] - x[i-1]) [ni]
    dxpeinv: Array       # 1/(x[i+1] - x[i]) [ni]
    dypsinv: Array       # 1/(y[j] - y[j-1]) [nj]
    dypninv: Array       # 1/(y[j+1] - y[j]) [nj]
    dzpbinv: Array       # 1/(z[k] - z[k-1]) [nk]
    dzptinv: Array       # 1/(z[k+1] - z[k]) [nk]
    
    # Cell volumes and face areas
    vol: Array           # Cell volumes [ni, nj, nk]
    areaij: Array        # xy-face areas [ni, nj]
    areaik: Array        # xz-face areas [ni, nk]
    areajk: Array        # yz-face areas [nj, nk]
    
    # Grid dimensions (convenience)
    ni: int
    nj: int
    nk: int


class State(NamedTuple):
    """Primary flow field variables (updated each timestep)"""
    uVel: Array          # x-velocity [ni, nj, nk]
    vVel: Array          # y-velocity [ni, nj, nk]
    wVel: Array          # z-velocity [ni, nj, nk]
    pressure: Array      # Pressure [ni, nj, nk]
    pp: Array            # Pressure correction [ni, nj, nk]
    enthalpy: Array      # Enthalpy [ni, nj, nk]
    temp: Array          # Temperature [ni, nj, nk]
    fracl: Array         # Liquid fraction [ni, nj, nk]


class StatePrev(NamedTuple):
    """Previous timestep values for transient terms"""
    unot: Array          # Previous u-velocity [ni, nj, nk]
    vnot: Array          # Previous v-velocity [ni, nj, nk]
    wnot: Array          # Previous w-velocity [ni, nj, nk]
    hnot: Array          # Previous enthalpy [ni, nj, nk]


class MaterialProps(NamedTuple):
    """Spatially-varying material properties"""
    vis: Array           # Viscosity [ni, nj, nk]
    diff: Array          # Thermal diffusivity [ni, nj, nk]
    den: Array           # Density [ni, nj, nk]
    tcond: Array         # Thermal conductivity [ni, nj, nk]


class DiscretCoeffs(NamedTuple):
    """FVM discretization coefficients (transient, rebuilt each solve)"""
    ap: Array            # Center coefficient [ni, nj, nk]
    ae: Array            # East neighbor [ni, nj, nk]
    aw: Array            # West neighbor [ni, nj, nk]
    an: Array            # North neighbor [ni, nj, nk]
    as_: Array           # South neighbor [ni, nj, nk] (as_ to avoid keyword)
    at: Array            # Top neighbor [ni, nj, nk]
    ab: Array            # Bottom neighbor [ni, nj, nk]
    su: Array            # Source term [ni, nj, nk]
    sp: Array            # Linearized source coefficient [ni, nj, nk]


class LaserState(NamedTuple):
    """Laser and toolpath state"""
    beam_x: float        # Current beam x-position [m]
    beam_y: float        # Current beam y-position [m]
    beam_z: float        # Current beam z-position [m]
    heatin: Array        # Heat flux at surface [ni, nj]
    laser_on: bool       # Laser on/off flag
    scanvel_x: float     # Scan velocity in x [m/s]
    scanvel_y: float     # Scan velocity in y [m/s]
    current_segment: int # Current toolpath segment index
    heat_total: float    # Total heat input this timestep [W]


class ToolPath(NamedTuple):
    """Toolpath data loaded from .crs file"""
    time: Array          # Time points [n_segments]
    x: Array             # x-positions [n_segments]
    y: Array             # y-positions [n_segments]
    z: Array             # z-positions [n_segments]
    laser_on: Array      # Laser on/off [n_segments] (0 or 1)
    n_segments: int      # Number of segments


class TimeState(NamedTuple):
    """Time stepping information"""
    timet: float         # Current simulation time [s]
    iteration: int       # Current iteration within timestep
    timestep: int        # Current timestep number
    converged: bool      # Convergence flag


class PoolDimensions(NamedTuple):
    """Melt pool dimensions"""
    length: float        # Pool length in x [m]
    width: float         # Pool width in y [m]
    depth: float         # Pool depth in z [m]
    volume: float        # Pool volume [m³]
    max_temp: float      # Maximum temperature [K]


class ConvergenceState(NamedTuple):
    """Convergence monitoring"""
    residual_h: float    # Enthalpy residual
    residual_u: float    # u-velocity residual
    residual_v: float    # v-velocity residual
    residual_w: float    # w-velocity residual
    residual_p: float    # Pressure residual
    heat_ratio: float    # Heat balance ratio (in/out)
    max_residual: float  # Maximum of all residuals


class OutputConfig(NamedTuple):
    """Output configuration"""
    output_interval: int     # Timesteps between outputs
    tecplot_interval: int    # Timesteps between Tecplot files
    output_dir: str          # Output directory path
    case_name: str           # Case name for file naming
