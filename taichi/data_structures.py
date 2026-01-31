"""
AM-CFD Taichi Implementation - Type Definitions

Converted from Fortran modules: mod_const.f90, mod_param.f90, and variable declarations
across all modules.

NOTE: Taichi uses dataclasses with ti.field() for GPU arrays instead of NamedTuples.
Parameter containers use dataclasses for Taichi compatibility.
"""

from dataclasses import dataclass
from typing import Optional
import taichi as ti
import numpy as np


@dataclass(unsafe_hash=True)
class PhysicsParams:
    """Physical constants and material properties (from mod_const.f90, mod_param.f90)"""
    # Specific heat coefficients (solid: Cp = acpa*T + acpb)
    acpa: float = 0.0           # Solid Cp coefficient a [J/(kg·K²)]
    acpb: float = 500.0         # Solid Cp coefficient b [J/(kg·K)]
    acpl: float = 800.0         # Liquid specific heat [J/(kg·K)]
    
    # Phase change temperatures
    tsolid: float = 1563.0      # Solidus temperature [K]
    tliquid: float = 1623.0     # Liquidus temperature [K]
    tpreheat: float = 300.0     # Preheat temperature [K]
    tvapor: float = 3000.0      # Vaporization temperature [K]
    
    # Enthalpy at phase boundaries (computed from T)
    hsmelt: float = 0.0         # Enthalpy at solidus [J/kg]
    hlcal: float = 0.0          # Enthalpy at liquidus [J/kg]
    hpreheat: float = 0.0       # Enthalpy at preheat [J/kg]
    
    # Latent heat
    hlatent: float = 2.7e5      # Latent heat of fusion [J/kg]
    
    # Density
    rho: float = 7800.0         # Reference density [kg/m³]
    rholiq: float = 7800.0      # Liquid density [kg/m³]
    
    # Thermal properties
    tcond: float = 25.0         # Thermal conductivity [W/(m·K)]
    tcondl: float = 25.0        # Liquid thermal conductivity [W/(m·K)]
    
    # Surface tension
    dgdt: float = -4.3e-4       # Surface tension gradient dγ/dT [N/(m·K)] (Marangoni)
    sigma_surf: float = 1.8     # Surface tension [N/m]
    
    # Radiation/convection
    emiss: float = 0.4          # Emissivity [-]
    sigma: float = 5.67e-8      # Stefan-Boltzmann constant [W/(m²·K⁴)]
    hconv: float = 10.0         # Convection coefficient [W/(m²·K)]
    tenv: float = 300.0         # Environment temperature [K]
    
    # Viscosity
    vis0: float = 6.0e-3        # Reference viscosity [Pa·s]
    
    # Gravity
    grav: float = 9.81          # Gravitational acceleration [m/s²]
    beta: float = 1.0e-4        # Thermal expansion coefficient [1/K]
    
    # Mushy zone
    darcy_coeff: float = 1.0e6  # Darcy coefficient for mushy zone damping


@dataclass(unsafe_hash=True)
class SimulationParams:
    """Numerical/simulation parameters (from mod_param.f90)"""
    # Time stepping
    delt: float = 1.0e-6        # Timestep [s]
    timax: float = 1.0e-3       # Maximum simulation time [s]
    
    # Under-relaxation factors
    urf_vel: float = 0.5        # Velocity under-relaxation
    urf_p: float = 0.3          # Pressure under-relaxation
    urf_h: float = 0.8          # Enthalpy under-relaxation
    
    # Convergence criteria
    max_iter: int = 1000        # Maximum iterations per timestep
    conv_tol: float = 5.0e-4    # Convergence tolerance
    
    # Grid dimensions
    ni: int = 100               # Number of cells in x
    nj: int = 100               # Number of cells in y
    nk: int = 50                # Number of cells in z
    
    # Domain size
    xlen: float = 1.0e-3        # Domain length in x [m]
    ylen: float = 1.0e-3        # Domain length in y [m]
    zlen: float = 0.5e-3        # Domain length in z [m]
    
    # Grid stretching parameters
    stretch_x: float = 1.0      # Power-law exponent for x
    stretch_y: float = 1.0      # Power-law exponent for y
    stretch_z: float = 1.0      # Power-law exponent for z
    
    # Derived index bounds (computed in __post_init__)
    istat: int = 0              # Start index in i (typically 0 in Python)
    jstat: int = 0              # Start index in j
    kstat: int = 0              # Start index in k
    iend: int = 0               # End index in i (ni)
    jend: int = 0               # End index in j (nj)
    kend: int = 0               # End index in k (nk)
    istatp1: int = 0            # istat + 1
    jstatp1: int = 0            # jstat + 1
    kstatp1: int = 0            # kstat + 1
    iendm1: int = 0             # iend - 1
    jendm1: int = 0             # jend - 1
    kendm1: int = 0             # kend - 1
    nim1: int = 0               # ni - 1
    njm1: int = 0               # nj - 1
    nkm1: int = 0               # nk - 1
    
    def __post_init__(self):
        """Compute derived index bounds."""
        # Start indices (0-based Python)
        self.istat = 0
        self.jstat = 0
        self.kstat = 0
        
        # End indices
        self.iend = self.ni
        self.jend = self.nj
        self.kend = self.nk
        
        # Start + 1
        self.istatp1 = self.istat + 1
        self.jstatp1 = self.jstat + 1
        self.kstatp1 = self.kstat + 1
        
        # End - 1
        self.iendm1 = self.iend - 1
        self.jendm1 = self.jend - 1
        self.kendm1 = self.kend - 1
        
        # Dimension - 1
        self.nim1 = self.ni - 1
        self.njm1 = self.nj - 1
        self.nkm1 = self.nk - 1


@dataclass
class LaserParams:
    """Laser parameters (from mod_param.f90, mod_laser.f90)"""
    power: float = 200.0        # Laser power [W]
    radius: float = 50.0e-6     # Beam radius (1/e²) [m]
    absorptivity: float = 0.35  # Absorptivity [-]
    efficiency: float = 1.0     # Laser efficiency [-]
    
    # Derived: peak heat flux
    peak_flux: float = 0.0      # Peak heat flux [W/m²]
    
    def __post_init__(self):
        """Compute peak flux after initialization."""
        if self.peak_flux == 0.0:
            import math
            self.peak_flux = 2.0 * self.power * self.absorptivity * self.efficiency / (
                math.pi * self.radius * self.radius)


@dataclass
class OutputConfig:
    """Output configuration"""
    output_interval: int = 100      # Timesteps between outputs
    tecplot_interval: int = 1000    # Timesteps between Tecplot files
    output_dir: str = './output'    # Output directory path
    case_name: str = 'amcfd'        # Case name for file naming


@dataclass
class ToolPath:
    """Toolpath data loaded from .crs file"""
    time: np.ndarray = None         # Time points [n_segments]
    x: np.ndarray = None            # x-positions [n_segments]
    y: np.ndarray = None            # y-positions [n_segments]
    z: np.ndarray = None            # z-positions [n_segments]
    laser_on: np.ndarray = None     # Laser on/off [n_segments] (0 or 1)
    n_segments: int = 0             # Number of segments
    
    def __post_init__(self):
        """Initialize empty arrays if not provided."""
        if self.time is None:
            self.time = np.array([0.0])
            self.x = np.array([0.0])
            self.y = np.array([0.0])
            self.z = np.array([0.0])
            self.laser_on = np.array([0], dtype=np.int32)
            self.n_segments = 1


@ti.data_oriented
class GridParams:
    """
    Computational grid using Taichi fields.
    
    This class holds immutable grid geometry after initialization.
    Uses ti.field() for GPU-compatible 1D and 3D arrays.
    """
    
    def __init__(self, ni: int, nj: int, nk: int):
        """Initialize grid with dimensions."""
        self.ni = ni
        self.nj = nj
        self.nk = nk
        
        # Cell centers (1D fields)
        self.x = ti.field(dtype=ti.f64, shape=(ni,))
        self.y = ti.field(dtype=ti.f64, shape=(nj,))
        self.z = ti.field(dtype=ti.f64, shape=(nk,))
        
        # Velocity face locations (staggered grid)
        self.xu = ti.field(dtype=ti.f64, shape=(ni,))
        self.yv = ti.field(dtype=ti.f64, shape=(nj,))
        self.zw = ti.field(dtype=ti.f64, shape=(nk,))
        
        # Cell dimensions
        self.dx = ti.field(dtype=ti.f64, shape=(ni,))
        self.dy = ti.field(dtype=ti.f64, shape=(nj,))
        self.dz = ti.field(dtype=ti.f64, shape=(nk,))
        
        # Inverse distances for discretization (efficiency)
        self.dxpwinv = ti.field(dtype=ti.f64, shape=(ni,))
        self.dxpeinv = ti.field(dtype=ti.f64, shape=(ni,))
        self.dypsinv = ti.field(dtype=ti.f64, shape=(nj,))
        self.dypninv = ti.field(dtype=ti.f64, shape=(nj,))
        self.dzpbinv = ti.field(dtype=ti.f64, shape=(nk,))
        self.dzptinv = ti.field(dtype=ti.f64, shape=(nk,))
        
        # Cell volumes and face areas
        self.vol = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.areaij = ti.field(dtype=ti.f64, shape=(ni, nj))
        self.areaik = ti.field(dtype=ti.f64, shape=(ni, nk))
        self.areajk = ti.field(dtype=ti.f64, shape=(nj, nk))


@ti.data_oriented
class State:
    """
    Primary flow field variables (updated each timestep).
    
    Uses Taichi fields for GPU computation.
    """
    
    def __init__(self, ni: int, nj: int, nk: int):
        """Initialize state arrays with dimensions."""
        self.ni = ni
        self.nj = nj
        self.nk = nk
        
        # Velocity components
        self.uVel = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.vVel = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.wVel = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        
        # Pressure fields
        self.pressure = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.pp = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        
        # Thermal fields
        self.enthalpy = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.temp = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.fracl = ti.field(dtype=ti.f64, shape=(ni, nj, nk))


@ti.data_oriented
class StatePrev:
    """Previous timestep values for transient terms."""
    
    def __init__(self, ni: int, nj: int, nk: int):
        """Initialize previous state arrays."""
        self.unot = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.vnot = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.wnot = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.tnot = ti.field(dtype=ti.f64, shape=(ni, nj, nk))   # Previous temperature
        self.hnot = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.fraclnot = ti.field(dtype=ti.f64, shape=(ni, nj, nk))  # Previous liquid fraction


@ti.data_oriented
class MaterialProps:
    """Spatially-varying material properties."""
    
    def __init__(self, ni: int, nj: int, nk: int):
        """Initialize material property arrays."""
        self.vis = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.diff = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.den = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.tcond = ti.field(dtype=ti.f64, shape=(ni, nj, nk))


@ti.data_oriented
class DiscretCoeffs:
    """FVM discretization coefficients (transient, rebuilt each solve)."""
    
    def __init__(self, ni: int, nj: int, nk: int):
        """Initialize discretization coefficient arrays."""
        self.ap = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.ae = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.aw = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.an = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.as_ = ti.field(dtype=ti.f64, shape=(ni, nj, nk))  # as_ to avoid keyword
        self.at = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.ab = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.su = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.sp = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        
        # Velocity correction coefficients (for SIMPLE algorithm)
        self.dux = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.dvy = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.dwz = ti.field(dtype=ti.f64, shape=(ni, nj, nk))


@ti.data_oriented
class LaserState:
    """Laser and toolpath state."""
    
    def __init__(self, ni: int, nj: int):
        """Initialize laser state."""
        self.beam_x = 0.0           # Current beam x-position [m]
        self.beam_y = 0.0           # Current beam y-position [m]
        self.beam_z = 0.0           # Current beam z-position [m]
        self.heatin = ti.field(dtype=ti.f64, shape=(ni, nj))  # Heat flux at surface
        self.laser_on = False       # Laser on/off flag
        self.scanvel_x = 0.0        # Scan velocity in x [m/s]
        self.scanvel_y = 0.0        # Scan velocity in y [m/s]
        self.current_segment = 0    # Current toolpath segment index
        self.heat_total = 0.0       # Total heat input this timestep [W]


@dataclass
class TimeState:
    """Time stepping information."""
    timet: float = 0.0          # Current simulation time [s]
    iteration: int = 0          # Current iteration within timestep
    timestep: int = 0           # Current timestep number
    converged: bool = False     # Convergence flag


@dataclass
class PoolDimensions:
    """Melt pool dimensions."""
    length: float = 0.0         # Pool length in x [m]
    width: float = 0.0          # Pool width in y [m]
    depth: float = 0.0          # Pool depth in z [m]
    volume: float = 0.0         # Pool volume [m³]
    max_temp: float = 0.0       # Maximum temperature [K]


@dataclass
class ConvergenceState:
    """Convergence monitoring."""
    residual_h: float = 0.0     # Enthalpy residual
    residual_u: float = 0.0     # u-velocity residual
    residual_v: float = 0.0     # v-velocity residual
    residual_w: float = 0.0     # w-velocity residual
    residual_p: float = 0.0     # Pressure residual
    heat_ratio: float = 0.0     # Heat balance ratio (in/out)
    max_residual: float = 0.0   # Maximum of all residuals
