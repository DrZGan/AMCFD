"""
AM-CFD JAX Implementation

Powder Bed Fusion CFD simulation for Additive Manufacturing.
Converted from Fortran to JAX for GPU acceleration.
"""

from .data_structures import (
    PhysicsParams,
    SimulationParams,
    LaserParams,
    GridParams,
    State,
    StatePrev,
    MaterialProps,
    DiscretCoeffs,
    LaserState,
    ToolPath,
    TimeState,
    PoolDimensions,
    ConvergenceState,
    OutputConfig,
)

from .param import (
    parse_input,
    load_toolpath,
    write_output_header,
    write_output_line,
    write_tecplot,
)

__version__ = "0.1.0"
__all__ = [
    # Types
    "PhysicsParams",
    "SimulationParams", 
    "LaserParams",
    "GridParams",
    "State",
    "StatePrev",
    "MaterialProps",
    "DiscretCoeffs",
    "LaserState",
    "ToolPath",
    "TimeState",
    "PoolDimensions",
    "ConvergenceState",
    "OutputConfig",
    # IO functions
    "parse_input",
    "load_toolpath",
    "write_output_header",
    "write_output_line",
    "write_tecplot",
]
