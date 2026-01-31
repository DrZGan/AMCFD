"""
AM-CFD Taichi Implementation - Toolpath Module

Translated from Fortran mod_toolpath.f90
Handles toolpath reading, coordinate tracking, and heat history.
"""

import numpy as np
from pathlib import Path
from typing import Optional

from data_structures import ToolPath, LaserState, TimeState, GridParams


def load_toolpath(filepath: str) -> ToolPath:
    """
    Load toolpath from .crs file. (Corresponds to read_toolpath in Fortran)
    
    Format: 5 columns - time, x, y, z, laser_on
    The toolmatrix in Fortran stores: time, x, y, z, power_flag
    
    Args:
        filepath: Path to .crs toolpath file
        
    Returns:
        ToolPath dataclass with arrays for each column
    """
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


def get_current_segment(timet: float, toolpath: ToolPath) -> int:
    """
    Find the current toolpath segment index for given time.
    
    Args:
        timet: Current simulation time [s]
        toolpath: ToolPath data
        
    Returns:
        Index of current segment (PathNum in Fortran)
    """
    # Find the segment where time[i] <= timet < time[i+1]
    for i in range(toolpath.n_segments - 1):
        if toolpath.time[i] <= timet < toolpath.time[i + 1]:
            return i
    
    # If past all segments, return last one
    return toolpath.n_segments - 1


def interpolate_position(timet: float, toolpath: ToolPath) -> tuple[float, float, float, float, float]:
    """
    Interpolate beam position and velocity from toolpath.
    
    Args:
        timet: Current simulation time [s]
        toolpath: ToolPath data
        
    Returns:
        Tuple of (beam_x, beam_y, beam_z, scanvel_x, scanvel_y)
    """
    seg = get_current_segment(timet, toolpath)
    
    if seg >= toolpath.n_segments - 1:
        # At or past end of toolpath
        return (
            float(toolpath.x[-1]),
            float(toolpath.y[-1]),
            float(toolpath.z[-1]),
            0.0,
            0.0
        )
    
    # Get segment endpoints
    t0 = toolpath.time[seg]
    t1 = toolpath.time[seg + 1]
    
    x0, x1 = toolpath.x[seg], toolpath.x[seg + 1]
    y0, y1 = toolpath.y[seg], toolpath.y[seg + 1]
    z0, z1 = toolpath.z[seg], toolpath.z[seg + 1]
    
    # Time fraction within segment
    dt = t1 - t0
    if dt > 0:
        frac = (timet - t0) / dt
    else:
        frac = 0.0
    
    # Interpolate position
    beam_x = x0 + frac * (x1 - x0)
    beam_y = y0 + frac * (y1 - y0)
    beam_z = z0 + frac * (z1 - z0)
    
    # Compute velocity (distance/time)
    if dt > 0:
        scanvel_x = (x1 - x0) / dt
        scanvel_y = (y1 - y0) / dt
    else:
        scanvel_x = 0.0
        scanvel_y = 0.0
    
    return beam_x, beam_y, beam_z, scanvel_x, scanvel_y


def read_coordinates(time_state: TimeState, toolpath: ToolPath,
                     laser_state: LaserState) -> LaserState:
    """
    Read/interpolate coordinates from toolpath. (From mod_toolpath.f90)
    
    Updates laser_state with current beam position, velocity, and laser on/off status.
    Also maintains coordinate history for RHF (Reheating Factor) calculation.
    
    Args:
        time_state: Current time information
        toolpath: Loaded toolpath data
        laser_state: Laser state to update
        
    Returns:
        Updated LaserState
    """
    timet = time_state.timet
    
    # Get current segment
    seg = get_current_segment(timet, toolpath)
    laser_state.current_segment = seg
    
    # Interpolate position and velocity
    beam_x, beam_y, beam_z, scanvel_x, scanvel_y = interpolate_position(timet, toolpath)
    
    laser_state.beam_x = beam_x
    laser_state.beam_y = beam_y
    laser_state.beam_z = beam_z
    laser_state.scanvel_x = scanvel_x
    laser_state.scanvel_y = scanvel_y
    
    # Determine if laser is on (toolmatrix(PathNum,5) >= 0.5 in Fortran)
    laser_state.laser_on = bool(toolpath.laser_on[seg] >= 1)
    
    return laser_state


def calc_rhf(coord_history: np.ndarray, current_pos: tuple[float, float],
             current_time: float, power: float,
             R: float = 0.2e-3, T: float = 2e-3, P0: float = 300.0) -> float:
    """
    Calculate Reheating Factor (RHF). (From mod_toolpath.f90 calcRHF subroutine)
    
    RHF measures the cumulative thermal effect from previous beam positions,
    accounting for spatial and temporal decay.
    
    Args:
        coord_history: Array of [time, x, y, z, power, ...] history
        current_pos: Current beam (x, y) position
        current_time: Current simulation time
        power: Current laser power
        R: Spatial influence radius [m]
        T: Temporal influence window [s]
        P0: Reference power [W]
        
    Returns:
        RHF value (normalized reheating factor)
    """
    RHFc = 3.176  # Normalization constant from Fortran
    RHF_sum = 0.0
    
    x_curr, y_curr = current_pos
    
    for i in range(len(coord_history)):
        if coord_history[i, 0] < -0.5:
            # Invalid/empty entry
            break
        
        t_hist = coord_history[i, 0]
        x_hist = coord_history[i, 1]
        y_hist = coord_history[i, 2]
        P_hist = coord_history[i, 4] if coord_history.shape[1] > 4 else power
        
        # Distance from historical position to current position
        disk = np.sqrt((x_hist - x_curr)**2 + (y_hist - y_curr)**2)
        
        # Time since historical point
        tk = current_time - t_hist
        
        # Apply spatial and temporal weighting
        if disk <= R and tk <= T and tk >= 0:
            RHFk = ((R - disk)**2 / R**2) * ((T - tk) / T) * (P_hist / P0)
            RHF_sum += RHFk
    
    return RHF_sum / RHFc


class CoordinateHistory:
    """
    Maintains coordinate history for RHF calculation.
    
    Corresponds to coordhistory array in Fortran.
    Columns: time, x, y, z, power, scan_speed, scanvel_x, scanvel_y
    """
    
    def __init__(self, max_entries: int = 1000):
        """Initialize coordinate history buffer."""
        self.max_entries = max_entries
        self.n_cols = 8
        # Initialize with -1 to indicate empty entries
        self.history = np.full((max_entries, self.n_cols), -1.0)
        self.count = 0
    
    def add_entry(self, timet: float, beam_x: float, beam_y: float, beam_z: float,
                  power: float, scan_speed: float, scanvel_x: float, scanvel_y: float) -> None:
        """
        Add a new entry to the history, shifting older entries down.
        
        This matches the Fortran logic in read_coordinates that shifts
        coordhistory entries to make room for the new one.
        """
        # Shift existing entries down
        if self.count > 0:
            if self.count < self.max_entries:
                # Shift from count down to 1
                self.history[1:self.count+1, :] = self.history[0:self.count, :]
            else:
                # Buffer full, shift everything
                self.history[1:, :] = self.history[:-1, :]
        
        # Insert new entry at position 0
        self.history[0, :] = [timet, beam_x, beam_y, beam_z,
                              power, scan_speed, scanvel_x, scanvel_y]
        
        if self.count < self.max_entries:
            self.count += 1
    
    def get_rhf(self, current_x: float, current_y: float, current_time: float,
                R: float = 0.2e-3, T: float = 2e-3, P0: float = 300.0) -> float:
        """Calculate RHF using stored history."""
        return calc_rhf(self.history, (current_x, current_y), current_time,
                        self.history[0, 4] if self.count > 0 else 0.0, R, T, P0)
