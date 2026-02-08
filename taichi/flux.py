"""
AM-CFD Taichi Implementation - Heat Flux Balance

Converted from Fortran module: mod_flux.f90

Computes heat fluxes across all domain boundaries and calculates
the energy conservation ratio for convergence monitoring.
"""

import taichi as ti
from data_structures import (
    State, StatePrev, GridParams, MaterialProps, LaserState,
    PhysicsParams, SimulationParams, ConvergenceState
)


@ti.data_oriented
class HeatFluxCalculator:
    """
    Computes boundary heat fluxes and energy balance.
    
    Corresponds to subroutine heat_fluxes in mod_flux.f90.
    """
    
    def __init__(self, ni: int, nj: int, nk: int):
        """Initialize flux accumulator fields."""
        self.ni = ni
        self.nj = nj
        self.nk = nk
        
        # Scalar accumulators for boundary fluxes
        self.flux_west = ti.field(dtype=ti.f64, shape=())
        self.flux_east = ti.field(dtype=ti.f64, shape=())
        self.flux_south = ti.field(dtype=ti.f64, shape=())
        self.flux_north = ti.field(dtype=ti.f64, shape=())
        self.flux_bottom = ti.field(dtype=ti.f64, shape=())
        self.flux_top = ti.field(dtype=ti.f64, shape=())
        
        # Heat accumulation and volumetric source
        self.accul = ti.field(dtype=ti.f64, shape=())
        self.heatvol = ti.field(dtype=ti.f64, shape=())
        
        # Total heat out and ratio
        self.heatout = ti.field(dtype=ti.f64, shape=())
        self.ratio = ti.field(dtype=ti.f64, shape=())
    
    @ti.kernel
    def compute_west_east_flux(self,
                               enthalpy: ti.template(),
                               diff: ti.template(),
                               dxpwinv: ti.template(),
                               areajk: ti.template()):
        """
        Compute heat fluxes through west (i=1) and east (i=ni) boundaries.
        
        Fortran equivalent:
            flux_west = sum over j,k of: diff(2,j,k) * (h(1,j,k) - h(2,j,k)) * dxpwinv(2) * areajk(j,k)
            flux_east = sum over j,k of: diff(nim1,j,k) * (h(ni,j,k) - h(nim1,j,k)) * dxpwinv(ni) * areajk(j,k)
        """
        self.flux_west[None] = 0.0
        self.flux_east[None] = 0.0
        
        ni = self.ni
        nj = self.nj
        nk = self.nk
        
        for j, k in ti.ndrange((1, nj - 1), (1, nk - 1)):
            # West boundary flux (i=0 is boundary, i=1 is first interior)
            fluxi1 = diff[1, j, k] * (enthalpy[0, j, k] - enthalpy[1, j, k]) * dxpwinv[1]
            self.flux_west[None] += areajk[j, k] * fluxi1
            
            # East boundary flux (i=ni-1 is last interior, i=ni-2 for nim1 in 0-based)
            fluxl1 = diff[ni - 2, j, k] * (enthalpy[ni - 1, j, k] - enthalpy[ni - 2, j, k]) * dxpwinv[ni - 1]
            self.flux_east[None] += areajk[j, k] * fluxl1
    
    @ti.kernel
    def compute_bottom_top_flux(self,
                                enthalpy: ti.template(),
                                diff: ti.template(),
                                dzpbinv: ti.template(),
                                areaij: ti.template()):
        """
        Compute heat fluxes through bottom (k=1) and top (k=nk) boundaries.
        """
        self.flux_bottom[None] = 0.0
        self.flux_top[None] = 0.0
        
        ni = self.ni
        nj = self.nj
        nk = self.nk
        
        for i, j in ti.ndrange((1, ni - 1), (1, nj - 1)):
            # Bottom boundary flux
            fluxk1 = diff[i, j, 1] * (enthalpy[i, j, 0] - enthalpy[i, j, 1]) * dzpbinv[1]
            self.flux_bottom[None] += areaij[i, j] * fluxk1
            
            # Top boundary flux
            fluxn1 = diff[i, j, nk - 2] * (enthalpy[i, j, nk - 1] - enthalpy[i, j, nk - 2]) * dzpbinv[nk - 1]
            self.flux_top[None] += areaij[i, j] * fluxn1
    
    @ti.kernel
    def compute_south_north_flux(self,
                                 enthalpy: ti.template(),
                                 diff: ti.template(),
                                 dypsinv: ti.template(),
                                 areaik: ti.template()):
        """
        Compute heat fluxes through south (j=1) and north (j=nj) boundaries.
        """
        self.flux_south[None] = 0.0
        self.flux_north[None] = 0.0
        
        ni = self.ni
        nj = self.nj
        nk = self.nk
        
        for i, k in ti.ndrange((1, ni - 1), (1, nk - 1)):
            # South boundary flux
            fluxj1 = diff[i, 1, k] * (enthalpy[i, 0, k] - enthalpy[i, 1, k]) * dypsinv[1]
            self.flux_south[None] += areaik[i, k] * fluxj1
            
            # North boundary flux
            fluxm1 = diff[i, nj - 2, k] * (enthalpy[i, nj - 1, k] - enthalpy[i, nj - 2, k]) * dypsinv[nj - 1]
            self.flux_north[None] += areaik[i, k] * fluxm1
    
    @ti.kernel
    def compute_accumulation(self,
                             enthalpy: ti.template(),
                             hnot: ti.template(),
                             fracl: ti.template(),
                             fraclnot: ti.template(),
                             den: ti.template(),
                             volume: ti.template(),
                             sourceinput: ti.template(),
                             hlatnt: ti.f64,
                             delt: ti.f64):
        """
        Compute heat accumulation and volumetric heat source.
        
        accul = sum of: volume * density * (dh + dfracl * hlatent) / delt
        heatvol = sum of: sourceinput * volume
        """
        self.accul[None] = 0.0
        self.heatvol[None] = 0.0
        
        ni = self.ni
        nj = self.nj
        nk = self.nk
        
        for i, j, k in ti.ndrange((1, ni - 1), (1, nj - 1), (1, nk - 1)):
            # Enthalpy change including latent heat contribution
            dh1 = (enthalpy[i, j, k] - hnot[i, j, k] + 
                   (fracl[i, j, k] - fraclnot[i, j, k]) * hlatnt)
            
            self.accul[None] += volume[i, j, k] * den[i, j, k] * dh1 / delt
            self.heatvol[None] += sourceinput[i, j, k] * volume[i, j, k]
    
    def compute_heat_balance(self, ahtoploss: float, heatinLaser: float) -> float:
        """
        Compute total heat balance and energy conservation ratio.
        
        Args:
            ahtoploss: Heat loss from top surface (convection + radiation)
            heatinLaser: Total laser heat input
            
        Returns:
            ratio: Energy conservation ratio (should be ~1.0 for good convergence)
        """
        # Total heat out through all boundaries
        heatout = (self.flux_north[None] + self.flux_bottom[None] + 
                   self.flux_west[None] + self.flux_east[None] + 
                   self.flux_south[None] - ahtoploss)
        
        self.heatout[None] = heatout
        
        # Energy balance ratio
        accul = self.accul[None]
        heatvol = self.heatvol[None]
        
        denominator = accul - heatout
        if abs(denominator) > 1.0e-20:
            ratio = (heatvol + heatinLaser) / denominator
        else:
            ratio = 1.0
        
        self.ratio[None] = ratio
        return ratio


def heat_fluxes(state: State,
                state_prev: StatePrev,
                mat_props: MaterialProps,
                grid: GridParams,
                physics: PhysicsParams,
                sim: SimulationParams,
                laser_state: LaserState,
                sourceinput: ti.template(),
                ahtoploss: float,
                conv: ConvergenceState,
                flux_calc: HeatFluxCalculator) -> ConvergenceState:
    """
    Compute heat balance ratio for convergence monitoring.
    
    This is the main interface function corresponding to subroutine heat_fluxes
    in mod_flux.f90.
    
    Args:
        state: Current flow field state
        state_prev: Previous timestep state (hnot, fraclnot)
        mat_props: Material properties (diff, den)
        grid: Computational grid (volumes, areas, inverse distances)
        physics: Physical parameters (hlatent)
        sim: Simulation parameters (delt)
        laser_state: Laser state (heat_total = heatinLaser)
        sourceinput: Volumetric heat source field
        ahtoploss: Heat loss from top surface
        conv: Convergence state to update
        flux_calc: Pre-allocated HeatFluxCalculator instance
        
    Returns:
        conv: Updated convergence state with heat_ratio
    """
    # Compute boundary fluxes
    flux_calc.compute_west_east_flux(
        state.enthalpy, mat_props.diff, grid.dxpwinv, grid.areajk)
    
    flux_calc.compute_bottom_top_flux(
        state.enthalpy, mat_props.diff, grid.dzpbinv, grid.areaij)
    
    flux_calc.compute_south_north_flux(
        state.enthalpy, mat_props.diff, grid.dypsinv, grid.areaik)
    
    # Compute heat accumulation
    flux_calc.compute_accumulation(
        state.enthalpy, state_prev.hnot,
        state.fracl, state_prev.fraclnot,
        mat_props.den, grid.vol,
        sourceinput,
        physics.hlatent, sim.delt)
    
    # Compute energy balance ratio
    conv.heat_ratio = flux_calc.compute_heat_balance(
        ahtoploss, laser_state.heat_total)
    
    return conv
