!______________________________________________________________________________
!
module cfd_utils
!______________________________________________________________________________
!
! Utility functions for CFD calculations
! - Harmonic mean for interface properties
! - Power-law coefficient for convection-diffusion
! - Darcy resistance for mushy zone penalty
! - temp_to_enthalpy for boundary initial conditions
!
    use constant
    use precision
    implicit none
    private

    public :: harmonic_mean, power_law_coeff, darcy_resistance, temp_to_enthalpy

contains

!------------------------------------------------------------------------------
! Temperature to enthalpy (solid: acpa*T^2/2+acpb*T; liquid: hlcal+acpl*(T-tliquid); mushy: linear)
!------------------------------------------------------------------------------
pure real(wp) function temp_to_enthalpy(T, acpa, acpb, acpl, tsolid, tliquid, hsmelt, hlcal, deltemp) result(h)
    real(wp), intent(in) :: T, acpa, acpb, acpl, tsolid, tliquid, hsmelt, hlcal, deltemp
    if (T >= tliquid) then
        h = hlcal + acpl * (T - tliquid)
    else if (T <= tsolid) then
        h = acpa * T**2 / 2.0_wp + acpb * T
    else
        h = hsmelt + (hlcal - hsmelt) * (T - tsolid) / deltemp
    endif
end function temp_to_enthalpy

!------------------------------------------------------------------------------
! Harmonic mean of two values with weighting fraction
! Used for interface properties (viscosity, diffusivity) between cells
!------------------------------------------------------------------------------
pure elemental real function harmonic_mean(val1, val2, frac) result(res)
    real, intent(in) :: val1, val2, frac
    res = val1 * val2 / (frac * val2 + (1.0 - frac) * val1)
end function harmonic_mean

!------------------------------------------------------------------------------
! Power-law coefficient for convection-diffusion scheme
! Returns: diff * max(0, (1 - 0.1*|Pe|)^5) + max(0, -flux)
! where Pe = flux/diff (Peclet number)
!------------------------------------------------------------------------------
pure elemental real function power_law_coeff(diff, flux) result(res)
    real, intent(in) :: diff, flux
    res = diff * max(0.0, (1.0 - 0.1*abs(flux/diff))**5) + max(0.0, -flux)
end function power_law_coeff

!------------------------------------------------------------------------------
! Darcy resistance coefficient for mushy zone momentum sink
! Implements Carman-Kozeny equation for flow through porous media
! term = 180 * viscosity / permeability * (1-fl)^2 / (fl + epsilon)
!------------------------------------------------------------------------------
pure elemental real function darcy_resistance(viscos, fracl) result(res)
    real, intent(in) :: viscos, fracl
    real, parameter :: perm_const = 1.0e-10  ! permeability constant (1e-5)^2
    real, parameter :: eps = 1.0e-3          ! small constant to avoid division by zero
    res = 180.0 * viscos / perm_const * (1.0 - fracl)**2 / (fracl + eps)
end function darcy_resistance

end module cfd_utils
