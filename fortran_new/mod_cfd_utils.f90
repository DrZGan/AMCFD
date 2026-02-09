!______________________________________________________________________________
!
module cfd_utils
!______________________________________________________________________________
!
! Utility functions for CFD calculations
! - Harmonic mean for interface properties
! - Power-law coefficient for convection-diffusion
! - Darcy resistance for mushy zone penalty
!
    use constant
    implicit none
    private

    public :: harmonic_mean, power_law_coeff, darcy_resistance, temp_to_enthalpy

contains

!------------------------------------------------------------------------------
! Harmonic mean of two values with weighting fraction
! Used for interface properties (viscosity, diffusivity) between cells
!------------------------------------------------------------------------------
pure elemental real(wp) function harmonic_mean(val1, val2, frac) result(res)
    real(wp), intent(in) :: val1, val2, frac
    res = val1 * val2 / (frac * val2 + (1.0_wp - frac) * val1)
end function harmonic_mean

!------------------------------------------------------------------------------
! Power-law coefficient for convection-diffusion scheme
! Returns: diff * max(0, (1 - 0.1*|Pe|)^5) + max(0, -flux)
! where Pe = flux/diff (Peclet number)
!------------------------------------------------------------------------------
pure elemental real(wp) function power_law_coeff(diff, flux) result(res)
    real(wp), intent(in) :: diff, flux
    res = diff * max(0.0_wp, (1.0_wp - 0.1_wp*abs(flux/diff))**5) + max(0.0_wp, -flux)
end function power_law_coeff

!------------------------------------------------------------------------------
! Darcy resistance coefficient for mushy zone momentum sink
! Implements Carman-Kozeny equation for flow through porous media
! term = 180 * viscosity / permeability * (1-fl)^2 / (fl + epsilon)
!------------------------------------------------------------------------------
pure elemental real(wp) function darcy_resistance(viscos, fracl) result(res)
    real(wp), intent(in) :: viscos, fracl
    real(wp), parameter :: perm_const = 1.0e-10_wp  ! permeability constant (1e-5)^2
    real(wp), parameter :: eps = 1.0e-3_wp           ! small constant to avoid division by zero
    res = 180.0_wp * viscos / perm_const * (1.0_wp - fracl)**2 / (fracl + eps)
end function darcy_resistance

!------------------------------------------------------------------------------
! Convert temperature to enthalpy using quadratic heat capacity fit
!------------------------------------------------------------------------------
pure elemental real(wp) function temp_to_enthalpy(temp, acpa, acpb) result(res)
    real(wp), intent(in) :: temp, acpa, acpb
    res = 0.5_wp * acpa * temp**2 + acpb * temp
end function temp_to_enthalpy

end module cfd_utils
