!______________________________________________________________________________
!
module constant
!______________________________________________________________________________
!
	implicit none

	! Physical constants
	real, parameter :: g = 9.8
	real, parameter :: pi = 3.1415926
	real, parameter :: sigm = 5.67e-8
	real, parameter :: great = 1.0e20
	real, parameter :: small = 1.0e-06

	! Maximum array dimensions
	integer, parameter :: nx = 1200, ny = 1200, nz = 180, nvar = 4
	integer, parameter :: nx1 = 7, ny1 = 7, nz1 = 7, ng = 5
	integer, parameter :: TOOLLINES = 1000, COORDLINES = 5000

	! Numerical parameters (previously magic numbers)
	real, parameter :: KOZENY_CARMAN_CONST = 180.0        ! Kozeny-Carman constant for mushy zone
	real, parameter :: DENDRITE_ARM_SQ = (1.0e-5)**2      ! Square of dendrite arm spacing
	real, parameter :: DARCY_SMALL = 1.0e-3               ! Small constant to avoid division by zero in Darcy term
	real, parameter :: SOLID_VIS = 1.0e10                  ! Artificial large viscosity for solid phase
	real, parameter :: TDMA_ZERO_GUARD = 1.0e-12          ! TDMA denominator zero guard threshold
	real, parameter :: TDMA_ZERO_OFFSET = 1.0e-13         ! TDMA denominator zero guard offset
	real, parameter :: CONV_HEAT_THRESHOLD = 5.0e-4       ! Convergence criterion for heating stage
	real, parameter :: CONV_COOL_THRESHOLD = 5.0e-7       ! Convergence criterion for cooling stage
	real, parameter :: ENERGY_RATIO_LOW = 0.99             ! Lower bound of energy conservation ratio
	real, parameter :: ENERGY_RATIO_HIGH = 1.01            ! Upper bound of energy conservation ratio

contains

	!--------------------------------------------------------------------------
	! Power-law convection-diffusion scheme coefficient
	! Computes: d * max(0, (1 - 0.1*|f|/d)^5) + max(0, -f)
	!--------------------------------------------------------------------------
	pure real function power_law_coeff(d, f) result(c)
		real, intent(in) :: d, f
		c = d * max(0.0, (1.0 - 0.1 * abs(f) / d)**5) + max(0.0, -f)
	end function power_law_coeff

	!--------------------------------------------------------------------------
	! Power-law coefficient for the "positive flow" direction
	! Computes: d * max(0, (1 - 0.1*|f|/d)^5) + max(0, f)
	!--------------------------------------------------------------------------
	pure real function power_law_coeff_pos(d, f) result(c)
		real, intent(in) :: d, f
		c = d * max(0.0, (1.0 - 0.1 * abs(f) / d)**5) + max(0.0, f)
	end function power_law_coeff_pos

	!--------------------------------------------------------------------------
	! Darcy resistance term for mushy zone momentum damping
	!--------------------------------------------------------------------------
	pure real function darcy_resistance(viscos, fl) result(term)
		real, intent(in) :: viscos, fl
		term = KOZENY_CARMAN_CONST * viscos / DENDRITE_ARM_SQ * (1.0 - fl)**2 / (fl + DARCY_SMALL)
	end function darcy_resistance

	!--------------------------------------------------------------------------
	! Zero out all coefficients for solid region (enforce zero velocity)
	!--------------------------------------------------------------------------
	subroutine zero_velocity_coeff(su_ijk, sp_ijk, ap_ijk, an_ijk, as_ijk, ae_ijk, aw_ijk, at_ijk, ab_ijk)
		real, intent(out) :: su_ijk, ap_ijk, an_ijk, as_ijk, ae_ijk, aw_ijk, at_ijk, ab_ijk
		real, intent(inout) :: sp_ijk
		su_ijk = 0.0
		an_ijk = 0.0
		as_ijk = 0.0
		ae_ijk = 0.0
		aw_ijk = 0.0
		at_ijk = 0.0
		ab_ijk = 0.0
		ap_ijk = great
	end subroutine zero_velocity_coeff

end module constant
