!______________________________________________________________________________
!
module sim_state
!______________________________________________________________________________
! Simulation state: derived constants, residuals, boundary values, laser state
!
	use precision
	use constant
	implicit none

	! Derived constants (set in initialize)
	real(wp) dgdt, deltemp, cpavg, hlcal, hlatnt, boufac

	! Residual scalars
	real(wp) resorm, refmom, ahtoploss

	! Boundary enthalpies
	real(wp) enthalpyWest, enthalpyEast, enthalpyNorth, enthalpyBottom, enthalpyPreheat

	! Laser/toolpath state
	integer TrackNum, PathNum
	real(wp) beam_pos, beam_posy
	real(wp) toolmatrix(TOOLLINES,5)
	real(wp) coordhistory(COORDLINES,8)
	real(wp) RHF

end module sim_state
