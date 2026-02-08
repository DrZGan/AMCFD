!______________________________________________________________________________
!
module simulation_state
!______________________________________________________________________________
! Simulation state variables: equation index, beam position, toolpath data, etc.
! Split from the original 'initialization' God Module for better encapsulation.
!
	use constant
	implicit none

	! Current equation index (1=u, 2=v, 3=w, 4=pp, 5=enthalpy)
	integer ivar

	! Toolpath tracking
	integer TrackNum, PathNum
	real toolmatrix(TOOLLINES, 5)
	real coordhistory(COORDLINES, 8)

	! Beam position
	real beam_pos, beam_posy

	! Solidification tracking field
	real solidfield(nx,ny,nz)

	! Repeated heating factor
	real RHF

	! Residual and reference values
	real resorm, refmom, ahtoploss

	! Boundary enthalpies (derived from boundary temperatures)
	real enthalpyWest, enthalpyEast, enthalpyNorth, enthalpyBottom, enthalpyPreheat

end module simulation_state
