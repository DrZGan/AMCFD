!______________________________________________________________________________
!
module initialization
!______________________________________________________________________________
! Re-exports all sub-modules so downstream code using 'use initialization'
! continues to work without changes. Contains derived parameters and the
! initialize subroutine.
!
	use geometry
	use constant
	use parameters
	use field_vars
	use coefficients
	use material_fields
	use simulation_state

	implicit none

	! Derived material parameters (computed in initialize)
	real dgdt, deltemp, cpavg, hlcal, hlatnt, boufac

	contains

subroutine initialize

	integer i, j, k

	! Derived material parameters
	dgdt    = dgdtp
	deltemp = tliquid - tsolid
	cpavg   = (acpa * tsolid + acpb + acpl) * 0.5
	hlcal   = hsmelt + cpavg * deltemp
	hlatnt  = hlfriz - hlcal
	boufac  = denl * g * beta

	! Convert boundary temperatures to enthalpies using common function
	enthalpyPreheat = temp_to_enthalpy(tempPreheat)
	enthalpyWest    = temp_to_enthalpy(tempWest)
	enthalpyEast    = temp_to_enthalpy(tempEast)
	enthalpyNorth   = temp_to_enthalpy(tempNorth)
	enthalpyBottom  = temp_to_enthalpy(tempBottom)

	! Initialize all field arrays
	do k = 1, nk
	do j = 1, nj
	do i = 1, ni
		vis(i,j,k) = viscos
		den(i,j,k) = denl
		uVel(i,j,k) = 0.0
		unot(i,j,k) = 0.0
		vVel(i,j,k) = 0.0
		vnot(i,j,k) = 0.0
		wVel(i,j,k) = 0.0
		wnot(i,j,k) = 0.0
		pressure(i,j,k) = 0.0
		pp(i,j,k) = 0.0
		enthalpy(i,j,k) = enthalpyPreheat
		hnot(i,j,k) = enthalpyPreheat
		temp(i,j,k) = tempPreheat
		tnot(i,j,k) = tempPreheat
		dux(i,j,k) = 0.0
		dvy(i,j,k) = 0.0
		dwz(i,j,k) = 0.0
		su(i,j,k) = 0.0
		sp(i,j,k) = 0.0
		fracl(i,j,k) = 0.0
		fraclnot(i,j,k) = 0.0
		diff(i,j,k) = (thconsa * tempPreheat + thconsb) / (acpa * tempPreheat / 2 + acpb)
		solidfield(i,j,k) = 0
	enddo
	enddo
	enddo

	coordhistory(1:COORDLINES, 1:8) = -1
	TrackNum = 0
	PathNum = 2

	beam_pos  = toolmatrix(PathNum, 2)
	beam_posy = toolmatrix(PathNum, 3)

!---- Enthalpy boundary conditions ----
	! i=1 plane (West)
	do j = 1, nj
	do k = 1, nk
		enthalpy(1, j, k) = enthalpyWest
	enddo
	enddo

	! i=ni plane (East)
	do j = 1, nj
	do k = 1, nk
		enthalpy(ni, j, k) = enthalpyEast
	enddo
	enddo

	! k=1 plane (Bottom)
	do i = 1, ni
	do j = 1, nj
		enthalpy(i, j, 1) = enthalpyBottom
	enddo
	enddo

	! j=nj plane (North)
	do i = 1, ni
	do k = 1, nk
		enthalpy(i, nj, k) = enthalpyNorth
	enddo
	enddo

	return

end subroutine initialize

end module initialization
