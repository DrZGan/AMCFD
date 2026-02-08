
!****************************************************************************
!
!  AM-CFD
!
!  Version 2.0 (Refactored)
!
!  July, 2020, With English Annotation, Northwestern University
!  Feb, 2026, Refactored for maintainability
!
!****************************************************************************

program main
	use geometry
	use initialization
	use discretization
	use dimensions
	use boundary
	use source
	use residue
	use solver
	use fluxes
	use parameters
	use printing
	use constant
	use entotemp
	use convergence
	use property
	use revision
	use laserinput
	use toolpath

	implicit none
	integer i, j, k
	real amaxres

	call read_data
	call read_toolpath

	call generate_grid
	call OpenFiles
	call initialize

	call StartTime

	itertot = 0
	timet = small

!------ Time stepping loop ------
10	timet = timet + delt

	niter = 0

	call laser_beam
	call read_coordinates
	call calcRHF

!------ Inner iteration loop ------
30	niter = niter + 1
	itertot = itertot + 1

!----- Solve energy equation (ivar=5) -----
	ivar = 5

	call properties
	call bound_condition
	call discretize
	call source_term
	call residual

	call enhance_converge_speed
	call solution_enthalpy

	call enthalpy_to_temp

	call pool_size
	if (tpeak .le. tsolid) goto 41
	call cleanuvw

!----- Solve momentum equations (ivar=1..4) -----
	do ivar = 1, 4
		call bound_condition
		call discretize
		call source_term
		call residual
		call solution_uvw
		call revision_p
	enddo
41	continue

!----- Convergence check -----
	call heat_fluxes
	amaxres = max(resorm, resoru, resorv, resorw)

	if (toolmatrix(PathNum, 5) .ge. 0.5) then
		if (amaxres .lt. CONV_HEAT_THRESHOLD .and. ratio .le. ENERGY_RATIO_HIGH .and. ratio .ge. ENERGY_RATIO_LOW) goto 50
	else
		if (resorh .lt. CONV_COOL_THRESHOLD) goto 50
	endif

	if (niter .lt. maxit) goto 30
50	continue

	call CalTime
	call outputres

	do k = 1, nk
	do j = 1, nj
	do i = 1, ni
		if (temp(i,j,k) .le. tsolid) then
			uVel(i,j,k) = 0.0
			vVel(i,j,k) = 0.0
			wVel(i,j,k) = 0.0
		endif
		unot(i,j,k)     = uVel(i,j,k)
		vnot(i,j,k)     = vVel(i,j,k)
		wnot(i,j,k)     = wVel(i,j,k)
		tnot(i,j,k)     = temp(i,j,k)
		hnot(i,j,k)     = enthalpy(i,j,k)
		fraclnot(i,j,k) = fracl(i,j,k)
	enddo
	enddo
	enddo
	call Cust_Out

	if (timet .lt. timax) goto 10

	call EndTime

	stop
	end
