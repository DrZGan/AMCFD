
!****************************************************************************
!
!  AM-CFD
!
!  Version 1.0
!
!  July, 2020, With English Annotation, Northwestern University
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
	integer i,j,k
	real amaxres

	call read_data
	call read_toolpath

	call generate_grid
	call OpenFiles
	call initialize

	call StartTime

	itertot=0
	timet=small

!------time stepping loop------------------------------------
	time_loop: do while (timet.lt.timax)
		timet=timet+delt
		niter=0

!-------Move laser and calculate freesurface-----------------
		call laser_beam
		call read_coordinates
		call calcRHF

!-----iteration loop within each time step----------------
		iter_loop: do while (niter.lt.maxit)
			niter=niter+1
			itertot=itertot+1

!-----ivar=5------------solve energy equation
			ivar=5

			call properties
			call bound_condition
			call discretize
			call source_term
			call residual

			call enhance_converge_speed
			call solution_enthalpy

			call enthalpy_to_temp

			call pool_size
			if(tpeak.le.tsolid) exit iter_loop
			call cleanuvw

!-----ivar=1,4------------solve momentum equations
			do ivar=1,4
				call bound_condition
				call discretize
				call source_term
				call residual
				call solution_uvw
				call revision_p
			enddo

!-----convergence criterion------------
			call heat_fluxes
			amaxres=max(resorm, resoru,resorv,resorw)

			if(toolmatrix(PathNum,5) .ge. 0.5)then
				! Laser on: transient-state criteria (heating stage)
				if(amaxres.lt.5.0e-4 .and. ratio.le.1.01.and.ratio.ge.0.99) exit iter_loop
			else
				! Laser off: transient-state criteria (cooling stage)
				if(resorh.lt.5.0e-7) exit iter_loop
			endif

		end do iter_loop

		call CalTime
		call outputres

		do k=1,nk
		do j=1,nj
		do i=1,ni
			if(temp(i,j,k).le.tsolid) then
				uVel(i,j,k)=0.0
				vVel(i,j,k)=0.0
				wVel(i,j,k)=0.0
			endif
			unot(i,j,k)=uVel(i,j,k)
			vnot(i,j,k)=vVel(i,j,k)
			wnot(i,j,k)=wVel(i,j,k)
			tnot(i,j,k)=temp(i,j,k)
			hnot(i,j,k)=enthalpy(i,j,k)
			fraclnot(i,j,k)=fracl(i,j,k)
		enddo
		enddo
		enddo
		call Cust_Out

	end do time_loop

	call EndTime

	stop
	end
