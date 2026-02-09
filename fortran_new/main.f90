
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
	real(wp) amaxres

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

!-----solve energy equation (formerly ivar=5)-----
			call properties
			call bound_enthalpy
			call discretize_enthalpy
			call source_enthalpy
			call calc_enthalpy_residual

			call enhance_converge_speed
			call solution_enthalpy

			call enthalpy_to_temp

			call pool_size
			if(tpeak.gt.tsolid) then
				call cleanuvw

!-----solve u-momentum (formerly ivar=1)-----
				call bound_uv(1)
				call discretize_momentum(1)
				call source_momentum(1)
				call calc_momentum_residual(uVel, resoru, .true.)
				call solution_uvw(uVel)

!-----solve v-momentum (formerly ivar=2)-----
				call bound_uv(2)
				call discretize_momentum(2)
				call source_momentum(2)
				call calc_momentum_residual(vVel, resorv, .false.)
				call solution_uvw(vVel)

!-----solve w-momentum (formerly ivar=3)-----
				call bound_w
				call discretize_momentum(3)
				call source_momentum(3)
				call calc_momentum_residual(wVel, resorw, .false.)
				call solution_uvw(wVel)

!-----solve pressure correction (formerly ivar=4)-----
				call bound_pp
				call discretize_pp
				call source_pp
				call calc_pressure_residual
				call solution_uvw(pp)
				call revision_p
			endif

!-----convergence criterion------------
			call heat_fluxes
			amaxres=max(resorm, resoru,resorv,resorw)

			if(toolmatrix(PathNum,5) .ge. laser_on_threshold)then
				! Laser on: transient-state criteria (heating stage)
				if(amaxres.lt.conv_res_heat .and. ratio.le.ratio_upper .and. ratio.ge.ratio_lower) exit iter_loop
			else
				! Laser off: transient-state criteria (cooling stage)
				if(resorh.lt.conv_res_cool) exit iter_loop
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
