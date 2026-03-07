
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
	use timing
	use omp_lib
	use local_enthalpy

	implicit none
	integer i,j,k
	integer step_idx
	integer ilo, ihi, jlo, jhi, klo, khi
	logical is_local
	real(wp) amaxres
	real(wp) t0, t1
	real(wp) t_step_wall

	call read_data
	call read_toolpath
	call generate_grid          ! sets ni,nj,nk and allocates geometry
	call allocate_fields(ni, nj, nk)
	call allocate_source(ni, nj, nk)
	call allocate_print(ni, nj, nk)
	call allocate_laser(ni, nj)
	call OpenFiles
	call initialize
	call init_thermal_history

	call StartTime

	itertot=0
	step_idx=0
	timet=small

!------time stepping loop------------------------------------
	time_loop: do while (timet.lt.timax)
		timet=timet+delt
		step_idx = step_idx + 1
		niter=0
		call StartStepTime

!-------Move laser and calculate freesurface-----------------
		call cpu_time(t0)
		call laser_beam
		call read_coordinates
!		call calcRHF   ! RHF disabled
		call get_enthalpy_region(step_idx, is_local, ilo, ihi, jlo, jhi, klo, khi)
		if (toolmatrix(PathNum,5) .lt. laser_on_threshold) then
			is_local = .false.
			ilo = 2; ihi = nim1
			jlo = 2; jhi = njm1
			klo = 2; khi = nkm1
		endif
		call update_localfield(ilo, ihi, jlo, jhi, klo, khi)
		call cpu_time(t1)
		t_laser = t_laser + (t1 - t0)

		if (is_local) then
			write(*,'(A,I6,A)')  '  Step', step_idx, ' => LOCAL enthalpy solve'
			write(9,'(A,I6,A)') '  Step', step_idx, ' => LOCAL enthalpy solve'
		else
			write(*,'(A,I6,A)')  '  Step', step_idx, ' => GLOBAL enthalpy solve'
			write(9,'(A,I6,A)') '  Step', step_idx, ' => GLOBAL enthalpy solve'
		endif

!-----iteration loop within each time step----------------
		iter_loop: do while (niter.lt.maxit)
			niter=niter+1
			itertot=itertot+1

!-----solve energy equation (formerly ivar=5)-----
			call cpu_time(t0)
			call properties(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_prop = t_prop + (t1 - t0)

			call cpu_time(t0)
			call bound_enthalpy(ilo, ihi, jlo, jhi, klo, khi, is_local)
			call cpu_time(t1)
			t_bound = t_bound + (t1 - t0)

			call cpu_time(t0)
			call discretize_enthalpy(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_discret = t_discret + (t1 - t0)
			t_discret_enthalpy = t_discret_enthalpy + (t1 - t0)

			call cpu_time(t0)
			call source_enthalpy(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_sour = t_sour + (t1 - t0)
			t_sour_enthalpy = t_sour_enthalpy + (t1 - t0)

			call cpu_time(t0)
			call calc_enthalpy_residual(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_resid = t_resid + (t1 - t0)

			call cpu_time(t0)
			call enhance_converge_speed(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_converge = t_converge + (t1 - t0)

			call cpu_time(t0)
			call solution_enthalpy(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_solve = t_solve + (t1 - t0)
			t_solve_enthalpy = t_solve_enthalpy + (t1 - t0)

			call cpu_time(t0)
			call enthalpy_to_temp(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_entot = t_entot + (t1 - t0)

			call cpu_time(t0)
			call pool_size
			call cpu_time(t1)
			t_dimen = t_dimen + (t1 - t0)

			if(tpeak.gt.tsolid) then
				call cpu_time(t0)
				call cleanuvw
				call cpu_time(t1)
				t_solve = t_solve + (t1 - t0)
				t_solve_cleanuvw = t_solve_cleanuvw + (t1 - t0)

!-----solve u-momentum (formerly ivar=1)-----
				call cpu_time(t0)
				call bound_uv(1)
				call cpu_time(t1)
				t_bound = t_bound + (t1 - t0)
				call cpu_time(t0)
				call discretize_momentum(1)
				call cpu_time(t1)
				t_discret = t_discret + (t1 - t0)
				t_discret_momentum = t_discret_momentum + (t1 - t0)
				call cpu_time(t0)
				call source_momentum(1)
				call cpu_time(t1)
				t_sour = t_sour + (t1 - t0)
				t_sour_momentum = t_sour_momentum + (t1 - t0)
				call cpu_time(t0)
				call calc_momentum_residual(uVel, resoru, .true.)
				call cpu_time(t1)
				t_resid = t_resid + (t1 - t0)
				call cpu_time(t0)
				call solution_uvw(uVel)
				call cpu_time(t1)
				t_solve = t_solve + (t1 - t0)
				t_solve_uvw = t_solve_uvw + (t1 - t0)

!-----solve v-momentum (formerly ivar=2)-----
				call cpu_time(t0)
				call bound_uv(2)
				call cpu_time(t1)
				t_bound = t_bound + (t1 - t0)
				call cpu_time(t0)
				call discretize_momentum(2)
				call cpu_time(t1)
				t_discret = t_discret + (t1 - t0)
				t_discret_momentum = t_discret_momentum + (t1 - t0)
				call cpu_time(t0)
				call source_momentum(2)
				call cpu_time(t1)
				t_sour = t_sour + (t1 - t0)
				t_sour_momentum = t_sour_momentum + (t1 - t0)
				call cpu_time(t0)
				call calc_momentum_residual(vVel, resorv, .false.)
				call cpu_time(t1)
				t_resid = t_resid + (t1 - t0)
				call cpu_time(t0)
				call solution_uvw(vVel)
				call cpu_time(t1)
				t_solve = t_solve + (t1 - t0)
				t_solve_uvw = t_solve_uvw + (t1 - t0)

!-----solve w-momentum (formerly ivar=3)-----
				call cpu_time(t0)
				call bound_w
				call cpu_time(t1)
				t_bound = t_bound + (t1 - t0)
				call cpu_time(t0)
				call discretize_momentum(3)
				call cpu_time(t1)
				t_discret = t_discret + (t1 - t0)
				t_discret_momentum = t_discret_momentum + (t1 - t0)
				call cpu_time(t0)
				call source_momentum(3)
				call cpu_time(t1)
				t_sour = t_sour + (t1 - t0)
				t_sour_momentum = t_sour_momentum + (t1 - t0)
				call cpu_time(t0)
				call calc_momentum_residual(wVel, resorw, .false.)
				call cpu_time(t1)
				t_resid = t_resid + (t1 - t0)
				call cpu_time(t0)
				call solution_uvw(wVel)
				call cpu_time(t1)
				t_solve = t_solve + (t1 - t0)
				t_solve_uvw = t_solve_uvw + (t1 - t0)

!-----solve pressure correction (formerly ivar=4)-----
				call cpu_time(t0)
				call bound_pp
				call cpu_time(t1)
				t_bound = t_bound + (t1 - t0)
				call cpu_time(t0)
				call discretize_pp
				call cpu_time(t1)
				t_discret = t_discret + (t1 - t0)
				t_discret_pp = t_discret_pp + (t1 - t0)
				call cpu_time(t0)
				call source_pp
				call cpu_time(t1)
				t_sour = t_sour + (t1 - t0)
				t_sour_pp = t_sour_pp + (t1 - t0)
				call cpu_time(t0)
				call calc_pressure_residual
				call cpu_time(t1)
				t_resid = t_resid + (t1 - t0)
				call cpu_time(t0)
				call solution_uvw(pp)
				call cpu_time(t1)
				t_solve = t_solve + (t1 - t0)
				t_solve_uvw = t_solve_uvw + (t1 - t0)
				call cpu_time(t0)
				call revision_p
				call cpu_time(t1)
				t_revise = t_revise + (t1 - t0)
			endif

!-----convergence criterion------------
			call cpu_time(t0)
			call heat_fluxes
			call cpu_time(t1)
			t_flux = t_flux + (t1 - t0)
			amaxres=max(resorm, resoru,resorv,resorw)

			if(toolmatrix(PathNum,5) .ge. laser_on_threshold)then
				! Laser on: transient-state criteria (heating stage)
				if (is_local) then
					if(resorh.lt.conv_res_heat) exit iter_loop
				else
					if(resorh.lt.conv_res_heat .and. ratio.le.ratio_upper .and. ratio.ge.ratio_lower) exit iter_loop
				endif
			else
				! Laser off: transient-state criteria (cooling stage)
				if(resorh.lt.conv_res_cool) exit iter_loop
			endif

		end do iter_loop

		call cpu_time(t0)
		call CalTime
		call outputres
		call cpu_time(t1)
		t_print = t_print + (t1 - t0)

		! Accumulate heating/cooling and local/global wall-clock time
		t_step_wall = omp_get_wtime() - t_step_start
		if (toolmatrix(PathNum,5) .ge. laser_on_threshold) then
			t_heating = t_heating + t_step_wall
			n_heating = n_heating + 1
		else
			t_cooling = t_cooling + t_step_wall
			n_cooling = n_cooling + 1
		endif
		if (is_local) then
			t_local_step  = t_local_step  + t_step_wall
			n_local_step  = n_local_step  + 1
		else
			t_global_step = t_global_step + t_step_wall
			n_global_step = n_global_step + 1
		endif

		call cpu_time(t0)
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
		call write_thermal_history(timet)
		call cpu_time(t1)
		t_other = t_other + (t1 - t0)

	end do time_loop

	call EndTime
	call finalize_thermal_history

	call write_timing_report(itertot, timet)
	stop
	end
