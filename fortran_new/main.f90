
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
	use defect_field

	implicit none
	integer i,j,k
	integer step_idx
	integer ilo, ihi, jlo, jhi, klo, khi
	logical is_local
	real(wp) amaxres
	real(wp) t0, t1
	real(wp) t_step_wall
	real(wp) wall_start, wall_elapsed

	call read_data
	call read_toolpath
	call generate_grid          ! sets ni,nj,nk and allocates geometry
	call allocate_fields(ni, nj, nk)
	call allocate_source(ni, nj, nk)
	call allocate_print(ni, nj, nk)
	call allocate_laser(ni, nj)
	call allocate_skipped(ni, nj, nk)
	call allocate_defect(ni, nj, nk)
	call OpenFiles
	call initialize
	call init_thermal_history

	call StartTime
	wall_start = omp_get_wtime()

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
		call update_localfield(ilo, ihi, jlo, jhi, klo, khi)
		call cpu_time(t1)
		t_laser = t_laser + (t1 - t0)

		call cpu_time(t0)
		call compute_delt_eff()
		call cpu_time(t1)
		t_skipped_mgmt = t_skipped_mgmt + (t1 - t0)

		if (is_local) then
			write(9,'(A,I6,A)') '  Step', step_idx, ' => LOCAL enthalpy solve'
		else
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
			call pool_size(ilo, ihi, jlo, jhi, klo, khi)
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
			if (.not. is_local) then
				call cpu_time(t0)
				call heat_fluxes
				call cpu_time(t1)
				t_flux = t_flux + (t1 - t0)
			else
				! Zero flux diagnostics (heat_fluxes skipped during local steps)
				flux_west=0.0_wp; flux_east=0.0_wp
				flux_top=0.0_wp;  flux_bottom=0.0_wp
				flux_north=0.0_wp; flux_south=0.0_wp
				heatout=0.0_wp; accul=0.0_wp; heatvol=0.0_wp; ratio=0.0_wp
			endif
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
		call update_skipped(ilo, ihi, jlo, jhi, klo, khi, is_local)
		call cpu_time(t1)
		t_skipped_mgmt = t_skipped_mgmt + (t1 - t0)

		call cpu_time(t0)
		call update_max_temp()
		call cpu_time(t1)
		t_defect = t_defect + (t1 - t0)

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
		if (is_local) then
			! Local step: velocity update for all cells
			!$OMP PARALLEL DO PRIVATE(i,j,k)
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
			enddo
			enddo
			enddo
			!$OMP END PARALLEL DO
			! Only update hnot/tnot/fraclnot for cells inside the local region
			!$OMP PARALLEL DO PRIVATE(i,j,k)
			do k=klo,khi
			do j=jlo,jhi
			do i=ilo,ihi
				tnot(i,j,k)=temp(i,j,k)
				hnot(i,j,k)=enthalpy(i,j,k)
				fraclnot(i,j,k)=fracl(i,j,k)
			enddo
			enddo
			enddo
			!$OMP END PARALLEL DO
			! Restore enthalpy/temp/fracl outside local region from hnot
			!$OMP PARALLEL DO PRIVATE(i,j,k)
			do k=1,nk
			do j=1,nj
			do i=1,ni
				if (i < ilo .or. i > ihi .or. j < jlo .or. j > jhi .or. k < klo .or. k > khi) then
					enthalpy(i,j,k)=hnot(i,j,k)
					temp(i,j,k)=tnot(i,j,k)
					fracl(i,j,k)=fraclnot(i,j,k)
				endif
			enddo
			enddo
			enddo
			!$OMP END PARALLEL DO
		else
			! Global step: update all cells
			!$OMP PARALLEL DO PRIVATE(i,j,k)
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
			!$OMP END PARALLEL DO
		endif
		call Cust_Out
		call write_thermal_history(timet)
		call cpu_time(t1)
		t_other = t_other + (t1 - t0)

	end do time_loop

	! Post-simulation defect analysis (before EndTime closes output file)
	call compute_defect_determ()
	call write_defect_report()

	call EndTime
	call finalize_thermal_history

	wall_elapsed = omp_get_wtime() - wall_start

	call write_timing_report(itertot, timet, wall_elapsed, file_prefix)
	stop
	end
