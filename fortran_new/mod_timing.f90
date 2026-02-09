!______________________________________________________________________________
!
module timing
!______________________________________________________________________________
! Accumulates CPU time per module and writes timing report to result/timing_report.txt
!
	use precision
	implicit none
	real(wp), save :: t_prop     = 0.0_wp   ! mod_prop (properties)
	real(wp), save :: t_bound   = 0.0_wp   ! mod_bound (bound_*)
	real(wp), save :: t_discret = 0.0_wp   ! mod_discret (discretize_*)
	real(wp), save :: t_sour    = 0.0_wp   ! mod_sour (source_*)
	real(wp), save :: t_resid   = 0.0_wp   ! mod_resid (calc_*_residual)
	real(wp), save :: t_converge= 0.0_wp   ! mod_converge (enhance_converge_speed)
	real(wp), save :: t_solve   = 0.0_wp   ! mod_solve (solution_*, cleanuvw)
	real(wp), save :: t_entot   = 0.0_wp   ! mod_entot (enthalpy_to_temp)
	real(wp), save :: t_dimen   = 0.0_wp   ! mod_dimen (pool_size)
	real(wp), save :: t_flux    = 0.0_wp   ! mod_flux (heat_fluxes)
	real(wp), save :: t_revise  = 0.0_wp   ! mod_revise (revision_p)
	real(wp), save :: t_print   = 0.0_wp   ! mod_print (CalTime, outputres, Cust_Out)
	real(wp), save :: t_laser   = 0.0_wp   ! laser/toolpath (laser_beam, read_coordinates, calcRHF)
	real(wp), save :: t_other   = 0.0_wp   ! copy loop and misc

	! Subroutine-level (for modules >10%): mod_sour, mod_discret, mod_solve
	real(wp), save :: t_sour_momentum  = 0.0_wp   ! source_momentum(1+2+3)
	real(wp), save :: t_sour_pp       = 0.0_wp   ! source_pp
	real(wp), save :: t_sour_enthalpy  = 0.0_wp   ! source_enthalpy
	real(wp), save :: t_discret_enthalpy = 0.0_wp ! discretize_enthalpy
	real(wp), save :: t_discret_momentum = 0.0_wp ! discretize_momentum(1+2+3)
	real(wp), save :: t_discret_pp    = 0.0_wp   ! discretize_pp
	real(wp), save :: t_solve_enthalpy  = 0.0_wp ! solution_enthalpy
	real(wp), save :: t_solve_uvw     = 0.0_wp   ! solution_uvw (u,v,w,pp)
	real(wp), save :: t_solve_cleanuvw = 0.0_wp  ! cleanuvw

	contains

subroutine write_timing_report(itertot, timet_end)
	integer, intent(in) :: itertot
	real(wp), intent(in) :: timet_end
	integer, parameter :: lun = 88, nmod = 14
	real(wp) :: t_total, t_sum
	real(wp) :: pct, pct_glob
	real(wp) :: t_list(nmod)
	character(len=20) :: name_list(nmod)
	integer :: idx(nmod), i, j, itmp

	t_total = t_prop + t_bound + t_discret + t_sour + t_resid + t_converge + &
	          t_solve + t_entot + t_dimen + t_flux + t_revise + t_print + t_laser + t_other
	if (t_total <= 0.0_wp) t_total = 1.0_wp

	t_list(1) = t_prop;    name_list(1) = 'mod_prop'
	t_list(2) = t_bound;   name_list(2) = 'mod_bound'
	t_list(3) = t_discret; name_list(3) = 'mod_discret'
	t_list(4) = t_sour;    name_list(4) = 'mod_sour'
	t_list(5) = t_resid;   name_list(5) = 'mod_resid'
	t_list(6) = t_converge; name_list(6) = 'mod_converge'
	t_list(7) = t_solve;   name_list(7) = 'mod_solve'
	t_list(8) = t_entot;   name_list(8) = 'mod_entot'
	t_list(9) = t_dimen;   name_list(9) = 'mod_dimen'
	t_list(10) = t_flux;   name_list(10) = 'mod_flux'
	t_list(11) = t_revise; name_list(11) = 'mod_revise'
	t_list(12) = t_print;  name_list(12) = 'mod_print'
	t_list(13) = t_laser;  name_list(13) = 'laser/toolpath'
	t_list(14) = t_other;  name_list(14) = 'other (copy/misc)'
	do i = 1, nmod
		idx(i) = i
	enddo
	! sort idx by t_list descending (bubble)
	do i = 1, nmod - 1
		do j = i + 1, nmod
			if (t_list(idx(j)) > t_list(idx(i))) then
				itmp = idx(i); idx(i) = idx(j); idx(j) = itmp
			endif
		enddo
	enddo

	open(unit=lun, file='result/timing_report.txt', action='write', status='replace')
	write(lun,'(a)') '============================================'
	write(lun,'(a)') '  AM-CFD Module Timing Report'
	write(lun,'(a)') '  (CPU time accumulated in main loop)'
	write(lun,'(a)') '============================================'
	write(lun,'(a,i0)') '  Total iterations (itertot): ', itertot
	write(lun,'(a,es15.6)') '  Simulation time reached: ', timet_end
	write(lun,'(a,f12.3,a)') '  Total CPU time: ', t_total, ' s'
	write(lun,'(a)') '--------------------------------------------'
	write(lun,'(a)') '  Module              |   Time(s)   |   Ratio(%)'
	write(lun,'(a)') '--------------------------------------------'
	do i = 1, nmod
		j = idx(i)
		pct = 100.0_wp * t_list(j) / t_total
		write(lun,'(a,a20,a,f12.3,a,f8.2,a)') '  ', trim(name_list(j)), '|', t_list(j), '  |', pct, '%'
	enddo
	write(lun,'(a)') '--------------------------------------------'
	t_sum = t_prop + t_bound + t_discret + t_sour + t_resid + t_converge + &
	        t_solve + t_entot + t_dimen + t_flux + t_revise + t_print + t_laser + t_other
	write(lun,'(a,f12.3)') '  Sum (check):         ', t_sum
	write(lun,'(a)') '============================================'
	write(lun,'(a)') ''
	write(lun,'(a)') '  Subroutine breakdown (modules with ratio >10%)'
	write(lun,'(a)') '  Time(s) = CPU time; Ratio(%) = within module; Glob%(%) = of total'
	write(lun,'(a)') '--------------------------------------------'
	if (t_sour > 0.0_wp) then
		write(lun,'(a)') '  mod_sour:'
		pct = 100.0_wp * t_sour_momentum / t_sour
		pct_glob = 100.0_wp * t_sour_momentum / t_total
		write(lun,'(a,f12.3,a,f8.2,a,f6.2,a)') '    source_momentum     ', t_sour_momentum, '  |', pct, '% mod |', pct_glob, '% glob'
		pct = 100.0_wp * t_sour_pp / t_sour
		pct_glob = 100.0_wp * t_sour_pp / t_total
		write(lun,'(a,f12.3,a,f8.2,a,f6.2,a)') '    source_pp           ', t_sour_pp, '  |', pct, '% mod |', pct_glob, '% glob'
		pct = 100.0_wp * t_sour_enthalpy / t_sour
		pct_glob = 100.0_wp * t_sour_enthalpy / t_total
		write(lun,'(a,f12.3,a,f8.2,a,f6.2,a)') '    source_enthalpy     ', t_sour_enthalpy, '  |', pct, '% mod |', pct_glob, '% glob'
		write(lun,'(a,f12.3)') '    (mod_sour sum)      ', t_sour_momentum + t_sour_pp + t_sour_enthalpy
		write(lun,'(a)') '--------------------------------------------'
	endif
	if (t_discret > 0.0_wp) then
		write(lun,'(a)') '  mod_discret:'
		pct = 100.0_wp * t_discret_enthalpy / t_discret
		pct_glob = 100.0_wp * t_discret_enthalpy / t_total
		write(lun,'(a,f12.3,a,f8.2,a,f6.2,a)') '    discretize_enthalpy ', t_discret_enthalpy, '  |', pct, '% mod |', pct_glob, '% glob'
		pct = 100.0_wp * t_discret_momentum / t_discret
		pct_glob = 100.0_wp * t_discret_momentum / t_total
		write(lun,'(a,f12.3,a,f8.2,a,f6.2,a)') '    discretize_momentum ', t_discret_momentum, '  |', pct, '% mod |', pct_glob, '% glob'
		pct = 100.0_wp * t_discret_pp / t_discret
		pct_glob = 100.0_wp * t_discret_pp / t_total
		write(lun,'(a,f12.3,a,f8.2,a,f6.2,a)') '    discretize_pp       ', t_discret_pp, '  |', pct, '% mod |', pct_glob, '% glob'
		write(lun,'(a,f12.3)') '    (mod_discret sum)   ', t_discret_enthalpy + t_discret_momentum + t_discret_pp
		write(lun,'(a)') '--------------------------------------------'
	endif
	if (t_solve > 0.0_wp) then
		write(lun,'(a)') '  mod_solve:'
		pct = 100.0_wp * t_solve_enthalpy / t_solve
		pct_glob = 100.0_wp * t_solve_enthalpy / t_total
		write(lun,'(a,f12.3,a,f8.2,a,f6.2,a)') '    solution_enthalpy   ', t_solve_enthalpy, '  |', pct, '% mod |', pct_glob, '% glob'
		pct = 100.0_wp * t_solve_uvw / t_solve
		pct_glob = 100.0_wp * t_solve_uvw / t_total
		write(lun,'(a,f12.3,a,f8.2,a,f6.2,a)') '    solution_uvw (4x)   ', t_solve_uvw, '  |', pct, '% mod |', pct_glob, '% glob'
		pct = 100.0_wp * t_solve_cleanuvw / t_solve
		pct_glob = 100.0_wp * t_solve_cleanuvw / t_total
		write(lun,'(a,f12.3,a,f8.2,a,f6.2,a)') '    cleanuvw            ', t_solve_cleanuvw, '  |', pct, '% mod |', pct_glob, '% glob'
		write(lun,'(a,f12.3)') '    (mod_solve sum)     ', t_solve_enthalpy + t_solve_uvw + t_solve_cleanuvw
		write(lun,'(a)') '--------------------------------------------'
	endif
	write(lun,'(a)') '============================================'
	close(lun)
end subroutine write_timing_report

end module timing
