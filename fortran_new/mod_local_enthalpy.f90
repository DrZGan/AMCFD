!______________________________________________________________________________
!
module local_enthalpy
!______________________________________________________________________________
! Time-step-level local/global enthalpy scheduling and local index management.
!
	use precision
	use geometry
	use parameters
	use sim_state
	use field_data
	implicit none

	contains

!********************************************************************
subroutine get_enthalpy_region(step_idx, is_local, ilo, ihi, jlo, jhi, klo, khi)
	integer, intent(in) :: step_idx
	logical, intent(out) :: is_local
	integer, intent(out) :: ilo, ihi, jlo, jhi, klo, khi

	integer :: period, phase

	! Default: global region
	ilo = 2; ihi = nim1
	jlo = 2; jhi = njm1
	klo = 2; khi = nkm1

	if (localnum <= 0) then
		is_local = .false.
		return
	endif

	period = localnum + 1
	phase = mod(step_idx - 1, period) + 1
	is_local = (phase <= localnum)

	if (.not. is_local) return

	call compute_local_region(ilo, ihi, jlo, jhi, klo, khi)
end subroutine get_enthalpy_region

!********************************************************************
subroutine compute_local_region(ilo, ihi, jlo, jhi, klo, khi)
	integer, intent(out) :: ilo, ihi, jlo, jhi, klo, khi

	integer :: i, j, k
	integer :: best_i, best_j
	real(wp) :: xmin, xmax, ymin, ymax, dmin
	logical :: found_i, found_j, found_k

	xmin = beam_pos  - local_half_x
	xmax = beam_pos  + local_half_x
	ymin = beam_posy - local_half_y
	ymax = beam_posy + local_half_y

	found_i = .false.
	ilo = nim1
	ihi = 2
	do i = 2, nim1
		if (x(i) >= xmin .and. x(i) <= xmax) then
			found_i = .true.
			ilo = min(ilo, i)
			ihi = max(ihi, i)
		endif
	enddo
	if (.not. found_i) then
		best_i = 2
		dmin = abs(x(2) - beam_pos)
		do i = 3, nim1
			if (abs(x(i) - beam_pos) < dmin) then
				dmin = abs(x(i) - beam_pos)
				best_i = i
			endif
		enddo
		ilo = best_i
		ihi = best_i
	endif

	found_j = .false.
	jlo = njm1
	jhi = 2
	do j = 2, njm1
		if (y(j) >= ymin .and. y(j) <= ymax) then
			found_j = .true.
			jlo = min(jlo, j)
			jhi = max(jhi, j)
		endif
	enddo
	if (.not. found_j) then
		best_j = 2
		dmin = abs(y(2) - beam_posy)
		do j = 3, njm1
			if (abs(y(j) - beam_posy) < dmin) then
				dmin = abs(y(j) - beam_posy)
				best_j = j
			endif
		enddo
		jlo = best_j
		jhi = best_j
	endif

	! Beam center z is top surface (z(nk)); local_depth_z extends downward.
	found_k = .false.
	klo = nkm1
	khi = nkm1
	do k = 2, nkm1
		if (z(nk) - z(k) <= local_depth_z) then
			klo = min(klo, k)
			found_k = .true.
		endif
	enddo
	if (.not. found_k) klo = nkm1

	ilo = max(2, min(ilo, nim1))
	ihi = max(2, min(ihi, nim1))
	jlo = max(2, min(jlo, njm1))
	jhi = max(2, min(jhi, njm1))
	klo = max(2, min(klo, nkm1))
	khi = max(2, min(khi, nkm1))

	if (ilo > ihi) then
		i = ilo; ilo = ihi; ihi = i
	endif
	if (jlo > jhi) then
		j = jlo; jlo = jhi; jhi = j
	endif
	if (klo > khi) then
		k = klo; klo = khi; khi = k
	endif
end subroutine compute_local_region

!********************************************************************
! Sets localfield for VTK: 0 = outside local cube, 1 = inside local cube,
! 2 = laser heating zone (override). Uses geometric local cube, not current
! step solve region, so 0/1/2 are visible regardless of local/global step.
subroutine update_localfield(ilo, ihi, jlo, jhi, klo, khi)
	integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
	integer :: i, j, k
	integer :: lilo, lihi, ljlo, ljhi, lklo, lkhi

	if (.not. allocated(localfield)) return

	localfield = 0.0_wp

	! 1 in geometric local cube (for visualization), independent of localnum.
	call compute_local_region(lilo, lihi, ljlo, ljhi, lklo, lkhi)
	!$OMP PARALLEL DO
	do k = lklo, lkhi
	do j = ljlo, ljhi
	do i = lilo, lihi
		localfield(i,j,k) = 1.0_wp
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

	! 2 in laser heating region (override 0 or 1)
	!$OMP PARALLEL DO
	do k = 2, nkm1
	do j = 2, njm1
	do i = 2, nim1
		if (toolmatrix(PathNum,5) > laser_on_threshold) then
			if (z(nk)-z(k) <= sourcedepth .and. abs(x(i)-beam_pos) <= sourcerad .and. abs(y(j)-beam_posy) <= sourcerad) then
				localfield(i,j,k) = 2.0_wp
			endif
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine update_localfield

end module local_enthalpy
