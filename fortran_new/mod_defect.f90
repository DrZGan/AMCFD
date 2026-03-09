!______________________________________________________________________________
!
module defect_field
!______________________________________________________________________________
! Defect detection module: computes max_temp and defect arrays over one layer.
! Supports multiple methods (currently: maxtemp_determ, maxtemp_stochas placeholder).
!
	use precision
	use geometry
	use initialization
	use parameters
	use sim_state
	use constant

	implicit none

	! --- maxtemp_determ parameters ---
	real(wp), parameter :: k_lof = 0.9_wp    ! lack-of-fusion calibration (range 0-1)
	real(wp), parameter :: k_kep = 1.0_wp    ! keyhole porosity scaling

	! --- Arrays (full X-Y, limited Z range) ---
	real(wp), allocatable :: max_temp(:,:,:)  ! (ni, nj, nk) - only k_def_lo:k_def_hi active
	real(wp), allocatable :: defect_arr(:,:,:)

	! --- Z range for defect computation ---
	integer :: k_def_lo, k_def_hi

	! --- Laser scanning range (X-Y) ---
	real(wp) :: x_scan_min, x_scan_max, y_scan_min, y_scan_max

	contains

!********************************************************************
subroutine allocate_defect(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk
	integer :: k

	allocate(max_temp(nni, nnj, nnk))
	allocate(defect_arr(nni, nnj, nnk))
	max_temp = tempPreheat
	defect_arr = 0.0_wp

	! Determine k range: from top (nkm1) down to z(nk) - layerheight
	k_def_hi = nkm1
	k_def_lo = nkm1
	do k = nkm1, 2, -1
		if (z(k) < z(nk) - layerheight) exit
		k_def_lo = k
	enddo

end subroutine allocate_defect

!********************************************************************
subroutine update_max_temp()
! Called every time step: update max_temp where current temp exceeds it.
! Only processes the single-layer Z range.
	integer :: i, j, k

	!$OMP PARALLEL DO PRIVATE(i, j, k)
	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		if (temp(i,j,k) > max_temp(i,j,k)) then
			max_temp(i,j,k) = temp(i,j,k)
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

end subroutine update_max_temp

!********************************************************************
subroutine compute_defect_determ()
! Called after simulation completes: compute defect array from max_temp.
	integer :: i, j, k

	! Step 1: Compute defect values
	!$OMP PARALLEL DO PRIVATE(i, j, k)
	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		if (max_temp(i,j,k) < tsolid) then
			defect_arr(i,j,k) = -k_lof
		else if (max_temp(i,j,k) > tboiling) then
			defect_arr(i,j,k) = k_kep * (max_temp(i,j,k) - tboiling) / tboiling
		else
			defect_arr(i,j,k) = 0.0_wp
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

	! Step 2: Determine max laser scanning range from toolpath
	call compute_scan_range()

	! Step 3: Clean up — zero defect outside scanning range
	!$OMP PARALLEL DO PRIVATE(i, j, k)
	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		if (x(i) < x_scan_min .or. x(i) > x_scan_max .or. &
		    y(j) < y_scan_min .or. y(j) > y_scan_max) then
			defect_arr(i,j,k) = 0.0_wp
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

end subroutine compute_defect_determ

!********************************************************************
subroutine maxtemp_stochas()
! Placeholder for future stochastic defect method.
! Currently does nothing.
end subroutine maxtemp_stochas

!********************************************************************
subroutine compute_scan_range()
! Read toolpath to find max laser scanning range in X-Y plane.
	integer :: n

	x_scan_min =  1.0e30_wp
	x_scan_max = -1.0e30_wp
	y_scan_min =  1.0e30_wp
	y_scan_max = -1.0e30_wp

	do n = 1, TOOLLINES
		if (toolmatrix(n,1) < -0.5_wp) exit   ! end of toolpath data
		if (toolmatrix(n,5) >= laser_on_threshold) then
			x_scan_min = min(x_scan_min, toolmatrix(n,2))
			x_scan_max = max(x_scan_max, toolmatrix(n,2))
			y_scan_min = min(y_scan_min, toolmatrix(n,3))
			y_scan_max = max(y_scan_max, toolmatrix(n,3))
		endif
	enddo

end subroutine compute_scan_range

!********************************************************************
subroutine write_defect_report()
! Compute defect metrics and write defect_report.txt and VTK output.
	integer :: i, j, k
	integer, parameter :: lun = 89
	real(wp) :: v_cell, v_total, v_defect, v_lof, v_kep
	real(wp) :: frac_defect, frac_lof, frac_kep
	real(kind=4) :: val4
	integer :: gridx, gridy, gridk, npts
	character(len=3) :: cTemp

	! --- Compute metrics over scan range × layer height ---
	v_total  = 0.0_wp
	v_defect = 0.0_wp
	v_lof    = 0.0_wp
	v_kep    = 0.0_wp

	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		if (x(i) < x_scan_min .or. x(i) > x_scan_max .or. &
		    y(j) < y_scan_min .or. y(j) > y_scan_max) cycle
		v_cell = volume(i,j,k)
		v_total = v_total + v_cell
		if (defect_arr(i,j,k) < 0.0_wp) then
			v_lof = v_lof + abs(defect_arr(i,j,k)) * v_cell
		else if (defect_arr(i,j,k) > 0.0_wp .and. defect_arr(i,j,k) <= 1.0_wp) then
			v_kep = v_kep + defect_arr(i,j,k) * v_cell
		endif
	enddo
	enddo
	enddo

	v_defect = v_lof + v_kep
	if (v_total > 0.0_wp) then
		frac_defect = v_defect / v_total
		frac_lof    = v_lof / v_total
		frac_kep    = v_kep / v_total
	else
		frac_defect = 0.0_wp
		frac_lof    = 0.0_wp
		frac_kep    = 0.0_wp
	endif

	! --- Write defect_report.txt ---
	open(unit=lun, file=trim(file_prefix)//'defect_report.txt', action='write', status='replace')
	write(lun,'(a)') '============================================'
	write(lun,'(a)') '  AM-CFD Defect Report'
	write(lun,'(a)') '  Method: maxtemp_determ'
	write(lun,'(a)') '============================================'
	write(lun,'(a)')
	write(lun,'(a,f10.6,a)') '  Defect fraction:           ', frac_defect * 100.0_wp, ' %'
	write(lun,'(a,f10.6,a)') '  Lack-of-fusion fraction:   ', frac_lof * 100.0_wp, ' %'
	write(lun,'(a,f10.6,a)') '  Keyhole porosity fraction:  ', frac_kep * 100.0_wp, ' %'
	write(lun,'(a)')
	write(lun,'(a,es12.5,a)') '  Defect volume:             ', v_defect, ' m^3'
	write(lun,'(a,es12.5,a)') '  Lack-of-fusion volume:     ', v_lof, ' m^3'
	write(lun,'(a,es12.5,a)') '  Keyhole porosity volume:    ', v_kep, ' m^3'
	write(lun,'(a,es12.5,a)') '  Total reference volume:     ', v_total, ' m^3'
	write(lun,'(a)')
	write(lun,'(a)') '  Max laser scanning range (X-Y plane):'
	write(lun,'(a,es12.5,a,es12.5,a)') '    X: [', x_scan_min, ', ', x_scan_max, '] m'
	write(lun,'(a,es12.5,a,es12.5,a)') '    Y: [', y_scan_min, ', ', y_scan_max, '] m'
	write(lun,'(a)')
	write(lun,'(a)') '  Layer Z range:'
	write(lun,'(a,i4,a,i4)') '    k indices: ', k_def_lo, ' to ', k_def_hi
	write(lun,'(a,es12.5,a,es12.5,a)') '    Z: [', z(k_def_lo), ', ', z(k_def_hi), '] m'
	write(lun,'(a,es12.5,a)') '    Layer height: ', layerheight, ' m'
	write(lun,'(a)')
	write(lun,'(a)') '  Parameters:'
	write(lun,'(a,f8.4)') '    k_lof = ', k_lof
	write(lun,'(a,f8.4)') '    k_kep = ', k_kep
	write(lun,'(a,f10.2,a)') '    T_solid   = ', tsolid, ' K'
	write(lun,'(a,f10.2,a)') '    T_boiling = ', tboiling, ' K'
	write(lun,'(a)') '============================================'
	close(lun)

	! Print summary to output file
	write(9,'(a)') ''
	write(9,'(a)') '  === Defect Analysis (maxtemp_determ) ==='
	write(9,'(a,f10.6,a)') '  Defect fraction:         ', frac_defect * 100.0_wp, ' %'
	write(9,'(a,f10.6,a)') '  Lack-of-fusion fraction: ', frac_lof * 100.0_wp, ' %'
	write(9,'(a,f10.6,a)') '  Keyhole porosity fraction:', frac_kep * 100.0_wp, ' %'

	! --- Write VTK files for max_temp and defect ---
	call write_defect_vtk('maxtemp', max_temp)
	call write_defect_vtk('defect', defect_arr)

end subroutine write_defect_report

!********************************************************************
subroutine write_defect_vtk(fieldname, field)
	character(len=*), intent(in) :: fieldname
	real(wp), intent(in) :: field(:,:,:)
	integer :: i, j, k, npts
	integer :: gridx, gridy, gridz
	real(kind=4) :: val4
	integer, parameter :: lun = 90

	gridx = nim1 - 2 + 1   ! 2:nim1
	gridy = njm1 - 2 + 1
	gridz = k_def_hi - k_def_lo + 1
	npts = gridx * gridy * gridz

	! ASCII header
	open(unit=lun, file=trim(file_prefix)//trim(fieldname)//'.vtk')
	write(lun,'(A)') '# vtk DataFile Version 3.0'
	write(lun,'(A)') 'AMCFD '//trim(fieldname)//' field'
	write(lun,'(A)') 'BINARY'
	write(lun,'(A)') 'DATASET STRUCTURED_GRID'
	write(lun,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', gridx, ' ', gridy, ' ', gridz
	write(lun,'(A,I0,A)') 'POINTS ', npts, ' float'
	close(lun)

	! Binary coordinates
	open(unit=lun, file=trim(file_prefix)//trim(fieldname)//'.vtk', &
	     access='stream', form='unformatted', position='append', convert='big_endian')
	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		val4 = real(x(i), 4); write(lun) val4
		val4 = real(y(j), 4); write(lun) val4
		val4 = real(z(k), 4); write(lun) val4
	enddo
	enddo
	enddo
	close(lun)

	! POINT_DATA header
	open(unit=lun, file=trim(file_prefix)//trim(fieldname)//'.vtk', position='append')
	write(lun,'(A,I0)') 'POINT_DATA ', npts
	close(lun)

	! Binary field data
	open(unit=lun, file=trim(file_prefix)//trim(fieldname)//'.vtk', &
	     access='stream', form='unformatted', position='append', convert='big_endian')
	write(lun) char(10), 'SCALARS '//trim(fieldname)//' float 1', char(10), 'LOOKUP_TABLE default', char(10)
	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		val4 = real(field(i,j,k), 4)
		write(lun) val4
	enddo
	enddo
	enddo
	close(lun)

end subroutine write_defect_vtk

end module defect_field
