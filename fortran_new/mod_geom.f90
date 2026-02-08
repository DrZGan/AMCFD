!______________________________________________________________________________
!
module geometry
!______________________________________________________________________________
!
	use constant
	use parameters
	implicit none

	real dimx, dimy, dimz
	!maximum coordinate values

	real x(nx), y(ny), z(nz), xu(nx), yv(ny), zw(nz), dxpwinv(nx), dypsinv(ny), dzpbinv(nz)

	real volume_u(nx,ny,nz), volume_v(nx,ny,nz), volume_w(nx,ny,nz), volume(nx,ny,nz)

	real areaij(nx,ny), areajk(ny,nz), areaik(nx,nz)

	real areauij(nx,ny), areauik(nx,nz), areavjk(ny,nz), areavij(nx,ny), areawik(nx,nz), areawjk(ny,nz)

	real fracx(nx), fracy(ny), fracz(nz)

	integer ni, nim1, nj, njm1, nk, nkm1

	contains

!--------------------------------------------------------------------------
! Generate 1D non-uniform grid for one direction (reusable for x, y, z)
!   nzone     : number of zones
!   zones     : length of each zone
!   ncv       : number of CVs in each zone
!   powr      : power-law exponent for each zone
!   face_coord: output face coordinates (xu/yv/zw)
!   node_coord: output node coordinates (x/y/z)
!   n_total   : output total node count
!   nm1       : output n_total - 1
!--------------------------------------------------------------------------
subroutine generate_1d_grid(nzone, zones, ncv, powr, face_coord, node_coord, n_total, nm1)
	integer, intent(in) :: nzone
	real, intent(in) :: zones(:), powr(:)
	integer, intent(in) :: ncv(:)
	real, intent(out) :: face_coord(:), node_coord(:)
	integer, intent(out) :: n_total, nm1

	integer :: i, j, ist
	real :: statloc, term

	face_coord(1:2) = 0.0
	ist = 2
	statloc = 0.0
	n_total = 2

	do i = 1, nzone
		n_total = n_total + ncv(i)
		do j = 1, ncv(i)
			if (powr(i) .ge. 0.0) then
				term = (real(j) / real(ncv(i)))**powr(i)
			else
				term = 1.0 - (1.0 - real(j) / real(ncv(i)))**(-powr(i))
			endif
			face_coord(j + ist) = statloc + zones(i) * term
		enddo
		ist = ist + ncv(i)
		statloc = statloc + zones(i)
	enddo

	nm1 = n_total - 1

	do i = 1, nm1
		node_coord(i) = (face_coord(i+1) + face_coord(i)) * 0.5
	enddo
	node_coord(n_total) = face_coord(n_total)

end subroutine generate_1d_grid

!--------------------------------------------------------------------------
! Generate the full 3D staggered grid
!--------------------------------------------------------------------------
subroutine generate_grid

	integer i, j, k

!---- Generate 1D grids for each direction using common subroutine ----
	call generate_1d_grid(nzx, xzone, ncvx, powrx, xu, x, ni, nim1)
	call generate_1d_grid(nzy, yzone, ncvy, powry, yv, y, nj, njm1)
	call generate_1d_grid(nzz, zzone, ncvz, powrz, zw, z, nk, nkm1)

!---- Inverse distance operators ----
	do i = 2, ni
		dxpwinv(i) = 1.0 / (x(i) - x(i-1))
	enddo

	do j = 2, nj
		dypsinv(j) = 1.0 / (y(j) - y(j-1))
	enddo

	do k = 2, nk
		dzpbinv(k) = 1.0 / (z(k) - z(k-1))
	enddo

!---- Interpolation fractions ----
	do i = 1, nim1
		fracx(i) = (x(i+1) - xu(i+1)) / (x(i+1) - x(i))
	enddo

	do j = 1, njm1
		fracy(j) = (y(j+1) - yv(j+1)) / (y(j+1) - y(j))
	enddo

	do k = 1, nkm1
		fracz(k) = (z(k+1) - zw(k+1)) / (z(k+1) - z(k))
	enddo

!---- Volumes ----
	do k = 2, nkm1
	do j = 2, njm1
	do i = 2, nim1
		volume(i,j,k)   = (xu(i+1) - xu(i)) * (yv(j+1) - yv(j)) * (zw(k+1) - zw(k))
		volume_u(i,j,k) = (x(i) - x(i-1))   * (yv(j+1) - yv(j)) * (zw(k+1) - zw(k))
		volume_v(i,j,k) = (xu(i+1) - xu(i))  * (y(j) - y(j-1))   * (zw(k+1) - zw(k))
		volume_w(i,j,k) = (xu(i+1) - xu(i))  * (yv(j+1) - yv(j)) * (z(k) - z(k-1))
	enddo
	enddo
	enddo

!---- Areas ----
	do j = 2, njm1
	do i = 2, nim1
		areaij(i,j)  = (xu(i+1) - xu(i)) * (yv(j+1) - yv(j))
		areauij(i,j) = (x(i) - x(i-1))   * (yv(j+1) - yv(j))
		areavij(i,j) = (xu(i+1) - xu(i))  * (y(j) - y(j-1))
	enddo
	enddo

	do k = 2, nkm1
	do i = 2, nim1
		areaik(i,k)  = (xu(i+1) - xu(i)) * (zw(k+1) - zw(k))
		areawik(i,k) = (xu(i+1) - xu(i)) * (z(k) - z(k-1))
		areauik(i,k) = (x(i) - x(i-1))   * (zw(k+1) - zw(k))
	enddo
	enddo

	do k = 2, nkm1
	do j = 2, njm1
		areajk(j,k)  = (yv(j+1) - yv(j)) * (zw(k+1) - zw(k))
		areavjk(j,k) = (y(j) - y(j-1))   * (zw(k+1) - zw(k))
		areawjk(j,k) = (yv(j+1) - yv(j)) * (z(k) - z(k-1))
	enddo
	enddo

	dimx = x(ni)
	dimy = y(nj)
	dimz = z(nk)

	return
end subroutine generate_grid

end module geometry
