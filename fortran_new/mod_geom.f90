!______________________________________________________________________________
!
module geometry
!______________________________________________________________________________
!
	use constant
	use parameters
	implicit none

	real(wp) dimx,dimy,dimz	
	!maximum coordinate values
	
	real(wp) x(nx),y(ny),z(nz),xu(nx),yv(ny),zw(nz),dxpwinv(nx),dypsinv(ny),dzpbinv(nz)
	!    x-value, y-value, z-value, u-value, v-value, z-value, 3 operaters

	real(wp) volume_u(nx,ny,nz),volume_v(nx,ny,nz),volume_w(nx,ny,nz),volume(nx,ny,nz)
	!volume of CV

	real(wp) areaij(nx,ny),areajk(ny,nz),areaik(nx,nz)
	!area of surface of scalar CV
		
	real(wp) areauij(nx,ny),areauik(nx,nz),areavjk(ny,nz),areavij(nx,ny),areawik(nx,nz),areawjk(ny,nz)
	!area of surface of velocity CV
	
	real(wp) fracx(nx),fracy(ny),fracz(nz)
	!operaters to calcu diffusion coefficient
		
	integer ni,nim1,nj,njm1,nk,nkm1 
	!node num in x-dirction, ni-1,node num in y-dirction, nj-1,node num in z-dirction, nk-1
			

	contains 

!----------------------------------------------------------------------
! Generate a 1D non-uniform grid for one coordinate direction.
! vel_grid: staggered velocity node positions (xu, yv, or zw)
! scalar_grid: scalar node positions (x, y, or z)
! n, nm1: total node count and n-1
!----------------------------------------------------------------------
subroutine generate_1d_grid(nzones, zones, ncv, powr, vel_grid, scalar_grid, n, nm1)
	integer, intent(in) :: nzones
	real(wp), intent(in) :: zones(:), powr(:)
	integer, intent(in) :: ncv(:)
	real(wp), intent(inout) :: vel_grid(:), scalar_grid(:)
	integer, intent(out) :: n, nm1
	integer :: i, j, ist
	real(wp) :: statloc, term

	vel_grid(1:2) = 0.0
	ist = 2
	statloc = 0.0
	n = 2
	do i = 1, nzones
		n = n + ncv(i)
		do j = 1, ncv(i)
			if (powr(i) .ge. 0.0) then
				term = (real(j,wp)/real(ncv(i),wp))**powr(i)
			else
				term = 1.0 - (1.0 - real(j,wp)/real(ncv(i),wp))**(-powr(i))
			endif
			vel_grid(j+ist) = statloc + zones(i)*term
		enddo
		ist = ist + ncv(i)
		statloc = statloc + zones(i)
	enddo
	nm1 = n - 1

	do i = 1, nm1
		scalar_grid(i) = (vel_grid(i+1) + vel_grid(i)) * 0.5
	enddo
	scalar_grid(n) = vel_grid(n)
end subroutine generate_1d_grid

subroutine generate_grid

	integer i,j,k

!----x grid---------------------------------------
	call generate_1d_grid(nzx, xzone, ncvx, powrx, xu, x, ni, nim1)

!-------y grids----------------------------
	call generate_1d_grid(nzy, yzone, ncvy, powry, yv, y, nj, njm1)

!-----------z grids----------------------------
	call generate_1d_grid(nzz, zzone, ncvz, powrz, zw, z, nk, nkm1)

!********************************************************************

	do i=2,ni    ! loop of uVel nodes
		dxpwinv(i)=1.0/(x(i)-x(i-1))  !reciprocal of nodes distance
	enddo

	do j=2,nj
		dypsinv(j)=1.0/(y(j)-y(j-1))
	enddo

	do k=2,nk
		dzpbinv(k)=1.0/(z(k)-z(k-1))
	enddo

!	interpolation
	do i=1,nim1
		fracx(i)=(x(i+1)-xu(i+1))/(x(i+1)-x(i))   !distance between node and interface/distance between two nodes
	enddo

	do j=1,njm1
		fracy(j)=(y(j+1)-yv(j+1))/(y(j+1)-y(j))
	enddo

	do k=1,nkm1
		fracz(k)=(z(k+1)-zw(k+1))/(z(k+1)-z(k))
	enddo

!---volumes-------------------------calcu volume of CV
	do k=2,nkm1
	do j=2,njm1
	do i=2,nim1
		volume(i,j,k)=(xu(i+1)-xu(i))*(yv(j+1)-yv(j))*(zw(k+1)-zw(k))
		volume_u(i,j,k)=(x(i)-x(i-1))*(yv(j+1)-yv(j))*(zw(k+1)-zw(k))
		volume_v(i,j,k)=(xu(i+1)-xu(i))*(y(j)-y(j-1))*(zw(k+1)-zw(k))
		volume_w(i,j,k)=(xu(i+1)-xu(i))*(yv(j+1)-yv(j))*(z(k)-z(k-1))
	enddo
	enddo
	enddo

!----------areas-----------------------calcu areas in three direction of CV
	do j=2,njm1
	do i=2,nim1
		areaij(i,j)=(xu(i+1)-xu(i))*(yv(j+1)-yv(j))
		areauij(i,j)=(x(i)-x(i-1))*(yv(j+1)-yv(j))
		areavij(i,j)=(xu(i+1)-xu(i))*(y(j)-y(j-1))
	enddo
	enddo
	
	do k=2,nkm1
	do i=2,nim1
		areaik(i,k)=(xu(i+1)-xu(i))*(zw(k+1)-zw(k))
		areawik(i,k)=(xu(i+1)-xu(i))*(z(k)-z(k-1))
		areauik(i,k)=(x(i)-x(i-1))*(zw(k+1)-zw(k))
	enddo
	enddo

	do k=2,nkm1
	do j=2,njm1
		areajk(j,k)=(yv(j+1)-yv(j))*(zw(k+1)-zw(k))
		areavjk(j,k)=(y(j)-y(j-1))*(zw(k+1)-zw(k))
		areawjk(j,k)=(yv(j+1)-yv(j))*(z(k)-z(k-1))
	enddo
	enddo

	dimx=x(ni)   !maximum of x-value
	dimy=y(nj)   !maximum of y-value
	dimz=z(nk)   !maximum of z-value

	return
end subroutine generate_grid
end module geometry
