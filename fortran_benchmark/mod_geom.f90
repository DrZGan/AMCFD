!______________________________________________________________________________
!
module geometry
!______________________________________________________________________________
!
	use constant
	use parameters
	implicit none

	!--------------------------------------------------------------
	! Derived type: mesh_t (encapsulates all geometry data)
	!--------------------------------------------------------------
	type :: mesh_t
		integer :: ni, nj, nk, nim1, njm1, nkm1
		real(wp) :: dimx, dimy, dimz
		! 1D grid arrays
		real(wp), allocatable :: x(:), y(:), z(:)
		real(wp), allocatable :: xu(:), yv(:), zw(:)
		real(wp), allocatable :: dxpwinv(:), dypsinv(:), dzpbinv(:)
		real(wp), allocatable :: fracx(:), fracy(:), fracz(:)
		! 3D control volumes
		real(wp), allocatable :: volume(:,:,:)
		real(wp), allocatable :: volume_u(:,:,:), volume_v(:,:,:), volume_w(:,:,:)
		! 2D face areas
		real(wp), allocatable :: areaij(:,:), areajk(:,:), areaik(:,:)
		real(wp), allocatable :: areauij(:,:), areauik(:,:)
		real(wp), allocatable :: areavjk(:,:), areavij(:,:)
		real(wp), allocatable :: areawik(:,:), areawjk(:,:)
	contains
		procedure :: init => mesh_init
		procedure :: destroy => mesh_destroy
	end type mesh_t

	!--------------------------------------------------------------
	! Module-level variables (backward compatibility — existing code uses these)
	!--------------------------------------------------------------
	real(wp) dimx,dimy,dimz

	! 1D grid arrays (allocatable, sized to actual grid)
	real(wp), allocatable :: x(:),y(:),z(:),xu(:),yv(:),zw(:)
	real(wp), allocatable :: dxpwinv(:),dypsinv(:),dzpbinv(:)

	! 3D control volumes (allocatable)
	real(wp), allocatable :: volume_u(:,:,:),volume_v(:,:,:),volume_w(:,:,:),volume(:,:,:)

	! 2D face areas (allocatable)
	real(wp), allocatable :: areaij(:,:),areajk(:,:),areaik(:,:)
	real(wp), allocatable :: areauij(:,:),areauik(:,:),areavjk(:,:),areavij(:,:),areawik(:,:),areawjk(:,:)

	! 1D interpolation fractions
	real(wp), allocatable :: fracx(:),fracy(:),fracz(:)

	integer ni,nim1,nj,njm1,nk,nkm1
	!node num in x-direction, ni-1, y-direction, nj-1, z-direction, nk-1
			

	contains

!********************************************************************
subroutine allocate_geometry(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk
	! 1D arrays
	allocate(x(nni), y(nnj), z(nnk))
	allocate(xu(nni), yv(nnj), zw(nnk))
	allocate(dxpwinv(nni), dypsinv(nnj), dzpbinv(nnk))
	allocate(fracx(nni), fracy(nnj), fracz(nnk))
	! 3D volumes
	allocate(volume(nni,nnj,nnk), volume_u(nni,nnj,nnk))
	allocate(volume_v(nni,nnj,nnk), volume_w(nni,nnj,nnk))
	! 2D areas
	allocate(areaij(nni,nnj), areajk(nnj,nnk), areaik(nni,nnk))
	allocate(areauij(nni,nnj), areauik(nni,nnk))
	allocate(areavjk(nnj,nnk), areavij(nni,nnj))
	allocate(areawik(nni,nnk), areawjk(nnj,nnk))
end subroutine allocate_geometry

!********************************************************************
subroutine deallocate_geometry()
	if (allocated(x)) deallocate(x,y,z,xu,yv,zw,dxpwinv,dypsinv,dzpbinv,fracx,fracy,fracz)
	if (allocated(volume)) deallocate(volume,volume_u,volume_v,volume_w)
	if (allocated(areaij)) deallocate(areaij,areajk,areaik,areauij,areauik,areavjk,areavij,areawik,areawjk)
end subroutine deallocate_geometry

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

!-----compute grid sizes first (ni,nj,nk)-----
	ni = 2; do i=1,nzx; ni = ni + ncvx(i); enddo; nim1 = ni - 1
	nj = 2; do i=1,nzy; nj = nj + ncvy(i); enddo; njm1 = nj - 1
	nk = 2; do i=1,nzz; nk = nk + ncvz(i); enddo; nkm1 = nk - 1

!-----allocate geometry arrays-----
	call allocate_geometry(ni, nj, nk)

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

!********************************************************************
subroutine mesh_init(self, nni, nnj, nnk)
	class(mesh_t), intent(inout) :: self
	integer, intent(in) :: nni, nnj, nnk
	self%ni = nni; self%nj = nnj; self%nk = nnk
	self%nim1 = nni-1; self%njm1 = nnj-1; self%nkm1 = nnk-1
	allocate(self%x(nni), self%y(nnj), self%z(nnk))
	allocate(self%xu(nni), self%yv(nnj), self%zw(nnk))
	allocate(self%dxpwinv(nni), self%dypsinv(nnj), self%dzpbinv(nnk))
	allocate(self%fracx(nni), self%fracy(nnj), self%fracz(nnk))
	allocate(self%volume(nni,nnj,nnk), self%volume_u(nni,nnj,nnk))
	allocate(self%volume_v(nni,nnj,nnk), self%volume_w(nni,nnj,nnk))
	allocate(self%areaij(nni,nnj), self%areajk(nnj,nnk), self%areaik(nni,nnk))
	allocate(self%areauij(nni,nnj), self%areauik(nni,nnk))
	allocate(self%areavjk(nnj,nnk), self%areavij(nni,nnj))
	allocate(self%areawik(nni,nnk), self%areawjk(nnj,nnk))
end subroutine mesh_init

!********************************************************************
subroutine mesh_destroy(self)
	class(mesh_t), intent(inout) :: self
	if (allocated(self%x)) deallocate(self%x,self%y,self%z)
	if (allocated(self%xu)) deallocate(self%xu,self%yv,self%zw)
	if (allocated(self%dxpwinv)) deallocate(self%dxpwinv,self%dypsinv,self%dzpbinv)
	if (allocated(self%fracx)) deallocate(self%fracx,self%fracy,self%fracz)
	if (allocated(self%volume)) deallocate(self%volume,self%volume_u,self%volume_v,self%volume_w)
	if (allocated(self%areaij)) deallocate(self%areaij,self%areajk,self%areaik)
	if (allocated(self%areauij)) deallocate(self%areauij,self%areauik)
	if (allocated(self%areavjk)) deallocate(self%areavjk,self%areavij)
	if (allocated(self%areawik)) deallocate(self%areawik,self%areawjk)
end subroutine mesh_destroy

end module geometry
