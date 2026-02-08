!______________________________________________________________________________
!
module geometry
!______________________________________________________________________________
!
	use constant
	use parameters
	implicit none

	real dimx,dimy,dimz	
	!maximum coordinate values
	
	real x(nx),y(ny),z(nz),xu(nx),yv(ny),zw(nz),dxpwinv(nx),dypsinv(ny),dzpbinv(nz)
	!    x-value, y-value, z-value, u-value, v-value, z-value, 3 operaters

	real volume_u(nx,ny,nz),volume_v(nx,ny,nz),volume_w(nx,ny,nz),volume(nx,ny,nz)
	!volume of CV

	real areaij(nx,ny),areajk(ny,nz),areaik(nx,nz)
	!area of surface of scalar CV
		
	real areauij(nx,ny),areauik(nx,nz),areavjk(ny,nz),areavij(nx,ny),areawik(nx,nz),areawjk(ny,nz)
	!area of surface of velocity CV
	
	real fracx(nx),fracy(ny),fracz(nz)
	!operaters to calcu diffusion coefficient
		
	integer ni,nim1,nj,njm1,nk,nkm1 
	!node num in x-dirction, ni-1,node num in y-dirction, nj-1,node num in z-dirction, nk-1
			

	contains 

subroutine generate_grid

	integer i,j,k,ist
	real statloc,term

!----x grid---------------------------------------
	xu(1:2) = 0.0        !
	ist = 2
	statloc = 0
	ni=2
	do i=1, nzx           ! loop in sub-region
		ni=ni+ncvx(i)  !calcu node num in each sub-region
		do j=1, ncvx(i)
			if(powrx(i).ge.0.0)then   !prow>=0
				term=(real(j)/real(ncvx(i)))**powrx(i)    !grids transit from fine to coarse
			else
				term=1.0-(1.0-real(j)/real(ncvx(i)))**(-powrx(i)) !grids transit from coarse to fine
			endif
			xu(j+ist) = statloc + xzone(i)*term
		enddo
		ist = ist + ncvx(i)
		statLoc = statLoc + xzone(i)
	enddo
	nim1=ni-1          !coordinate value of nim1 scalar nodes can be interpolately got from that of ni uVel nodes

	do i=1,nim1
		x(i)=(xu(i+1)+xu(i))*0.5  ! central interpolation
	enddo
	x(ni)=xu(ni)     ! last coordinate value of scalar node equals to that of uVel nodes

!-------y grids----------------------------
	yv(1:2) = 0.0
	ist = 2
	statloc = 0.0
	nj=2
	do i=1, nzy
		nj=nj+ncvy(i)
		do j=1, ncvy(i)
			if(powry(i).ge.0.0)then
				term=(real(j)/real(ncvy(i)))**powry(i)
			else
				term=1.0-(1.0-real(j)/real(ncvy(i)))**(-powry(i))
			endif
			yv(j+ist) = statloc + yzone(i)*term
		enddo
		ist = ist + ncvy(i)
		statLoc = statLoc + yzone(i)
	enddo
	njm1=nj-1
	
	do i=1,njm1
		y(i)=(yv(i+1)+yv(i))*0.5
	enddo
	y(nj)=yv(nj)

!-----------z grids----------------------------
	zw(1:2) = 0.0
	ist = 2
	statloc = 0.0
	nk=2
	do i = 1,nzz
		nk=nk+ncvz(i)
		do j = 1, ncvz(i)
			if(powrz(i).ge.0.0)then
				term=(real(j)/real(ncvz(i)))**powrz(i)
			else
				term=1.0-(1.0-real(j)/real(ncvz(i)))**(-powrz(i))
			endif
			zw(j+ist) = statloc + zzone(i)*term
		enddo
		ist = ist + ncvz(i)
		statLoc = statLoc + zzone(i)
	enddo
	nkm1=nk-1

	do i=1,nkm1
	z(i)=(zw(i+1)+zw(i))*0.5
	enddo
	z(nk)=zw(nk)

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
