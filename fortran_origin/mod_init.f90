!______________________________________________________________________________
!
module initialization
!______________________________________________________________________________
!
	use geometry
	use constant
	use parameters

	implicit none

 
	real dgdt, deltemp,cpavg,hlcal,hlatnt,boufac
 	!temperature coefficient of surface tension,
        !temperature difference between solidus and liquidus,
        !average capacity,
        !calucation enthalpy,
        !latent enthalpy,
        !diffusion of solid,
        !diffusion of liquid,
        !bouyancy/temperature rise,
        !density*scanning speed

	real vis(nx,ny,nz),diff(nx,ny,nz),den(nx,ny,nz)
	!viscosity matrix, diffusion matrix, density matrix

	real uVel(nx,ny,nz),vVel(nx,ny,nz),wVel(nx,ny,nz),unot(nx,ny,nz),vnot(nx,ny,nz),wnot(nx,ny,nz)
	!uvw velocity matrix, previous time velocity matrix
		
	real pressure(nx,ny,nz),pp(nx,ny,nz),enthalpy(nx,ny,nz),hnot(nx,ny,nz),temp(nx,ny,nz),tnot(nx,ny,nz)
	!pressure, enthalpy, temperature, and previous pressure, enthalpy, temperature.

	real dux(nx,ny,nz),dvy(nx,ny,nz),dwz(nx,ny,nz),su(nx,ny,nz),sp(nx,ny,nz)
	!velocity rise matrix, source term matrix su and sp

	integer ivar	!main   index of variable equation

	real phi(nx,ny,nz,nvar)  ! four diemnsion matrix, four varialbe at each space point: u,v,w,p
	equivalence (phi(1,1,1,1),uVel(1,1,1)),(phi(1,1,1,2),vVel(1,1,1)),(phi(1,1,1,3),wVel(1,1,1)),	&
		(phi(1,1,1,4),pp(1,1,1))

	real fracl(nx,ny,nz),fraclnot(nx,ny,nz)	
	!volume fraction of liquid matrix, and previous matrix
	

	real resorm,refmom,ahtoploss	
	!  pressure residual error, reference residual error and total heat loss at top surface
	
	real ap(nx,ny,nz),an(nx,ny,nz),as(nx,ny,nz),ae(nx,ny,nz),aw(nx,ny,nz),at(nx,ny,nz),ab(nx,ny,nz), &
		apnot(nx,ny,nz)	

	real enthalpyWest,enthalpyEast,enthalpyNorth,enthalpyBottom,enthalpyPreheat

	integer TrackNum, PathNum
	real solidfield(nx,ny,nz)
	real beam_pos, beam_posy
	real toolmatrix(TOOLLINES,5)
	real coordhistory(COORDLINES,8)
	real RHF
	
	contains

subroutine initialize
		
	integer i,j,k
!********************************************************************
	
	dgdt=dgdtp                    !temperature coefficient of surface tension
	deltemp = tliquid - tsolid    !difference between solidus and liquidus
	cpavg = (acpa*tsolid+acpb+acpl)*0.5        !average heat capacity
	hlcal = hsmelt+cpavg*deltemp  !calcu enthalpy at liquidus
	hlatnt = hlfriz - hlcal       !calcu latent heat
	!write(*,*) hlatnt
	boufac = denl*g*beta         !bouyancy factor=density/temperature rise

	enthalpyPreheat = 0.5*acpa*tempPreheat**2+acpb*tempPreheat!(tempPreheat-tsolid)*(acpa*tempPreheat/2+acpb)+hsmelt  !translate preheat temp to enthalpy at the boundary
	enthalpyWest =  0.5*acpa*tempWest**2+acpb*tempWest!(tempWest-tsolid)*(acpa*tempWest/2+acpb)+hsmelt
	enthalpyEast =  0.5*acpa*tempEast**2+acpb*tempEast!(tempEast-tsolid)*(acpa*tempEast/2+acpb)+hsmelt
	enthalpyNorth =  0.5*acpa*tempNorth**2+acpb*tempNorth!(tempNorth-tsolid)*(acpa*tempNorth/2+acpb)+hsmelt
	enthalpyBottom =  0.5*acpa*tempBottom**2+acpb*tempBottom!(tempBottom-tsolid)*(acpa*tempBottom/2+acpb)+hsmelt

	do k=1,nk
	do j=1,nj
	do i=1,ni
		vis(i,j,k)= viscos         !inital vis to be viscos
		den(i,j,k)=denl
		uVel(i,j,k)=0.0
		unot(i,j,k)=0.0
		vVel(i,j,k)=0.0
		vnot(i,j,k)=0.0
		wVel(i,j,k)=0.0
		wnot(i,j,k)=0.0
		pressure(i,j,k)=0.0
		pp(i,j,k)=0.0
		enthalpy(i,j,k)=enthalpyPreheat !inital enthalpy to be preheat enthalpy
		hnot(i,j,k)=enthalpyPreheat
		temp(i,j,k)=tempPreheat
		tnot(i,j,k)=tempPreheat
		dux(i,j,k)=0.0
		dvy(i,j,k)=0.0
		dwz(i,j,k)=0.0
		su(i,j,k)=0.0
		sp(i,j,k)=0.0

		fracl(i,j,k)=0.0               !!inital volume fraction of liquid
		fraclnot(i,j,k)=0.0            !!inital previous vol fraction of liq

		diff(i,j,k)=(thconsa*tempPreheat+thconsb)/(acpa*tempPreheat/2+acpb)          !diffusion coefficient= thermal conductivity * capacity

		solidfield(i,j,k)=0

	enddo
	enddo
	enddo	
  
  	coordhistory(1:COORDLINES,1:8)=-1
	TrackNum=0
	PathNum=2
	
	beam_pos=toolmatrix(PathNum,2)
	beam_posy=toolmatrix(PathNum,3)

!------------enthalpy BC: i=1 plane------------inital enthalpy at boundary 
	do j=1,nj
	do k=1,nk
		enthalpy(1,j,k)=enthalpyWest
	enddo
	enddo

!-----i=ni plane------------------
	do j=1,nj
	do k=1,nk
		enthalpy(ni,j,k)=enthalpyEast
	enddo
	enddo

!-----k=1 plane--------------
	do i=1,ni
	do j=1,nj
		enthalpy(i,j,1)=enthalpyBottom
	enddo
	enddo

!-----j=nj plane---------------
	do i=1,ni
	do k=1,nk
		enthalpy(i,nj,k)=enthalpyNorth
	enddo
	enddo


	return

end subroutine initialize

end module initialization

