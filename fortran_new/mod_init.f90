!______________________________________________________________________________
!
module initialization
!______________________________________________________________________________
!
	use geometry
	use constant
	use parameters
	use cfd_utils

	implicit none

 
	real(wp) dgdt, deltemp,cpavg,hlcal,hlatnt,boufac
 	!temperature coefficient of surface tension,
        !temperature difference between solidus and liquidus,
        !average capacity,
        !calucation enthalpy,
        !latent enthalpy,
        !diffusion of solid,
        !diffusion of liquid,
        !bouyancy/temperature rise,
        !density*scanning speed

	real(wp) vis(nx,ny,nz),diff(nx,ny,nz),den(nx,ny,nz)
	!viscosity matrix, diffusion matrix, density matrix

	real(wp) uVel(nx,ny,nz),vVel(nx,ny,nz),wVel(nx,ny,nz),unot(nx,ny,nz),vnot(nx,ny,nz),wnot(nx,ny,nz)
	!uvw velocity matrix, previous time velocity matrix
		
	real(wp) pressure(nx,ny,nz),pp(nx,ny,nz),enthalpy(nx,ny,nz),hnot(nx,ny,nz),temp(nx,ny,nz),tnot(nx,ny,nz)
	!pressure, enthalpy, temperature, and previous pressure, enthalpy, temperature.

	real(wp) dux(nx,ny,nz),dvy(nx,ny,nz),dwz(nx,ny,nz),su(nx,ny,nz),sp(nx,ny,nz)
	!velocity rise matrix, source term matrix su and sp

	! ivar removed: specific subroutines are called directly from main.f90
	! phi/equivalence removed: field arrays are passed explicitly to the TDMA solver

	real(wp) fracl(nx,ny,nz),fraclnot(nx,ny,nz)	
	!volume fraction of liquid matrix, and previous matrix
	

	real(wp) resorm,refmom,ahtoploss	
	!  pressure residual error, reference residual error and total heat loss at top surface
	
	real(wp) ap(nx,ny,nz),an(nx,ny,nz),as(nx,ny,nz),ae(nx,ny,nz),aw(nx,ny,nz),at(nx,ny,nz),ab(nx,ny,nz), &
		apnot(nx,ny,nz)	

	real(wp) enthalpyWest,enthalpyEast,enthalpyNorth,enthalpyBottom,enthalpyPreheat

	integer TrackNum, PathNum
	real(wp) solidfield(nx,ny,nz)
	real(wp) beam_pos, beam_posy
	real(wp) toolmatrix(TOOLLINES,5)
	real(wp) coordhistory(COORDLINES,8)
	real(wp) RHF
	
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

	enthalpyPreheat = temp_to_enthalpy(tempPreheat, acpa, acpb)  !translate preheat temp to enthalpy at the boundary
	enthalpyWest = temp_to_enthalpy(tempWest, acpa, acpb)
	enthalpyEast = temp_to_enthalpy(tempEast, acpa, acpb)
	enthalpyNorth = temp_to_enthalpy(tempNorth, acpa, acpb)
	enthalpyBottom = temp_to_enthalpy(tempBottom, acpa, acpb)

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

