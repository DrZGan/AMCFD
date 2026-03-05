!______________________________________________________________________________
!
module initialization
!______________________________________________________________________________
! Thin wrapper: re-exports field_data, coeff_data, sim_state for backward compat.
! Contains allocate_fields, deallocate_fields, and initialize.
!
	use geometry
	use constant
	use parameters
	use cfd_utils
	use field_data
	use coeff_data
	use sim_state

	implicit none

	contains

!********************************************************************
subroutine allocate_fields(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk
	call allocate_field_data(nni, nnj, nnk)
	call allocate_coeff_data(nni, nnj, nnk)
end subroutine allocate_fields

!********************************************************************
subroutine deallocate_fields()
	call deallocate_field_data()
	call deallocate_coeff_data()
end subroutine deallocate_fields

!********************************************************************
subroutine initialize

	integer i,j,k
!********************************************************************

	dgdt=dgdtp                    !temperature coefficient of surface tension
	deltemp = tliquid - tsolid    !difference between solidus and liquidus
	cpavg = (acpa*tsolid+acpb+acpl)*0.5        !average heat capacity
	hlcal = hsmelt+cpavg*deltemp  !calcu enthalpy at liquidus
	hlatnt = hlfriz - hlcal       !calcu latent heat
	boufac = denl*g*beta         !bouyancy factor=density/temperature rise

	enthalpyPreheat = temp_to_enthalpy(tempPreheat, acpa, acpb, acpl, tsolid, tliquid, hsmelt, hlcal, deltemp)
	enthalpyWest = temp_to_enthalpy(tempWest, acpa, acpb, acpl, tsolid, tliquid, hsmelt, hlcal, deltemp)
	enthalpyEast = temp_to_enthalpy(tempEast, acpa, acpb, acpl, tsolid, tliquid, hsmelt, hlcal, deltemp)
	enthalpyNorth = temp_to_enthalpy(tempNorth, acpa, acpb, acpl, tsolid, tliquid, hsmelt, hlcal, deltemp)
	enthalpyBottom = temp_to_enthalpy(tempBottom, acpa, acpb, acpl, tsolid, tliquid, hsmelt, hlcal, deltemp)

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
	RHF=1.0_wp   ! RHF disabled: use 1.0 so source terms use sourcedepth_rhf = sourcedepth, etc.

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
