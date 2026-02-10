!______________________________________________________________________________
!
module source
!______________________________________________________________________________
!
	use geometry
	use initialization
	use parameters
	use laserinput
	use dimensions
	use boundary
	use cfd_utils
	implicit none

	real(wp), allocatable :: sourceinput(:,:,:)

	contains

!********************************************************************
subroutine allocate_source(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk
	if (allocated(sourceinput)) deallocate(sourceinput)
	allocate(sourceinput(nni, nnj, nnk))
end subroutine allocate_source

!*********************************************************************
subroutine zero_solid_coefficients(i,j,k)
	integer, intent(in) :: i,j,k

	su(i,j,k)=0.0
	an(i,j,k)=0.0
	as(i,j,k)=0.0
	ae(i,j,k)=0.0
	aw(i,j,k)=0.0
	at(i,j,k)=0.0
	ab(i,j,k)=0.0
	ap(i,j,k)=great
end subroutine zero_solid_coefficients

!*********************************************************************
! Generic momentum source term (replaces source_u, source_v, source_w)
! idir: 1=u, 2=v, 3=w
!*********************************************************************
subroutine source_momentum(idir)
	integer, intent(in) :: idir
	integer i,j,k
	real(wp) fracl_stag,term,tw,tlc
	real(wp), parameter :: perm_const = 1.0e-10_wp, eps_darcy = 1.0e-3_wp
	real(wp) darcy_c0

	darcy_c0 = 180.0_wp * viscos / perm_const
!-----Darcy resistance in mushy zone-----
	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(fracl_stag, term, tw)
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		select case(idir)
		case(1); fracl_stag=fracl(i,j,k)*(1.0_wp-fracx(i-1))+fracl(i-1,j,k)*fracx(i-1)
		case(2); fracl_stag=fracl(i,j,k)*(1.0_wp-fracy(j-1))+fracl(i,j-1,k)*fracy(j-1)
		case(3); fracl_stag=fracl(i,j,k)*(1.0_wp-fracz(k-1))+fracl(i,j,k-1)*fracz(k-1)
		end select
		if(fracl_stag.gt.0.0_wp) then
!------mushy zone (Darcy resistance, inlined)--------
			term = darcy_c0 * (1.0_wp - fracl_stag)**2 / (fracl_stag + eps_darcy)
			select case(idir)
			case(1); sp(i,j,k)=sp(i,j,k)-term*volume_u(i,j,k)
			case(2); sp(i,j,k)=sp(i,j,k)-term*volume_v(i,j,k)
			case(3)
				sp(i,j,k)=sp(i,j,k)-term*volume_w(i,j,k)
!-----buoyancy (w only)----------------
				tw=temp(i,j,k)*(1.0-fracz(k-1))+temp(i,j,k-1)*fracz(k-1)
				su(i,j,k) = su(i,j,k)+boufac*volume_w(i,j,k)*(tw-tsolid)
			end select
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo

!-----top boundary treatment (u and v only)-----
	if(idir.eq.1 .or. idir.eq.2) then
		do j=jstat,jend
		do i=istatp1,iendm1
			select case(idir)
			case(1); su(i,j,nkm1)=su(i,j,nkm1)+at(i,j,nkm1)*uVel(i,j,nk)
			case(2); su(i,j,nkm1)=su(i,j,nkm1)+at(i,j,nkm1)*vVel(i,j,nk)
			end select
			sp(i,j,nkm1)=sp(i,j,nkm1)-at(i,j,nkm1)
			at(i,j,nkm1)=0.0
		enddo
		enddo
	endif

!-----assembly and under-relaxation-------------------------------------
	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(tlc)
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		ap(i,j,k)=an(i,j,k)+as(i,j,k)+ae(i,j,k)+aw(i,j,k)+at(i,j,k)+ab(i,j,k)+apnot(i,j,k)-sp(i,j,k)
		select case(idir)
		case(1)
			dux(i,j,k)=areajk(j,k)/ap(i,j,k)
			ap(i,j,k)=ap(i,j,k)/urfu
			su(i,j,k)=su(i,j,k)+(1.-urfu)*ap(i,j,k)*uVel(i,j,k)
			dux(i,j,k)=dux(i,j,k)*urfu
			tlc=min(temp(i,j,k),temp(i-1,j,k))
		case(2)
			dvy(i,j,k)=areaik(i,k)/ap(i,j,k)
			ap(i,j,k)=ap(i,j,k)/urfv
			su(i,j,k)=su(i,j,k)+(1.-urfv)*ap(i,j,k)*vVel(i,j,k)
			dvy(i,j,k)=dvy(i,j,k)*urfv
			tlc=min(temp(i,j,k),temp(i,j-1,k))
		case(3)
			dwz(i,j,k)=areaij(i,j)/ap(i,j,k)
			ap(i,j,k)=ap(i,j,k)/urfw
			su(i,j,k)=su(i,j,k)+(1.-urfw)*ap(i,j,k)*wVel(i,j,k)
			dwz(i,j,k)=dwz(i,j,k)*urfw
			tlc=min(temp(i,j,k),temp(i,j,k-1))
		end select
!------zero velocity in solid (inlined to avoid call overhead in hot loop)-------
		if(tlc.le.tsolid) then
			su(i,j,k)=0.0_wp
			an(i,j,k)=0.0_wp
			as(i,j,k)=0.0_wp
			ae(i,j,k)=0.0_wp
			aw(i,j,k)=0.0_wp
			at(i,j,k)=0.0_wp
			ab(i,j,k)=0.0_wp
			ap(i,j,k)=great
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine source_momentum

!********************************************************************
subroutine source_pp
	integer i,j,k

	do k=kstat,nkm1
!$OMP PARALLEL
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		ap(i,j,k)=an(i,j,k)+as(i,j,k)+ae(i,j,k)+aw(i,j,k)+at(i,j,k)+ab(i,j,k)-sp(i,j,k)
		if(temp(i,j,k).le.tsolid)then
			su(i,j,k)=0.0_wp
			an(i,j,k)=0.0_wp
			as(i,j,k)=0.0_wp
			ae(i,j,k)=0.0_wp
			aw(i,j,k)=0.0_wp
			at(i,j,k)=0.0_wp
			ab(i,j,k)=0.0_wp
			ap(i,j,k)=great
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine source_pp

!********************************************************************
subroutine source_enthalpy(ilo, ihi, jlo, jhi, klo, khi)
	integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
	integer i,j,k
	real(wp) volht,flew,flns,fltb

!-----source term + latent heat (merged into one pass for better cache use)------
	do k=klo,khi
!$OMP PARALLEL PRIVATE(volht, flew, flns, fltb)
!$OMP DO
	do j=jlo,jhi
	do i=ilo,ihi
		if(toolmatrix(PathNum,5) .gt. laser_on_threshold)then
			if(z(nk)-z(k).le.sourcedepth) then
				sourceinput(i,j,k)=alaspowvol*alasfact/pi/sourcerad**2/sourcedepth*alasetavol*&
				exp(-alasfact/sourcerad**2*((beam_pos-x(i))**2+(beam_posy-y(j))**2))
			else
				sourceinput(i,j,k)=0.0_wp
			endif
		else
			sourceinput(i,j,k)=0.0_wp
		endif
		su(i,j,k)=su(i,j,k)+volume(i,j,k)*sourceinput(i,j,k)
		volht=volume(i,j,k)*hlatnt*den(i,j,k)/delt
		su(i,j,k)=su(i,j,k)-volht*(fracl(i,j,k)-fraclnot(i,j,k))
		flew=areajk(j,k)*(max(uVel(i,j,k),0.0_wp)*fracl(i-1,j,k)-max(-uVel(i,j,k),0.0_wp)*fracl(i,j,k) &
			+max(-uVel(i+1,j,k),0.0_wp)*fracl(i+1,j,k)-max(uVel(i+1,j,k),0.0_wp)*fracl(i,j,k))
		flns=areaik(i,k)*(max(vVel(i,j,k),0.0_wp)*fracl(i,j-1,k)-max(-vVel(i,j,k),0.0_wp)*fracl(i,j,k)  &
			+max(-vVel(i,j+1,k),0.0_wp)*fracl(i,j+1,k)-max(vVel(i,j+1,k),0.0_wp)*fracl(i,j,k))
		fltb=areaij(i,j)*(max(wVel(i,j,k),0.0_wp)*fracl(i,j,k-1)-max(-wVel(i,j,k),0.0_wp)*fracl(i,j,k) &
			+max(-wVel(i,j,k+1),0.0_wp)*fracl(i,j,k+1)-max(wVel(i,j,k+1),0.0_wp)*fracl(i,j,k))
		su(i,j,k)=su(i,j,k)+den(i,j,k)*hlatnt*(flew+flns+fltb)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

!----- boundary condition transfers ------
!----- k=nk & k=1 ------
	do j=jlo,jhi
	do i=ilo,ihi
		su(i,j,klo)=su(i,j,klo)+ab(i,j,klo)*enthalpy(i,j,klo-1)
		sp(i,j,klo)=sp(i,j,klo)-ab(i,j,klo)
		ab(i,j,klo)=0.0
		su(i,j,khi)=su(i,j,khi)+at(i,j,khi)*enthalpy(i,j,khi+1)
		sp(i,j,khi)=sp(i,j,khi)-at(i,j,khi)
		at(i,j,khi)=0.0
	enddo
	enddo

!----- j=1 & j=nj ------
	do k=klo,khi
	do i=ilo,ihi
		su(i,jlo,k)=su(i,jlo,k)+as(i,jlo,k)*enthalpy(i,jlo-1,k)
		sp(i,jlo,k)=sp(i,jlo,k)-as(i,jlo,k)
		as(i,jlo,k)=0.0
		su(i,jhi,k)=su(i,jhi,k)+an(i,jhi,k)*enthalpy(i,jhi+1,k)
		sp(i,jhi,k)=sp(i,jhi,k)-an(i,jhi,k)
		an(i,jhi,k)=0.0
	enddo
	enddo

!----- i=1 & i=ni---
	do k=klo,khi
	do j=jlo,jhi
		su(ilo,j,k)=su(ilo,j,k)+aw(ilo,j,k)*enthalpy(ilo-1,j,k)
		sp(ilo,j,k)=sp(ilo,j,k)-aw(ilo,j,k)
		aw(ilo,j,k)=0.0
		su(ihi,j,k)=su(ihi,j,k)+ae(ihi,j,k)*enthalpy(ihi+1,j,k)
		sp(ihi,j,k)=sp(ihi,j,k)-ae(ihi,j,k)
		ae(ihi,j,k)=0.0
	enddo
	enddo

!-----assembly and under-relaxation------
	do k=klo,khi
!$OMP PARALLEL
!$OMP DO
	do j=jlo,jhi
	do i=ilo,ihi
		ap(i,j,k)=an(i,j,k)+as(i,j,k)+ae(i,j,k)+aw(i,j,k)+at(i,j,k)+ab(i,j,k)+apnot(i,j,k)-sp(i,j,k)
!-----under-relaxation---
		ap(i,j,k)=ap(i,j,k)/urfh
		su(i,j,k)=su(i,j,k)+(1.-urfh)*ap(i,j,k)*enthalpy(i,j,k)
	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine source_enthalpy

end module source
