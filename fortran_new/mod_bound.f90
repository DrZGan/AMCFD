!______________________________________________________________________________
!
module boundary
!______________________________________________________________________________
!
	use initialization
	use parameters
	use laserinput
	use dimensions
	use cfd_utils

	implicit none

	contains

!********************************************************************
! Generic u/v boundary condition (replaces bound_u, bound_v)
! idir: 1=u, 2=v
!********************************************************************
subroutine bound_uv(idir)
	integer, intent(in) :: idir
	integer i,j
	real(wp) dtd,fracl_stag,vis1,term1

!-----k=nk (Marangoni stress on top surface)
	do j=jstat,jend
	do i=istatp1,iendm1
		select case(idir)
		case(1)
			dtd=(temp(i,j,nk)-temp(i-1,j,nk))*dxpwinv(i)
			fracl_stag=fracl(i,j,nk)*(1.0-fracx(i-1))+fracl(i-1,j,nk)*fracx(i-1)
			vis1 = harmonic_mean(vis(i,j,nkm1), vis(i-1,j,nkm1), 1.0-fracx(i-1))
			term1=fracl_stag*dgdt*dtd/(vis1*dzpbinv(nk))
			uVel(i,j,nk)=uVel(i,j,nkm1)+term1
		case(2)
			dtd=(temp(i,j,nk)-temp(i,j-1,nk))*dypsinv(j)
			fracl_stag=fracl(i,j,nk)*(1.0-fracy(j-1))+fracl(i,j-1,nk)*fracy(j-1)
			vis1 = harmonic_mean(vis(i,j,nkm1), vis(i,j-1,nkm1), 1.0-fracy(j-1))
			term1=fracl_stag*dgdt*dtd/(vis1*dzpbinv(nk))
			vVel(i,j,nk)=vVel(i,j,nkm1)+term1
		end select
	enddo
	enddo

!----- in solid (boundary conditions)
	select case(idir)
	case(1)
		uVel(istat,jstat:jend,kstat:nkm1)=0.0
		uVel(iend,jstat:jend,kstat:nkm1)=0.0
	case(2)
		vVel(istat,jstat:jend,kstat:nkm1)=0.0
		vVel(iend,jstat:jend,kstat:nkm1)=0.0
	end select
end subroutine bound_uv

!********************************************************************
subroutine bound_w
!-----in solid
	wVel(istat,jstat:jend,kstat:nkm1)=0.0
	wVel(iend,jstat:jend,kstat:nkm1)=0.0
end subroutine bound_w

!********************************************************************
subroutine bound_pp
!----- pp velocities in solid
	pp(istat,jstat:jend,kstat:nkm1)=0.0
	pp(iend,jstat:jend,kstat:nkm1)=0.0
end subroutine bound_pp

!********************************************************************
subroutine bound_enthalpy
	integer i,j,k
	real(wp) ctmp1,hlossradia,hlossconvec

!-----k=nk top surface
	ahtoploss=0.0
	do j=2,njm1
	do i=2,nim1
		hlossradia=emiss*sigm*(temp(i,j,nk)**4-tempAmb**4)
		hlossconvec=htckn*(temp(i,j,nk)-tempAmb)
		ctmp1=diff(i,j,nkm1)*dzpbinv(nk)
		enthalpy(i,j,nk)=enthalpy(i,j,nkm1)+(heatin(i,j)-hlossradia-hlossradia)/ctmp1
		ahtoploss=ahtoploss+(hlossradia+hlossradia)*areaij(i,j)
	enddo
	enddo

!-----k=1 bottom surface
	do j=2,njm1
	do i=2,nim1
		hlossconvec=htck1*(temp(i,j,1)-tempAmb)+emiss*sigm*(temp(i,j,1)**4-tempAmb**4)
		ctmp1=diff(i,j,2)*dzpbinv(2)
		enthalpy(i,j,1)=enthalpy(i,j,2)-hlossconvec/ctmp1
	enddo
	enddo

!-----west and east
	do j=2,njm1
	do k=2,nkm1
		hlossconvec=htci*(temp(1,j,k)-tempAmb)
		ctmp1=diff(2,j,k)*dxpwinv(2)
		enthalpy(1,j,k)=enthalpy(2,j,k)-hlossconvec/ctmp1

		hlossconvec=htci*(temp(ni,j,k)-tempAmb)
		ctmp1=diff(nim1,j,k)*dxpwinv(nim1)
		enthalpy(ni,j,k)=enthalpy(nim1,j,k)-hlossconvec/ctmp1
	enddo
	enddo

!-----north and south
	do i=2,nim1
	do k=2,nkm1
		hlossconvec=htcj*(temp(i,1,k)-tempAmb)
		ctmp1=diff(i,2,k)*dypsinv(2)
		enthalpy(i,1,k)=enthalpy(i,2,k)-hlossconvec/ctmp1

		hlossconvec=htcj*(temp(i,nj,k)-tempAmb)
		ctmp1=diff(i,njm1,k)*dypsinv(njm1)
		enthalpy(i,nj,k)=enthalpy(i,njm1,k)-hlossconvec/ctmp1
	enddo
	enddo
end subroutine bound_enthalpy

end module boundary
