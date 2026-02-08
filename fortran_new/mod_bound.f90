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

subroutine bound_condition

!********************************************************************

	integer i,j,k
	real dtdx,dtdy,term1,ctmp1
	real fraclu,fraclv
	real hlossradia,hlossconvec
	real visu1,visv1

!********************************************************************
	select case(ivar)
	case(1)
		call bound_u
	case(2)
		call bound_v
	case(3)
		call bound_w
	case(4)
		call bound_pp
	case(5)
		call bound_enthalpy
	end select

	return
end subroutine bound_condition

!********************************************************************
subroutine bound_u
	integer i,j,k
	real dtdx,fraclu,visu1,term1

!-----k=nk (Marangoni stress on top surface)
	do j=jstat,jend
	do i=istatp1,iendm1
		dtdx=(temp(i,j,nk)-temp(i-1,j,nk))*dxpwinv(i)
		fraclu=fracl(i,j,nk)*(1.0-fracx(i-1))+fracl(i-1,j,nk)*fracx(i-1)
		visu1 = harmonic_mean(vis(i,j,nkm1), vis(i-1,j,nkm1), 1.0-fracx(i-1))
		term1=fraclu*dgdt*dtdx/(visu1*dzpbinv(nk))
		uVel(i,j,nk)=uVel(i,j,nkm1)+term1
	enddo
	enddo

!----- in solid (boundary conditions)
	uVel(istat,jstat:jend,kstat:nkm1)=0.0
	uVel(iend,jstat:jend,kstat:nkm1)=0.0
end subroutine bound_u

!********************************************************************
subroutine bound_v
	integer i,j,k
	real dtdy,fraclv,visv1,term1

!-----k=nk (Marangoni stress on top surface)
	do j=jstat,jend
	do i=istatp1,iendm1
		dtdy=(temp(i,j,nk)-temp(i,j-1,nk))*dypsinv(j)
		fraclv=fracl(i,j,nk)*(1.0-fracy(j-1))+fracl(i,j-1,nk)*fracy(j-1)
		visv1 = harmonic_mean(vis(i,j,nkm1), vis(i,j-1,nkm1), 1.0-fracy(j-1))
		term1=fraclv*dgdt*dtdy/(visv1*dzpbinv(nk))
		vVel(i,j,nk)=vVel(i,j,nkm1)+term1
	enddo
	enddo

!-----in solid
	vVel(istat,jstat:jend,kstat:nkm1)=0.0
	vVel(iend,jstat:jend,kstat:nkm1)=0.0
end subroutine bound_v

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
	real ctmp1,hlossradia,hlossconvec

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
