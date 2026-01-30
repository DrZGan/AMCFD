!______________________________________________________________________________
!
module boundary
!______________________________________________________________________________
!
	use initialization
	use parameters
	use laserinput
	use dimensions

	implicit none

	contains 

subroutine bound_condition
	
!********************************************************************
	
	integer i,j,k
	real dtdx,dtdy,term1,acptp,ctmp1,dtemp,htcbtm,tmp
	real fraclu,fraclv
	real hlossradia,hlossconvec		       
	real visu1,visv1



	
!********************************************************************
	goto (100,200,300,400,500)ivar

!********************************************************************
100	continue

!-----k =nk                                                  !Marangoni stress on the top suface
	do j=jstat,jend
	do i=istatp1,iendm1
		
		dtdx=(temp(i,j,nk)-temp(i-1,j,nk))*dxpwinv(i)

		fraclu=fracl(i,j,nk)*(1.0-fracx(i-1))+fracl(i-1,j,nk)*fracx(i-1)
		visu1=vis(i,j,nkm1)*vis(i-1,j,nkm1)/(vis(i-1,j,nkm1)*(1.0-fracx(i-1))+vis(i,j,nkm1)*fracx(i-1))
		term1=fraclu*dgdt*dtdx/(visu1*dzpbinv(nk))

		uVel(i,j,nk)=uVel(i,j,nkm1)+term1

	enddo
	enddo

	
!----- in solid
	uVel(istat,jstat:jend,kstat:nkm1)=0.0                 !at left boundary of i direction, u=0
	uVel(iend,jstat:jend,kstat:nkm1)=0.0 		      !at right boundary of i direction, u=0


	
	return

!********************************************************************
200	continue
!-----k=nk 
	do j=jstat,jend
	do i=istatp1,iendm1
		
		dtdy=(temp(i,j,nk)-temp(i,j-1,nk))*dypsinv(j)
	
		fraclv=fracl(i,j,nk)*(1.0-fracy(j-1))+fracl(i,j-1,nk)*fracy(j-1)
		visv1=vis(i,j,nkm1)*vis(i,j-1,nkm1)/(vis(i,j-1,nkm1)*(1.0-fracy(j-1))+vis(i,j,nkm1)*fracy(j-1))
		term1=fraclv*dgdt*dtdy/(visv1*dzpbinv(nk))

		vVel(i,j,nk)=vVel(i,j,nkm1)+term1


	enddo
	enddo

!-----in solid
	vVel(istat,jstat:jend,kstat:nkm1)=0.0     !at boundary of i direction, v=0
	vVel(iend,jstat:jend,kstat:nkm1)=0.0
	return

!********************************************************************
300	continue
!-----in solid
	wVel(istat,jstat:jend,kstat:nkm1)=0.0     !at boundary of i direction, w=0
	wVel(iend,jstat:jend,kstat:nkm1)=0.0


	return

!********************************************************************
400	continue
!----- pp velocities in solid
	pp(istat,jstat:jend,kstat:nkm1)=0.0     !at boundary of i direction, p=0
	pp(iend,jstat:jend,kstat:nkm1)=0.0
	return

!********************************************************************
500	continue

!-----k=nk top surface
	
	ahtoploss=0.0

	do j=2,njm1   !

	do i=2,nim1
		hlossradia=emiss*sigm*(temp(i,j,nk)**4-tempAmb**4)     !heat loss of radiation
		hlossconvec=htckn*(temp(i,j,nk)-tempAmb)               !heat loss of convection
		ctmp1=diff(i,j,nkm1)*dzpbinv(nk)	
		enthalpy(i,j,nk)=enthalpy(i,j,nkm1)+(heatin(i,j)-hlossradia-hlossradia)/ctmp1	
		!discreted enthalpy boudary condition  
		!ctmp1=k(nkim1)/cp(nkim1)/(z(nk-z(nkm1))

		ahtoploss=ahtoploss+(hlossradia+hlossradia)*areaij(i,j) !total heat loss at top surface
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
!-----    west and east 
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
!-----    north and south 
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
	return
!********************************************************************

end subroutine bound_condition

end module boundary
