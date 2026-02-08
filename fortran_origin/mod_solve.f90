!______________________________________________________________________________
!
module solver
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions
	
	contains

	subroutine solution_uvw
	
	implicit none
	integer i,j,k,ksweep,jsweep,isweep
	real pr,qr,d,denom

	dimension pr(nx),qr(nx)

!********************************************************************
	do ksweep=1,2
	do k=nkm1,kstat,-1
	do jsweep=1,2
!$OMP PARALLEL PRIVATE(pr, qr, d, denom)
!$OMP DO
	do j=jstat,jend
		i=istat
		pr(i)=0.0
		qr(i)=phi(i,j,k,ivar)

		do i=istatp1,iendm1
			d = at(i,j,k)*phi(i,j,k+1,ivar)+ab(i,j,k)*phi(i,j,k-1,ivar)+an(i,j,k)*phi(i,j+1,k,ivar) &
				+as(i,j,k)*phi(i,j-1,k,ivar)+su(i,j,k)
			denom=ap(i,j,k)-aw(i,j,k)*pr(i-1)

		    	if(denom.le.1e-12 .and.  denom.ge.0) denom=denom+1e-13     !!avoid divide zero
		    	if(denom.ge.-1e-12 .and.  denom.lt.0) denom=denom-1e-13     !!avoid divide zero


			pr(i)=ae(i,j,k)/(denom)
			qr(i)=(d+aw(i,j,k)*qr(i-1))/(denom)
		enddo
!-----back---------------- 
		do i=iendm1, istatp1, -1
			phi(i,j,k,ivar)=pr(i)*phi(i+1,j,k,ivar)+qr(i)
		enddo
	enddo
!$OMP END PARALLEL
	enddo
	enddo
	enddo





	return
end subroutine solution_uvw
!********************************************************************
subroutine solution_enthalpy

	implicit none
	integer i,j,k,ksweep,jsweep,isweep
	real pr,qr,d,denom
	dimension pr(nx),qr(nx)

!********************************************************************
	goto (100,200,300,400,500)ivar

!********************************************************************
100	continue
200	continue
300	continue
400	continue
500	continue
!-----TDMA

	do ksweep=1,2                          !scan k direction twice

	do k=nkm1,2,-1
	do jsweep=1,2                          !scen j direction twice
!$OMP PARALLEL PRIVATE(pr, qr, d, denom)
!$OMP DO
	do j=2,njm1
		pr(1)=0.0
		qr(1)=enthalpy(1,j,k)

		do i=2,nim1
			d = at(i,j,k)*enthalpy(i,j,k+1)+ab(i,j,k)*enthalpy(i,j,k-1)+an(i,j,k)*enthalpy(i,j+1,k)+ &
				as(i,j,k)*enthalpy(i,j-1,k)+su(i,j,k)
			denom=ap(i,j,k)-aw(i,j,k)*pr(i-1)

		    	if(denom.le.1e-12 .and.  denom.ge.0) denom=denom+1e-13     !!avoid divide zero
		    	if(denom.ge.-1e-12 .and.  denom.lt.0) denom=denom-1e-13     !!avoid divide zero


			pr(i)=ae(i,j,k)/(denom)
			qr(i)=(d+aw(i,j,k)*qr(i-1))/(denom)
		enddo

!-----back 
		do i=nim1,2,-1
			enthalpy(i,j,k)=pr(i)*enthalpy(i+1,j,k)+qr(i)
		enddo

	enddo
!$OMP END PARALLEL
	enddo
	enddo

	enddo
	return
!********************************************************************

end subroutine solution_enthalpy


subroutine cleanuvw

	implicit none
	integer i,j,k
	real tulc,tvlc,twlc

	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(tulc, tvlc, twlc)
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		tulc=min(temp(i,j,k),temp(i+1,j,k))
		tvlc=min(temp(i,j,k),temp(i,j+1,k))
		twlc=min(temp(i,j,k),temp(i,j,k+1))
		if(tulc.le.tsolid) uVel(i+1,j,k)=0.0
		if(tvlc.le.tsolid) vVel(i,j+1,k)=0.0
		if(twlc.le.tsolid) wVel(i,j,k+1)=0.0

		if(temp(i,j,nk) .ge. tboiling)then
			uVel(i,j,k)=0.0
			vVel(i,j,k)=0.0
			wVel(i,j,k)=0.0
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	return

end subroutine cleanuvw

end module
