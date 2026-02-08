!______________________________________________________________________________
!
module revision
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions

	contains

	subroutine revision_p
	implicit none
	integer i,j,k
	real tulc,tvlc,twlc

	goto (500,500,500,400,500)ivar

!********************************************************************
400	continue
	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(tulc, tvlc, twlc)
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1

		tulc=min(temp(i,j,k),temp(i-1,j,k))
		if(tulc.gt.tsolid)	uVel(i,j,k)=uVel(i,j,k)+dux(i,j,k)*(pp(i-1,j,k)-pp(i,j,k))

		tvlc=min(temp(i,j,k),temp(i,j-1,k))
		if(tvlc.gt.tsolid)	vVel(i,j,k)=vVel(i,j,k)+dvy(i,j,k)*(pp(i,j-1,k)-pp(i,j,k))

		twlc=min(temp(i,j,k),temp(i,j,k-1))
		if(twlc.gt.tsolid)	wVel(i,j,k)=wVel(i,j,k)+dwz(i,j,k)*(pp(i,j,k-1)-pp(i,j,k))
	enddo
	enddo
!$OMP END PARALLEL
	enddo

!-------------------------
	do k=kstat,nkm1
!$OMP PARALLEL 
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		if(temp(i,j,k).gt.tsolid) then
			pressure(i,j,k)=pressure(i,j,k)+urfp*pp(i,j,k)
			pp(i,j,k)=0.0
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	return

!********************************************************************
500	return
	
end subroutine revision_p
end module revision
