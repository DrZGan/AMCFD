!______________________________________________________________________________
!
module residue
!______________________________________________________________________________
!
	use constant
	use geometry
	use initialization
	use dimensions
	use parameters
	
	implicit none

	real resoru,resorv,resorw,resorh  !residual error of u v w h
	contains

subroutine residual
	integer i,j,k
	real sumd,resor,abs,umaxt,denom,dtpvar,sumh 
	go to (100,200,300,400,500) ivar

!********************************************************************
100	continue
	sumd=0.0

	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(resor)
!$OMP DO REDUCTION(+: sumd)
	do j=jstat,jend
	do i=istatp1,iendm1

		resor=an(i,j,k)*uVel(i,j+1,k)+as(i,j,k)*uVel(i,j-1,k)+ae(i,j,k)*uVel(i+1,j,k)+aw(i,j,k)*uVel(i-1,j,k) &
			+at(i,j,k)*uVel(i,j,k+1)+ab(i,j,k)*uVel(i,j,k-1)+su(i,j,k)-ap(i,j,k)*uVel(i,j,k) 
	
	
		sumd=sumd+abs(resor)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	umaxt=maxval(abs(uVel(istatp1:iendm1,jstat:jend,nk)))

!---reference momentum---------
	refmom=0.25*pi*MIN(width,alen,depth)**2*denl*umaxt**2
	
!----normalized residual--------- 
	resoru=sumd/refmom
	return

!********************************************************************
200	continue
	sumd=0.0
	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(resor)
!$OMP DO REDUCTION(+: sumd)
	do j=jstat,jend
	do i=istatp1,iendm1
		resor=an(i,j,k)*vVel(i,j+1,k)+as(i,j,k)*vVel(i,j-1,k)+ae(i,j,k)*vVel(i+1,j,k)+aw(i,j,k)*vVel(i-1,j,k) &
			+at(i,j,k)*vVel(i,j,k+1)+ab(i,j,k)*vVel(i,j,k-1)+su(i,j,k)-ap(i,j,k)*vVel(i,j,k)

	
		sumd=sumd+abs(resor)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	resorv=sumd/refmom
	return

!********************************************************************
300	continue
	sumd=0.0
	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(resor)
!$OMP DO REDUCTION(+: sumd)
	do j=jstat,jend
	do i=istatp1,iendm1
		resor=an(i,j,k)*wVel(i,j+1,k)+as(i,j,k)*wVel(i,j-1,k)+ae(i,j,k)*wVel(i+1,j,k)+aw(i,j,k)*wVel(i-1,j,k) &
			+at(i,j,k)*wVel(i,j,k+1)+ab(i,j,k)*wVel(i,j,k-1)+su(i,j,k)-ap(i,j,k)*wVel(i,j,k)

		sumd=sumd+abs(resor)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	resorw=sumd/refmom
	return

!********************************************************************
400	continue
!----- normalized mass source
	denom=0.0
	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(dtpvar)
!$OMP DO REDUCTION(+: denom)
	do j=jstat,jend
	do i=istatp1,iendm1
		dtpvar=(abs(uVel(i,j,k))+abs(uVel(i+1,j,k)))*areajk(j,k)+(abs(vVel(i,j,k))+abs(vVel(i,j+1,k))) &
				*areaik(i,k)+(abs(wVel(i,j,k))+abs(wVel(i,j,k+1)))*areaij(i,j)
		denom=denom+0.5*abs(dtpvar)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	denom=denom*denl
	resorm=resorm/(denom+small)
	return

!********************************************************************
500	continue

	sumh=0.0
	sumd=0.0

	do k=2,nkm1
!$OMP PARALLEL PRIVATE(resor)
!$OMP DO REDUCTION(+: sumd, sumh)
	do j=2,njm1
	do i=2,nim1
		resor=(an(i,j,k)*enthalpy(i,j+1,k)+as(i,j,k)*enthalpy(i,j-1,k)+ae(i,j,k)*enthalpy(i+1,j,k)+ &
			aw(i,j,k)*enthalpy(i-1,j,k)+at(i,j,k)*enthalpy(i,j,k+1)+ab(i,j,k)*enthalpy(i,j,k-1)+ &
			su(i,j,k))/ap(i,j,k)-enthalpy(i,j,k)
		sumd=sumd+abs(resor)
		sumh=sumh+abs(enthalpy(i,j,k))
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	resorh=sumd/(sumh+small)
	return
!********************************************************************

end subroutine residual
end module residue
