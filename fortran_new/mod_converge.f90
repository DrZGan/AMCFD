!______________________________________________________________________________
!
module convergence
!______________________________________________________________________________
!
	use initialization
	use dimensions

	contains

subroutine enhance_converge_speed(ilo, ihi, jlo, jhi, klo, khi)
	implicit none
	integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
	integer i,j,k
	real(wp), allocatable :: bl(:),blp(:),blm(:),blc(:),delh(:),pib(:),qib(:)
	real(wp) denom

	allocate(bl(ni), blp(ni), blm(ni), blc(ni), delh(ni), pib(ni), qib(ni))




	bl(1:ni)=0.0
	blp(1:ni)=0.0
	blm(1:ni)=0.0
	blc(1:ni)=0.0

	do k=klo,khi
!$OMP PARALLEL 
!$OMP DO REDUCTION(+: bl, blp, blm, blc)
	do j=jlo,jhi
	do i=ilo,ihi
	
	bl(i)=bl(i)+ap(i,j,k)-an(i,j,k)-as(i,j,k)-at(i,j,k)-ab(i,j,k)
		blp(i)=blp(i)+ae(i,j,k)
		blm(i)=blm(i)+aw(i,j,k)
		blc(i)=blc(i)+ae(i,j,k)*enthalpy(i+1,j,k)+aw(i,j,k)*enthalpy(i-1,j,k)+an(i,j,k)*enthalpy(i,j+1,k)+ &
			as(i,j,k)*enthalpy(i,j-1,k)+at(i,j,k)*enthalpy(i,j,k+1)+ab(i,j,k)*enthalpy(i,j,k-1)+su(i,j,k) &
			-ap(i,j,k)*enthalpy(i,j,k)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	pib(ilo)=blp(ilo)/bl(ilo)
	qib(ilo)=blc(ilo)/bl(ilo)

	do i=ilo+1,ihi
		denom=bl(i)-blm(i)*pib(i-1)
		pib(i)=blp(i)/denom
		qib(i)=(blc(i)+blm(i)*qib(i-1))/denom
	enddo

	delh(ihi)=qib(ihi)

	do i=ihi-1,ilo,-1
		delh(i)=pib(i)*delh(i+1)+qib(i)
	enddo

	do k=klo,khi
!$OMP PARALLEL 
!$OMP DO
	do j=jlo,jhi
	do i=ilo,ihi
		enthalpy(i,j,k)=enthalpy(i,j,k)+delh(i)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	deallocate(bl, blp, blm, blc, delh, pib, qib)
	return

end subroutine enhance_converge_speed
end module convergence
