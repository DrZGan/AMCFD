!______________________________________________________________________________
!
module convergence
!______________________________________________________________________________
!
	use initialization
	use dimensions

	contains

subroutine enhance_converge_speed   
	implicit none
	integer i,j,k
	real bl,blp,blm,blc       
	real pib,qib,denom,delh      !Defination method from F77,  denom is real, others are arraies with the length of nx
		
	dimension bl(nx), blp(nx), blm(nx), blc(nx), delh(nx), pib(nx),qib(nx)




	bl(1:ni)=0.0
	blp(1:ni)=0.0
	blm(1:ni)=0.0
	blc(1:ni)=0.0

	do k=2,nkm1
!$OMP PARALLEL 
!$OMP DO REDUCTION(+: bl, blp, blm, blc)
	do j=2,njm1
	do i=2,nim1
	
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

	pib(2)=blp(2)/bl(2)
	qib(2)=blc(2)/bl(2)

	do i=3,nim1
		denom=bl(i)-blm(i)*pib(i-1)
		pib(i)=blp(i)/denom
		qib(i)=(blc(i)+blm(i)*qib(i-1))/denom
	enddo

	delh(nim1)=qib(nim1)

	do i=nim1-1,2,-1
		delh(i)=pib(i)*delh(i+1)+qib(i)
	enddo

	do k=2,nkm1
!$OMP PARALLEL 
!$OMP DO
	do j=2,njm1
	do i=2,nim1
		enthalpy(i,j,k)=enthalpy(i,j,k)+delh(i)
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	return

end subroutine enhance_converge_speed
end module convergence
