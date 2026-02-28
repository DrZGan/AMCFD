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

	real(wp) resoru,resorv,resorw,resorh  !residual error of u v w h
	contains

!********************************************************************
! Generic momentum residual calculation for u, v, w velocities
!********************************************************************
subroutine calc_momentum_residual(vel, resor_out, calc_refmom)
	real(wp), intent(in) :: vel(:,:,:)
	real(wp), intent(out) :: resor_out
	logical, intent(in) :: calc_refmom
	integer i,j,k
	real(wp) sumd,resor,umaxt

	sumd=0.0

	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(resor)
!$OMP DO REDUCTION(+: sumd)
	do j=jstat,jend
	do i=istatp1,iendm1
		resor=an(i,j,k)*vel(i,j+1,k)+as(i,j,k)*vel(i,j-1,k)+ae(i,j,k)*vel(i+1,j,k)+aw(i,j,k)*vel(i-1,j,k) &
			+at(i,j,k)*vel(i,j,k+1)+ab(i,j,k)*vel(i,j,k-1)+su(i,j,k)-ap(i,j,k)*vel(i,j,k)
		sumd=sumd+abs(resor)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	if(calc_refmom) then
		umaxt=maxval(abs(uVel(istatp1:iendm1,jstat:jend,nk)))
		refmom=0.25*pi*MIN(width,alen,depth)**2*denl*umaxt**2
		refmom=max(refmom, small)
	endif

	resor_out=sumd/refmom
end subroutine calc_momentum_residual

!********************************************************************
subroutine calc_pressure_residual
	integer i,j,k
	real(wp) denom,dtpvar

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
end subroutine calc_pressure_residual

!********************************************************************
subroutine calc_enthalpy_residual
	integer i,j,k
	real(wp) sumd,sumh,resor

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
end subroutine calc_enthalpy_residual

end module residue
