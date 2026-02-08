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

	real resoru, resorv, resorw, resorh

	contains

!--------------------------------------------------------------------------
! Generic momentum residual calculation (reusable for u, v, w)
!--------------------------------------------------------------------------
subroutine calc_momentum_residual(vel, sumd)
	real, intent(in) :: vel(nx,ny,nz)
	real, intent(out) :: sumd
	integer :: i, j, k
	real :: resor

	sumd = 0.0
	do k = kstat, nkm1
!$OMP PARALLEL PRIVATE(resor)
!$OMP DO REDUCTION(+: sumd)
	do j = jstat, jend
	do i = istatp1, iendm1
		resor = an(i,j,k) * vel(i,j+1,k) + as(i,j,k) * vel(i,j-1,k) &
			+ ae(i,j,k) * vel(i+1,j,k) + aw(i,j,k) * vel(i-1,j,k) &
			+ at(i,j,k) * vel(i,j,k+1) + ab(i,j,k) * vel(i,j,k-1) &
			+ su(i,j,k) - ap(i,j,k) * vel(i,j,k)
		sumd = sumd + abs(resor)
	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine calc_momentum_residual

!--------------------------------------------------------------------------
! Main residual dispatch
!--------------------------------------------------------------------------
subroutine residual
	integer i, j, k
	real sumd, resor, abs, umaxt, denom, dtpvar, sumh

	select case (ivar)

!============================================================================
	case (1)   ! u-residual
	call calc_momentum_residual(uVel, sumd)
	umaxt = maxval(abs(uVel(istatp1:iendm1, jstat:jend, nk)))
	refmom = 0.25 * pi * MIN(width, alen, depth)**2 * denl * umaxt**2
	resoru = sumd / refmom
	return

!============================================================================
	case (2)   ! v-residual
	call calc_momentum_residual(vVel, sumd)
	resorv = sumd / refmom
	return

!============================================================================
	case (3)   ! w-residual
	call calc_momentum_residual(wVel, sumd)
	resorw = sumd / refmom
	return

!============================================================================
	case (4)   ! pressure-correction (normalized mass source)
	denom = 0.0
	do k = kstat, nkm1
!$OMP PARALLEL PRIVATE(dtpvar)
!$OMP DO REDUCTION(+: denom)
	do j = jstat, jend
	do i = istatp1, iendm1
		dtpvar = (abs(uVel(i,j,k)) + abs(uVel(i+1,j,k))) * areajk(j,k) &
			+ (abs(vVel(i,j,k)) + abs(vVel(i,j+1,k))) * areaik(i,k) &
			+ (abs(wVel(i,j,k)) + abs(wVel(i,j,k+1))) * areaij(i,j)
		denom = denom + 0.5 * abs(dtpvar)
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	denom = denom * denl
	resorm = resorm / (denom + small)
	return

!============================================================================
	case (5)   ! enthalpy residual
	sumh = 0.0
	sumd = 0.0
	do k = 2, nkm1
!$OMP PARALLEL PRIVATE(resor)
!$OMP DO REDUCTION(+: sumd, sumh)
	do j = 2, njm1
	do i = 2, nim1
		resor = (an(i,j,k) * enthalpy(i,j+1,k) + as(i,j,k) * enthalpy(i,j-1,k) &
			+ ae(i,j,k) * enthalpy(i+1,j,k) + aw(i,j,k) * enthalpy(i-1,j,k) &
			+ at(i,j,k) * enthalpy(i,j,k+1) + ab(i,j,k) * enthalpy(i,j,k-1) &
			+ su(i,j,k)) / ap(i,j,k) - enthalpy(i,j,k)
		sumd = sumd + abs(resor)
		sumh = sumh + abs(enthalpy(i,j,k))
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	resorh = sumd / (sumh + small)
	return

	end select

end subroutine residual
end module residue
