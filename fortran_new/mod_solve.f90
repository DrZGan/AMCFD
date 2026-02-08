!______________________________________________________________________________
!
module solver
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions

	contains

!--------------------------------------------------------------------------
! Solve momentum equations (u, v, w, or pp) using line-by-line TDMA
!--------------------------------------------------------------------------
subroutine solution_uvw

	implicit none
	integer i, j, k, ksweep, jsweep
	real pr, qr, d, denom
	dimension pr(nx), qr(nx)

	do ksweep = 1, 2
	do k = nkm1, kstat, -1
	do jsweep = 1, 2
!$OMP PARALLEL PRIVATE(pr, qr, d, denom)
!$OMP DO
	do j = jstat, jend
		i = istat
		pr(i) = 0.0
		qr(i) = phi(i,j,k,ivar)

		do i = istatp1, iendm1
			d = at(i,j,k) * phi(i,j,k+1,ivar) + ab(i,j,k) * phi(i,j,k-1,ivar) &
				+ an(i,j,k) * phi(i,j+1,k,ivar) + as(i,j,k) * phi(i,j-1,k,ivar) + su(i,j,k)
			denom = ap(i,j,k) - aw(i,j,k) * pr(i-1)

			if (denom .le. TDMA_ZERO_GUARD .and. denom .ge. 0) denom = denom + TDMA_ZERO_OFFSET
			if (denom .ge. -TDMA_ZERO_GUARD .and. denom .lt. 0) denom = denom - TDMA_ZERO_OFFSET

			pr(i) = ae(i,j,k) / (denom)
			qr(i) = (d + aw(i,j,k) * qr(i-1)) / (denom)
		enddo
!----- back substitution -----
		do i = iendm1, istatp1, -1
			phi(i,j,k,ivar) = pr(i) * phi(i+1,j,k,ivar) + qr(i)
		enddo
	enddo
!$OMP END PARALLEL
	enddo
	enddo
	enddo

	return
end subroutine solution_uvw

!--------------------------------------------------------------------------
! Solve enthalpy equation using line-by-line TDMA
!--------------------------------------------------------------------------
subroutine solution_enthalpy

	implicit none
	integer i, j, k, ksweep, jsweep
	real pr, qr, d, denom
	dimension pr(nx), qr(nx)

	! ivar=5 always reaches here; no goto dispatch needed
	do ksweep = 1, 2
	do k = nkm1, 2, -1
	do jsweep = 1, 2
!$OMP PARALLEL PRIVATE(pr, qr, d, denom)
!$OMP DO
	do j = 2, njm1
		pr(1) = 0.0
		qr(1) = enthalpy(1, j, k)

		do i = 2, nim1
			d = at(i,j,k) * enthalpy(i,j,k+1) + ab(i,j,k) * enthalpy(i,j,k-1) &
				+ an(i,j,k) * enthalpy(i,j+1,k) + as(i,j,k) * enthalpy(i,j-1,k) + su(i,j,k)
			denom = ap(i,j,k) - aw(i,j,k) * pr(i-1)

			if (denom .le. TDMA_ZERO_GUARD .and. denom .ge. 0) denom = denom + TDMA_ZERO_OFFSET
			if (denom .ge. -TDMA_ZERO_GUARD .and. denom .lt. 0) denom = denom - TDMA_ZERO_OFFSET

			pr(i) = ae(i,j,k) / (denom)
			qr(i) = (d + aw(i,j,k) * qr(i-1)) / (denom)
		enddo

!----- back substitution -----
		do i = nim1, 2, -1
			enthalpy(i,j,k) = pr(i) * enthalpy(i+1,j,k) + qr(i)
		enddo
	enddo
!$OMP END PARALLEL
	enddo
	enddo
	enddo

	return

end subroutine solution_enthalpy


subroutine cleanuvw

	implicit none
	integer i, j, k
	real tulc, tvlc, twlc

	do k = kstat, nkm1
!$OMP PARALLEL PRIVATE(tulc, tvlc, twlc)
!$OMP DO
	do j = jstat, jend
	do i = istatp1, iendm1
		tulc = min(temp(i,j,k), temp(i+1,j,k))
		tvlc = min(temp(i,j,k), temp(i,j+1,k))
		twlc = min(temp(i,j,k), temp(i,j,k+1))
		if (tulc .le. tsolid) uVel(i+1,j,k) = 0.0
		if (tvlc .le. tsolid) vVel(i,j+1,k) = 0.0
		if (twlc .le. tsolid) wVel(i,j,k+1) = 0.0

		if (temp(i,j,nk) .ge. tboiling) then
			uVel(i,j,k) = 0.0
			vVel(i,j,k) = 0.0
			wVel(i,j,k) = 0.0
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	return

end subroutine cleanuvw

end module
