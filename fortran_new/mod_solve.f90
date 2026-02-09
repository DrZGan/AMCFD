!______________________________________________________________________________
!
module solver
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions
	
	contains

!********************************************************************
! Unified line-by-line TDMA solver for 3D fields.
! Solves the system defined by coefficient arrays (an,as,ae,aw,at,ab,ap,su)
! for the given field, using double k-sweep and j-sweep.
!
! field:  the 3D array to solve (uVel, vVel, wVel, pp, or enthalpy)
! klo,khi: k-direction loop bounds
! jlo,jhi: j-direction loop bounds
! ibc:     i-index of boundary cell (pr/qr initialization)
! ilo,ihi: i-direction solve bounds (forward sweep ilo..ihi)
!********************************************************************
subroutine tdma_solve_3d(field, klo, khi, jlo, jhi, ibc, ilo, ihi)
	real(wp), intent(inout) :: field(nx,ny,nz)
	integer, intent(in) :: klo, khi, jlo, jhi, ibc, ilo, ihi

	integer :: i, j, k, ksweep, jsweep
	real(wp) :: pr(nx), qr(nx), d, denom

	do ksweep=1,2
	do k=khi,klo,-1
	do jsweep=1,2
!$OMP PARALLEL PRIVATE(pr, qr, d, denom)
!$OMP DO
	do j=jlo,jhi
		pr(ibc)=0.0
		qr(ibc)=field(ibc,j,k)

		do i=ilo,ihi
			d = at(i,j,k)*field(i,j,k+1)+ab(i,j,k)*field(i,j,k-1)+an(i,j,k)*field(i,j+1,k) &
				+as(i,j,k)*field(i,j-1,k)+su(i,j,k)
			denom=ap(i,j,k)-aw(i,j,k)*pr(i-1)

		    	if(denom.le.1e-12_wp .and.  denom.ge.0) denom=denom+1e-13_wp     !!avoid divide zero
		    	if(denom.ge.-1e-12_wp .and.  denom.lt.0) denom=denom-1e-13_wp     !!avoid divide zero

			pr(i)=ae(i,j,k)/(denom)
			qr(i)=(d+aw(i,j,k)*qr(i-1))/(denom)
		enddo
!-----back substitution----------------
		do i=ihi, ilo, -1
			field(i,j,k)=pr(i)*field(i+1,j,k)+qr(i)
		enddo
	enddo
!$OMP END PARALLEL
	enddo
	enddo
	enddo
end subroutine tdma_solve_3d

!********************************************************************
! Convenience wrappers that call tdma_solve_3d with correct bounds
!********************************************************************
subroutine solution_uvw(field)
	real(wp), intent(inout) :: field(nx,ny,nz)
	call tdma_solve_3d(field, kstat, nkm1, jstat, jend, istat, istatp1, iendm1)
end subroutine solution_uvw

subroutine solution_enthalpy
	call tdma_solve_3d(enthalpy, 2, nkm1, 2, njm1, 1, 2, nim1)
end subroutine solution_enthalpy


!********************************************************************
subroutine cleanuvw

	implicit none
	integer i,j,k
	real(wp) tulc,tvlc,twlc

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
end subroutine cleanuvw

end module
