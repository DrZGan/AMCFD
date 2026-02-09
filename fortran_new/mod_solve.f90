!______________________________________________________________________________
!
module solver
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions
	
	! Set .true. to use tdma_solve_3d_2 for comparison; .false. = original tdma_solve_3d
	logical, parameter :: use_tdma2 = .true.
	
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



subroutine tdma_solve_3d_2(field, klo, khi, jlo, jhi, ibc, ilo, ihi)
	use omp_lib
	implicit none
	real(wp), intent(inout) :: field(nx,ny,nz)
	integer, intent(in) :: klo, khi, jlo, jhi, ibc, ilo, ihi
  
	integer :: i, j, k, ksweep
	integer :: nth, tid
	real(wp), allocatable :: prbuf(:,:), qrbuf(:,:)   ! per-thread buffers
	real(wp) :: d, denom
	integer :: nthreads
  
	! determine number of threads to allocate buffers for
	!$OMP PARALLEL
	!$OMP MASTER
	  nthreads = omp_get_num_threads()
	!$OMP END MASTER
	!$OMP END PARALLEL
  
	if (nthreads < 1) nthreads = 1
	allocate(prbuf(nx, nthreads))
	allocate(qrbuf(nx, nthreads))
  
	!$OMP PARALLEL DEFAULT(NONE) &
	!$OMP SHARED(field,prbuf,qrbuf,klo,khi,jlo,jhi,ibc,ilo,ihi,at,ab,an,as,ae,aw,ap,su) &
	!$OMP PRIVATE(tid,k,j,i,ksweep,d,denom)
	  tid = omp_get_thread_num() + 1  ! Fortran 1-based index for buffer
	  do ksweep = 1, 2
		do k = khi, klo, -1
		  !$OMP DO SCHEDULE(STATIC)
		  do j = jlo, jhi
			prbuf(ibc, tid) = 0.0_wp
			qrbuf(ibc, tid) = field(ibc, j, k)
  
			! forward sweep along i; make it vector-friendly
			do i = ilo, ihi
			  d = at(i,j,k)*field(i,j,k+1) + ab(i,j,k)*field(i,j,k-1) &
				  + an(i,j,k)*field(i,j+1,k) + as(i,j,k)*field(i,j-1,k) + su(i,j,k)
  
			  denom = ap(i,j,k) - aw(i,j,k)*prbuf(i-1, tid)
  
			  ! clamp denom away from zero (branch-free style)
			  if (abs(denom) < 1.0e-12_wp) then
				 denom = sign(1.0e-12_wp, denom)    ! preserves sign
			  end if
  
			  prbuf(i, tid) = ae(i,j,k) / denom
			  qrbuf(i, tid) = (d + aw(i,j,k) * qrbuf(i-1, tid)) / denom
			end do
  
			! back substitution (i decreasing)
			do i = ihi, ilo, -1
			  field(i,j,k) = prbuf(i, tid) * field(i+1,j,k) + qrbuf(i, tid)
			end do
		  end do   ! j
		  !$OMP END DO
		  !$OMP BARRIER
		end do   ! k
	  end do   ! ksweep
	!$OMP END PARALLEL
  
	deallocate(prbuf, qrbuf)
  end subroutine tdma_solve_3d_2





!********************************************************************
! Convenience wrappers: call tdma_solve_3d or tdma_solve_3d_2 per use_tdma2
!********************************************************************
subroutine solution_uvw(field)
	real(wp), intent(inout) :: field(nx,ny,nz)
	if (use_tdma2) then
		call tdma_solve_3d_2(field, kstat, nkm1, jstat, jend, istat, istatp1, iendm1)
	else
		call tdma_solve_3d(field, kstat, nkm1, jstat, jend, istat, istatp1, iendm1)
	endif
end subroutine solution_uvw

subroutine solution_enthalpy
	if (use_tdma2) then
		call tdma_solve_3d_2(enthalpy, 2, nkm1, 2, njm1, 1, 2, nim1)
	else
		call tdma_solve_3d(enthalpy, 2, nkm1, 2, njm1, 1, 2, nim1)
	endif
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
