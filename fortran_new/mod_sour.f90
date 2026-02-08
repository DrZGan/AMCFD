!______________________________________________________________________________
!
module source
!______________________________________________________________________________
!
	use geometry
	use initialization
	use parameters
	use laserinput
	use dimensions
	use boundary
	implicit none

	real sourceinput(nx,ny,nz)

	contains

subroutine source_term

	integer i, j, k
	real fraclu, fraclv, fraclw, tw, volht
	real term, term1, term3
	real tulc, tvlc, twlc
	real flew, flns, fltb
	real alasetavol_rhf, sourcerad_rhf, sourcedepth_rhf

	select case (ivar)

!============================================================================
	case (1)   ! u-momentum source terms

	do k = kstat, nkm1
!$OMP PARALLEL PRIVATE(fraclu, term, term1)
!$OMP DO
	do j = jstat, jend
	do i = istatp1, iendm1
		fraclu = fracl(i,j,k) * (1.0 - fracx(i-1)) + fracl(i-1,j,k) * fracx(i-1)
		if (fraclu .gt. 0) then
			term = darcy_resistance(viscos, fraclu)
			sp(i,j,k) = sp(i,j,k) - term * volume_u(i,j,k)
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo

!----- k=nk top boundary treatment -----
	do j = jstat, jend
	do i = istatp1, iendm1
		su(i,j,nkm1) = su(i,j,nkm1) + at(i,j,nkm1) * uVel(i,j,nk)
		sp(i,j,nkm1) = sp(i,j,nkm1) - at(i,j,nkm1)
		at(i,j,nkm1) = 0.0
	enddo
	enddo

!----- ap, under-relaxation, zero velocity -----
	do k = kstat, nkm1
!$OMP PARALLEL PRIVATE(tulc)
!$OMP DO
	do j = jstat, jend
	do i = istatp1, iendm1
		ap(i,j,k) = an(i,j,k) + as(i,j,k) + ae(i,j,k) + aw(i,j,k) + at(i,j,k) + ab(i,j,k) + apnot(i,j,k) - sp(i,j,k)
		dux(i,j,k) = areajk(j,k) / ap(i,j,k)
		ap(i,j,k)  = ap(i,j,k) / urfu
		su(i,j,k)  = su(i,j,k) + (1. - urfu) * ap(i,j,k) * uVel(i,j,k)
		dux(i,j,k) = dux(i,j,k) * urfu

		tulc = min(temp(i,j,k), temp(i-1,j,k))
		if (tulc .le. tsolid) then
			call zero_velocity_coeff(su(i,j,k), sp(i,j,k), ap(i,j,k), an(i,j,k), as(i,j,k), ae(i,j,k), aw(i,j,k), at(i,j,k), ab(i,j,k))
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	return

!============================================================================
	case (2)   ! v-momentum source terms

	do k = kstat, nkm1
!$OMP PARALLEL PRIVATE(fraclv, term, term1)
!$OMP DO
	do j = jstat, jend
	do i = istatp1, iendm1
		fraclv = fracl(i,j,k) * (1.0 - fracy(j-1)) + fracl(i,j-1,k) * fracy(j-1)
		if (fraclv .gt. 0) then
			term = darcy_resistance(viscos, fraclv)
			sp(i,j,k) = sp(i,j,k) - term * volume_v(i,j,k)
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo

!----- k=nk top boundary treatment -----
	do j = jstat, jend
	do i = istatp1, iendm1
		su(i,j,nkm1) = su(i,j,nkm1) + at(i,j,nkm1) * vVel(i,j,nk)
		sp(i,j,nkm1) = sp(i,j,nkm1) - at(i,j,nkm1)
		at(i,j,nkm1) = 0.0
	enddo
	enddo

	do k = kstat, nkm1
!$OMP PARALLEL PRIVATE(tvlc)
!$OMP DO
	do j = jstat, jend
	do i = istatp1, iendm1
		ap(i,j,k) = an(i,j,k) + as(i,j,k) + ae(i,j,k) + aw(i,j,k) + at(i,j,k) + ab(i,j,k) + apnot(i,j,k) - sp(i,j,k)
		dvy(i,j,k) = areaik(i,k) / ap(i,j,k)
		ap(i,j,k)  = ap(i,j,k) / urfv
		su(i,j,k)  = su(i,j,k) + (1. - urfv) * ap(i,j,k) * vVel(i,j,k)
		dvy(i,j,k) = dvy(i,j,k) * urfv

		tvlc = min(temp(i,j,k), temp(i,j-1,k))
		if (tvlc .le. tsolid) then
			call zero_velocity_coeff(su(i,j,k), sp(i,j,k), ap(i,j,k), an(i,j,k), as(i,j,k), ae(i,j,k), aw(i,j,k), at(i,j,k), ab(i,j,k))
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	return

!============================================================================
	case (3)   ! w-momentum source terms

	do k = kstat, nkm1
!$OMP PARALLEL PRIVATE(fraclw, term, term1, tw)
!$OMP DO
	do j = jstat, jend
	do i = istatp1, iendm1
		fraclw = fracl(i,j,k) * (1.0 - fracz(k-1)) + fracl(i,j,k-1) * fracz(k-1)
		if (fraclw .gt. 0) then
			term = darcy_resistance(viscos, fraclw)
			sp(i,j,k) = sp(i,j,k) - term * volume_w(i,j,k)
			! Buoyancy
			tw = temp(i,j,k) * (1.0 - fracz(k-1)) + temp(i,j,k-1) * fracz(k-1)
			su(i,j,k) = su(i,j,k) + boufac * volume_w(i,j,k) * (tw - tsolid)
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	do k = kstat, nkm1
!$OMP PARALLEL PRIVATE(twlc)
!$OMP DO
	do j = jstat, jend
	do i = istatp1, iendm1
		ap(i,j,k) = an(i,j,k) + as(i,j,k) + ae(i,j,k) + aw(i,j,k) + at(i,j,k) + ab(i,j,k) + apnot(i,j,k) - sp(i,j,k)
		dwz(i,j,k) = areaij(i,j) / ap(i,j,k)
		ap(i,j,k)  = ap(i,j,k) / urfw
		su(i,j,k)  = su(i,j,k) + (1. - urfw) * ap(i,j,k) * wVel(i,j,k)
		dwz(i,j,k) = dwz(i,j,k) * urfw

		twlc = min(temp(i,j,k), temp(i,j,k-1))
		if (twlc .le. tsolid) then
			call zero_velocity_coeff(su(i,j,k), sp(i,j,k), ap(i,j,k), an(i,j,k), as(i,j,k), ae(i,j,k), aw(i,j,k), at(i,j,k), ab(i,j,k))
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	return

!============================================================================
	case (4)   ! pressure correction source terms

	do k = kstat, nkm1
!$OMP PARALLEL
!$OMP DO
	do j = jstat, jend
	do i = istatp1, iendm1
		ap(i,j,k) = an(i,j,k) + as(i,j,k) + ae(i,j,k) + aw(i,j,k) + at(i,j,k) + ab(i,j,k) - sp(i,j,k)
		if (temp(i,j,k) .le. tsolid) then
			call zero_velocity_coeff(su(i,j,k), sp(i,j,k), ap(i,j,k), an(i,j,k), as(i,j,k), ae(i,j,k), aw(i,j,k), at(i,j,k), ab(i,j,k))
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	return

!============================================================================
	case (5)   ! enthalpy source terms

	sourcedepth_rhf = sourcedepth * (RHF)**2
	alasetavol_rhf  = sourcedepth_rhf * 3700
	sourcerad_rhf   = sourcedepth_rhf * 0.37

!----- volumetric heat source -----
	do k = 2, nkm1
!$OMP PARALLEL
!$OMP DO
	do j = 2, njm1
	do i = 2, nim1
		if (toolmatrix(PathNum,5) .gt. 0.5) then
			if (z(nk) - z(k) .le. sourcedepth_rhf) then
				sourceinput(i,j,k) = alaspowvol * alasfact / pi / sourcerad_rhf**2 / sourcedepth_rhf * alasetavol_rhf * &
					exp(-alasfact / sourcerad_rhf**2 * ((beam_pos - x(i))**2 + (beam_posy - y(j))**2))
			else
				sourceinput(i,j,k) = 0.0
			endif
		else
			sourceinput(i,j,k) = 0.0
		endif
		su(i,j,k) = su(i,j,k) + volume(i,j,k) * sourceinput(i,j,k)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

!----- latent heat source terms -----
	do k = 2, nkm1
!$OMP PARALLEL PRIVATE(volht, flew, flns, fltb)
!$OMP DO
	do j = 2, njm1
	do i = 2, nim1
		volht = volume(i,j,k) * hlatnt * den(i,j,k) / delt
		su(i,j,k) = su(i,j,k) - volht * (fracl(i,j,k) - fraclnot(i,j,k))
		flew = areajk(j,k) * (max(uVel(i,j,k), 0.0) * fracl(i-1,j,k) - max(-uVel(i,j,k), 0.0) * fracl(i,j,k) &
			+ max(-uVel(i+1,j,k), 0.0) * fracl(i+1,j,k) - max(uVel(i+1,j,k), 0.0) * fracl(i,j,k))
		flns = areaik(i,k) * (max(vVel(i,j,k), 0.0) * fracl(i,j-1,k) - max(-vVel(i,j,k), 0.0) * fracl(i,j,k) &
			+ max(-vVel(i,j+1,k), 0.0) * fracl(i,j+1,k) - max(vVel(i,j+1,k), 0.0) * fracl(i,j,k))
		fltb = areaij(i,j) * (max(wVel(i,j,k), 0.0) * fracl(i,j,k-1) - max(-wVel(i,j,k), 0.0) * fracl(i,j,k) &
			+ max(-wVel(i,j,k+1), 0.0) * fracl(i,j,k+1) - max(wVel(i,j,k+1), 0.0) * fracl(i,j,k))
		su(i,j,k) = su(i,j,k) + den(i,j,k) * hlatnt * (flew + flns + fltb)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

!----- Transfer boundary coefficients to source terms -----
	! k=nk & k=1
	do j = 2, njm1
	do i = 2, nim1
		su(i,j,2)    = su(i,j,2)    + ab(i,j,2)    * enthalpy(i,j,1)
		sp(i,j,2)    = sp(i,j,2)    - ab(i,j,2)
		ab(i,j,2)    = 0.0
		su(i,j,nkm1) = su(i,j,nkm1) + at(i,j,nkm1) * enthalpy(i,j,nk)
		sp(i,j,nkm1) = sp(i,j,nkm1) - at(i,j,nkm1)
		at(i,j,nkm1) = 0.0
	enddo
	enddo

	! j=1 & j=nj
	do k = 2, nkm1
	do i = 2, nim1
		su(i,2,k)    = su(i,2,k)    + as(i,2,k)    * enthalpy(i,1,k)
		sp(i,2,k)    = sp(i,2,k)    - as(i,2,k)
		as(i,2,k)    = 0.0
		su(i,njm1,k) = su(i,njm1,k) + an(i,njm1,k) * enthalpy(i,nj,k)
		sp(i,njm1,k) = sp(i,njm1,k) - an(i,njm1,k)
		an(i,njm1,k) = 0.0
	enddo
	enddo

	! i=1 & i=ni
	do k = 2, nkm1
	do j = 2, njm1
		su(2,j,k)    = su(2,j,k)    + aw(2,j,k)    * enthalpy(1,j,k)
		sp(2,j,k)    = sp(2,j,k)    - aw(2,j,k)
		aw(2,j,k)    = 0.0
		su(nim1,j,k) = su(nim1,j,k) + ae(nim1,j,k) * enthalpy(ni,j,k)
		sp(nim1,j,k) = sp(nim1,j,k) - ae(nim1,j,k)
		ae(nim1,j,k) = 0.0
	enddo
	enddo

	! ap and under-relaxation for enthalpy
	do k = 2, nkm1
!$OMP PARALLEL
!$OMP DO
	do j = 2, njm1
	do i = 2, nim1
		ap(i,j,k) = an(i,j,k) + as(i,j,k) + ae(i,j,k) + aw(i,j,k) + at(i,j,k) + ab(i,j,k) + apnot(i,j,k) - sp(i,j,k)
		ap(i,j,k) = ap(i,j,k) / urfh
		su(i,j,k) = su(i,j,k) + (1. - urfh) * ap(i,j,k) * enthalpy(i,j,k)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	return

	end select

end subroutine source_term
end module source
