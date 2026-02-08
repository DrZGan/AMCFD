!______________________________________________________________________________
!
module property
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions

	contains

subroutine properties
	implicit none
	integer i, j, k, jind
	real diffs, diffl
	real dL_mix, vMag, visT, diffT

	do k = 1, nk
!$OMP PARALLEL PRIVATE(jind, dL_mix, vMag, visT, diffT, diffs, diffl)
!$OMP DO
	do j = 1, nj
	do i = 1, ni
		visT = 0
		diffT = visT / 0.9

		diffs = (thconsa * temp(i,j,k) + thconsb) / (acpa * temp(i,j,k) + acpb)
		diffl = thconl / acpl
		vis(i,j,k) = (viscos + visT)
		diff(i,j,k) = (diffl + diffT)
		den(i,j,k) = denl

		if (temp(i,j,k) .ge. tliquid) cycle
			diff(i,j,k) = diffs
			vis(i,j,k) = SOLID_VIS
			den(i,j,k) = dens

			if (z(nk) - z(k) .le. layerheight .and. solidfield(i,j,k) .le. 0.5) then
				den(i,j,k) = pden
				vis(i,j,k) = SOLID_VIS
				diff(i,j,k) = (pthcona * temp(i,j,k) + pthconb) / (pcpa * temp(i,j,k) + pcpb)
			endif

		if (temp(i,j,k) .le. tsolid) cycle
			diff(i,j,k) = (fracl(i,j,k) * diffl + (1.0 - fracl(i,j,k)) * diffs)
			vis(i,j,k) = (viscos + visT)
			den(i,j,k) = (fracl(i,j,k) * denl + (1.0 - fracl(i,j,k)) * dens)
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	return

end subroutine properties
end module property
