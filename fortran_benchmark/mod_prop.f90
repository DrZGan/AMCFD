!______________________________________________________________________________
!
module property
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions

	contains

subroutine properties(ilo, ihi, jlo, jhi, klo, khi)
	implicit none
	integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
	integer i,j,k, jind
	real(wp) diffs, diffl
	real(wp) dL_mix, vMag, visT, diffT



	do k=klo,khi
!$OMP PARALLEL PRIVATE(jind, dL_mix, vMag, visT, diffT,diffs,diffl)
!$OMP DO
	do j=jlo,jhi
	do i=ilo,ihi

	
		visT=0 !dens*dL_mix*0.3*vMag
		diffT=visT/0.9

		diffs=(thconsa*temp(i,j,k)+thconsb)/(acpa*temp(i,j,k)+acpb)
		diffl=thconl/acpl
		vis(i,j,k)=(viscos+visT)                !temperature is larger than liquidus 
		diff(i,j,k)=(diffl+diffT)
		den(i,j,k)=denl

		if(temp(i,j,k).ge.tliquid) cycle          
			diff(i,j,k)=diffs          !temperature is less than solidus
			vis(i,j,k)=vis_solid
			den(i,j,k)=dens

			if(z(nk)-z(k) .le. layerheight .and. solidfield(i,j,k) .le. powder_threshold)then    !powder properties
				den(i,j,k)=pden
				vis(i,j,k)=vis_solid
				diff(i,j,k)=(pthcona*temp(i,j,k)+pthconb)/(pcpa*temp(i,j,k)+pcpb)
			endif

		if(temp(i,j,k).le.tsolid) cycle   !temperature is between liquidus and solidus           
			diff(i,j,k)=(fracl(i,j,k)*diffl+(1.0-fracl(i,j,k))*diffs)
			vis(i,j,k)=(viscos+visT)
			den(i,j,k)=(fracl(i,j,k)*denl+(1.0-fracl(i,j,k))*dens)

	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine properties	
end module property

