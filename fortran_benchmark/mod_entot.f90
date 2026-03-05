!______________________________________________________________________________
!
module entotemp
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions

	contains
subroutine enthalpy_to_temp(ilo, ihi, jlo, jhi, klo, khi)
	implicit none
	integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
	integer i,j,k



	do k=max(1,klo-1),min(nk,khi+1)
!$OMP PARALLEL
!$OMP DO
	do j=max(1,jlo-1),min(nj,jhi+1)
	do i=max(1,ilo-1),min(ni,ihi+1)
		if(enthalpy(i,j,k).ge.hlcal) then
			fracl(i,j,k)=1.0
			temp(i,j,k)=(enthalpy(i,j,k)-hlcal)/acpl+tliquid
		elseif(enthalpy(i,j,k).le.hsmelt) then
			fracl(i,j,k)=0.0
			!temp(i,j,k)=tsolid-(hsmelt-enthalpy(i,j,k))/acp
			temp(i,j,k)=(sqrt(acpb**2+2*acpa*enthalpy(i,j,k))-acpb)/acpa
		else
			fracl(i,j,k)=(enthalpy(i,j,k)-hsmelt)/(hlcal-hsmelt)
			temp(i,j,k)=deltemp*fracl(i,j,k)+tsolid
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	return
!********************************************************************

end subroutine enthalpy_to_temp
end module entotemp
