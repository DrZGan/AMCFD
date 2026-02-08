!______________________________________________________________________________
!
module entotemp
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions

	contains
	subroutine enthalpy_to_temp
	implicit none
	integer i,j,k



	do k=1,nk
!$OMP PARALLEL
!$OMP DO
	do j=1,nj
	do i=1,ni
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
