!______________________________________________________________________________
!
module fluxes
!______________________________________________________________________________
!
	use initialization
	use parameters
	use laserinput
	use source



	implicit none
	real heatout,flux_west,flux_east,flux_top,flux_bottom,flux_north,flux_south,accul,heatvol,ratio
	contains
	
subroutine heat_fluxes

	integer i,j,k
	real fluxi1,fluxl1,fluxk1,fluxn1,fluxm1,fluxj1,dh1
		

!-----i=1 & i=ni-----
	flux_west=0.0
	flux_east=0.0
	do k=2,nkm1
	do j=2,njm1
		fluxi1=diff(2,j,k)*(enthalpy(1,j,k)-enthalpy(2,j,k))*dxpwinv(2)           
		fluxl1=diff(nim1,j,k)*(enthalpy(ni,j,k)-enthalpy(nim1,j,k))*dxpwinv(ni)

		flux_west=flux_west+areajk(j,k)*(fluxi1)         
		flux_east=flux_east+areajk(j,k)*(fluxl1)
	enddo
	enddo

!********************************************************************
!-----k=nk and k=1------------
	flux_bottom=0.0
	flux_top=0.0
	do j=2,njm1
	do i=2,nim1
		fluxk1=diff(i,j,2)*(enthalpy(i,j,1)-enthalpy(i,j,2))*dzpbinv(2)
		fluxn1=diff(i,j,nkm1)*(enthalpy(i,j,nk)-enthalpy(i,j,nkm1))*dzpbinv(nk)
		flux_bottom=flux_bottom+areaij(i,j)*fluxk1
		flux_top=flux_top+areaij(i,j)*fluxn1
	enddo
	enddo

!********************************************************************
!-----j=1 and j=nj--------
	flux_north=0.0
	flux_south=0.0
	do k=2,nkm1
	do i=2,nim1
		fluxj1=diff(i,2,k)*(enthalpy(i,1,k)-enthalpy(i,2,k))*dypsinv(2)		
		fluxm1=diff(i,njm1,k)*(enthalpy(i,nj,k)-enthalpy(i,njm1,k))*dypsinv(nj)
		flux_south=flux_south+areaik(i,k)*fluxj1
		flux_north=flux_north+areaik(i,k)*fluxm1
	enddo
	enddo

!********************************************************************
!-----heat accumulation--------
		accul=0.0
		heatvol=0.0
		do k=2,nkm1
		do j=2,njm1
		do i=2,nim1
			dh1=enthalpy(i,j,k)-hnot(i,j,k)+(fracl(i,j,k)-fraclnot(i,j,k))*hlatnt
			accul=accul+volume(i,j,k)*den(i,j,k)*dh1/delt
			heatvol=heatvol+sourceinput(i,j,k)*volume(i,j,k)
		enddo
		enddo
		enddo


!********************************************************************
	heatout=flux_north+flux_bottom+flux_west+flux_east+flux_south-ahtoploss   !total heat loss 
	
	ratio=(heatvol+heatinLaser)/(accul-heatout)
	
	return



end subroutine heat_fluxes
end module fluxes
