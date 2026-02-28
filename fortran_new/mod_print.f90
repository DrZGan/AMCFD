!______________________________________________________________________________
!
module printing
!______________________________________________________________________________
!
	use dimensions
	use geometry
	use initialization
	use residue
	use fluxes
	use parameters
	use laserinput
	use boundary

	implicit none
	integer itertot,niter  !main   

	real(wp) aAveSec           !how many seconds are needed for each iteration
	real(wp), allocatable :: auvl(:,:,:),avvl(:,:,:),awvl(:,:,:)  ! velocity field at central nodes

	integer, private:: i,j,k,ist,gridx, gridy, gridz ! calcu how many grids should be output in different axis
	integer, private:: itimestart,itimeend
	dimension iTimeStart(8),iTimeEnd(8)      !start time and end time
	real(wp) data_and_time,adtdys,area


	contains

!********************************************************************
subroutine allocate_print(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk
	allocate(auvl(nni,nnj,nnk), avvl(nni,nnj,nnk), awvl(nni,nnj,nnk))
end subroutine allocate_print

!********************************************************************
subroutine StartTime
	call date_and_time(values = iTimeStart)     !fortran API, get system time
	write(6,800)iTimeStart(1:3),iTimeStart(5:7)
	write(9,800)iTimeStart(1:3),iTimeStart(5:7)
800	format(2x,'Date: ',I4,'-',I2,'-',I2,2x,'time: ',2(I2,' :'),I2,/)
end subroutine StartTime

!********************************************************************
subroutine outputres


	
	if(tpeak.gt.tsolid) then
		umax=maxval(abs(uVel(istatp1:iendm1,jstat:jend,kstat:nkm1)))
		vmax=maxval(abs(vVel(istatp1:iendm1,jstat:jend,kstat:nkm1)))
		wmax=maxval(abs(wVel(istatp1:iendm1,jstat:jend,kstat:nkm1)))


	else
		umax=0.0
		vmax=0.0
		wmax=0.0
		
	endif

	



		write(9,3)timet,niter,aAveSec,itertot,resorh,resorm,resoru,resorv,resorw
		write(9,5)tpeak,umax,vmax,wmax,alen,depth,width
		write(9,2)flux_north,flux_south,ahtoploss,ahtoploss,flux_bottom,flux_west,flux_east,heatout,accul,heatinLaser,heatvol,ratio
		write(6,3)timet,niter,aAveSec,itertot,resorh,resorm,resoru,resorv,resorw
		write(6,5)tpeak,umax,vmax,wmax,alen,depth,width
		write(6,2)flux_north,flux_south,ahtoploss,ahtoploss,flux_bottom,flux_west,flux_east,heatout,accul,heatinLaser,heatvol,ratio
		write(6,7)coordhistory(1,1),coordhistory(1,2),coordhistory(1,3),coordhistory(1,4), &
		coordhistory(1,5),coordhistory(1,6),coordhistory(1,7),coordhistory(1,8)
		write(9,7)coordhistory(1,1),coordhistory(1,2),coordhistory(1,3),coordhistory(1,4), &
			coordhistory(1,5),coordhistory(1,6),coordhistory(1,7),coordhistory(1,8)
3		format('  time  iter  time/iter  tot_iter  res_enth  res_mass     res_u   res_v   res_w',/, &
		es9.2,1x,i4,2x,f7.3,3x,i7,2x,es8.1,2x,es8.1,1x,(3(es8.1,1x)))
2		format('  north  south  top  toploss  bottom  west  east  hout  accu  hin heatvol ratio',/, &
		3(f7.1),4(f7.1),4(f6.1),f7.2)
5		format('  Tmax        umax       vmax         wmax       length       depth     width',/, &
		f9.2,2x,3(es9.2,3x),3(es9.2,3x))
7		format('  time    beam_posx  beam_posy  beam_posz  power  scanspeed   speedx   speedy',/, &
		es9.2,3(es11.3),f7.1,3(f9.3),/)
		
			
	
end subroutine outputres


!********************************************************************
subroutine Cust_Out

	character(len=10) outputfilename
	integer outputnum
	integer outputintervel
	character( len = 3 ) :: cTemp
	integer :: npts
	real(kind=4) :: val4


	outputintervel=5               ! reserve transient datas every XX time steps
	
	if(Mod(INT(timet/delt),outputintervel).ne.0)	return       

	outputnum=INT(timet/delt)/outputintervel
	write( cTemp,'(i3)' ) outputnum
	! Open file for ASCII header
	open(unit=41,file='./result/vtkmov' //trim(adjustl( cTemp ))// '.vtk')


	gridx=0
	gridy=0
	gridz=0


	do i=2,nim1
		gridx=gridx+1
	enddo

	do j=2,njm1
		gridy=gridy+1
	enddo

	do k=2,nkm1
		gridz=gridz+1
	enddo
	
	
	npts = gridx * gridy * gridz

	! Write VTK legacy format header (ASCII)
	write(41,'(A)') '# vtk DataFile Version 3.0'
	write(41,'(A)') 'AMCFD Simulation Output'
	write(41,'(A)') 'BINARY'
	write(41,'(A)') 'DATASET STRUCTURED_GRID'
	write(41,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', gridx, ' ', gridy, ' ', gridz
	write(41,'(A,I0,A)') 'POINTS ', npts, ' float'
	close(41)
	
	! Reopen for binary append to write coordinates
	open(unit=41,file='./result/vtkmov' //trim(adjustl( cTemp ))// '.vtk', &
	     access='stream', form='unformatted', position='append', convert='big_endian')
	
	! Write point coordinates in binary
	do k=2,nkm1
	do j=2,njm1
	do i=2,nim1
		val4 = real(x(i), 4)
		write(41) val4
		val4 = real(y(j), 4)
		write(41) val4
		val4 = real(z(k), 4)
		write(41) val4
	enddo
	enddo
	enddo
	
	! Write POINT_DATA header (need to reopen as formatted to write ASCII)
	close(41)
	open(unit=41,file='./result/vtkmov' //trim(adjustl( cTemp ))// '.vtk', &
	     position='append')
	write(41,'(A,I0)') 'POINT_DATA ', npts
	close(41)
	
	! Reopen for binary append for field data
	open(unit=41,file='./result/vtkmov' //trim(adjustl( cTemp ))// '.vtk', &
	     access='stream', form='unformatted', position='append', convert='big_endian')



	do i=2,nim1
!$OMP PARALLEL
!$OMP DO
	do j=2,njm1
	do k=2,nkm1
		auvl(i,j,k)=(uVel(i,j,k)+uVel(i+1,j,k))*0.5
		avvl(i,j,k)=(vVel(i,j,k)+vVel(i,j+1,k))*0.5
		awvl(i,j,k)=(wVel(i,j,k)+wVel(i,j,k+1))*0.5
		if(temp(i,j,k).le.tsolid) then
			auvl(i,j,k)=0.0
			avvl(i,j,k)=0.0
			awvl(i,j,k)=0.0
		endif

	enddo
	enddo
!$OMP END PARALLEL
	enddo
 
!-----top plane--------
	do i=2,nim1
	do j=2,njm1
		auvl(i,j,nk)=(uVel(i,j,nk)+uVel(i+1,j,nk))*0.5
		avvl(i,j,nk)=(vVel(i,j,nk)+vVel(i,j+1,nk))*0.5
		if(temp(i,j,nk).le.tsolid) then
			auvl(i,j,nk)=0.0
			avvl(i,j,nk)=0.0
		endif
	enddo
	enddo

!-----symmetry plane
	do i=2,nim1
	do k=2,nkm1
		auvl(i,1,k)=(uVel(i,1,k)+uVel(i+1,1,k))*0.5
		awvl(i,1,k)=(wVel(i,1,k)+wVel(i,1,k+1))*0.5
		if(temp(i,1,k).le.tsolid) then
			auvl(i,1,k)=0.0
			awvl(i,1,k)=0.0
		endif
	enddo
	enddo

!-----left plane---------
	do j=2,njm1
	do k=2,nkm1
		avvl(1,j,k)=(vVel(1,j,k)+vVel(1,j+1,k))*0.5
		awvl(1,j,k)=(wVel(1,j,k)+wVel(1,j,k+1))*0.5
		if(temp(1,j,k).le.tsolid) then
			avvl(1,j,k)=0.0
			awvl(1,j,k)=0.0
		endif
	enddo
	enddo

	
	call write_vtk_vector(41, 'Velocity', auvl, avvl, awvl)
	call write_vtk_scalar(41, 'T', temp)
	call write_vtk_scalar(41, 'vis', vis)
	call write_vtk_scalar(41, 'diff', diff)
	call write_vtk_scalar(41, 'den', den)
	call write_vtk_scalar(41, 'solidID', solidfield)

	close(41)


end subroutine Cust_Out

!********************************************************************
subroutine write_vtk_vector(unit, name, ufield, vfield, wfield)
	integer, intent(in) :: unit
	character(len=*), intent(in) :: name
	real(wp), intent(in) :: ufield(:,:,:), vfield(:,:,:), wfield(:,:,:)
	integer i,j,k
	real(kind=4) :: val4

	write(unit) char(10), 'VECTORS ' // trim(name) // ' float', char(10)
	do k=2,nkm1
	do j=2,njm1
	do i=2,nim1
		val4 = real(ufield(i,j,k), 4)
		write(unit) val4
		val4 = real(vfield(i,j,k), 4)
		write(unit) val4
		val4 = real(wfield(i,j,k), 4)
		write(unit) val4
	enddo
	enddo
	enddo
end subroutine write_vtk_vector

!********************************************************************
subroutine write_vtk_scalar(unit, name, field)
	integer, intent(in) :: unit
	character(len=*), intent(in) :: name
	real(wp), intent(in) :: field(:,:,:)
	integer i,j,k
	real(kind=4) :: val4

	write(unit) char(10), 'SCALARS ' // trim(name) // ' float 1', char(10), 'LOOKUP_TABLE default', char(10)
	do k=2,nkm1
	do j=2,njm1
	do i=2,nim1
		val4 = real(field(i,j,k), 4)
		write(unit) val4
	enddo
	enddo
	enddo
end subroutine write_vtk_scalar


!********************************************************************
subroutine CalTime
	integer isecused
	call date_and_time(values = iTimeEnd)
	iSecUsed=86400*(iTimeEnd(3)-iTimeStart(3))+3600*(iTimeEnd(5)-iTimeStart(5))+60* &    !calcu the time has been used
		(iTimeEnd(6)-iTimeStart(6))+iTimeEnd(7)-iTimeStart(7)
	aAveSec=real(iSecUsed,wp)/real(itertot,wp)     ! calcu how many seconds are needed for each iteration
end subroutine CalTime

!********************************************************************
subroutine EndTime
	call date_and_time(values = iTimeEnd)
	write(6,807)iTimeEnd(1:3),iTimeEnd(5:7)
	write(9,807)iTimeEnd(1:3),iTimeEnd(5:7)
807	format(2x,'Date: ',I4,'-',I2,'-',I2,2x,'time: ',2(I2,':'),I2,/)
	if(iTimeEnd(7).lt.iTimeStart(7)) then
		iTimeEnd(7)=iTimeEnd(7)+60
		iTimeEnd(6)=iTimeEnd(6)-1
	endif
	if(iTimeEnd(6).lt.iTimeStart(6)) then
		iTimeEnd(6)=iTimeEnd(6)+60
		iTimeEnd(5)=iTimeEnd(5)-1
	endif
	if(iTimeEnd(5).lt.iTimeStart(5))	 iTimeEnd(5)=iTimeEnd(5)+24
	write(6,808)(iTimeEnd(5)-iTimeStart(5)),(iTimeEnd(6)-iTimeStart(6)),(iTimeEnd(7)-iTimeStart(7))
	write(9,808)(iTimeEnd(5)-iTimeStart(5)),(iTimeEnd(6)-iTimeStart(6)),(iTimeEnd(7)-iTimeStart(7))
808	format(2x,'Total time used:',I6,2x,'hr',I6,2x,'m',I6,2x,'s',/)

!----- close output file------------
	close(9)

	


end subroutine EndTime


!********************************************************************
subroutine OpenFiles

	
	open(unit=9,file='./result/output.txt')


	write(41,*)'TITLE = "Thermo-Capillary Flow in Laser-Generated Melt Pool"'

end subroutine OpenFiles

end module printing
