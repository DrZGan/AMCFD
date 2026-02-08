!______________________________________________________________________________
!
module laserinput
!______________________________________________________________________________
!
	use geometry
	use constant
	use initialization

	implicit none

	integer istart,jstart,imin,imax,jmin,jmax,kmin,kmax
	!index of mid-length of melt pool, maximum and minimum indexes of melt pool dimension

	real heatin(nx,ny), heatinLaser, peakhin, timet
	!heat flux matrix, heat flux [J/s], peak laser density, time, x-value of central laser beam, y-value of central laser beam

	contains 

subroutine laser_beam


	integer i,j,k,iout,jout

	real xloc,rb2,varlas,xdist,ydist,dist2, xprime, yprime, zprime

	if(timet .gt. toolmatrix(PathNum,1) .and. toolmatrix(PathNum+1,1) .ge. -0.5)then
		PathNum = PathNum+1
		if(toolmatrix(PathNum,5) .ge. 0.5) TrackNum=TrackNum+1
	endif

	!print*, PathNum

	scanvelx=(toolmatrix(PathNum,2)-toolmatrix(PathNum-1,2))/(toolmatrix(PathNum,1)-toolmatrix(PathNum-1,1))
	scanvely=(toolmatrix(PathNum,3)-toolmatrix(PathNum-1,3))/(toolmatrix(PathNum,1)-toolmatrix(PathNum-1,1))
	beam_pos = beam_pos + delt*scanvelx 
	beam_posy = beam_posy+ delt*scanvely
	

	iout=0            !reserve intial index of x[]
	jout=0            !reserve intial index of y[]


	xloc=beam_pos     !inital x-coodinate value of laser beam center


	do i=2,nim1
		if (xloc.le.x(i)) go to 701   !x(i)>=xloc, goto701, iout=i-1 at this time
		iout=i
	enddo
701	if(abs(xloc-x(iout+1)).lt.abs(xloc-x(iout))) iout=iout+1    !judge which is more close to xloc, x(i) or x(i-1)
	istart=iout        !reserve index of x[], which is nearest from xloc

	!print*, istart
	!***!Find the index j of y(), which y(j)~beam_posy****

	do j=2,njm1
		if (beam_posy.le.y(j)) go to 702   !y(j)>=beam_posy, goto702, jout=j-1 at this time    
		jout=j
	enddo
702	if(abs(beam_posy-y(jout+1)).lt.abs(beam_posy-y(jout))) jout=jout+1    !judge which is more close to xloc, x(i) or x(i-1)
	jstart=jout        !reserve index of y[], which is nearest from beam_posy
	!print*, jstart
	!************************************************
     
	heatin=0.0         !2D matrix is unit matrix
	heatinLaser=0.0    
	rb2=alasrb**2                     !  rb2=square of beam radius

	if(toolmatrix(PathNum,5) .gt. 0.5)then
		varlas=alaspow*alaseta            !  varlas= effective laser power
	else
		varlas=0
	endif

	peakhin=alasfact*varlas/(pi*rb2)  !  peakhin= peak power density

	do i = 1,ni
		xdist=beam_pos-x(i)       ! distance between x-value and laser beam center
		do j=1,nj

			ydist=beam_posy-y(j)   !distance between y-value and pre-position of y

			!dist2 = xdist ** 2 + y(j) ** 2
			dist2=xdist**2+ydist**2

			heatin(i,j)=peakhin*exp(-alasfact*dist2/rb2)  ! Gaussian distribution
			heatinLaser=heatinLaser+areaij(i,j)*heatin(i,j) ! calcu total heat flux

		enddo
	enddo
			
!------------------initial dimension of the pool---------------
	imin=istart-2
	imax=istart+2

	jmin=jstart-2
	jmax=jstart+2

	kmin=nk-4
	kmax=nkm1



	return
end subroutine laser_beam

end module laserinput
