!______________________________________________________________________________
!
module laserinput
!______________________________________________________________________________
!
	use geometry
	use constant
	use initialization

	implicit none

	integer istart, jstart, imin, imax, jmin, jmax, kmin, kmax

	real heatin(nx,ny), heatinLaser, peakhin, timet

	contains

subroutine laser_beam

	integer i, j, k, iout, jout
	real xloc, rb2, varlas, xdist, ydist, dist2

	if (timet .gt. toolmatrix(PathNum,1) .and. toolmatrix(PathNum+1,1) .ge. -0.5) then
		PathNum = PathNum + 1
		if (toolmatrix(PathNum,5) .ge. 0.5) TrackNum = TrackNum + 1
	endif

	scanvelx = (toolmatrix(PathNum,2) - toolmatrix(PathNum-1,2)) / (toolmatrix(PathNum,1) - toolmatrix(PathNum-1,1))
	scanvely = (toolmatrix(PathNum,3) - toolmatrix(PathNum-1,3)) / (toolmatrix(PathNum,1) - toolmatrix(PathNum-1,1))
	beam_pos  = beam_pos  + delt * scanvelx
	beam_posy = beam_posy + delt * scanvely

	iout = 0
	jout = 0
	xloc = beam_pos

	! Find nearest i-index to beam center
	do i = 2, nim1
		if (xloc .le. x(i)) exit
		iout = i
	enddo
	if (abs(xloc - x(iout+1)) .lt. abs(xloc - x(iout))) iout = iout + 1
	istart = iout

	! Find nearest j-index to beam center
	do j = 2, njm1
		if (beam_posy .le. y(j)) exit
		jout = j
	enddo
	if (abs(beam_posy - y(jout+1)) .lt. abs(beam_posy - y(jout))) jout = jout + 1
	jstart = jout

	! Calculate heat flux distribution
	heatin = 0.0
	heatinLaser = 0.0
	rb2 = alasrb**2

	if (toolmatrix(PathNum,5) .gt. 0.5) then
		varlas = alaspow * alaseta
	else
		varlas = 0
	endif

	peakhin = alasfact * varlas / (pi * rb2)

	do i = 1, ni
		xdist = beam_pos - x(i)
		do j = 1, nj
			ydist = beam_posy - y(j)
			dist2 = xdist**2 + ydist**2
			heatin(i,j) = peakhin * exp(-alasfact * dist2 / rb2)
			heatinLaser = heatinLaser + areaij(i,j) * heatin(i,j)
		enddo
	enddo

	! Initial dimension estimates of the pool
	imin = istart - 2
	imax = istart + 2
	jmin = jstart - 2
	jmax = jstart + 2
	kmin = nk - 4
	kmax = nkm1

	return
end subroutine laser_beam

end module laserinput
