!______________________________________________________________________________
!
module toolpath
!______________________________________________________________________________
!
	use laserinput
	implicit none

	
	contains 

subroutine read_toolpath

	integer i,j,k,status

	open(unit=10,file='./ToolFiles/B26.crs',form='formatted')
	
	i=1
	do k = 1, TOOLLINES
		read(10,*, IOSTAT = status) (toolmatrix(i,j),j=1,5)
		if(status .eq. 0)then
			i=i+1
			cycle
		else if(status .lt. 0)then
			toolmatrix(i,1:5)=-1
			exit
		else
			print*, 'input toolpath error'
			exit
		endif
	end do
!********show toolmatrix*********
!	do i = 1, TOOLLINES
!
!		if(toolmatrix(i,1) .ne. -1)then
!			write(6,11) toolmatrix(i,1),toolmatrix(i,2),toolmatrix(i,3),toolmatrix(i,4),toolmatrix(i,5)
!11			format(5(es14.4))
!		else
!			write(6,11) toolmatrix(i,1),toolmatrix(i,2),toolmatrix(i,3),toolmatrix(i,4),toolmatrix(i,5)
!		 	exit
!		 endif
!	end do

	return
end subroutine read_toolpath

subroutine read_coordinates

	integer i,j,k, CoordNum

	CoordNum=1
	
	do i=1, COORDLINES
		if(coordhistory(i,1) .ge. 0) CoordNum = CoordNum+1 
	end do

!	print*, CoordNum

	if(CoordNum .eq. 1)then

	else if(CoordNum .gt. 1 .and. CoordNum .le. COORDLINES)then
		do i=CoordNum,2,-1
			do j=1,8
			coordhistory(i,j)=coordhistory(i-1,j)
			end do
		end do 

	else
		do i=COORDLINES,2,-1
			do j=1,8
			coordhistory(i,j)=coordhistory(i-1,j)
			end do
		end do
	endif

	coordhistory(1,1)=timet
	coordhistory(1,2)=beam_pos
	coordhistory(1,3)=beam_posy
	coordhistory(1,4)=z(nk)
	coordhistory(1,5)=toolmatrix(PathNum,5)*(alaspow+alaspowvol)
	coordhistory(1,6)=sqrt(scanvelx**2+scanvely**2)
	coordhistory(1,7)=scanvelx
	coordhistory(1,8)=scanvely


	do i=2,nim1   !!!! record solidID
!$OMP PARALLEL
!$OMP DO
	do j=2,njm1
	do k=2,nkm1

		if(temp(i,j,k).ge.tsolid) solidfield(i,j,k)=TrackNum
	enddo
	enddo
!$OMP END PARALLEL
	enddo
 

!*********show coordhistory*********
!	do i = 1, COORDLINES
!
!		if(coordhistory(i,1) .ne. -1)then
!			write(6,12) coordhistory(i,1),coordhistory(i,2),coordhistory(i,3),coordhistory(i,4),coordhistory(i,5),coordhistory(i,6),coordhistory(i,7),coordhistory(i,8)
!12			format(8(es14.4))
!		else
!			write(6,12) coordhistory(i,1),coordhistory(i,2),coordhistory(i,3),coordhistory(i,4),coordhistory(i,5),coordhistory(i,6),coordhistory(i,7),coordhistory(i,8)
!		 	exit
!		 endif
!	end do


	return
end subroutine read_coordinates


subroutine calcRHF

	integer i
	real disk, tk, Pk, R, T, P0, RHFk, RHFc

	RHFk=0
	RHF=0

	RHFc=3.176
	P0=300
	R=0.2e-3
	T=2e-3
	do i = 1, COORDLINES

		if(coordhistory(i,1) .ge. -0.5)then

			disk=sqrt((coordhistory(i,2)-coordhistory(1,2))**2+(coordhistory(i,3)-coordhistory(1,3))**2)
			tk=coordhistory(1,1)-coordhistory(i,1)
			if(disk .le. R .and. tk .le. T)then
				Pk=coordhistory(i,5)
				RHFk=(R-disk)**2/R**2*(T-tk)/T*Pk/P0

				RHF=RHF+RHFk
			endif
		else 
		 	exit
		 endif	
	end do

	RHF=RHF/RHFc

	write(6,20) RHF
	write(9,20) RHF
20	format('RHF= ',f8.3,/)


	return
end subroutine calcRHF
end module toolpath