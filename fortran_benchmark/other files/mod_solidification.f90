!______________________________________________________________________________
!
module solidification
!______________________________________________________________________________
!
	use initialization
	use dimensions

	contains

subroutine thermalparam

	implicit none
	integer i, j, k
	real gradx, grady, gradz
	
	if(tpeak .LE. tsolid)		return


!	Calculate temperature gradient and solidification growth rate


	print*, 'Calculating thermal parameters ...'

	do k=nkm1,1,-1
!$OMP PARALLEL PRIVATE(gradx, grady, gradz)
!$OMP DO
	do j=1,nj
	do i=1,ni

		if(temp(i,j,k) .ge. tliquid)then
			gradx=(temp(i,j,k)-temp(i-2,j,k))/(x(i)-x(i-1))/2
			grady=(temp(i,j+1,k)-temp(i,j-1,k))/(y(j+1)-y(i))/2
			gradz=(temp(i,j,k+1)-temp(i,j,k-1))/(z(k)-z(k-1))/2
			grad(i-1,j,k)=sqrt(gradx**2+grady**2+gradz**2)
			
			if(grad(i-1,j,k) .ge. 1e8) grad(i-1,j,k)=1e8
			rate(i-1,j,k)=scanvel*gradx/grad(i-1,j,k)

			exit
		endif 

 
	enddo
	enddo
!$OMP END PARALLEL
	enddo


	return

end subroutine  


end module


