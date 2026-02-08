!______________________________________________________________________________
!
module constant
!______________________________________________________________________________
!
	implicit none
	real g,pi,sigm,great,small
	parameter(g=9.8,pi=3.1415926,sigm=5.67e-8,great=1.0e20,small=1.0e-06)

	integer nx,ny,nz,nvar  		!maximum node number at xyz direction, maximum equation number
	integer nx1,ny1,nz1,ng 		!maximum region number
	integer TOOLLINES, COORDLINES
	parameter (nx=1200,ny=1200,nz=180,nvar=4)
	parameter (nx1=7,ny1=7,nz1=7,ng=5)
	parameter (TOOLLINES=1000, COORDLINES=5000)
	
	
end module constant
