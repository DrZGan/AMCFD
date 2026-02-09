!______________________________________________________________________________
!
module constant
!______________________________________________________________________________
!
	implicit none
	real g,pi,sigm,great,small
	parameter(g=9.8,pi=3.1415926,sigm=5.67e-8,great=1.0e20,small=1.0e-06)
	real, parameter :: conv_res_heat = 5.0e-4
	real, parameter :: conv_res_cool = 5.0e-7
	real, parameter :: ratio_upper = 1.01
	real, parameter :: ratio_lower = 0.99
	real, parameter :: laser_on_threshold = 0.5
	real, parameter :: powder_threshold = 0.5  ! solidfield threshold for powder region
	real, parameter :: vis_solid = 1.0e10      ! effective viscosity in solid region

	integer nx,ny,nz,nvar  		!maximum node number at xyz direction, maximum equation number
	integer nx1,ny1,nz1,ng 		!maximum region number
	integer TOOLLINES, COORDLINES
	parameter (nx=1200,ny=1200,nz=180,nvar=4)
	parameter (nx1=7,ny1=7,nz1=7,ng=5)
	parameter (TOOLLINES=1000, COORDLINES=5000)
	
	
end module constant
