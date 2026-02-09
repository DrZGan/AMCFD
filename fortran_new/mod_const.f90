!______________________________________________________________________________
!
module constant
!______________________________________________________________________________
!
	use precision
	implicit none
	real(wp) g,pi,sigm,great,small
	parameter(g=9.8_wp,pi=3.14159265358979323846_wp,sigm=5.670374419e-8_wp,great=1.0e20_wp,small=1.0e-06_wp)
	real(wp), parameter :: conv_res_heat = 5.0e-4_wp
	real(wp), parameter :: conv_res_cool = 5.0e-7_wp
	real(wp), parameter :: ratio_upper = 1.01_wp
	real(wp), parameter :: ratio_lower = 0.99_wp
	real(wp), parameter :: laser_on_threshold = 0.5_wp
	real(wp), parameter :: powder_threshold = 0.5_wp  ! solidfield threshold for powder region
	real(wp), parameter :: vis_solid = 1.0e10_wp      ! effective viscosity in solid region

	integer nx1,ny1,nz1,ng 		!maximum region number (for zone arrays in mod_param)
	integer TOOLLINES, COORDLINES
	parameter (nx1=7,ny1=7,nz1=7,ng=5)
	parameter (TOOLLINES=1000, COORDLINES=5000)
	
	
end module constant
