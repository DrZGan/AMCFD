!______________________________________________________________________________
!
module field_vars
!______________________________________________________________________________
! Field variable arrays: velocity, pressure, enthalpy, temperature
! Split from the original 'initialization' God Module for better encapsulation.
!
	use constant
	implicit none

	! Velocity fields and their time-level-n ("not") counterparts
	real uVel(nx,ny,nz), vVel(nx,ny,nz), wVel(nx,ny,nz)
	real unot(nx,ny,nz), vnot(nx,ny,nz), wnot(nx,ny,nz)

	! Pressure and pressure correction
	real pressure(nx,ny,nz), pp(nx,ny,nz)

	! Enthalpy and temperature fields
	real enthalpy(nx,ny,nz), hnot(nx,ny,nz)
	real temp(nx,ny,nz), tnot(nx,ny,nz)

	! 4D combined variable array for u,v,w,pp (used by TDMA solver)
	real phi(nx,ny,nz,nvar)
	equivalence (phi(1,1,1,1), uVel(1,1,1)), (phi(1,1,1,2), vVel(1,1,1)), &
		(phi(1,1,1,3), wVel(1,1,1)), (phi(1,1,1,4), pp(1,1,1))

	! Liquid fraction
	real fracl(nx,ny,nz), fraclnot(nx,ny,nz)

end module field_vars
