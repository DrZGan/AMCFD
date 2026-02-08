!______________________________________________________________________________
!
module coefficients
!______________________________________________________________________________
! Discretization coefficient arrays and source terms.
! Split from the original 'initialization' God Module for better encapsulation.
!
	use constant
	implicit none

	! Discretization coefficients (neighbor links)
	real ap(nx,ny,nz), an(nx,ny,nz), as(nx,ny,nz)
	real ae(nx,ny,nz), aw(nx,ny,nz)
	real at(nx,ny,nz), ab(nx,ny,nz)
	real apnot(nx,ny,nz)

	! Source term arrays
	real su(nx,ny,nz), sp(nx,ny,nz)

	! Velocity correction coefficients
	real dux(nx,ny,nz), dvy(nx,ny,nz), dwz(nx,ny,nz)

end module coefficients
