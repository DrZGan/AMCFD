!______________________________________________________________________________
!
module material_fields
!______________________________________________________________________________
! Material property field arrays: viscosity, diffusion, density.
! Split from the original 'initialization' God Module for better encapsulation.
!
	use constant
	implicit none

	real vis(nx,ny,nz), diff(nx,ny,nz), den(nx,ny,nz)

end module material_fields
