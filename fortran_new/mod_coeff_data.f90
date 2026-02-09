!______________________________________________________________________________
!
module coeff_data
!______________________________________________________________________________
! Discretization coefficients, material property arrays, source terms
!
	use precision
	implicit none

	! Material property arrays
	real(wp), allocatable :: vis(:,:,:), diff(:,:,:), den(:,:,:)

	! Discretization coefficients
	real(wp), allocatable :: ap(:,:,:), an(:,:,:), as(:,:,:), ae(:,:,:)
	real(wp), allocatable :: aw(:,:,:), at(:,:,:), ab(:,:,:), apnot(:,:,:)

	! Source terms
	real(wp), allocatable :: su(:,:,:), sp(:,:,:)

	! Velocity correction coefficients
	real(wp), allocatable :: dux(:,:,:), dvy(:,:,:), dwz(:,:,:)

	contains

subroutine allocate_coeff_data(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk
	allocate(vis(nni,nnj,nnk), diff(nni,nnj,nnk), den(nni,nnj,nnk))
	allocate(ap(nni,nnj,nnk), an(nni,nnj,nnk), as(nni,nnj,nnk), ae(nni,nnj,nnk))
	allocate(aw(nni,nnj,nnk), at(nni,nnj,nnk), ab(nni,nnj,nnk), apnot(nni,nnj,nnk))
	allocate(su(nni,nnj,nnk), sp(nni,nnj,nnk))
	allocate(dux(nni,nnj,nnk), dvy(nni,nnj,nnk), dwz(nni,nnj,nnk))
end subroutine allocate_coeff_data

subroutine deallocate_coeff_data()
	if (allocated(vis)) deallocate(vis,diff,den)
	if (allocated(ap)) deallocate(ap,an,as,ae,aw,at,ab,apnot)
	if (allocated(su)) deallocate(su,sp)
	if (allocated(dux)) deallocate(dux,dvy,dwz)
end subroutine deallocate_coeff_data

end module coeff_data
