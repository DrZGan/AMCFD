!______________________________________________________________________________
!
module field_data
!______________________________________________________________________________
! Physical field arrays (velocity, pressure, enthalpy, temperature, liquid fraction)
!
	use precision
	implicit none

	!--------------------------------------------------------------
	! Derived type: field_3d_t (generic 3D field with init/destroy)
	!--------------------------------------------------------------
	type :: field_3d_t
		real(wp), allocatable :: data(:,:,:)
	contains
		procedure :: init => field_init
		procedure :: destroy => field_destroy
	end type field_3d_t

	!--------------------------------------------------------------
	! Module-level variables (backward compatibility)
	!--------------------------------------------------------------
	! Velocity fields + previous time step
	real(wp), allocatable :: uVel(:,:,:), vVel(:,:,:), wVel(:,:,:)
	real(wp), allocatable :: unot(:,:,:), vnot(:,:,:), wnot(:,:,:)

	! Pressure fields
	real(wp), allocatable :: pressure(:,:,:), pp(:,:,:)

	! Enthalpy and temperature fields + previous
	real(wp), allocatable :: enthalpy(:,:,:), hnot(:,:,:)
	real(wp), allocatable :: temp(:,:,:), tnot(:,:,:)

	! Liquid fraction fields
	real(wp), allocatable :: fracl(:,:,:), fraclnot(:,:,:)

	! Solid field tracking
	real(wp), allocatable :: solidfield(:,:,:)
	real(wp), allocatable :: localfield(:,:,:)

	contains

subroutine allocate_field_data(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk
	allocate(uVel(nni,nnj,nnk), vVel(nni,nnj,nnk), wVel(nni,nnj,nnk))
	allocate(unot(nni,nnj,nnk), vnot(nni,nnj,nnk), wnot(nni,nnj,nnk))
	allocate(pressure(nni,nnj,nnk), pp(nni,nnj,nnk))
	allocate(enthalpy(nni,nnj,nnk), hnot(nni,nnj,nnk))
	allocate(temp(nni,nnj,nnk), tnot(nni,nnj,nnk))
	allocate(fracl(nni,nnj,nnk), fraclnot(nni,nnj,nnk))
	allocate(solidfield(nni,nnj,nnk))
	allocate(localfield(nni,nnj,nnk))
	localfield = 0.0_wp
end subroutine allocate_field_data

subroutine deallocate_field_data()
	if (allocated(uVel)) deallocate(uVel,vVel,wVel,unot,vnot,wnot)
	if (allocated(pressure)) deallocate(pressure,pp)
	if (allocated(enthalpy)) deallocate(enthalpy,hnot,temp,tnot)
	if (allocated(fracl)) deallocate(fracl,fraclnot)
	if (allocated(solidfield)) deallocate(solidfield)
	if (allocated(localfield)) deallocate(localfield)
end subroutine deallocate_field_data

!********************************************************************
subroutine field_init(self, nni, nnj, nnk)
	class(field_3d_t), intent(inout) :: self
	integer, intent(in) :: nni, nnj, nnk
	allocate(self%data(nni, nnj, nnk))
	self%data = 0.0_wp
end subroutine field_init

!********************************************************************
subroutine field_destroy(self)
	class(field_3d_t), intent(inout) :: self
	if (allocated(self%data)) deallocate(self%data)
end subroutine field_destroy

end module field_data
