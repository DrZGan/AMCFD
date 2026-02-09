!______________________________________________________________________________
!
module precision
!______________________________________________________________________________
!
! Defines working precision for the entire codebase.
! Double precision: selected_real_kind(15, 307)
! Single precision: selected_real_kind(6, 37)
!
	implicit none
	integer, parameter :: wp = selected_real_kind(6, 37)    ! single precision
end module precision
