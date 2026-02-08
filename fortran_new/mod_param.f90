!______________________________________________________________________________
!
module parameters
!______________________________________________________________________________
!
	use constant

	implicit none
	real alaspow, alaseta, alasrb, alasfact, scanvelx, scanvely
	real alaspowvol, alasetavol, sourcerad, sourcedepth
	real dens, denl, viscos, tsolid, tliquid, tboiling, hsmelt, hlfriz, acpa, acpb, acpl, thconsa, thconsb, thconl, beta, emiss, dgdtp
	real layerheight, pden, pcpa, pcpb, pthcona, pthconb
	real delt, timax, urfu, urfv, urfw, urfp, urfh
	real xzone(nx1), yzone(ny1), zzone(nz1), powrx(nx1), powry(ny1), powrz(nz1)
	real htci, htcj, htck1, htckn, tempWest, tempEast, tempNorth, tempBottom, tempPreheat, tempAmb
	integer nzx, nzy, nzz, maxit
	integer ncvx(nx1), ncvy(ny1), ncvz(nz1)

	namelist / process_parameters / alaspow, alaseta, alasrb, alasfact
	namelist / volumetric_parameters / alaspowvol, alasetavol, sourcerad, sourcedepth
	namelist / material_properties / dens, denl, viscos, tsolid, tliquid, tboiling, hsmelt, hlfriz, acpa, acpb, acpl, &
		thconsa, thconsb, thconl, beta, emiss, dgdtp
	namelist / powder_properties / layerheight, pden, pcpa, pcpb, pthcona, pthconb
	namelist / numerical_relax / maxit, delt, timax, urfu, urfv, urfw, urfp, urfh
	namelist / boundary_conditions / htci, htcj, htck1, htckn, tempWest, tempEast, tempNorth, &
		tempBottom, tempPreheat, tempAmb

	contains

!--------------------------------------------------------------------------
! Convert temperature to enthalpy (solid phase, quadratic Cp model)
! h(T) = 0.5 * acpa * T^2 + acpb * T
!--------------------------------------------------------------------------
pure real function temp_to_enthalpy(T) result(h)
	real, intent(in) :: T
	h = 0.5 * acpa * T**2 + acpb * T
end function temp_to_enthalpy

!--------------------------------------------------------------------------
! Read input parameters from file
!--------------------------------------------------------------------------
subroutine read_data

	integer i, j, k

	open(unit=10, file='./inputfile/input_param.txt', form='formatted')

!-----geometrical parameters-------------------------------

	read(10, *)
	read(10, *) nzx
	read(10, *) (xzone(i), i=1, nzx)
	read(10, *) (ncvx(i), i=1, nzx)
	read(10, *) (powrx(i), i=1, nzx)
	read(10, *) nzy
	read(10, *) (yzone(i), i=1, nzy)
	read(10, *) (ncvy(i), i=1, nzy)
	read(10, *) (powry(i), i=1, nzy)
	read(10, *) nzz
	read(10, *) (zzone(i), i=1, nzz)
	read(10, *) (ncvz(i), i=1, nzz)
	read(10, *) (powrz(i), i=1, nzz)

	READ (10, NML=process_parameters)
	READ (10, NML=volumetric_parameters)
	READ (10, NML=material_properties)
	READ (10, NML=powder_properties)
	READ (10, NML=numerical_relax)
	READ (10, NML=boundary_conditions)

	close(10)
	return
end subroutine read_data

end module parameters
