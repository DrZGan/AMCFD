
!****************************************************************************
!
!  AM-CFD
!
!  Version 1.0
!
!  July, 2020, With English Annotation, Northwestern University
!
!****************************************************************************

program main
	use geometry
	use initialization
	use discretization
	use dimensions
	use boundary
	use source
	use residue
	use solver
	use fluxes
	use parameters
	use printing
	use constant
	use entotemp
	use convergence
	use property
	use revision
	use laserinput
	use toolpath

	implicit none
	integer i,j,k   ! ipoweroff=1 means laser off (move into inti.f90)
	real amaxres              ! residual error
	

	call read_data            ! read input_param into variables in memory
	call read_toolpath

	call generate_grid        ! calculate the xyz locations of nodes (staggered grid), area and volume of CV and some operators
	call OpenFiles            ! create/open files and write down the title
	call initialize           ! initial parameters and varialbes 




	call StartTime                                  ! get system time, input to screen and output.txt

	itertot=0        ! total iteration number

	
	timet=small      ! set timet=small (1e-7)

!------iteration Start------------------------------------Start total iteration index of GOTO 10
	
10	timet=timet+delt    !delt: time step 


	niter=0          ! iteration number in each time step
!-------Move laser and calculate freesurface-----------------!

	call laser_beam
	call read_coordinates
	call calcRHF
	
!-----iteration loop----------------Start iteration in each time step, index of GOTO 30



30	niter=niter+1         !interation num in time step+1
	itertot=itertot+1     !total interation num+1

!-----ivar=5------------solve ivar=5 equation (solve energy equation)


	ivar=5     

    call properties               !update temperature-dependent thermo-properties (viscosity, diffusion coefficient)          
	call bound_condition          !thermal boudary condition
	call discretize               !calc coefficients of discreted enthalpy conservation equation
	call source_term              !calc source term coefficients to complete discrete equation
	call residual                 !calc residual error in the domain
	

	call enhance_converge_speed   !calc residual error in each x direction slice to enhance converge speed (TDMA)
	call solution_enthalpy        !solve enthalpy conservation equation (Line-by-line TDMA)


	call enthalpy_to_temp         !translate from enthalpy to temperature


	call pool_size                !get melt pool dimension, start and end index of i,j,k to detemine fluid region
	if(tpeak.le.tsolid) goto 41   !if peak temperature less than solidus, skipping solution process of momentum eqution
	call cleanuvw                 !give velocity outside the melt pool a ZERO

!-----ivar=1,4------------
	do ivar=1,4
		call bound_condition   !
		call discretize        !
		call source_term       !
		call residual	       !
		call solution_uvw      !
		call revision_p        !
	enddo
41	continue                           

!-----convergence criterion------------

	call heat_fluxes                       !calcu criteria of energy conservation (ratio)
	amaxres=max(resorm, resoru,resorv,resorw)     

	if(toolmatrix(PathNum,5) .ge. 0.5)then	  ! laser is on
		if(amaxres.lt.5.0e-4 .and. ratio.le.1.01.and.ratio.ge.0.99) goto 50  ! transient-state criteria (heating stage) 
	else
		if(resorh.lt.5.0e-7) goto 50  ! transient-state criteria (cooling stage) 
	endif

	if(niter.lt.maxit) goto 30
50	continue


	call CalTime				!calcu current time

	call outputres			! output parameters to file and screen

		do k=1,nk
		do j=1,nj
		do i=1,ni
			if(temp(i,j,k).le.tsolid) then    ! 
				uVel(i,j,k)=0.0
				vVel(i,j,k)=0.0
				wVel(i,j,k)=0.0
			endif
			unot(i,j,k)=uVel(i,j,k)           !reserve previous time velocity, temperature, enthalpy and volume fraction of liquid
			vnot(i,j,k)=vVel(i,j,k)
			wnot(i,j,k)=wVel(i,j,k)
			tnot(i,j,k)=temp(i,j,k)
			hnot(i,j,k)=enthalpy(i,j,k)
			fraclnot(i,j,k)=fracl(i,j,k)

		enddo
		enddo
		enddo
		call Cust_Out 			!output thermal cycle, velocity and temperature at surface and symmetry plane to file



	if(timet.lt.timax) goto 10
	
	call EndTime                     ! output end time, close files


	stop
	end


