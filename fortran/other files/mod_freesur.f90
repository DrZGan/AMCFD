!______________________________________________________________________________
!
module freesur
!______________________________________________________________________________
!
	use initialization
	use parameters
	use geometry
	use constant
	use boundary
	use laserinput
	implicit none
 	real  fi(nx,ny), zzm(nx,ny)  !surface defomation (Z minus is positive direction), previous surf defomation
	
    real  normx(nx,ny,nz), normy(nx,ny,nz), normz(nx,ny,nz), curv(nx,ny,nz)
	real  vd ,factor !volume addition
	real apdis(nx,ny)

	contains

subroutine freesurface
	
	real  alamopt, dx, xmid, alam1, alam2, fmid, f 
	integer i, j, k
	real  uVel_CT, vVel_CT, wVel_CT, surf_tens_grad

	if(Mod(INT(timet/delt),1).ne.0)	return 
	if(maxval(temp(1:ni,1:nj,nk)) .le. tsolid)	return 

!-----calculate pressure

	
	do j=1,nj
	do i=1,ni
		apdis(i,j)=7e4*exp(-5*((beam_pos-x(i))**2+(beam_posy-y(j))**2)/(alasrb**2))
	enddo
	enddo

!	Calculate optimum lagrangian parameter and free surface
	print*, 'Calculating free surface profile ...'

	vd=0

	alam1=0          ! inital value of lamda (Adjustiable)
	alam2=-1e5

	fmid=deltav(alam2)   !caclulate f and fmid using alam1 and alam2
	f=deltav(alam1)


	write(6,*) "f=",f,"fmid=",fmid

	if(f*fmid.ge.0 .and. f.ne.vd) print*, 'Warning: lamda should be bracked!'

	if(f.lt.0.)then       !search satisfied lamda by using Bisection Method
		alamopt=alam1
		dx=alam2-alam1
	else
		alamopt=alam2
		dx=alam1-alam2
	endif


	do j=1,50  !50 loops
		dx=dx*0.5
		xmid=alamopt+dx

		fmid=deltav(xmid)

		if(fmid.le.0.)alamopt=xmid

		if(abs(dx).lt.1.e-5.or.fmid.eq.0.) goto 500
	enddo


	print*, 'too many bisections'

500	continue
	
	
	write(6,1100) xmid,abs(fmid)
1100	format('Lamda: ', es14.4,2X,', Delta volume: ', es14.4)
!-----output free surface profile
	print*,"max value of fi= ", maxval(fi),"   min value of fi= ", minval(fi)


!-------------------------------


	do k=1,nkm1
!$OMP PARALLEL PRIVATE(factor)
!$OMP DO
	do j=1,nj
	do i=1,ni
		
		factor=fi(i,j)/z(nkm1)
			
		zr(i,j,k)=z(k)*(1.0-factor)
                
	enddo
	enddo
!$OMP END PARALLEL
	enddo


!$OMP PARALLEL
!$OMP DO
	do j=2,njm1     ! calculate normal vector and curvature
	do i=2,nim1
  



	enddo
	enddo
!$OMP END PARALLEL


	return

end subroutine  



real function deltav(alamda)
!	The following function in conjuction with energy balance
!	routine define the function (delta Volume) 
!	vd	-	volume of the additional droplet
	implicit none
	real alamda
	real delv	
	integer i, j, k
	
!	write(6,*) "alamda=",alamda

	call ebal(alamda)
	delv=0.0

	do j=2,njm1
	do i=2,nim1
		
		delv=delv+areaij(i,j)*fi(i,j)
		
	enddo
	enddo

	deltav=delv+vd

!	write(6,*) "deltav=",deltav

	return 

end function



subroutine ebal(alamda)
!	This subroutine define the free surface for a given
!	lamda (lagrangian parameter) 
!	gama - surface tension coefficient
	implicit none
	integer i, j, k

	real alamda
	real gama, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, diff1, difmax, &
	urfree, dely, dely1, delyp, delym, delx, delx1, delxp, delxm, dzdx, dzdy, dzdyp1, dzdym1, d2zdxdy,term1, term2
	integer niters
	gama=1
	c1=denl*g
	
!----- under-relaxation for free surface
	urfree=1.0
!	iteration starts here
	do niters=1,800  !  max niters=800

!$OMP PARALLEL PRIVATE(dely,dely1,delyp,delym,delx,delx1,delxp,delxm,c2,dzdx,dzdy,dzdyp1,dzdym1,d2zdxdy,c3,c4,c5,c6,c7,c8,c9,c10,c11,term1,term2)
!$OMP DO
!-----	calculate coefficients: c1 ... c11
		do j=2,njm1
			dely=yv(j+1)-yv(j-1)
			dely1=y(j+1)-y(j-1)
			delyp=y(j+1)-y(j)
			delym=y(j)-y(j-1)
		do i=2,nim1
			delx=xu(i+1)-xu(i-1)
			delx1=x(i+1)-x(i-1)
			delxp=x(i+1)-x(i)
			delxm=x(i)-x(i-1)
			if(temp(i,j,nk).lt.tsolid) cycle   !!if peak temp is less than solidus, NO calcu
			c2=apdis(i,j)+alamda
			dzdx=(zzm(i+1,j)-zzm(i-1,j))/delx1
			dzdy=(zzm(i,j+1)-zzm(i,j-1))/dely1
			dzdyp1=(zzm(i+1,j+1)-zzm(i+1,j-1))/dely1
			dzdym1=(zzm(i-1,j+1)-zzm(i-1,j-1))/dely1
			d2zdxdy=(dzdyp1-dzdym1)/delx1
			c3=(1.0+dzdx**2+dzdy**2)**1.5
			c4=1.0+dzdy**2
			c5=1.0+dzdx**2
			c6=2.0*dzdx*dzdy*d2zdxdy
			c7=zzm(i+1,j)/(delxp*delx)+zzm(i-1,j)/(delxm*delx)
			c8=zzm(i,j+1)/(delyp*dely)+zzm(i,j-1)/(delym*dely)
			c9=1.0/(delxp*delx)+1.0/(delxm*delx)
			c10=1.0/(delyp*dely)+1.0/(delym*dely)
			c11=gama/c3
			term1=c2+c11*(c4*c7+c5*c8-c6)
			term2=c1+c11*(c4*c9+c5*c10)
			fi(i,j)=max(-2e-3,urfree*term1/term2+(1.0-urfree)*zzm(i,j))   !!!!fi_min=-2e-3

			enddo
		enddo
!$OMP END PARALLEL

!-----	calculation error
		difmax=0.0

		do j=1,njm1
		do i=2,nim1
			diff1=abs((fi(i,j)-zzm(i,j))/zzm(i,j))
			difmax=amax1(diff1,difmax)
		enddo
		enddo
	!	write(6,*) "difmax=",difmax

		if(difmax.lt.5.0e-4) return
!----- set zzm equal to fi
		zzm=fi
	enddo
end subroutine ebal



end module


