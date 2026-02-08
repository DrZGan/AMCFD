!______________________________________________________________________________
!
module dimensions
!______________________________________________________________________________
!
	use initialization
	use parameters
	use laserinput
	implicit none

	real alen,depth,width,hpeak,tpeak,umax,vmax,wmax
        !   length, depth, half-width, peak enthalpy, peak temperature, maximum velocity at each direction

	integer	istat,jstat,kstat,iend,jend,istatp1,iendm1  !
	!assign in subroutine pool_size in mod_dimen. (object: solve momentum equations only in fluid region)

	contains

subroutine pool_size
	
	integer i,j,k
	real dtdxxinv,dtdzzinv,dtdyyinv,dep, wid
	real xxmax,xxmin,yymax,yymin

	tpeak=maxval(temp(1:ni,1:nj,1:nk))    !tpeak is the maximum temperature in the whole domain 
!	tpeak=maxval(temp(1:ni,1,nk))    
	if(tpeak.le.tsolid) then         !if tpeak less than solidus, melt pool dimension is zero 
		alen=0.0
		depth=0.0
		width=0.0
	return                           !return
	endif

!-----length--------------------- 
	imax=istart
	imin=istart
	alen=0.0

	do i=istart,nim1
		imax=i
		if (temp(i,jstart,nk).le.tsolid) exit
	end do
	dtdxxinv = (x(imax)-x(imax+1))/(temp(imax,jstart,nk)-temp(imax+1,jstart,nk))
	xxmax = x(imax) + (tsolid - temp(imax,jstart,nk))*dtdxxinv                !maximum position value at i direction

	do i=istart,2,-1
		imin=i
		if (temp(i,jstart,nk).lt.tsolid) exit
	end do
	dtdxxinv = (x(imin)-x(imin-1))/(temp(imin,jstart,nk)-temp(imin-1,jstart,nk))
	xxmin = x(imin)+(tsolid - temp(imin,jstart,nk))*dtdxxinv                  !minimum position value at i direction
	alen=xxmax-xxmin						     	  !melt pool length at i direction
	!print*, imax
!-----depth--------------------- 
	kmin = nkm1
	depth = 0.0

	do 20 i=2,nim1
	do 10 k=nkm1,2,-1 
	if (temp(i,jstart,k).lt.tsolid) go to 20
		kmin=min(kmin,k)

10	continue
20	continue

	kmin=kmin-1
	if (kmin.eq.1) then
		depth=z(nk)-z(1)
		go to 40  !!!

	endif
	do 30 i=2,nim1
		if (temp(i,jstart,kmin+1).lt.tsolid) go to 30
		dtdzzinv = (z(kmin)-z(kmin-1))/(temp(i,jstart,kmin)-temp(i,jstart,kmin-1))
		dep = z(nk)-z(kmin)+(temp(i,jstart,kmin)-tsolid)*dtdzzinv
		depth=max(dep,depth)
30	continue
40	continue
	kmax=nkm1

!-----width------------------------- 
	jmax=jstart
	jmin=jstart
	width=0.0

	do 60 i=2,nim1
	do 50 j=jstart,njm1
		if (temp(i,j,nk).lt.tsolid) go to 60
		jmax=max(jmax,j)
50	continue
60	continue

	jmax=jmax+1
	if(jmax.eq.jstart) go to 80
	do 70 i=2,nim1
		if (temp(i,jmax-1,nk).lt.tsolid) go to 70
		dtdyyinv = (y(jmax)-y(jmax+1))/(temp(i,jmax,nk)-temp(i,jmax+1,nk))
		wid = y(jmax)+(tsolid - temp(i,jmax,nk))*dtdyyinv
		yymax=max(wid,yymax)
70	continue
80	continue

	do 95 i=2,nim1
	do 85 j=jstart,2,-1
		if (temp(i,j,nk).lt.tsolid) go to 95
		jmin=min(jmin,j)
85	continue
95	continue

	yymin=y(jstart)
	jmin=jmin-1
	if(jmin.eq.jstart) go to 97
	do 96 i=2,nim1
		if (temp(i,jmin+1,nk).lt.tsolid) go to 96
		dtdyyinv = (y(jmin)-y(jmin-1))/(temp(i,jmin,nk)-temp(i,jmin-1,nk))
		wid = y(jmin)+(tsolid - temp(i,jmin,nk))*dtdyyinv
		yymin=min(wid,yymin)
96	continue
97	continue
	width = yymax-yymin



!----- define solution domain for momentum equations----------------------------
	istat=max(imin-3,2)
	iend=min(imax+3,nim1)

	jstat=max(jmin-3,2)   
	jend=min(jmax+2,njm1)

	kstat=max(kmin-2,3)
	!kstat=2
	istatp1=istat+1
	iendm1=iend-1

	return

end subroutine pool_size

end module dimensions
