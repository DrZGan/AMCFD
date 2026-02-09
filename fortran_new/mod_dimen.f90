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

	tpeak=maxval(temp(1:ni,1:nj,1:nk))
	if(tpeak.le.tsolid) then
		alen=0.0
		depth=0.0
		width=0.0
		return
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
	xxmax = x(imax) + (tsolid - temp(imax,jstart,nk))*dtdxxinv

	do i=istart,2,-1
		imin=i
		if (temp(i,jstart,nk).lt.tsolid) exit
	end do
	dtdxxinv = (x(imin)-x(imin-1))/(temp(imin,jstart,nk)-temp(imin-1,jstart,nk))
	xxmin = x(imin)+(tsolid - temp(imin,jstart,nk))*dtdxxinv
	alen=xxmax-xxmin

!-----depth---------------------
	kmin = nkm1
	depth = 0.0

	outer_depth: do i=2,nim1
		do k=nkm1,2,-1
			if (temp(i,jstart,k).lt.tsolid) cycle outer_depth
			kmin=min(kmin,k)
		end do
	end do outer_depth

	kmin=kmin-1
	if (kmin.eq.1) then
		depth=z(nk)-z(1)
	else
		do i=2,nim1
			if (temp(i,jstart,kmin+1).lt.tsolid) cycle
			dtdzzinv = (z(kmin)-z(kmin-1))/(temp(i,jstart,kmin)-temp(i,jstart,kmin-1))
			dep = z(nk)-z(kmin)+(temp(i,jstart,kmin)-tsolid)*dtdzzinv
			depth=max(dep,depth)
		end do
	endif
	kmax=nkm1

!-----width-------------------------
	jmax=jstart
	jmin=jstart
	width=0.0
	yymax=y(jstart)  ! Bug fix: initialize yymax before use
	yymin=y(jstart)

	outer_jmax: do i=2,nim1
		do j=jstart,njm1
			if (temp(i,j,nk).lt.tsolid) cycle outer_jmax
			jmax=max(jmax,j)
		end do
	end do outer_jmax

	jmax=jmax+1
	if(jmax.ne.jstart) then
		do i=2,nim1
			if (temp(i,jmax-1,nk).lt.tsolid) cycle
			dtdyyinv = (y(jmax)-y(jmax+1))/(temp(i,jmax,nk)-temp(i,jmax+1,nk))
			wid = y(jmax)+(tsolid - temp(i,jmax,nk))*dtdyyinv
			yymax=max(wid,yymax)
		end do
	endif

	outer_jmin: do i=2,nim1
		do j=jstart,2,-1
			if (temp(i,j,nk).lt.tsolid) cycle outer_jmin
			jmin=min(jmin,j)
		end do
	end do outer_jmin

	jmin=jmin-1
	if(jmin.ne.jstart) then
		do i=2,nim1
			if (temp(i,jmin+1,nk).lt.tsolid) cycle
			dtdyyinv = (y(jmin)-y(jmin-1))/(temp(i,jmin,nk)-temp(i,jmin-1,nk))
			wid = y(jmin)+(tsolid - temp(i,jmin,nk))*dtdyyinv
			yymin=min(wid,yymin)
		end do
	endif
	width = yymax-yymin

!----- define solution domain for momentum equations----------------------------
	istat=max(imin-3,2)
	iend=min(imax+3,nim1)

	jstat=max(jmin-3,2)
	jend=min(jmax+2,njm1)

	kstat=max(kmin-2,3)
	istatp1=istat+1
	iendm1=iend-1

end subroutine pool_size

end module dimensions
