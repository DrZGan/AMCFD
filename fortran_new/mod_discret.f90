!______________________________________________________________________________
!
module discretization
!______________________________________________________________________________
!
	use constant
	use geometry
	use initialization
	use parameters
	use dimensions
	use cfd_utils

	implicit none

	contains

!********************************************************************
subroutine discretize_u
	integer i,j,k
	real(wp) vn,vs,ue,uw,wt,wb,fn,fs,fe,fw,ft,fb
	real(wp) visu,visun,visus,visn,viss,vise,visw,visut,visub,vist,visb
	real(wp) ds,dn,de,dw,dt,db,delf,cp0,cp1
	real(wp) dudxp,dudxm,dvdxp,dvdxm,dwdxp,dwdxm

	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb, visu, visun, visus, visn, viss, vise, visw, visut, visub, &
!$OMP vist, visb, dn, ds, dw, de, dt, db, delf, cp0, cp1, dudxp, dudxm, dvdxp, dvdxm, dwdxp, dwdxm)
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
!-----velocities at cv faces------------------------
		vn=vVel(i,j+1,k)*(1.0-fracx(i-1))+vVel(i-1,j+1,k)*fracx(i-1)
		vs=vVel(i,j,k)*(1.0-fracx(i-1))+vVel(i-1,j,k)*fracx(i-1)
		ue=(uVel(i+1,j,k)+uVel(i,j,k))*0.5
		uw=(uVel(i-1,j,k)+uVel(i,j,k))*0.5
		wt=wVel(i,j,k+1)*(1.0-fracx(i-1))+wVel(i-1,j,k+1)*fracx(i-1)
		wb=wVel(i,j,k)*(1.0-fracx(i-1))+wVel(i-1,j,k)*fracx(i-1)
!-----convection coefficients---------------------------
		fn = vn*den(i,j,k)*areauik(i,k)
		fs = vs*den(i,j,k)*areauik(i,k)
		fe = ue*den(i,j,k)*areajk(j,k)
		fw = uw*den(i,j,k)*areajk(j,k)
		ft = wt*den(i,j,k)*areauij(i,j)
		fb = wb*den(i,j,k)*areauij(i,j)
!-----viscosity at cv faces (harmonic mean)------------------------------------
		visu = harmonic_mean(vis(i,j,k), vis(i-1,j,k), fracx(i-1))
		visun = harmonic_mean(vis(i,j+1,k), vis(i-1,j+1,k), fracx(i-1))
		visus = harmonic_mean(vis(i,j-1,k), vis(i-1,j-1,k), fracx(i-1))
		if(j.eq.njm1) then
			visn=visu
		else
			visn = harmonic_mean(visu, visun, fracy(j))
		endif
		if(j.eq.2) then
			viss=visu
		else
			viss = harmonic_mean(visu, visus, 1.0-fracy(j-1))
		endif

		vise=vis(i,j,k)
		visw=vis(i-1,j,k)
		visut = harmonic_mean(vis(i,j,k+1), vis(i-1,j,k+1), fracx(i-1))
		visub = harmonic_mean(vis(i,j,k-1), vis(i-1,j,k-1), fracx(i-1))

		if(k.eq.nkm1) then
			vist=visu
		else
			vist = harmonic_mean(visut, visu, 1.0-fracz(k))
		endif
		if(k.eq.2) then
			visb=visu
		else
			visb = harmonic_mean(visub, visu, fracz(k-1))
		endif
!-----diffusion coefficients----------------------------
		dn = visn*areauik(i,k)*dypsinv(j+1)
		ds = viss*areauik(i,k)*dypsinv(j)
		de = vise*areajk(j,k)/(xu(i+1)-xu(i))
		dw = visw*areajk(j,k)/(xu(i)-xu(i-1))
		dt = vist*areauij(i,j)*dzpbinv(k+1)
		db = visb*areauij(i,j)*dzpbinv(k)
!-----coefficients (power law scheme)------------------------
		an(i,j,k) = power_law_coeff(dn, fn)
		as(i,j,k) = power_law_coeff(ds, -fs)
		ae(i,j,k) = power_law_coeff(de, fe)
		aw(i,j,k) = power_law_coeff(dw, -fw)
		at(i,j,k) = power_law_coeff(dt, ft)
		ab(i,j,k) = power_law_coeff(db, -fb)
		apnot(i,j,k)=den(i,j,k)*volume_u(i,j,k)/delt
!-----su and sp--------------------------------
		delf=fn-fs+fe-fw+ft-fb
		cp0=max(0.0,delf)
		cp1=min(0.0,delf)
		su(i,j,k)=-cp1*uVel(i,j,k)
		su(i,j,k)=su(i,j,k)+areajk(j,k)*(pressure(i-1,j,k)-pressure(i,j,k))
		sp(i,j,k)=-cp0
		su(i,j,k)=su(i,j,k)+apnot(i,j,k)*unot(i,j,k)

		dudxp =(uVel(i+1,j,k)-uVel(i,j,k))/(xu(i+1)-xu(i))
		dudxm =(uVel(i,j,k)-uVel(i-1,j,k))/(xu(i)-xu(i-1))
		su(i,j,k) =su(i,j,k)+(vise*dudxp-visw*dudxm)*areajk(j,k)
		dvdxp =(vVel(i,j+1,k)-vVel(i-1,j+1,k))*dxpwinv(i)
		dvdxm =(vVel(i,j,k)-vVel(i-1,j,k))*dxpwinv(i)
		su(i,j,k) =su(i,j,k)+(visn*dvdxp-viss*dvdxm)*areauik(i,k)
		dwdxp=(wVel(i,j,k+1)-wVel(i-1,j,k+1))*dxpwinv(i)
		dwdxm=(wVel(i,j,k)-wVel(i-1,j,k))*dxpwinv(i)
		su(i,j,k)=su(i,j,k)+(vist*dwdxp-visb*dwdxm)*areauij(i,j)
	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine discretize_u

!********************************************************************
subroutine discretize_v
	integer i,j,k
	real(wp) vn,vs,ue,uw,wt,wb,fn,fs,fe,fw,ft,fb
	real(wp) visn,viss,visv,visve,visvw,vise,visw,visvt,visvb,vist,visb
	real(wp) ds,dn,de,dw,dt,db,delf,cp0,cp1
	real(wp) dudyp,dudym,dvdyp,dvdym,dwdyp,dwdym

	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb, visn, viss, visv, visve, visvw, vise, visw, visvt, visvb, &
!$OMP vist, visb, dn, ds, dw, de, dt, db, delf, cp0, cp1, dudyp, dudym, dvdyp, dvdym, dwdyp, dwdym)
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		vn=(vVel(i,j,k)+vVel(i,j+1,k))*0.5
		vs=(vVel(i,j,k)+vVel(i,j-1,k))*0.5
		ue=uVel(i+1,j,k)*(1.0-fracy(j-1))+uVel(i+1,j-1,k)*fracy(j-1)
		uw=uVel(i,j,k)*(1.0-fracy(j-1))+uVel(i,j-1,k)*fracy(j-1)
		wt=wVel(i,j,k+1)*(1.0-fracy(j-1))+wVel(i,j-1,k+1)*fracy(j-1)
		wb=wVel(i,j,k)*(1.0-fracy(j-1))+wVel(i,j-1,k)*fracy(j-1)
!-----convection coefficients----------------------------
		fn = vn*den(i,j,k)*areaik(i,k)
		fs = vs*den(i,j,k)*areaik(i,k)
		fe = ue*den(i,j,k)*areavjk(j,k)
		fw = uw*den(i,j,k)*areavjk(j,k)
		ft = wt*den(i,j,k)*areavij(i,j)
		fb = wb*den(i,j,k)*areavij(i,j)
!-----viscosity at cv faces (harmonic mean)------------------------------------
		visn=vis(i,j,k)
		viss=vis(i,j-1,k)
		visv = harmonic_mean(vis(i,j,k), vis(i,j-1,k), fracy(j-1))
		visve = harmonic_mean(vis(i+1,j,k), vis(i+1,j-1,k), fracy(j-1))
		visvw = harmonic_mean(vis(i-1,j,k), vis(i-1,j-1,k), fracy(j-1))
		if(i.eq.nim1) then
			vise=visv
		else
			vise = harmonic_mean(visv, visve, fracx(i))
		endif
		if(i.eq.2) then
			visw=visv
		else
			visw = harmonic_mean(visv, visvw, 1.0-fracx(i-1))
		endif

		visvt = harmonic_mean(vis(i,j,k+1), vis(i,j-1,k+1), fracy(j-1))
		visvb = harmonic_mean(vis(i,j,k-1), vis(i,j-1,k-1), fracy(j-1))
		if(k.eq.nkm1) then
			vist=visv
		else
			vist = harmonic_mean(visv, visvt, fracz(k))
		endif
		if(k.eq.2) then
			visb=visv
		else
			visb = harmonic_mean(visv, visvb, 1.0-fracz(k-1))
		endif
!-----diffusion coefficients----------------------------
		dn=visn*areaik(i,k)/(yv(j+1)-yv(j))
		ds=viss*areaik(i,k)/(yv(j)-yv(j-1))
		de=vise*areavjk(j,k)*dxpwinv(i+1)
		dw=visw*areavjk(j,k)*dxpwinv(i)
		dt=vist*areavij(i,j)*dzpbinv(k+1)
		db=visb*areavij(i,j)*dzpbinv(k)
!-----coefficients (power law scheme)--------------------------------
		an(i,j,k) = power_law_coeff(dn, fn)
		as(i,j,k) = power_law_coeff(ds, -fs)
		ae(i,j,k) = power_law_coeff(de, fe)
		aw(i,j,k) = power_law_coeff(dw, -fw)
		at(i,j,k) = power_law_coeff(dt, ft)
		ab(i,j,k) = power_law_coeff(db, -fb)
		apnot(i,j,k)=den(i,j,k)*volume_v(i,j,k)/delt
!-----su and sp----------------------------------------
		delf=fn-fs+fe-fw+ft-fb
		cp0=max(0.0,delf)
		cp1=min(0.0,delf)
		su(i,j,k)=-cp1*vVel(i,j,k)
		su(i,j,k)=su(i,j,k)+areaik(i,k)*(pressure(i,j-1,k)-pressure(i,j,k))
		sp(i,j,k)=-cp0
		su(i,j,k)=su(i,j,k)+apnot(i,j,k)*vnot(i,j,k)

		dudyp =(uVel(i+1,j,k)-uVel(i+1,j-1,k))*dypsinv(j)
		dudym =(uVel(i,j,k)-uVel(i,j-1,k))*dypsinv(j)
		su(i,j,k) =su(i,j,k)+(vise*dudyp-visw*dudym)*areavjk(j,k)
		dvdyp=(vVel(i,j+1,k)-vVel(i,j,k))/(yv(j+1)-yv(j))
		dvdym=(vVel(i,j,k)-vVel(i,j-1,k))/(yv(j)-yv(j-1))
		su(i,j,k) =su(i,j,k)+(visn*dvdyp-viss*dvdym)*areaik(i,k)
		dwdyp=(wVel(i,j,k+1)-wVel(i,j-1,k+1))*dypsinv(j)
		dwdym=(wVel(i,j,k)-wVel(i,j-1,k))*dypsinv(j)
		su(i,j,k)=su(i,j,k)+(vist*dwdyp-visb*dwdym)*areavij(i,j)

	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine discretize_v

!********************************************************************
subroutine discretize_w
	integer i,j,k
	real(wp) vn,vs,ue,uw,wt,wb,fn,fs,fe,fw,ft,fb
	real(wp) vis_w,viswn,visws,visn,viss,viswe,visww,vise,visw,vist,visb
	real(wp) ds,dn,de,dw,dt,db,delf,cp0,cp1
	real(wp) dudzp,dudzm,dvdzp,dvdzm,dwdzp,dwdzm

	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft,fb, vis_w, viswn, visws, visn, viss, viswe, visww, vise,visw, &
!$OMP vist, visb, dn, ds, de, dw, dt, db, delf, cp0, cp1, dudzp, dudzm, dvdzp, dvdzm, dwdzp, dwdzm)
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		vn=vVel(i,j+1,k)*(1.0-fracz(k-1))+vVel(i,j+1,k-1)*fracz(k-1)
		vs=vVel(i,j,k)*(1.0-fracz(k-1))+vVel(i,j,k-1)*fracz(k-1)
		ue=uVel(i+1,j,k)*(1.0-fracz(k-1))+uVel(i+1,j,k-1)*fracz(k-1)
		uw=uVel(i,j,k)*(1.0-fracz(k-1))+uVel(i,j,k-1)*fracz(k-1)
		wt=(wVel(i,j,k)+wVel(i,j,k+1))*0.5
		wb=(wVel(i,j,k)+wVel(i,j,k-1))*0.5
!-----calculate convection coefficients------------------------------------------------
		fn = vn*den(i,j,k)*areawik(i,k)
		fs = vs*den(i,j,k)*areawik(i,k)
		fe = ue*den(i,j,k)*areawjk(j,k)
		fw = uw*den(i,j,k)*areawjk(j,k)
		ft = wt*den(i,j,k)*areaij(i,j)
		fb = wb*den(i,j,k)*areaij(i,j)
!-----viscosity at cv faces (harmonic mean)--------------------------------------------
		vis_w = harmonic_mean(vis(i,j,k), vis(i,j,k-1), fracz(k-1))
		viswn = harmonic_mean(vis(i,j+1,k), vis(i,j+1,k-1), fracz(k-1))
		visws = harmonic_mean(vis(i,j-1,k), vis(i,j-1,k-1), fracz(k-1))
		if(j.eq.njm1) then
			visn=vis_w
		else
			visn = harmonic_mean(vis_w, viswn, fracy(j))
		endif
		if(j.eq.2) then
			viss=vis_w
		else
			viss = harmonic_mean(vis_w, visws, 1.0-fracy(j-1))
		endif
		viswe = harmonic_mean(vis(i+1,j,k), vis(i+1,j,k-1), fracz(k-1))
		visww = harmonic_mean(vis(i-1,j,k), vis(i-1,j,k-1), fracz(k-1))
		if(i.eq.nim1) then
			vise=vis_w
		else
			vise = harmonic_mean(vis_w, viswe, fracx(i))
		endif
		if(i.eq.2) then
			visw=vis_w
		else
			visw = harmonic_mean(vis_w, visww, 1.0-fracx(i-1))
		endif
		vist=vis(i,j,k)
		visb=vis(i,j,k-1)
!-----diffusion coefficients------------------------
		dn=visn*areawik(i,k)*dypsinv(j+1)
		ds=viss*areawik(i,k)*dypsinv(j)
		de=vise*areawjk(j,k)*dxpwinv(i+1)
		dw=visw*areawjk(j,k)*dxpwinv(i)
		dt=vist*areaij(i,j)/(zw(k+1)-zw(k))
		db=visb*areaij(i,j)/(zw(k)-zw(k-1))
!-----coefficients (power law scheme)------------------------------------------------
		an(i,j,k) = power_law_coeff(dn, fn)
		as(i,j,k) = power_law_coeff(ds, -fs)
		ae(i,j,k) = power_law_coeff(de, fe)
		aw(i,j,k) = power_law_coeff(dw, -fw)
		at(i,j,k) = power_law_coeff(dt, ft)
		ab(i,j,k) = power_law_coeff(db, -fb)
		apnot(i,j,k)=den(i,j,k)*volume_w(i,j,k)/delt
!-----su and sp------------------------------------
		delf=fn-fs+fe-fw+ft-fb
		cp0=max(0.0,delf)
		cp1=min(0.0,delf)
		su(i,j,k)=-cp1*wVel(i,j,k)
		su(i,j,k)=su(i,j,k)+areaij(i,j)*(pressure(i,j,k-1)-pressure(i,j,k))
		sp(i,j,k)=-cp0
		su(i,j,k)=su(i,j,k)+apnot(i,j,k)*wnot(i,j,k)

		dudzp  =(uVel(i+1,j,k)-uVel(i+1,j,k-1))*dzpbinv(k)
		dudzm  =(uVel(i,j,k)-uVel(i,j,k-1))*dzpbinv(k)
		su(i,j,k) =su(i,j,k)+(vise*dudzp-visw*dudzm)*areawjk(j,k)
		dvdzp=(vVel(i,j+1,k)-vVel(i,j+1,k-1))*dzpbinv(k)
		dvdzm=(vVel(i,j,k)-vVel(i,j,k-1))*dzpbinv(k)
		su(i,j,k) =su(i,j,k)+(visn*dvdzp-viss*dvdzm)*areawik(i,k)
		dwdzp=(wVel(i,j,k+1)-wVel(i,j,k))/(zw(k+1)-zw(k))
		dwdzm=(wVel(i,j,k)-wVel(i,j,k-1))/(zw(k)-zw(k-1))
		su(i,j,k)=su(i,j,k)+(vist*dwdzp-visb*dwdzm)*areaij(i,j)

	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine discretize_w

!********************************************************************
subroutine discretize_pp
	integer i,j,k
	real(wp) vn,vs,ue,uw,wt,wb,fn,fs,fe,fw,ft,fb,delf

	resorm=0.0
	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb, delf)
!$OMP DO REDUCTION(+: resorm)
	do j=jstat,jend
	do i=istatp1,iendm1
!-----main coefficients------------------------------------
		an(i,j,k) = areaik(i,k)*dvy(i,j+1,k)*den(i,j,k)
		as(i,j,k) = areaik(i,k)*dvy(i,j,k)*den(i,j,k)
		ae(i,j,k) = areajk(j,k)*dux(i+1,j,k)*den(i,j,k)
		aw(i,j,k) = areajk(j,k)*dux(i,j,k)*den(i,j,k)
		at(i,j,k) = areaij(i,j)*dwz(i,j,k+1)*den(i,j,k)
		ab(i,j,k) = areaij(i,j)*dwz(i,j,k)*den(i,j,k)
!-----velocities at cv faces----------------------------
		vn = vVel(i,j+1,k)
		vs = vVel(i,j,k)
		ue = uVel(i+1,j,k)
		uw = uVel(i,j,k)
		wt = wVel(i,j,k+1)
		wb = wVel(i,j,k)

		fn = vn*areaik(i,k)*den(i,j,k)
		fs = vs*areaik(i,k)*den(i,j,k)
		fe = ue*areajk(j,k)*den(i,j,k)
		fw = uw*areajk(j,k)*den(i,j,k)
		ft = wt*areaij(i,j)*den(i,j,k)
		fb = wb*areaij(i,j)*den(i,j,k)

		delf=fn-fs+fe-fw+ft-fb
		sp(i,j,k)=0.0
		su(i,j,k)=-delf

		resorm=resorm+abs(delf)

	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine discretize_pp

!********************************************************************
subroutine discretize_enthalpy
	integer i,j,k
	real(wp) vn,vs,ue,uw,wt,wb,fn,fs,fe,fw,ft,fb
	real(wp) difn,dife,dift,dn,de,dt,ds,dw,db,tmp1

	do k=2,nkm1
!$OMP PARALLEL PRIVATE(vn, ue, wt, fn, fe, ft, difn, dife, dift, dn, de, dt, tmp1)
!$OMP DO
	do j=2,njm1
	do i=2,nim1
		vn = vVel(i,j+1,k)
		ue = uVel(i+1,j,k)
		wt = wVel(i,j,k+1)
!-----convection coefficients--------------------------------------------------------
		fn = den(i,j,k)*vn*areaik(i,k)
		fe = den(i,j,k)*ue*areajk(j,k)
		ft = den(i,j,k)*wt*areaij(i,j)
!-----diffusion coefficients (harmonic mean)----------------------------------------
		if(j.eq.njm1) then
			difn=diff(i,j,k)
		else
			difn = harmonic_mean(diff(i,j,k), diff(i,j+1,k), fracy(j))
		endif
		if(i.eq.nim1) then
			dife=diff(i,j,k)
		else
			dife = harmonic_mean(diff(i,j,k), diff(i+1,j,k), fracx(i))
		endif
		if(k.eq.nkm1) then
			dift=diff(i,j,k)
		else
			dift = harmonic_mean(diff(i,j,k), diff(i,j,k+1), fracz(k))
		endif

		dn = difn*areaik(i,k)*dypsinv(j+1)
		de = dife*areajk(j,k)*dxpwinv(i+1)
		dt = dift*areaij(i,j)*dzpbinv(k+1)
!-----coefficients (power law scheme)----------------------------------------
		tmp1 = dn*max(0.0,(1.0-0.1*(abs(fn)/dn))**5)
		an(i,j,k) = tmp1+max(0.0,-fn)
		as(i,j+1,k) = tmp1+max(0.0,fn)

		tmp1 = de*max(0.0,(1.0-0.1*(abs(fe)/de))**5)
		ae(i,j,k) = tmp1+max(0.0,-fe)
		aw(i+1,j,k) = tmp1+max(0.0,fe)

		tmp1 = dt*max(0.0,(1.0-0.1*(abs(ft)/dt))**5)
		at(i,j,k) = tmp1+max(0.0,-ft)
		ab(i,j,k+1) = tmp1+max(0.0,ft)

		apnot(i,j,k)=den(i,j,k)/delt*volume(i,j,k)

		sp(i,j,k)=0.0
		su(i,j,k)=apnot(i,j,k)*hnot(i,j,k)

	enddo
	enddo
!$OMP END PARALLEL
	enddo

	j=2
	do k=2,nkm1
		do i=2,nim1
			vs = vVel(i,j,k)
			fs = den(i,j,k)*vs*areaik(i,k)
			ds = diff(i,j,k)*areaik(i,k)*dypsinv(j)
			as(i,j,k) = ds*max(0.0,(1.0-0.1*(abs(fs)/ds))**5)+max(0.0,fs)
		enddo
	enddo

	i=2
	do k=2,nkm1
		do j=2,njm1
			uw = uVel(i,j,k)
			fw = den(i,j,k)*uw*areajk(j,k)
			dw = diff(i,j,k)*areajk(j,k)*dxpwinv(i)
			aw(i,j,k) = dw*max(0.0,(1.0-0.1*(abs(fw)/dw))**5)+max(0.0,fw)
		enddo
	enddo

	k=2
	do j=2,njm1
		do i=2,nim1
			wb = wVel(i,j,k)
			fb = den(i,j,k)*wb*areaij(i,j)
			db = diff(i,j,k)*areaij(i,j)*dzpbinv(k)
			ab(i,j,k) = db*max(0.0,(1.0-0.1*(abs(fb)/db))**5)+max(0.0,fb)
		enddo
	enddo
end subroutine discretize_enthalpy

end module discretization
