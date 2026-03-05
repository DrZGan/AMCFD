module species_transport
	use geometry
	use initialization
	use discretization
	use dimensions
	use source
	use residue
	use parameters
	use constant
	implicit none


	contains
 
	
	subroutine species

	integer i,j,k
	real(8) vn,vs,ue,uw,wt,wb,fn,fs,fe,fw,ft,fb

	real(8) ds,dn,de,dw,dt,db
	real(8) delf,cp0,cp1
	real(8) densbydt,mass 

    	real(8) massdiffusivity(nx,ny,nz) ! = rou*D=7500*7e-7
  

	do k=1,nz
!$OMP PARALLEL
!$OMP DO
	do j=1,ny
	do i=1,nx
    		massdiffusivity(i,j,k)=5.25e-3
    		if(temp(i,j,k).le.tsolid)massdiffusivity(i,j,k)=1e-11
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	do k=2,nkm1
!$OMP PARALLEL PRIVATE(vn,ue,wt,fn,fe,ft,dn,de,dt,mass,densbydt)
!$OMP DO
	do j=2,njm1
	do i=2,nim1
		vn = vVel(i,j+1,k)*2
		ue = uVel(i+1,j,k)*2
		wt = wVel(i,j,k+1)*2
		
		
		
!-----convection coefficients--------------------------------------------------------
		fn = dens*vn*areaik(i,k)
		fe = dens*ue*areajk(j,k)
		ft = dens*wt*areaij(i,j)
!-----diffusion coefficients----------------------------------------
     
	 	dn = massdiffusivity(i,j,k)*areaik(i,k)*dypsinv(j+1)
		de = massdiffusivity(i,j,k)*areajk(j,k)*dxpwinv(i+1)
		dt = massdiffusivity(i,j,k)*areaij(i,j)*dzpbinv(k+1)
!-----coefficients (power law scheme)----------------------------------------
		mass = dn*max(0.0,(1.0-0.1*(abs(fn)/dn))**5)
		an(i,j,k) = mass+max(0.0,-fn)
		as(i,j+1,k) = mass+max(0.0,fn)

		mass = de*max(0.0,(1.0-0.1*(abs(fe)/de))**5)
		ae(i,j,k) = mass+max(0.0,-fe)
		aw(i+1,j,k) = mass+max(0.0,fe)

		mass = dt*max(0.0,(1.0-0.1*(abs(ft)/dt))**5)
		at(i,j,k) = mass+max(0.0,-ft)
		ab(i,j,k+1) = mass+max(0.0,ft)
	    	densbydt = dens/delt
		acpnot(i,j,k)=densbydt*volume(i,j,k)

		sp(i,j,k)=0.0
		su(i,j,k)=acpnot(i,j,k)*concentrationnot(i,j,k)
	
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	j=2
	do k=2,nkm1
		do i=2,nim1

			vs = vVel(i,j,k)*2
			
			fs = dens*vs*areaik(i,k)
			ds = massdiffusivity(i,j,k)*areaik(i,k)*dypsinv(j) 
			as(i,j,k) = ds*max(0.0,(1.0-0.1*(abs(fs)/ds))**5)+max(0.0,fs)
		enddo
	enddo

	i=2
	do k=2,nkm1
		do j=2,njm1

			uw = uVel(i,j,k)*2
			
			fw = dens*uw*areajk(j,k)
			dw = massdiffusivity(i,j,k)*areajk(j,k)*dxpwinv(i)
			aw(i,j,k) = dw*max(0.0,(1.0-0.1*(abs(fw)/dw))**5)+max(0.0,fw)
		enddo
	enddo

	k=2
	do j=2,njm1
		do i=2,nim1

			wb = wVel(i,j,k)*2
			
			fb = dens*wb*areaij(i,j)
			db = massdiffusivity(i,j,k)*areaij(i,j)*dzpbinv(k)
			ab(i,j,k) = db*max(0.0,(1.0-0.1*(abs(fb)/db))**5)+max(0.0,fb)
		enddo
	enddo

	end subroutine  species



	
subroutine solution_species

!-----TDMA

	implicit none
	integer i,j,k,ksweep,jsweep
	real(8) pr,qr,d,denom
	dimension pr(nx),qr(nx)

	do ksweep=1,2
	do k=nkm1,2,-1
	do jsweep=1,2
!$OMP PARALLEL PRIVATE(pr,qr,d,denom)
!$OMP DO
	do j=2,njm1
		pr(1)=0.0
		qr(1)=concentration(1,j,k)
		do i=2,nim1
			d = at(i,j,k)*concentration(i,j,k+1)+ab(i,j,k)*concentration(i,j,k-1)+an(i,j,k)*concentration(i,j+1,k)+ &
				as(i,j,k)*concentration(i,j-1,k)+su(i,j,k)
			denom=ap(i,j,k)-aw(i,j,k)*pr(i-1)
			pr(i)=ae(i,j,k)/denom
			qr(i)=(d+aw(i,j,k)*qr(i-1))/denom
		enddo

!-----back 
		do i=nim1,2,-1
			concentration(i,j,k)=pr(i)*concentration(i+1,j,k)+qr(i)
		enddo
	enddo
!$OMP END PARALLEL
	enddo
	enddo
	enddo

	return
end subroutine solution_species



subroutine boundary_species

	
	integer i,j,k
!-----k=nk
	do j=2,njm1
	do i=2,nim1
		if(temp(i,j,nk) .gt. tliquid) then
			concentration(i,j,nk)=2*(13-concentration(i,j,nk-1))/dzpbinv(nk)/5.25e-3+concentration(i,j,nk-1)
			!concentration(i,j,nk)=1

		else
			concentration(i,j,nk)=concentration(i,j,nk-1)
		endif
		!-concentration(i,j,nk-1)
	enddo
	enddo

!-----k=1 

	do j=2,njm1
	do i=2,nim1
		concentration(i,j,1)=0.0 
	enddo
	enddo

!-----j

	concentration(:,1,:)=concentration(:,2,:)
	concentration(:,nj,:)=concentration(:,nj-1,:)
!-----i
	concentration(1,:,:)=concentration(2,:,:)
	concentration(ni,:,:)=concentration(ni-1,:,:)
end subroutine boundary_species




subroutine enhance_species_speed
	implicit none
	integer i,j,k
	real(8) bl,blp,blm,blc       
	real(8) pib,qib,denom,delh
		
	dimension bl(nx), blp(nx), blm(nx), blc(nx), delh(nx), pib(nx),qib(nx)

	bl(1:ni)=0.0
	blp(1:ni)=0.0
	blm(1:ni)=0.0
	blc(1:ni)=0.0

	do k=2,nkm1
!$OMP PARALLEL 
!$OMP DO REDUCTION(+:bl,blp,blm,blc)
	do j=2,njm1
	do i=2,nim1
		bl(i)=bl(i)+ap(i,j,k)-an(i,j,k)-as(i,j,k)-at(i,j,k)-ab(i,j,k)
		blp(i)=blp(i)+ae(i,j,k)
		blm(i)=blm(i)+aw(i,j,k)
		blc(i)=blc(i)+ae(i,j,k)*concentration(i+1,j,k)+aw(i,j,k)*concentration(i-1,j,k)+an(i,j,k)*concentration(i,j+1,k)+ &
			as(i,j,k)*concentration(i,j-1,k)+at(i,j,k)*concentration(i,j,k+1)+ab(i,j,k)*concentration(i,j,k-1)+su(i,j,k) &
			-ap(i,j,k)*concentration(i,j,k)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	pib(2)=blp(2)/bl(2)
	qib(2)=blc(2)/bl(2)

	do i=3,nim1
		denom=bl(i)-blm(i)*pib(i-1)
		pib(i)=blp(i)/denom
		qib(i)=(blc(i)+blm(i)*qib(i-1))/denom
	enddo
	delh(nim1)=qib(nim1)

	do i=nim1-1,2,-1
		delh(i)=pib(i)*delh(i+1)+qib(i)
	enddo

	do k=2,nkm1
!$OMP PARALLEL
!$OMP DO 
	do j=2,njm1
	do i=2,nim1
		concentration(i,j,k)=concentration(i,j,k)+delh(i)
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	
	return

end subroutine enhance_species_speed

subroutine source_species

	integer i,j,k
	real(8) fraclu,fraclv,fraclw,variable1,volht
	real(8) term,term1,term3
	real(8) tulc,tvlc,twlc
	real(8) flew,flns,fltb

   	kp=1     !!!!!

!	if(scanvel.gt.small)then
!		do k=2,nkm1
!		do j=2,njm1
!			term1=areajk(j,k)*rhoscan
!			do i=2,nim1
!				su(i,j,k)=su(i,j,k)+concentration(i-1,j,k)*term1
!				sp(i,j,k)=sp(i,j,k)-term1
!			enddo
!		enddo
!		enddo
!	endif
	
	
	if(.not.steady) then
		variable1=dens/delt
	do k=2,nkm1
!$OMP PARALLEL PRIVATE(volht)
!$OMP DO
	do j=2,njm1
	do i=2,nim1
!		volht=volume(i,j,k)*variable1
		if(fracl(i,j,k).gt.0.0.and.(fraclnot(i,j,k)-fracl(i,j,k)).gt.0)then
		!su(i,j,k)=su(i,j,k)-kp*concentration(i,j,k)*variable1*volume(i,j,k)*(fraclnot(i,j,k)-fracl(i,j,k)) &
                 !  +variable1*volume(i,j,k)*((1-fracl(i,j,k))*concentration(i,j,k)-(1-fraclnot(i,j,k))*concentrationnot(i,j,k))
       		!sp(i,j,k)=volht*(1-fracl(i,j,k))
		!su(i,j,k)=su(i,j,k)-variable1*((1-fracl(i,j,k))*concentration(i,j,k)-(1-fraclnot(i,j,k))*concentrationnot(i,j,k))
		!concentrationaa(i,j,k)=concentration(i,j,k)*fracl(i,j,k)+(1-fracl(i,j,k))*concentration(i,j,k)*kp
	    	endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	endif

!----- k=nk & k=1 ------
	do j=2,njm1
	do i=2,nim1
		su(i,j,2)=su(i,j,2)+ab(i,j,2)*concentration(i,j,1)
		sp(i,j,2)=sp(i,j,2)-ab(i,j,2)
		ab(i,j,2)=0.0
		su(i,j,nkm1)=su(i,j,nkm1)+at(i,j,nkm1)*concentration(i,j,nk)
		sp(i,j,nkm1)=sp(i,j,nkm1)-at(i,j,nkm1)
		at(i,j,nkm1)=0.0
	enddo
	enddo

!----- j=1 & j=nj ------
	do k=2,nkm1
	do i=2,nim1
		su(i,2,k)=su(i,2,k)+as(i,2,k)*concentration(i,1,k)
		sp(i,2,k)=sp(i,2,k)-as(i,2,k)
		as(i,2,k)=0.0
		su(i,njm1,k)=su(i,njm1,k)+an(i,njm1,k)*concentration(i,nj,k)
		sp(i,njm1,k)=sp(i,njm1,k)-an(i,njm1,k)
		an(i,njm1,k)=0.0
	enddo
	enddo

!----- i=1 & i=ni---
	do k=2,nkm1
	do j=2,njm1
		su(2,j,k)=su(2,j,k)+aw(2,j,k)*concentration(1,j,k)
		sp(2,j,k)=sp(2,j,k)-aw(2,j,k)
		aw(2,j,k)=0.0
		su(nim1,j,k)=su(nim1,j,k)+ae(nim1,j,k)*concentration(ni,j,k)
		sp(nim1,j,k)=sp(nim1,j,k)-ae(nim1,j,k)
		ae(nim1,j,k)=0.0
	enddo
	enddo

	do k=2,nkm1
!$OMP PARALLEL
!$OMP DO
	do j=2,njm1
	do i=2,nim1
		ap(i,j,k)=an(i,j,k)+as(i,j,k)+ae(i,j,k)+aw(i,j,k)+at(i,j,k)+ab(i,j,k)+acpnot(i,j,k)-sp(i,j,k) 

!-----under-relaxation---
		ap(i,j,k)=ap(i,j,k)/urfh             ! 
		su(i,j,k)=su(i,j,k)+(1.-urfh)*ap(i,j,k)*concentration(i,j,k)
	enddo
	enddo
!$OMP END PARALLEL
	enddo


end subroutine source_species

subroutine residual_species
	integer i,j,k
	real(8) sumc,resor,sumd

	sumc=0.0
	sumd=0.0

	do k=2,nkm1
!$OMP PARALLEL
!$OMP DO REDUCTION(+:resor,sumd,sumc)
	do j=2,njm1
	do i=2,nim1
		resor=(an(i,j,k)*concentration(i,j+1,k)+as(i,j,k)*concentration(i,j-1,k)+ae(i,j,k)*concentration(i+1,j,k)+ &
			aw(i,j,k)*concentration(i-1,j,k)+at(i,j,k)*concentration(i,j,k+1)+ab(i,j,k)*concentration(i,j,k-1)+ &
			su(i,j,k))/ap(i,j,k)-concentration(i,j,k)
		sumd=sumd+abs(resor)
		sumc=sumc+abs(concentration(i,j,k))
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	resorc=sumd/(sumc+small)
	return
	
end subroutine residual_species

subroutine update_species
	integer i,j,k
	do k=1,nz
!$OMP PARALLEL
!$OMP DO
	do j=1,ny
	do i=1,nx
    		!concentration(i,j,k)=concentrationaa(i,j,k)
		!if(fracl(i,j,k).eq.0.0 .and. concentration(i,j,k).ne.1.0) then
		    !concentration(i,j,k)=concentration(i,j,k)*kp
		!endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	return
	
end subroutine update_species
end module
