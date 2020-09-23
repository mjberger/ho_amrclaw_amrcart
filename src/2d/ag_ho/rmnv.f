c
c ---------------------------------------------------------------------
c
      subroutine rmnv(ql,qr,qgdnv)
      implicit double precision (a-h,o-z)
      dimension ql(4), qr(4), qgdnv(4)
      data      itno/2/, small/1.d-6/
c
      gamma = 1.4
      gamma1 = gamma - 1.
      gp1g2i = (gamma+1.d0)/(2.d0*gamma)
c
c     # Riemann solver in x-direction
c     # Given states ql and qr, compute flux at interfaces.

c     # change variables from q1 = rho, q2 = rho*u, ...  to
c     # rho, u, p:
c
         rhol = ql(1)
         rhor = qr(1)
         ul = ql(2)/rhol
         ur = qr(2)/rhor
         utl = ql(3)/rhol
         utr = qr(3)/rhor
         pl = gamma1* (ql(4)- .5d0* (ql(2)*ql(2)/ql(1)  
     &      + ql(3)* ql(3)/ql(1)))
         pr = gamma1* (qr(4)- .5d0* (qr(2)*qr(2)/qr(1)  
     &      + qr(3)* qr(3)/qr(1)))
c
c     # iterate to find pstar, ustar:
c
      clsl = gamma*rhol*pl
      wm = dsqrt(clsl)
         clsr = gamma*rhor*pr
         wp = wm
         wm = dsqrt(clsr)
         pstar = (wp*pr+ wm*pl- wm*wp
     &      *(ur- ul))/(wp+ wm)
         pstar = dmax1(pstar,small)
c
      do 30 iter = 1, itno
           wsqp = clsl*(1.d0+gp1g2i*(pstar/pl-1.d0))
           wp   = dsqrt(wsqp)
           wsqm = clsr*(1.d0+gp1g2i*(pstar/pr-1.d0))
           wm   = dsqrt(wsqm)
c 
           zp = 2.d0*wp*wsqp/(wsqp+clsl)
           zm = 2.d0*wm*wsqm/(wsqm+clsr)
c 
           ustarm = ur-(pr-pstar)/wm
           ustarp = ul+(pl-pstar)/wp
           zsum = zp+zm
           epsp = zp*zm*(ustarm-ustarp)/zsum
           pstar = pstar-epsp
           pstar = dmax1(pstar,small)
  30       continue
c
         ustar = (zp*ustarm+zm*ustarp)/zsum
         sein  = -dsign(1.d0,ustar)
	 if (sein .ge. 0.0) then
            uo = ur
            po = pr
            ro = rhor
            utgdnv = utr
	 else
            uo = ul
            po = pl
            ro = rhol
            utgdnv = utl
	 endif
         co = dsqrt(po*gamma/ro)

         dummy = pstar-po
         wo = gamma*ro*po*(1.d0+gp1g2i*dummy/po)
         rstar = ro/(1.d0-ro*dummy/wo)
         wo = dsqrt(wo)
         cstar = dsqrt(gamma*pstar/rstar)
         wso = sein*uo+co
         wsi = sein*ustar+cstar
         ushok = sein*ustar+wo/rstar
	 if (dummy .ge. 0.d0) then
            wsi = ushok
            wso = ushok
	 else
            wsi = wsi
            wso = wso
	 endif
   
         t4 = (wso+wsi)/dmax1(wso-wsi,wso+wsi,small)
         t4 = .5d0*(t4+1.d0)
         t5  = 1.d0-t4
         t3  = t4*cstar+t5*co
         rgdnv = t4*rstar+t5*ro
         ugdnv = t4*ustar+t5*uo
         pgdnv = t3*t3*rgdnv/gamma
   
	 if (wso .ge. 0.d0) then
            pgdnv = pgdnv
            rgdnv = rgdnv
            ugdnv = ugdnv
	 else
            pgdnv = po
            rgdnv = ro
            ugdnv = uo
	 endif
   
	 if (wsi .lt. 0.d0) then
            pgdnv = pgdnv
            rgdnv = rgdnv
            ugdnv = ugdnv
	 else
            pgdnv = pstar
            rgdnv = rstar
            ugdnv = ustar
	 endif

         qgdnv(1) = rgdnv
         qgdnv(2) = rgdnv*ugdnv
         qgdnv(3) = rgdnv*utgdnv
         qgdnv(4) = pgdnv/gamma1+.5*rgdnv*(ugdnv**2+utgdnv**2)
c
      return
      end
