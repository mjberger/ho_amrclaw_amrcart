c
c -------------------------------------------------------------
c
      subroutine prm(qr,ql,rx,ixmin,ixmax,iflip,msize)
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      parameter ( m1size = 1033)
      dimension qr(msize,4),ql(msize,4)
      dimension rx(msize,4)
      dimension rhor(m1size),ur(m1size),pr(m1size),utr(m1size)
      dimension rhol(m1size),ul(m1size),pl(m1size),utl(m1size)
      dimension clsqp(m1size),clsqm(m1size),pstar(m1size)
      dimension zp(m1size),zm(m1size),ustarm(m1size),ustarp(m1size)
      dimension zsum(m1size)
      data      itno/2/, small/1.d-6/
c
      gp1g2i = (gamma+1.d0)/(2.d0*gamma)
c
c     # vector Riemann solver
c     # Given primitive variables, compute
c     # conserved variables at interface by solving riemann problem
c     # only solve neighboring vertical interfaces
c
c      iflip = 1 indicates the normal velocity is in
c                qr(k,2) and the tangential in qr(k,3)
c            = 0 indicates the reverse.
c   done so a single riemann solver can be used for both columns and rows.
c
      inorm = 3 - iflip
      itan  = 2 + iflip

      do 10 k = ixmin, ixmax
         rhor(k) = qr(k,1)
         ur(k)   = qr(k,inorm)
         utr(k)  = qr(k,itan)
         pr(k)   = qr(k,4)
         rhol(k) = ql(k,1)
         ul(k)   = ql(k,inorm)
         utl(k)  = ql(k,itan)
         pl(k)   = ql(k,4)
  10     continue
c
c     # iterate to find pstar, ustar:
c
      do 20 k = ixmin, ixmax
         clsqp(k) = gamma*rhol(k)*pl(k)
         clsqm(k) = gamma*rhor(k)*pr(k)
         wp = dsqrt(clsqp(k))
         wm = dsqrt(clsqm(k))
         pstar(k) = (wp*pr(k)+ wm*pl(k)- wm*wp
     &      *(ur(k)- ul(k)))/(wp+ wm)
         pstar(k) = dmax1(pstar(k),small)
  20     continue
c
      do 30 iter = 1, itno
	do 30 k = ixmin, ixmax
           wsqp = clsqp(k)*(1.d0+gp1g2i*(pstar(k)/pl(k)-1.d0))
           wp   = dsqrt(wsqp)
           wsqm = clsqm(k)*(1.d0+gp1g2i*(pstar(k)/pr(k)-1.d0))
           wm   = dsqrt(wsqm)
c 
           zp(k) = 2.d0*wp*wsqp/(wsqp+clsqp(k))
           zm(k) = 2.d0*wm*wsqm/(wsqm+clsqm(k))
c 
           ustarm(k) = ur(k)-(pr(k)-pstar(k))/wm
           ustarp(k) = ul(k)+(pl(k)-pstar(k))/wp
           zsum(k) = zp(k)+zm(k)
           epsp = zp(k)*zm(k)*(ustarm(k)-ustarp(k))/zsum(k)
           pstar(k) = pstar(k)-epsp
           pstar(k) = dmax1(pstar(k),small)
  30       continue
c
      do 40 k = ixmin, ixmax
         ustar = (zp(k)*ustarm(k)+zm(k)*ustarp(k))/zsum(k)
         sein  = -dsign(1.d0,ustar)
	 if (sein .ge. 0.d0) then
            uo = ur(k)
            po = pr(k)
            ro = rhor(k)
            utgdnv = utr(k)
	 else
	    uo = ul(k)
            po = pl(k)
            ro = rhol(k)
            utgdnv = utl(k)
	 endif
         co = dsqrt(po*gamma/ro)

         dummy = pstar(k)-po
         wo = gamma*ro*po*(1.d0+gp1g2i*dummy/po)
         rstar = ro/(1.d0-ro*dummy/wo)
         wo = dsqrt(wo)
         cstar = dsqrt(gamma*pstar(k)/rstar)
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
            pgdnv = pstar(k)
            rgdnv = rstar
            ugdnv = ustar
	 endif

         rx(k,1) = rgdnv*ugdnv
         rx(k,inorm) = rx(k,1)*ugdnv + pgdnv
         rx(k,itan) = rx(k,1)*utgdnv
         rx(k,4) = rx(k,1)*(pgdnv*gamma/(rgdnv*gamma1)+
     &                .5d0*(ugdnv**2 + utgdnv**2))
  40     continue
c
c
       return
       end
