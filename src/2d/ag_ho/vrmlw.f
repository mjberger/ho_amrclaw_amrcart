c
c -------------------------------------------------------------
c
      subroutine vrm(qr,ql,rx,ixmin,ixmax,iflip,msize)
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      common /comlw/ dt,h
      parameter ( m1size = 4033)
      dimension qr(msize,4),ql(msize,4)
      dimension rx(msize,4)
      dimension rhor(m1size),ur(m1size),pr(m1size),utr(m1size)
      dimension rhol(m1size),ul(m1size),pl(m1size),utl(m1size)

      data      itno/2/, small/1.d-6/
c
c
c     # law wendroff vector Riemann solver
c     # Given primitive variables, compute
c     # conserved variables at interface by solving riemann problem
c     # only solve neighboring vertical interfaces
c
c      iflip = 1 indicates the normal velocity is in
c                qr(k,2) and the tangential in qr(k,3)
c            = 0 indicates the reverse.
c   done so a single riemann solver can be used for both columns and rows.
c
      if (ixmax .gt. m1size) then
	 write(6,*)" need to increase m1size in vrm to ", ixmax
         stop
      endif

      inorm = 3 - iflip
      itan  = 2 + iflip
      dtoh = dt/h

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
      do 20 k = ixmin, ixmax
          el = pl(k)/gamma1 + .5*rhol(k)*(ul(k)**2+utl(k)**2)
          er = pr(k)/gamma1 + .5*rhor(k)*(ur(k)**2+utr(k)**2)
          frhor = rhor(k)*ur(k)
          frhol = rhol(k)*ul(k)
          fmomr = rhor(k)*ur(k)**2 + pr(k)
          fmoml = rhol(k)*ul(k)**2 + pl(k)
          ftanr = rhor(k)*ur(k)*utr(k)
          ftanl = rhol(k)*ul(k)*utl(k)
          fengr = ur(k)*(er + pr(k))
          fengl = ul(k)*(el + pl(k))

           rhohalf     = .5*(rhor(k)+rhol(k) - dtoh*(frhor- frhol))
           rmomhalf     = .5*(rhor(k)*ur(k)+rhol(k)*ul(k) -
     .                       dtoh*(fmomr - fmoml))
           tanhalf     = .5*(rhor(k)*utr(k)+rhol(k)*utl(k) -
     .                       dtoh*(ftanr-ftanl))
           enhalf      = .5*(er+el - dtoh*(fengr-fengl))
           phalf       = gamma1 * (enhalf - .5*rmomhalf**2/rhohalf 
     .                                    - .5*tanhalf**2/rhohalf)

           rx(k,1)     =  rmomhalf
           rx(k,inorm) =  rmomhalf**2/rhohalf + phalf
           rx(k,itan)  =  rmomhalf*tanhalf/rhohalf
           rx(k,4)     =  rmomhalf*(enhalf+phalf)/rhohalf
 20   continue

c
       return
       end
