c
c -------------------------------------------------------------
c
      subroutine vrm(qr,ql,rx,ixmin,ixmax,iflip,msize)
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      parameter ( m1size = 4033)
      dimension qr(msize,4),ql(msize,4)
      dimension rx(msize,4)
      dimension rhor(m1size),ur(m1size),pr(m1size),utr(m1size)
      dimension rhol(m1size),ul(m1size),pl(m1size),utl(m1size)

      data      itno/2/, small/1.d-6/
c
c
c     # van Leer vector Riemann solver
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

        orhor = 1.d0/rhor(k)
        unr   = ur(k)
        utanr = utr(k)

        orhol = 1.d0/rhol(k)
        unl   = ul(k)
        utanl = utl(k)

        cl2 = gamma*pl(k)*orhol
        cr2 = gamma*pr(k)*orhor
        cl = sqrt(cl2)
        cr = sqrt(cr2)
   
        odenom = .5d0/(gamma*gamma-1.d0)

        if (unl .ge. cl) then   ! completely supersonic in left state
          fl_rho     = rhol(k)*unl
          fl_thisMom = rhol(k)*unl*unl + pl(k)
          fl_tanMom  = rhol(k)*unl*utanl
          etot       = .5d0*rhol(k)*(unl*unl+utanl*utanl) + pl(k)/gamma1
          fl_rhoE    = unl*(etot + pl(k))
        else if (unl > -cl) then  ! subsonic left state (either in or out)
          f1l = gamma1*unl + 2.d0*cl
          fplus  = rhol(k) * .25d0 * (unl+cl)*(unl+cl)/cl
          fl_rho = fplus
          fl_thisMom = fplus * f1l /gamma
          fl_tanMom  = fplus * utanl
          fl_rhoE    = fplus * (utanl*utanl*.5d0 +  f1l*f1l*odenom)
        else 
          fl_rho     = 0.d0
          fl_thisMom = 0.d0
          fl_tanMom  = 0.d0
          fl_rhoE    = 0.d0
        endif

        if (unr .le. -cr) then   ! completely supersonic in right state
          fr_rho     = rhor(k)*unr
          fr_thisMom = rhor(k)*unr*unr
          fr_tanMom  = rhor(k)*unr*utanr
          etot       = .5d0*rhor(k)*(unr*unr+utanr*utanr) + pr(k)/gamma1
          fr_rhoE    = unr*(etot + pr(k))
        else if (unr < cr) then  ! subsonoic right state (either in or out)
          f1r = gamma1*unr - 2.d0*cr
          fminus = -rhor(k)*.25d0 * (unr-cr)*(unr-cr)/cr
          fr_rho  = fminus
          fr_thisMom = fminus*f1r /gamma
          fr_tanMom  = fminus * utanr
          fr_rhoE    = fminus * (utanr*utanr*.5d0 + f1r*f1r*odenom)
        else
          fr_rho     = 0.d0
          fr_thisMom = 0.d0
          fr_tanMom  = 0.d0
          fr_rhoE    = 0.d0
        endif


         rx(k,1)     = fl_rho + fr_rho
         rx(k,inorm) = fl_thisMom + fr_thisMom
         rx(k,itan)  = fl_tanMom  + fr_tanMom
         rx(k,4)     = fl_rhoE + fr_rhoE

 20   continue

c
       return
       end
