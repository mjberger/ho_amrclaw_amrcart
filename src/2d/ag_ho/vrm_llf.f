c
c -------------------------------------------------------------
c
      subroutine vrm_llf(qr,ql,rx,ixmin,ixmax,iflip,msize)

      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold


      parameter ( m1size = 5000 )

      dimension qr(msize,4),ql(msize,4)
      dimension ur(4), ul(4)
      dimension rx(msize,4),fl(4),fr(4)


      dimension clsqp(m1size),clsqm(m1size),pstar(m1size)
      dimension zp(m1size),zm(m1size),ustarm(m1size),ustarp(m1size)
      dimension zsum(m1size)
      data      itno/2/, small/1.d-6/
c
      if (ixmax .gt. m1size) then
         write(6,*)" need to increase m1size in vrm to ", ixmax
         stop
      endif

c
c     # vector Riemann solver
c     # Given conserved variables, compute
c     # conserved variables at interface by solving riemann problem
c     # only solve neighboring vertical interfaces
c     # as in Van Leer MUSCL paper # 5 (1979), JCP
c
c      iflip = 1 indicates the normal velocity is in
c                qr(k,2) and the tangential in qr(k,3)
c            = 0 indicates the reverse.
c
      inorm = 3 - iflip
      itan  = 2 + iflip
      do 10 k = ixmin, ixmax


         rhol = ql(k,1)
         rhoul   = ql(k,inorm)
         rhoutl  = ql(k,itan)
         El   = ql(k,4)
         pl = gamma1*(El-.5d0*(rhoul*rhoul+rhoutl*rhoutl)/rhol)
         cl = dsqrt(gamma * pl / rhol)
         sl = max(dabs(rhoul/rhol + cl), dabs(rhoul/rhol - cl) )

         rhor = qr(k,1)
         rhour   = qr(k,inorm)
         rhoutr  = qr(k,itan)
         Er   = qr(k,4)
         pr = gamma1*(Er-.5d0*(rhour*rhour+rhoutr*rhoutr)/rhor)
         cr = dsqrt(gamma * pr / rhor)
         sr = max(dabs(rhour/rhor + cr), dabs(rhour/rhor - cr) )


         s = max(sl,sr)

         ul(1) = rhol
         ul(inorm) = rhoul
         ul(itan) = rhoutl
         ul(4) = El

         ur(1) = rhor
         ur(inorm) = rhour
         ur(itan) = rhoutr
         ur(4) = Er

         fl(1) = rhoul
         fl(inorm) = rhoul*rhoul/rhol + pl
         fl(itan) = rhoul*rhoutl/rhol
         fl(4) = (El + pl) * rhoul/rhol

         fr(1) = rhour
         fr(inorm) = rhour*rhour/rhor + pr
         fr(itan) = rhour*rhoutr/rhor
         fr(4) = (Er + pr) * rhour/rhor

         rx(k,:) = 0.5d0*(fl(:)+fr(:)) + 0.5d0*s*(ul(:)-ur(:))

  10     continue






      return
      end
