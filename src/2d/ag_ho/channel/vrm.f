c
c -------------------------------------------------------------
c     LLF ON PHYSICAL VARIABLES

      subroutine vrm(qr,ql,rx,ixmin,ixmax,iflip,msize)

      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold

      common /order2/ ssw, quad, nolimiter
      logical quad
      parameter ( m1size = 5000 )

      dimension qr(msize,4),ql(msize,4)
      dimension uconl(4), uconr(4)
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
c
c      iflip = 1 indicates the normal velocity is in
c                qr(k,2) and the tangential in qr(k,3)
c            = 0 indicates the reverse.
c
      inorm = 3 - iflip
      itan  = 2 + iflip
      do 10 k = ixmin, ixmax

         if(quad) then ! everything is in primitive variables
         rhol = ql(k,1)
         ul   = ql(k,inorm)
         utl  = ql(k,itan)
         pl   = ql(k,4)
         El   = pl / (gamma-1.d0) + 0.5d0 * rhol * (ul**2 + utl**2)
         cl   = dsqrt(gamma * pl / rhol)
         sl   = max(dabs(ul + cl), dabs(ul - cl) )

         rhor = qr(k,1)
         ur   = qr(k,inorm)
         utr  = qr(k,itan)
         pr   = qr(k,4)
         Er   = pr / (gamma-1.d0) + 0.5d0 * rhor * (ur**2 + utr**2)
         cr   = dsqrt(gamma * pr / rhor)
         sr   = max(dabs(ur + cr), dabs(ur - cr) )

         s = max(sl,sr)

         uconl(1) = rhol
         uconl(inorm) = rhol*ul
         uconl(itan) = rhol*utl
         uconl(4) = El

         uconr(1) = rhor
         uconr(inorm) = rhor*ur
         uconr(itan) = rhor*utr
         uconr(4) = Er

         fl(1) = rhol * ul
         fl(inorm) = rhol*ul*ul + pl
         fl(itan) = rhol*ul*utl
         fl(4) = (El + pl) * ul

         fr(1) = rhor * ur
         fr(inorm) = rhor * ur * ur + pr
         fr(itan) = rhor * ur * utr
         fr(4) = (Er + pr) * ur
         else ! things are in conserved variables
         rhol = ql(k,1)
         ul   = ql(k,inorm)/ql(k,1)
         pl   = (gamma - 1.d0) * (ql(k,4)
     .               -0.5d0*(ql(k,inorm)**2 + ql(k,itan)**2)/ ql(k,1) )
         cl   = dsqrt(gamma * pl / rhol)
         sl   = max(dabs(ul + cl), dabs(ul - cl) )

         rhor = qr(k,1)
         ur   = qr(k,inorm)/qr(k,1)
         pr   = (gamma - 1.d0) * (qr(k,4)
     .                -0.5d0*(qr(k,inorm)**2 + qr(k,itan)**2)/ qr(k,1) )
         cr   = dsqrt(gamma * pr / rhor)
         sr   = max(dabs(ur + cr), dabs(ur - cr) )

         s = max(sl,sr)

         uconl(1) = ql(k,1)
         uconl(inorm) = ql(k,inorm)
         uconl(itan) = ql(k,itan)
         uconl(4) = ql(k,4)

         uconr(1) = qr(k,1)
         uconr(inorm) = qr(k,inorm)
         uconr(itan) = qr(k,itan)
         uconr(4) = qr(k,4)

         fl(1) = ql(k,inorm)
         fl(inorm) = (ql(k,inorm)**2)/ql(k,1) + pl
         fl(itan) = ql(k,inorm)*ql(k,itan)/ql(k,1)
         fl(4) = (ql(k,4) + pl) *  ql(k,inorm) /  ql(k,1)

         fr(1) = qr(k,inorm)
         fr(inorm) = (qr(k,inorm)**2)/qr(k,1) + pr
         fr(itan) = qr(k,inorm)*qr(k,itan)/qr(k,1)
         fr(4) = (qr(k,4) + pr) *  qr(k,inorm) /  qr(k,1)
         endif

         rx(k,:) = 0.5d0*(fl(:)+fr(:)) + 0.5d0*s*(uconl(:)-uconr(:))

  10     continue






      return
      end
