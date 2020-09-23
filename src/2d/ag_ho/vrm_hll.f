c
c -------------------------------------------------------------
c
      subroutine vrm(qr,ql,rx,ixmin,ixmax,iflip,msize)
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      parameter ( m1size = 4033)
      dimension qr(msize,4),ql(msize,4)
      dimension rx(msize,4), rlx(4), rrx(4), stateR(4), stateL(4)
      dimension rhor(m1size),ur(m1size),pr(m1size),utr(m1size)
      dimension rhol(m1size),ul(m1size),pl(m1size),utl(m1size)

      data      itno/2/, small/1.d-6/
c
c
c     # HLL vector Riemann solver
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
        aL  = dsqrt(gamma*pl(k)/rhol(k))
        aR  = dsqrt(gamma*pr(k)/rhor(k))

c       sL = ul(k) - aL; 
c       sR = ur(k) + aR; 
        sL = min((ul(k) - aL), (ur(k) - aR)) ! fancier option  
        sR = max((ul(k) + aL), (ur(k) + aR))

        if (0. .lt. sL) then   ! use left flux

         rx(k,1)     = rhol(k)*ul(k)
         rx(k,inorm) = rhol(k)*ul(k)*ul(k) + pl(k)
         rx(k,itan)  = rhol(k)*ul(k)*utl(k)
         etot        = .5*rhol(k)*(ul(k)*ul(k)+utl(k)*utl(k)) + 
     .                  pl(k)/gamma1
         rx(k,4)     = ul(k)*(etot + pl(k))

        else if (0. .gt. sR) then   ! flux is right flux

         rx(k,1)     = rhor(k)*ur(k)
         rx(k,inorm) = rhor(k)*ur(k)*ur(k) + pr(k)
         rx(k,itan)  = rhor(k)*ur(k)*utr(k)
         etot        = .5*rhor(k)*(ur(k)*ur(k)+utr(k)*utr(k)) + 
     .                  pr(k)/gamma1
         rx(k,4)     = ur(k)*(etot + pr(k))

        else
        ! middle flux
         rlx(1)     = rhol(k)*ul(k)
         rlx(inorm) = rhol(k)*ul(k)*ul(k) + pl(k)
         rlx(itan)  = rhol(k)*ul(k)*utl(k)
         etotl        = .5*rhol(k)*(ul(k)*ul(k)+utl(k)*utl(k)) +
     .                  pl(k)/gamma1
         rlx(4)     = ul(k)*(etotl + pl(k))

         rrx(1)     = rhor(k)*ur(k)
         rrx(inorm) = rhor(k)*ur(k)*ur(k) + pr(k)
         rrx(itan)  = rhor(k)*ur(k)*utr(k)
         etotr        = .5*rhor(k)*(ur(k)*ur(k)+utr(k)*utr(k)) +
     .                  pr(k)/gamma1
         rrx(4)     = ur(k)*(etotr + pr(k))

         stateL(1)     = rhol(k)
         stateL(inorm) = rhol(k)*ul(k)
         stateL(itan)  = rhol(k)*utl(k)
         stateL(4)     = etotl

         stateR(1)     = rhor(k)
         stateR(inorm) = rhor(k)*ur(k)
         stateR(itan)  = rhor(k)*utr(k)
         stateR(4)     = etotr

         rx(k,:)=(sR*rlx(:)-sL*rrx(:)
     .                       +sL*sR*( stateR(:)-stateL(:) ) )/ (sR-sL)
        endif
 20   continue

c
       return
       end
