c
c -------------------------------------------------------------
c
      subroutine vrm(qr,ql,rx,ixmin,ixmax,iflip,msize)

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
      dimension dir(2), dn(2)
      data      itno/2/, small/1.d-6/
c
      pi = 3.14159265358979d0

      if (ixmax .gt. m1size) then
         write(6,*)" need to increase m1size in vrm to ", ixmax
         stop
      endif


      call getdir(dir)

c      iflip = 1 indicates the normal velocity is in
c                qr(k,2) and the tangential in qr(k,3)
c            = 0 indicates the reverse.

      inorm = 3 - iflip
      itan  = 2 + iflip
      dn = 0.d0
      if(iflip .eq. 1) then
        dn(1) = 1.d0
      elseif(iflip .eq. 0) then
        dn(2) = 1.d0
      endif

      dot = dir(1)*dn(1) + dir(2)*dn(2)

      rx = 0.d0
      do 10 k = ixmin, ixmax
         rx(k,:) = 0.5d0*dot * (ql(k,1)+qr(k,1))
     .            +0.5d0*abs(dot)*( ql(k,1)-qr(k,1) )
         rx(k,2:4) = 0.d0

  10     continue






      return
      end


      subroutine getdir(dir)
      implicit double precision (a-h,o-z)
      dimension dir(2)
      pi = 3.14159265358979d0
!      dangle = pi/3.d0 + pi/2.d0
!      dangle = pi/3.d0
      dangle = pi/3.d0
      dir(1) = cos(dangle)
      dir(2) = sin(dangle)
      return
      end
