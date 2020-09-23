c
c
c
c
      function fbody(x,y)
      implicit double precision (a-h,o-z)
      dimension xuc(10),yuc(10),xlc(10),ylc(10)
      dimension xupw(500),yupw(500)
      logical first
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
       common /skewll/ y0skew,dely
       common /exwavey/ np,ipad,R1,R2,rp,dtshift
       data first /.true./
       pi = 3.14159265358979d0
c
c  negative inside the body (exterior to the domain), positive otherwise.
c
c     fbody = 1.
c     return
c
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     .22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32)iprob
c
    1 continue
c     # quarter cylinder
      r = .5001d0
      fbody = (x-1.0d0)**2 + (y-0.d0)**2 - r**2
      if (fbody.eq.0.d0) fbody = 1d-5
      return
c
    2 continue
      fbody = 1.d0
!      if (y .lt. .00001) fbody = -1.   ! test horizontal wall

      return

!   22 continue
!      fbody = 1.d0
!      if (y .lt. (0.1d0-.0001) ) fbody = -1.   ! test horizontal wall
!      return

   22 continue
      fbody = 1.d0
      if (y .lt. .0001 ) fbody = -1.   ! test horizontal wall
      return

!   23 continue
!      fbody = 1.d0
!      dyi = 0.00999d0
!      if ( (y .lt. dyi .and. x > 0.5d0) .or.
!     .     (y .lt. -(x-0.5d0)/sqrt(3.d0) + dyi) .and. x <= 0.5d0) then
!       fbody = -1.   ! test horizontal wall
!      endif
!      return

   23 continue
      fbody = 1.d0
      if ( (y .lt. 0.00001d0 .and. x > 0.5d0) .or.
     .     (y .lt. -(x-0.5d0)/sqrt(3.d0)) .and. x <= 0.5d0) then
       fbody = -1.   ! test horizontal wall
      endif
      return



   24 continue
c     # channel
      if ((y .ge. .1d0*x+.11d0) .and. (y .le. .1d0*x+.81d0)) then
        fbody = 1.0d0
      else
        fbody = -1.d0
      endif

      return
   25 continue

      xint = 2.d0 - 0.00001d0
      yint = 0.d0

c      dangle = 0.321750554d0
      dangle = pi/5.d0
      f1 = (y-yint) - tan(dangle) * (x-xint)
      f2 = (y-xint) - tan(dangle) * (x-yint)

      if ( f1 > 0.d0 .and. f2 < 0.d0 ) then
          fbody = 1.d0
      else
          fbody = -1.d0
      endif


      return
   26 continue
c     # bunch of overlapping circles
      fbody = 1.d0


      xstart = 1.d0
      xend = 3.d0
      numcircles = 4
      dx = (xend-xstart)/(numcircles-1)
      xc = xstart
      yc = 2.d0
      r = 0.5d0

      do n=1,numcircles
        if( (x-xc)**2 /(r/1.25)**2+(y-yc)**2 / r**2- 1. <0) then
            fbody = -1.d0
        endif
        xc = xc+dx
      end do

      return

  27  continue
      fbody = 1.d0
      R1 = 0.751d0
      R2 = 1.251d0
      if(   (x-1.5d0)**2 +(y-1.5d0)**2  < R1**2
     . .or. (x-1.5d0)**2 +(y-1.5d0)**2  > R2**2 ) then
            fbody = -1.d0
      endif




      return

 28   continue
      fbody = 1.d0
      R1 = 1.00d0
      R2 = 1.384d0
      if(   (x**2 +y**2  < R1**2)
     . .or. (x**2 +y**2  > R2**2) ) then
            fbody = -1.d0
      endif
      return

 29   continue
      fbody = 1.d0
      ! best so far
      R1 = 0.759d0
      R2 = 1.251d0
      rp = 0.25d0
      np = 10
!      dtshift = .0d0

      ! np = 10 works for 25,50,100 and has a nice overlapping example
      ! when np = 100
!      R1 = 0.751d0
!      R2 = 1.251d0
!      rp = 0.25d0

!      np = 7
!      R1 = 0.751d0
!      R2 = 1.251d0
!      rp = 0.25d0

      theta = atan2(y-1.5d0, x-1.5d0)+dtshift
      R = dsqrt( (x-1.5d0)**2 + (y-1.5d0)**2 )
      if(  (x-1.5d0)**2 +(y-1.5d0)**2  < (R1+rp*sin(np*theta))**2
     ..or. (x-1.5d0)**2+(y-1.5d0)**2>(R2+rp*sin(np*theta))**2 ) then
            fbody = -1.d0
      endif

      return

 30   continue
c     #NACA 0012 airfoil
      fbody = 1.d0

      ddt = 1.d-1


      if( x < 0.25 .or. x > 0.75d0) then
      return
      endif

      xp = (x-0.25)/0.5d0
      yp = (y-0.25)/0.5d0

      daf = 5.d0 * ddt*(0.2969 * dsqrt(xp) -0.1260 * xp
     . - 0.3516 *xp**2
     .              +0.2843 * xp**3 -0.1036*xp**4)
      if( dabs(yp) < daf) then
        fbody = -1.d0
      endif

      return
c
c again to match rotbox on zero origin domain
c
    3 continue
c     # nozzle
      fbody = 1.d0
      pi2 = 2.d0*datan(1.d0)
c      fbody = (w1-w2)*(datan(alf*(x-.5d0))/pi2)**2 + w2 - abs(y-.5d0)
      return
c
    4 continue
      fbody = 1.d0
!c     # piecewise linear channel (e.g. ramp)
!c     # (xlc(1),ylc(1)), ..., (xlc(nlc),ylc(nlc)) defines lower wall
!c     # (xuc(1),yuc(1)), ..., (xuc(nuc),yuc(nuc)) defines upper wall
!c
!c     # 40 degree ramp:
!      xlc(1) = 0.d0
!      ylc(1) = .8142d0
!      xlc(2) = .85d0
!      ylc(2) = .101d0
!      xlc(3) = 1.d0
!      ylc(3) = .101d0
!      nlc = 3
!      nuc = 0
!c
!c     # 20 degree ramp aligned with grid:
!      xlc(1) = 0.d0
!      ylc(1) = .04649702d0
!      xlc(2) = .1d0
!      ylc(2) = .0101d0
!      xlc(3) = 1.d0
!      ylc(3) = .0101d0
!      nlc = 3
!      nuc = 0
!c
!c     # 20 degree ramp:
!      xlc(1) = 0.d0
!      ylc(1) = .0101d0
!      xlc(2) = .1d0
!      ylc(2) = .0101d0
!      xlc(3) = 1.d0
!      ylc(3) = .3376732d0
!      nlc = 3
!      nuc = 0
!c
!      fbody = 1.
!      if (nuc.eq.0) go to 43
!      do 41 i=1,nuc
!   41    if (x.lt.xuc(i)) go to 42
!   42 if (i.eq.1) i=2
!      if (i.gt.nuc) i = nuc
!      ybdry = yuc(i-1)+(x-xuc(i-1))*(yuc(i)-yuc(i-1))/(xuc(i)-xuc(i-1))
!      if (y.gt.ybdry) fbody = -1.
!   43 if (nlc.eq.0. .or. fbody.eq.-1.) go to 46
!      do 44 i=1,nlc
!   44    if (x.lt.xlc(i)) go to 45
!   45 if (i.eq.1) i = 2
!      if (i.gt.nlc) i = nlc
!      ybdry = ylc(i-1)+(x-xlc(i-1))*(ylc(i)-ylc(i-1))/(xlc(i)-xlc(i-1))
!      if (y.lt.ybdry) fbody = -1.
!   46 continue
!c      if (r.eq.0.) return
!c      fbody = fbody * ((x-cx1)**2 + (y-cy1)**2 - r**2)
      return
c
    5 continue
      fbody = 1.d0
c     # channel
c      f1 = (x-1.0d0)*beta - (y-cy2)*alf
c      f2 = -x*beta + (y - cy2 + (beta+w1)/alf)*alf
c      f3 = (x-cx1)*alf + (y-cy1)*beta
c      if (f1.gt.0d0 .and. f2.gt.0d0 .and. f3.gt.0d0) then
c          fbody = 1.d0
c        else
c          fbody = -1.d0
c        endif
c      if (r.eq.0.) return
c      fbody = fbody * ((x-.5d0)**2 + (y-.5d0)**2 - r**2)
      return
c
    6 continue
c     #NACA 0012 airfoil
c     theta = (8.d0 * datan(1.d0)*4.d0)/180.d0
      theta = 0.d0
      alf = dcos(theta)
      beta = dsin(theta)
c trailing edge location (cx1,cy1). w1 = length of chord.
c     cx1 = 03.050d0
c     cy1 = 03.0d0
      cx1 = 06.050d0
      cy1 = 06.0d0
      w1 = 0.1d0
      xi = 1.0089d0*(w1 + (x-cx1)*alf + (y-cy1)*beta)/w1
      eta = 1.0089d0*(-(x-cx1)*beta + (y-cy1)*alf)/w1
      fbody = dabs(eta) - .6d0*(0.2969d0*dsqrt(dabs(xi)) 
     &        + (((-0.1015d0*xi
     &        + .2843d0)*xi - .3516d0)*xi - .126d0)*xi) 
     &        + 100.d0*(1.-dsign(1.d0,xi))
      return
c
    7 continue
c     # multiple cylinders
      fbody = 1.
c     if ((x-.4)**2+(y-.507)**2 - (.103)**2 .lt. 0.) fbody = -1.
c     if ((x-.5)**2+(y-.227)**2 - (.096)**2 .lt. 0.) fbody = -1.
c     if ((x-.5)**2+(y-.787)**2 - (.096)**2 .lt. 0.) fbody = -1.
      if ((x-.3)**2+(y-.30)**2 - (.12)**2 .lt. 0.) fbody = -1.
      if ((x-.4)**2+(y-.70)**2 - (.15)**2 .lt. 0.) fbody = -1.
c     if (y.lt. .03) fbody = -1.
      if (y.gt. .97) fbody = -1.
      return
c
    8 continue
c     # wall at 45 degrees
c     if (y.gt. 0.95d0-x) then
c       fbody = 1.0d0
c     else
c       fbody = -1.0d0
c     endif
c     # wall at small angle
      y0skew = 0.04999d0
c     dely = 0.2d0
      dely = 0.0d0
      if (y.gt. y0skew + dely*x) then
        fbody = 1.0d0
      else
        fbody = -1.0d0
      endif
      return
c
    9 continue
c     # kumar's nozzle
      fbody = 1.d0
c     if (y.gt. .29999d0) fbody = -1.d0
      if (y.gt. .39999d0) fbody = -1.d0
      if (x.lt. .1743d0 .or. x.gt. .68252d0) return
      if (y.gt. .10001d0) return                  
      if (x.lt. .35d0) then
          if (y.gt. .1d0-.1659d0*(x-.1743d0)) fbody = -1.d0
        else
          if (y.gt. .07085d0+.08766d0*(x-.35d0)) fbody = -1.d0
        endif
      return
c
   10 continue
c     # 2 cylinder run
      fbody1 = (x-.43)**2 + (y-.3)**2 - .0201
      fbody2 = (x-.5)**2 + (y-.75)**2 - .0201
      if ((fbody1 .gt. 0.d0) .and. (fbody2 .gt. 0.d0)) fbody = 1.d0
      if ((fbody1 .gt. 0.d0) .and. (fbody2 .lt. 0.d0)) fbody = fbody2
      if ((fbody1 .lt. 0.d0) .and. (fbody2 .gt. 0.d0)) fbody = fbody1
      return
c
   11 continue
      fbody = 1.d0
c     # single cylinder
      temp = (x-0.5d0)**2 + (y-0.5d0)**2 - .15d0**2
      if(temp < 0.d0) fbody = -1.d0
c     fbody = (x-.5)**2 + (y-.5)**2 - .037
c     fbody = (x-.5)**2 + (y-.5)**2 - .143**2
c     fbody = (x-.5)**2 + (y-.5)**2 - .235**2
      return
c
   12 continue
c     # 90 degree bend
      x0 = .1d0
      y0 = .9d0
      r1 = .59d0
      r2 = .89d0
      fbody = -1.
      if (x.le.x0 .and. y.le.y0-r1 .and. y.ge.y0-r2) fbody = 1.
      if (y.ge.y0 .and. x.ge.x0+r1 .and. x.le.x0+r2) fbody = 1.
      if (x.gt.x0 .and. y.lt.y0) then
         r = dsqrt((x-x0)**2 + (y-y0)**2)
         if (r.ge.r1 .and. r.le.r2) fbody = 1.
      endif
      return

   13 continue
c     # hemicircle
c     r = .25001d0
c     fbody = (x-.5d0)**2 + (y-0.d0)**2 - r**2
c     if (fbody.eq.0.d0) fbody = 1d-5
      r = .5001d0
      fbody = (x-0d0)**2 + (y-0.2d0)**2 - r**2
      if (fbody.eq.0.d0) fbody = 1d-5
      return
c
   14 continue
c    # no body
      fbody = 1.d0
      return
c
   15 continue
c     if (y .lt. .39d0 .and. y .gt. .00000000001d0) then
c     if (y .gt. .025001d0) then
      alpha = .11d0
      alpha = .333333333d0
      if (y .lt. .1749d0+alpha*x .and. y .gt. .025000d0+alpha*x) then
       fbody = 1.0d0
      else
       fbody = -1.d0
      endif
      return
c
   16 continue
      if ((y .ge. .1d0*x+.11d0) .and. (y .le. .1d0*x+.81d0)) then
c     if ((y .ge. -.13999999+x) .and. (y .le. x+.1399999999)) then
c     if ((y .ge. -.2+x) .and. (y .le. x+.2)) then
       fbody = 1.0d0
      else
       fbody = -1.d0
      endif
      return
c
   17 continue
c     # single cylinder for swirling flow
      fbody = (x-2.0d0)**2 + (y-2.0d0)**2 - 1.000567d0**2
c     fbody = (x-.5)**2 + (y-.5)**2 - .1923538406**2
      return
c 
   18 continue
c     # simple wave around curved wall, from Whitham, p. 204
c
c     # lower wall:
c
      if ((x.le.0.1d0 .and. y.ge.0.302d0) .or. (x.gt.0.1d0 .and.
     &     x.le.0.7d0 .and. y.ge.0.302d0-0.3d0*(x-0.1d0)**2).or.(x.gt.0.7d0
     &     .and. y.ge..194d0-0.36d0*(x-0.7d0))) then
         fbody = 1.d0
        else
         fbody = -1.d0
        endif
c     return
c
c     # upper wall:
c
      if (first) then
c     #  generate upper boundary points:
         nupw = 300
         call p18(xupw,yupw,nupw)
         first = .false.
      endif
c     
      if (fbody.gt. 0.) then
         do 1810 i=1,nupw
            if (x .lt. xupw(i)) go to 1820
 1810    continue
         write(6,*) 'error in fbody...  upper wall'
         write(6,*) '   xupw doesn"t go far enough.'
         write(6,*) '   try increasing h in p18'
 1820    continue
         if (i.eq.1) i = 2
c     # use linear interp between points i-1 and i to get y on wall:
         yw = yupw(i-1) + (x-xupw(i-1)) * (yupw(i)-yupw(i-1))/
     &        (xupw(i)-xupw(i-1))
         if (y.gt.yw) fbody = -1.d0
      endif
      return
c
 19   continue
c
c  supersonic vortex  in circular channel.  inner radius r1, outer r2.
c  center of circles at the origin
c
      r1 = 1.d0
c     r2 = 1.50
      r2 = 1.3820
      fbody = -1.d0
c     r = dsqrt(x**2 + y**2)
      rsq = (x**2 + y**2)
c     if (x .lt. 0.d0) then
c        if (y .gt. r1 .and. y .lt. r2) fbody = 1.0d0
c      else if (y .lt. 0.d0) then
c        if (x .gt. r1 .and. x .lt. r2) fbody = 1.0d0
c      else if (r.ge.r1 .and. r.le.r2) then
c        fbody = 1.d0
c      endif
c     if (r .ge. r1 .and. r .le. r2) fbody = 1.d0
      if (rsq .ge. r1*r1 .and. rsq .le. r2*r2) fbody = 1.d0

      return

 20   continue
c
c  annulus, as in 1d irregular grid paper. need to shift since 
c  code assumes lower left corner of domain is origin
c  but paper has center of circle at the origin
c
      xshift = 1.50
      yshift = 1.50  ! in case want to break symmetry
      r1 = .75d0
      r2 = 1.25d0
      fbody = -1.d0
      xuse = x - xshift
      yuse = y - yshift
      rsq = (xuse**2 + yuse**2)
      if (rsq .ge. r1*r1 .and. rsq .le. r2*r2) fbody = 1.d0
c
c stick problem 16 into prob 20 but with fixed velocity and 1 variable
c      if ((y .ge. .1*x+.11) .and. (y .le. .1*x+.81)) then
c       fbody = 1.0d0
c      else
c       fbody = -1.d0
c      endif

      return


 21   continue
c     # 30 degree ramp: starts ta x=.5 has BAD CELL BLOWING UP
      xwall = 1.d0/6.d0
      ywall = (3.-xwall)*tan(30*pi/180)
      if (x .lt. xwall) then 
         if (y .gt.0.d0) then
             fbody = 1.d0
         else
             fbody = -1.d0
         endif
      else 
        slope = (ywall - 0.d0)/(3.00d0-xwall)
        ybndry = (x-xwall)*slope
        if (y .ge. ybndry) then
           fbody = 1.d0
        else
           fbody = -1.d0
        endif
      endif
      return

 31   continue
      fbody = 1.d0
!      R1 = 1.00d0
!      R2 = 1.384d0
!      if(   (x**2 +y**2  < R1**2)
!     . .or. (x**2 +y**2  > R2**2) ) then
!            fbody = -1.d0
!      endif
      return
 32   continue
      fbody = 1.d0
      ! best so far
      R1 = 0.759d0
      R2 = 1.251d0
      rp = 0.25d0
      np = 10

      theta = atan2(y,x)
      R = dsqrt( x**2 + y**2 )
      if(   x**2 + y**2  < (R1+rp*sin(np*theta))**2
     . .or. x**2 + y**2  > (R2+rp*sin(np*theta))**2
     . .or. x+y < 0.25d0) then
            fbody = -1.d0
      endif

      return



      return

      end
