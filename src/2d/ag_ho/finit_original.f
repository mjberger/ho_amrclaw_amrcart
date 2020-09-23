c
c -------------------------------------------------------------
c
        subroutine finit(val,nvar,maxip1,maxjp1,
     1               corn1,corn2,hx,hy,irr,mitot,mjtot,lstgrd)
c
       implicit double precision (a-h,o-z)
       dimension val(maxip1,maxjp1,nvar),irr(mitot,mjtot)
       dimension r(4), s(4), ww(4)
       common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                 ismp,gradThreshold
       common /skewll/ y0skew,dely
      include "cirr.i"
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
c
       pi = 3.14159265358979d0

       r = (/-5.7735026918962562d-01, 5.7735026918962595d-01,
     .        5.7735026918962562d-01, -5.7735026918962595d-01 /)
       s = (/ -5.7735026918962595d-01, -5.7735026918962562d-01,
     .         5.7735026918962595d-01, 5.7735026918962562d-01 /)
       ww = (/ 1.0000000000000000d+00, 1.0000000000000000d+00,
     .         1.0000000000000000d+00, 1.0000000000000000d+00 /)




c
c      # default values:  (Mach 2.81 shock at x=.27)
c      # sound speeds : left: 1.57, right: 1
c      # shock speed: 2.81
c      # ul + cl = 3.63
       gamma1 = .4d0
       sloc = 0.27d0
       rhol = 5.1432d0
       ul = 2.04511d0
       vl = 0.d0
c      #  for 45 degree channel
c      ul = 2.04511d0/dsqrt(2.d0)
c      vl = 2.04511d0/dsqrt(2.d0)
c
       pl = 9.04545d0
c
       rhor = 1.4d0
       ur = 0.0d0
c       ur = 1.d0
       vr = 0.0d0
       pr = 1.0d0
c
       go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
     .        17,18,19,20,21,22,23,24,25,26) iprob
    1  continue
c      # Mach 8 flow into quarter cylinder:
       sloc = 10000.d0
       rhol = 1.4d0
       ul = 8.0d0
       vl = 0.0d0
       pl = 1.0d0

       go to 40
    2  continue
       alpha = 0.d0
       alf = alpha*pi/180.d0
c      rnx = -sin(alf)
c      rny = cos(alf)
c      tx = rny
c      ty = -rnx

       Rmach = .2

       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 93 i = 1, maxip1
       do 93 j = 1, maxjp1
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          kirr = irr(iadj,jadj)
          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
             xcen = xcirr(kirr)
             ycen = ycirr(kirr)
          else
             xcen = corn1 + (dfloat(i)-1.5d0)*hx 
             ycen = corn2 + (dfloat(j)-1.5d0)*hy 
             kirr = lstgrd
          endif

c  test code using linear input 
c             rn = -sin(alf)*xcen + cos(alf)*ycen
c             t  =  cos(alf)*xcen + sin(alf)*ycen
c             veln = 0.
c             velt = t
c             ut = cos(alf)*velt - sin(alf)*veln
c             vt = sin(alf)*velt + cos(alf)*veln

          dangle = atan(  ycen / (xcen - 1.d0/6.d0) )
          if(dangle > 1.0471975512d0 .or. xcen < 1.d0/6.d0) then
              ut =  dcos(pi/6.d0)*8.25d0
              vt = -dsin(pi/6.d0)*8.25d0
              rhot = 8.d0
              pt = 116.5d0
          else
              ut = 0.d0
              vt = 0.d0
              rhot = gamma
              pt = 1.d0/(gamma-1.d0)
          endif

          val(i,j,1) = rhot
          val(i,j,2) = rhot*ut
          val(i,j,3) = rhot*vt
          val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)
 93    continue

       return
  22  continue
       alpha = 0.d0
       alf = alpha*pi/180.d0
c      rnx = -sin(alf)
c      rny = cos(alf)
c      tx = rny
c      ty = -rnx

       Rmach = .2

       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 913 i = 1, maxip1
       do 913 j = 1, maxjp1
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          kirr = irr(iadj,jadj)
          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
             xcen = xcirr(kirr)
             ycen = ycirr(kirr)
          else
             xcen = corn1 + (dfloat(i)-1.5d0)*hx
             ycen = corn2 + (dfloat(j)-1.5d0)*hy
             kirr = lstgrd
          endif


          if(xcen < 1.d0) then
              ut = 5.d0/4.d0
              vt = 0.d0
              rhot = 56.d0/15.d0
              pt = 9.d0/2.d0
          else
              ut = 0.d0
              vt = 0.d0
              rhot = gamma
              pt = 1.d0/(gamma-1.d0)
          endif

          val(i,j,1) = rhot
          val(i,j,2) = rhot*ut
          val(i,j,3) = rhot*vt
          val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)
 913    continue

       return
  23  continue
       alpha = 0.d0
       alf = alpha*pi/180.d0

       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 923 i = 1, maxip1
       do 923 j = 1, maxjp1
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          kirr = irr(iadj,jadj)
          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
             xcen = xcirr(kirr)
             ycen = ycirr(kirr)
          else
             xcen = corn1 + (dfloat(i)-1.5d0)*hx
             ycen = corn2 + (dfloat(j)-1.5d0)*hy
             kirr = lstgrd
          endif



             dspos = 0.5d0 + ycen/dsqrt(3.d0)


             if(xcen < dspos) then
                 ut =  dcos(pi/6.d0)*8.25
                 vt = -dsin(pi/6.d0)*8.25
                 rhot = 8.
                 pt = 116.5
             else
                 ut = 0.d0
                 vt = 0.d0
                 rhot = gamma
                 pt = 1.d0
             endif

          val(i,j,1) = rhot
          val(i,j,2) = rhot*ut
          val(i,j,3) = rhot*vt
          val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)
 923    continue

       return
  24  continue

       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 9123 i = 1, maxip1
       do 9123 j = 1, maxjp1
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          kirr = irr(iadj,jadj)
          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
             xcen = xcirr(kirr)
             ycen = ycirr(kirr)
          else
             xcen = corn1 + (dfloat(i)-1.5d0)*hx
             ycen = corn2 + (dfloat(j)-1.5d0)*hy
             kirr = lstgrd
          endif

c          dangle = pi/6.d0
c          dangle = pi
c          dangle = pi/4.d0
          dangle = pi/5.d0

          dnx = cos(dangle)
          dny = sin(dangle)
          xi = dnx*(xcen+0.25d0) + dny*(ycen+0.25d0)

c          if(i .eq. 5 .and. j .eq. 4) then
c            print *,"here"
c          endif

          rhol = 8.d0
          ul = 8.25d0
          pl = 116.5d0
          cl = dsqrt( gamma * pl / rhol )

       solumag = ( (gamma - 1.d0)*ul + 2.d0*(cl + xi) ) / (gamma + 1.d0)
          ut   = solumag * dnx
          vt   = solumag * dny
          rhot =( ( rhol**gamma * (solumag - xi)**2.d0 ) / (gamma*pl) )
     .**( 1.d0/(gamma-1.d0) )
          pt   = ( pl/(rhol**gamma) ) * rhot**gamma

          val(i,j,1) = rhot
          val(i,j,2) = rhot*ut
          val(i,j,3) = rhot*vt
          val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)
c            if(iadj .eq. 53 .and. jadj .eq. 4) then
c            print *, iadj,jadj
c            endif



 9123 continue

       return


  25  continue

       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 8123 i = 1, maxip1
       do 8123 j = 1, maxjp1
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          kirr = irr(iadj,jadj)
          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
             xcen = xcirr(kirr)
             ycen = ycirr(kirr)
          else
             xcen = corn1 + (dfloat(i)-1.5d0)*hx
             ycen = corn2 + (dfloat(j)-1.5d0)*hy
             kirr = lstgrd
          endif



          pi = 3.14159265358979d0
          dangle = pi/5.d0
          dnx = cos(dangle)
          dny = sin(dangle)


          rhot = sin(dnx*xcen + dny*ycen)
!           rhot = dnx*xcen + dny*ycen




c          rhot = dnx*xcen + dny*ycen

          ut = 0.d0
          vt = 0.d0
          pt = 1.d0


          val(i,j,1) = rhot
          val(i,j,2) = rhot*ut
          val(i,j,3) = rhot*vt
          val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)


 8123 continue




!       xlow = corn1 - lwidth*hx
!       ylow = corn2 - lwidth*hy
!       do 8123 i = 1, maxip1
!       do 8123 j = 1, maxjp1
!          iadj = i + lwidth - 1
!          jadj = j + lwidth - 1
!          kirr = irr(iadj,jadj)
!          if (kirr .ne. lstgrd) then ! only full cells for now
!             continue
!          endif
!
!
!!
!!          pi = 3.14159265358979d0
!!          dangle = pi/5.d0
!!          dnx = cos(dangle)
!!          dny = sin(dangle)
!!
!!          do 1212 q = 1, 4
!!            rhot = sin(dnx*xcen + dny*ycen)
!!1212      continue
!
!          ut = 0.d0
!          vt = 0.d0
!          pt = 1.d0
!
!
!          val(i,j,1) = rhot
!          val(i,j,2) = rhot*ut
!          val(i,j,3) = rhot*vt
!          val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)
!
!
! 8123 continue




      return


  26  continue
       alpha = 0.d0
       alf = alpha*pi/180.d0

       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 44 i = 1, maxip1
       do 44 j = 1, maxjp1
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          kirr = irr(iadj,jadj)
          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
             xcen = xcirr(kirr)
             ycen = ycirr(kirr)
          else
             xcen = corn1 + (dfloat(i)-1.5d0)*hx
             ycen = corn2 + (dfloat(j)-1.5d0)*hy
             kirr = lstgrd
          endif



             dspos = 0.5d0


             if(xcen < dspos) then
                 ut = 8.25d0
                 vt = 0.d0
                 rhot = 8.
                 pt = 116.5
             else
                 ut = 0.d0
                 vt = 0.d0
                 rhot = gamma
                 pt = 1.d0
             endif

          val(i,j,1) = rhot
          val(i,j,2) = rhot*ut
          val(i,j,3) = rhot*vt
          val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)
 44    continue

      return








    3  continue
       go to 40
    4  continue
c

c      # Mach 3 flow into ramp:
       sloc = 10000.d0
       rhol = 1.4d0
       ul = 3.0d0
       vl = 0.0d0
       pl = 1.0d0
c
c      # rotated by 20 degrees, for ramp aligned with grid:
c      ul = 2.8819078d0
c      vl = -1.02606d0
       go to 40
    5  continue
       go to 40
    6  continue
       sloc = 10000.d0
       rhol = 1.0d0
       pl = 1.0d0
       rmach = .8d0
c      rmach = .302
c      alpha = 1.25* datan(1.0d0)/45.
c      alpha = 0.
       c = dsqrt(gamma*pl/rhol)
       ul = rmach * c * dcos(alpha)
       vl = rmach * c * dsin(alpha)
       go to 40
    7  continue
c      sloc = .42
       sloc = .15
       go to 40
    8  continue
       sloc = -1.
       ur = 1.
       vr = -1.
c      go to 40
c
c      # flow parallel to slightly skewed wall:
       alf = 1.d0/dsqrt(1.d0 + dely**2)
       beta = dely/dsqrt(1.d0 + dely**2)
       rslope = dely
       b0 = .2*rslope+.2/rslope+y0skew
       do 45 i = 1, maxip1
       do 45 j = 1, maxjp1
             x = corn1 + (i-2.5)*hx
             y = corn2 + (j-2.5)*hy
             if (y .gt. -1.d0/rslope*x+b0) then
                rho = 1.4d0
                utot = 0.d0
                p = 1.d0
             else
                rho = 5.1432d0
                utot = 2.04511
                p = 9.04545d0
             endif
             z = alf*x + beta*y
c            utot = 2.d0*(1.d0-z)**2
c            utot = 1.d0
             val(i,j,1) = rho
             val(i,j,2) = alf*utot*val(i,j,1)
             val(i,j,3) = beta*utot*val(i,j,1)
             val(i,j,4) = p/gamma1 + 0.5d0*(val(i,j,2)**2 + 
     &                     val(i,j,3)**2)/val(i,j,1)
   45        continue
       return
       go to 40
c
    9  continue
       sloc = -10000.
       ur = 1.9d0
       go to 40
c
   10  continue
       go to 40
   11  continue
c      # one  cell before cylinder
c  change initial conditions for Mach 2 shock as in Helzel et al
c       pl = pr * 4.5
c       rhol = rhor * 2.6666666
c       cl = sqrt(gamma*pl/rhol)
c       ul = 2. * cl
c      sloc = .299d0
c       rhol =  56.d0/15.d0
c       ul =  4.d0/5.d0
c       pl = 9.d0/2.d0
c       sloc = .199d0
       rhol =  8.d0
       ul =  8.25d0
       vl = 0.d0
       pl = 116.5d0
       sloc = 1.d0

       go to 40
   12  continue
c      # Mach 2.2 shock
       sloc = .10e0
       rhol = 4.13171d0
       ul = 1.45455d0
       pl = 5.48d0
       go to 40
   13  continue
       sloc = .23d0
       go to 40
   14  continue
       sloc =35.2d0
       ul = .20d0
       vl = 0.d0
       rhol = 1.d0 
c      pl = 1.d0   change normalization to match cart3d
       pl = 1.d0/1.4d0   ! where sound speed=1, p = 1/gamma
       rhor = .426d0
       ur = 1.427d0
       pr = .303d0
       go to 40
c
c     ----------------------------
c      # expanding shock:
       do 225 i = 1, maxip1
       do 225 j = 1, maxjp1
             x = corn1 + (i-1.5d0)*hx - 0.5d0
             y = corn2 + (j-1.5d0)*hy - 0.5d0
             r2 = x**2 + y**2
             rho = 1.4d0
             u = 0.d0
             v = 0.d0
             p = 1.d0
             if (r2.lt..04) then
                rho = 1.4d0
                p = 100.d0
                endif
             val(i,j,1) = rho
             val(i,j,2) = rho*u
             val(i,j,3) = rho*v
             val(i,j,4) = p/gamma1 + 0.5d0*rho*(u**2 + v**2)
  225        continue
       return
c     ----------------------------
c
       go to 40
   15  continue
       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 193 i = 1, maxip1
          do 193 j = 1, maxjp1
             iadj = i + lwidth - 1
             jadj = j + lwidth - 1
             kirr = irr(iadj,jadj)
             if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
                xcen = xcirr(kirr)
                ycen = ycirr(kirr)
             else
                xcen = corn1 + (dfloat(i)-1.5d0)*hx 
                ycen = corn2 + (dfloat(j)-1.5d0)*hy 
                kirr = lstgrd
                call makep(poly(1,1,kirr),iadj,jadj,xlow,ylow,hx,hy)
             endif
             call p15tru(xcen,ycen,rhot,ut,vt,pt,poly(1,1,kirr))
             val(i,j,1) = rhot
             val(i,j,2) = rhot*ut
             val(i,j,3) = rhot*vt
             val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)
 193      continue
          return
          sloc =1.2d0
          rhol = 1.d0
          ul = 2.d0
          vl = 0.d0
          pl = 1.0d0
c     rhor = .426d0
          go to 40
c
c rotated channel
 16    continue
       pl = 1.0d0
       ul = .1d0
       vl = .01d0
c      ul = .0d0
c      vl = .00d0
       do 1610 i = 1, maxip1
       do 1610 j = 1, maxjp1
          iuse = i + lwidth - 1
          juse = j + lwidth - 1
          kuse = irr(iuse,juse)
          if (kuse .ne. -1 .and. kuse .ne. lstgrd) then
             xcen = xcirr(kuse)
             ycen = ycirr(kuse)
          else
             xcen = (i-1.5d0)*hx + corn1
             ycen = (j-1.5d0)*hy + corn2
          endif
          rhol = ycen - .1d0*xcen + .5d0
c         rhol = 1.d0
          val(i,j,1) = rhol
          val(i,j,2) = ul*rhol
          val(i,j,3) = vl*rhol
          val(i,j,4) = pl/gamma1 + 0.5d0*rhol*(ul**2 + vl**2)
 1610     continue
       return
c
 17    continue
c     # flow around a cylinder
       go to 177
       rho = 30.d0
       do 125 i = 1, maxip1
       do 125 j = 1, maxjp1
             k = irr(i+lwidth-1,j+lwidth-1)
             if (k .eq. lstgrd) then
                x = corn1 + (i-1.5d0)*hx - 0.5d0
                y = corn2 + (j-1.5d0)*hy - 0.5d0
             else
                x = xcirr(k)
                y = ycirr(k)
             endif
             r2 = x**2 + y**2
             if (r2.lt. 0.04d0) r2 = 0.04d0
             u = y/r2
             v = -x/r2
             p = 20.d0 - 0.5d0/r2
c     
c     ----------------------------
c     # expanding shock:
             rho = 1.4d0
             u = 0.d0
             v = 0.d0
             p = 1.d0
             if (r2.lt..04) then
                rho = 1.4d0
                p = 100.d0
             endif
c     ----------------------------
c     
             val(i,j,1) = rho
             val(i,j,2) = rho*u
             val(i,j,3) = rho*v
             val(i,j,4) = p/gamma1 + 0.5d0*rho*(u**2 + v**2)
 125       continue
          return
 177      continue
          do 216 j = 1, maxjp1
             do 216 i = 1, maxip1
                iuse = i + lwidth - 1
                juse = j + lwidth - 1
                kuse = irr(iuse,juse)
                if (kuse .eq. lstgrd .or. kuse .eq. -1) then
                   x = corn1 + (i-1.5d0)*hx 
                   y = corn2 + (j-1.5d0)*hy
                else
                   x = xcirr(kuse) 
                   y = ycirr(kuse)
                endif

                call p17tru(x,y,rho,u,v,pr)

                if (pr .le. 0.d0) then
c     
c     point inside cylinder, set fake value
                   val(i,j,1) = 1.d0
                   val(i,j,2) = .0d0
                   val(i,j,3) = 0.d0
                   val(i,j,4) = 1.d0/.4d0
                else
                   val(i,j,1) = rho
                   val(i,j,2) = rho * u
                   val(i,j,3) = rho * v
                   val(i,j,4) = pr/(gamma1) + .5d0*rho*(u*u+v*v)
                endif
                if (kuse .ne. lstgrd .and. kuse .ne. -1) then
                   xdif = x - 2.d0
                   ydif = y - 2.d0
                   angle = datan2(ydif,xdif)*180.d0/3.1415926535D0
                   cp = (pr-1.)/(1.4/2.*1.*.05*.05)
                   ii = i-1
                   jj = j-1
                   write(34,345) angle,cp,xdif,ydif,ii,jj
 345               format(4f20.10,2i3)
                endif
 216          continue
             return
c
   18  continue
c      sloc = 0.1d0
c      sloc = 0.06d0
c
c      # initialize to exact solution (rarefaction wave around bend) using
c      # the construction from Whitham:
c
       sloc = -1.d0
       do 1810 i=1,maxip1
          do 1810 j=1,maxjp1
             iadj = i + lwidth - 1
             jadj = j + lwidth - 1
             kirr = irr(iadj,jadj)
             if (kirr.ne.-1 .and. kirr.ne.lstgrd) then
c     # use centroid of irregular cell:
                xcen = xcirr(kirr)
                ycen = ycirr(kirr)
             else
                xcen = (i-1.5d0)*hx + corn1
                ycen = (j-1.5d0)*hy + corn2
             endif
             call p18tru(xcen,ycen,rhot,ut,vt,pt)
             val(i,j,1) = rhot
             val(i,j,2) = rhot*ut
             val(i,j,3) = rhot*vt
             val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)
 1810     continue
          return
c     
 19    continue
       rhoSolid = 1.d0
       uSolid = 0.d0
       vSolid = 0.d0
       pSolid = 1.d0
       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 199 i = 1, maxip1
       do 199 j = 1, maxjp1
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          kirr = irr(iadj,jadj)
          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
             xcen = xcirr(kirr)
             ycen = ycirr(kirr)
          else if (kirr .eq. lstgrd) then
             xcen = corn1 + (dfloat(i)-1.5d0)*hx 
             ycen = corn2 + (dfloat(j)-1.5d0)*hy 
             call makep(poly(1,1,kirr),iadj,jadj,xlow,ylow,hx,hy)
          endif
          if (kirr .ne. -1) then
          call p19tru(xcen,ycen,rhot,ut,vt,pt,poly(1,1,kirr),kirr)
          else
            rhot = rhoSolid
            ut = uSolid
            vt = vSolid
            pt = pSolid
          endif
c        call p19old(xcen,ycen,rhot,ut,vt,pt)
          val(i,j,1) = rhot
          val(i,j,2) = rhot*ut
          val(i,j,3) = rhot*vt
          val(i,j,4) = pt/gamma1 + 0.5d0*rhot*(ut**2 + vt**2)
 199      continue
          return

   40  continue
c
       eps = 1.d-4
       ishock =  (sloc-corn1)/hx + 2
       jshock =  (sloc-corn2)/hy + 2
       do 50 i = 1, maxip1
       do 50 j = 1, maxjp1
          xcen = (i-1.5)*hx + corn1
          ycen = (j-1.5)*hy + corn2
          if (i .lt. ishock) then
             val(i,j,1) = rhol
             val(i,j,2) = rhol*ul
             val(i,j,3) = rhol*vl
             val(i,j,4) = pl/gamma1 + 0.5d0*rhol*(ul**2+vl**2)
          else
             val(i,j,1) = rhor
             val(i,j,2) = rhor*ur
             val(i,j,3) = rhor*vr
             val(i,j,4) = pr/gamma1 + 0.5d0*rhor*(ur**2 + vr**2)
c        else
c        val(i,j,1) = (rhol+rhor)/2.
c        val(i,j,2) = (rhol*ul+rhor*ur)/2.
c        val(i,j,3) = (rhol*vl+rhor*vr)/2.
c        val(i,j,4) =  (pr/gamma1 + 0.5d0*rhor*(ur**2 + vr**2) +
c     .                 pl/gamma1 + 0.5d0*rhol*(ul**2 + vl**2))/2.
          endif
 50    continue
       go to 99
c
 20   continue
c
c annulus test problem from 1d irregular paper, on shifted domain so lower left is origin
c
      xshift = 1.5d0
      yshift = 1.5d0
      pi = 3.14159265357989d0

      do 60 i = 1, maxip1
      do 60 j = 1, maxjp1
         iuse = i + lwidth - 1
         juse = j + lwidth - 1
         kuse = irr(iuse,juse)
         if (kuse .ne. -1 .and. kuse .ne. lstgrd) then
            xcen = xcirr(kuse) - xshift
            ycen = ycirr(kuse) - yshift
         else
            xcen = (i-1.5d0)*hx + corn1 - xshift
            ycen = (j-1.5d0)*hy + corn2 - yshift
         endif
         theta = atan2(ycen,xcen)
         w = .5d0*( erf((pi/6.d0-theta)/sqrt(4.d0/100.d0)) + 
     .              erf((pi/6.d0+theta)/sqrt(4.d0/100.d0)) )
         val(i,j,1) = w
c
c try constant w to test solid body rotation velocities
c         val(i,j,1) = 1.d0

c
c stick channel into problem 20 with fixed velocity, 1 unknown, put shift back in for prob 16
c         val(i,j,1) = ycen+yshift - .1*(xcen+xshift) + .5
 60   continue
      return
c
 21   continue
c     ## Mach 10 flow 30 degree ramp
       sloc = .5d0
       rhor = 1.4d0
       ur = 0.0d0
       vr = 0.0d0
c      ur = cos(pi/6.d0)
c      vr = sin(pi/6.d0)
       pr = 1.0d0
       rhol = 8.d0
c      ul  = 8.25d0*cos(pi/6.d0)
c      vl  =-8.25d0*sin(pi/6.d0)
       ul  = 8.25d0
       vl  = 0.d0
       pl  = 116.5d0
       go to 40

 99   return
      end
