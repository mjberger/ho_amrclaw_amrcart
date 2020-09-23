c
c---------------------------------------------------------------------
c
c

      subroutine irregflux_gauss(q,firreg,irr,mitot,mjtot,
     &                 dx,dy,lstgrd,xlow,ylow,mptr,
     &                 qx,qy,qxx,qxy,qyy,qxxx,qxxy,qxyy,qyyy,
     &                 lwidth, msize, time)
      implicit double precision (a-h,o-z)
      parameter ( mtsize=9033)

      include "../cirr.i"
      include "../quadrature.i"

      dimension q(mitot,mjtot,4), irr(mitot,mjtot)
      dimension qx(mitot,mjtot,4),qy(mitot,mjtot,4)
      dimension qxx(mitot,mjtot,4),qxy(mitot,mjtot,4),qyy(mitot,mjtot,4)
      dimension qxxx(mitot,mjtot,4),qxxy(mitot,mjtot,4)
      dimension qxyy(mitot,mjtot,4), qyyy(mitot,mjtot,4)
      dimension firreg(-1:irrsize,4)

c
c     # Compute irregular flux in cut cells. Solve Riemann problems
c     # assume there are less than mtsize irregular cell sides on
c     # any particular grid.
c
c ###
      dimension
     &   xnrml(irrsize),  ynrml(irrsize),   rleng(irrsize),
     &   fl(4),fr(4),
     &   gl(4),gr(4), fllf(4),
     &   statel(4),stater(4)

      logical  debug,quad, nolimiter
      integer  xrp, yrp
      integer  icell,jcell

      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      common  /order2/  ssw, quad, nolimiter

      data       xrp/1/, yrp/0/
      data  debug/.false./
      dimension Uout(4)




      if(ihob .eq. 0) then
      call irregflux_gauss_lo(q,firreg,irr,mitot,mjtot,
     &                 dx,dy,lstgrd,xlow,ylow,mptr,
     &                 qx,qy,qxx,qxy,qyy,qxxx,qxxy,qxyy,qyyy,
     &                 lwidth, msize, time)
      else
      call irregflux_gauss_ho(q,firreg,irr,mitot,mjtot,
     &                 dx,dy,lstgrd,xlow,ylow,mptr,
     &                 qx,qy,qxx,qxy,qyy,qxxx,qxxy,qxyy,qyyy,
     &                 lwidth, msize, time)
      endif

      return
      end subroutine



      subroutine irregflux_gauss_lo(q,firreg,irr,mitot,mjtot,
     &                 dx,dy,lstgrd,xlow,ylow,mptr,qx,qy,qxx,qxy,qyy,
     &                 qxxx,qxxy,qxyy,qyyy,
     &                 lwidth, msize, time)
      implicit double precision (a-h,o-z)
      parameter ( mtsize=9033)

      include "../cirr.i"
      include "../quadrature.i"

      dimension q(mitot,mjtot,4), irr(mitot,mjtot)
      dimension qx(mitot,mjtot,4),qy(mitot,mjtot,4)
      dimension qxx(mitot,mjtot,4),qxy(mitot,mjtot,4),qyy(mitot,mjtot,4)
      dimension qxxx(mitot,mjtot,4),qxxy(mitot,mjtot,4)
      dimension qxyy(mitot,mjtot,4)
      dimension qyyy(mitot,mjtot,4)
      dimension firreg(-1:irrsize,4), val(4)
      dimension dir(2), uout(4)
c
c     # Compute irregular flux in cut cells. Solve Riemann problems
c     # assume there are less than mtsize irregular cell sides on 
c     # any particular grid.
c
c ###
      dimension 
     &   xnrml(irrsize),  ynrml(irrsize),   rleng(irrsize),
     &   qlnrml(mtsize,4),qrnrml(mtsize,4), fnrml(mtsize,4)
     
      logical  debug,quad, nolimiter
      integer  xrp, yrp
      integer  icell,jcell

      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      common  /order2/  ssw, quad, nolimiter

      data       xrp/1/, yrp/0/
      data  debug/.false./
c
c     icell= 28
c     jcell= 12
c
c
c       ========================================================
c
c     # loop around boundary: compute outside boxes
c
c       =======================================================
c
      k = lstgrd
      do 100 ik=1,mtsize
         k = iabs(nxtirr(k))
         if (k.eq.0) go to 120
c
c     # march around the sides of this cell to find the irregular side
c
         do 20 kside=1,6
            if (poly(kside+2,1,k).eq.-11.) then
               x1 = poly(kside,1,k)
               y1 = poly(kside,2,k)
               x2 = poly(kside+1,1,k)
               y2 = poly(kside+1,2,k)
               go to 25 
            endif
 20      continue
 25      continue
c
c        # boundary segment face:
c        ------------------------
c
         ix0 = ix(k)
         iy0 = iy(k)
c     # compute (alf,beta) = unit normal to boundary pointing in.
         hsx1 = x2
         hsx2 = x1
         hsy1 = y2
         hsy2 = y1
         rlen = dsqrt((hsy1-hsy2)**2 + (hsx1-hsx2)**2)
         alf = (hsy1-hsy2)/rlen
         beta = (hsx2-hsx1)/rlen
         xnrml(ik) = alf
         ynrml(ik) = beta
         rleng(ik) = rlen
c     
c     # compute boundary data  - q already has primitive variables 
c     


         bxpt = .5d0*(hsx1+hsx2)
         bypt = .5d0*(hsy1+hsy2)
         x0   = xcirr(k)
         y0   = ycirr(k)

         do 33 nn = 1,nlinequad
            bxpt = 0.5d0*(1+rline(nn))*hsx1+0.5d0*(1-rline(nn))*hsx2
            bypt = 0.5d0*(1+rline(nn))*hsy1+0.5d0*(1-rline(nn))*hsy2


!         xdif = (bxpt-x0)
!         ydif = (bypt-y0)
!
!        call evalU(Uout,xdif, ydif, dx, dy,
!     .                 q, qx, qy,
!     .                 qxx, qxy, qyy,
!     .                 qxxx, qxxy, qxyy, qyyy,
!     .                 ix0, iy0, k, mitot, mjtot, 4, 1)
!
!
!      rho = Uout(1)
!
!
!        xright = xlow + mitot*dx
!        ytop   = ylow + mjtot*dy
!
!      call f(val, bxpt,bypt, time, iprob)
!      call getdir(dir,bxpt,bypt)
!
!      dot = alf*dir(1) +beta*dir(2)
!
!
!         firreg(k,1) = firreg(k,1) + ( 0.5d0*dot * ( val(1) + rho )
!     .            +0.5d0*abs(dot)*(  val(1) -rho  ) )*rlen*0.5d0
!         firreg(k,2) = 0.d0
!         firreg(k,3) = 0.d0
!         firreg(k,4) = 0.d0

         firreg(k,:) = 0.d0

 33   continue




 100  continue
      write(6,*) '>>> ERROR in irreg3... ik > mtsize'
      stop
 120  continue


      ikmax = ik-1


      return
      end


      subroutine irregflux_gauss_ho(q,firreg,irr,mitot,mjtot,
     &                 dx,dy,lstgrd,xlow,ylow,mptr,
     &                 qx,qy,qxx,qxy,qyy,qxxx,qxxy,qxyy,qyyy,
     &                 lwidth, msize, time)
      implicit double precision (a-h,o-z)
      parameter ( mtsize=9033)

      include "../cirr.i"
      include "../quadrature.i"

      dimension q(mitot,mjtot,4), irr(mitot,mjtot)
      dimension qx(mitot,mjtot,4),qy(mitot,mjtot,4)
      dimension qxx(mitot,mjtot,4),qxy(mitot,mjtot,4),qyy(mitot,mjtot,4)
      dimension qxxx(mitot,mjtot,4),qxxy(mitot,mjtot,4)
      dimension qxyy(mitot,mjtot,4), qyyy(mitot,mjtot,4)
      dimension firreg(-1:irrsize,4)

      dimension val(4), dir(2)
c
c     # Compute irregular flux in cut cells. Solve Riemann problems
c     # assume there are less than mtsize irregular cell sides on
c     # any particular grid.
c
c ###
      dimension
     &   xnrml(irrsize),  ynrml(irrsize),   rleng(irrsize),
     &   fl(4),fr(4),
     &   gl(4),gr(4), fllf(4),
     &   statel(4),stater(4)

      logical  debug,quad, nolimiter
      integer  xrp, yrp
      integer  icell,jcell

      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      common  /order2/  ssw, quad, nolimiter

      data       xrp/1/, yrp/0/
      data  debug/.false./
      dimension Uout(4)
      dimension err(4)
c
c     icell= 28
c     jcell= 12
c
c
c       ========================================================
c
c     # loop around boundary: compute outside boxes
c
c       =======================================================
c

      err = 0.d0

      k = lstgrd
      do 100 ik=1,mtsize
         k = iabs(nxtirr(k))
         if (k.eq.0) go to 120

       ivert = 1
       do 213 while (poly(ivert+1,1,k) .ne. -11.)
         ivert = ivert + 1
  213  continue


      if(ssw .eq. 2  .or. ssw .eq. -2) then
      x1 = poly(ivert-2,1,k)
      y1 = poly(ivert-2,2,k)

      x2 = bdry(1,1,k)
      y2 = bdry(1,2,k)

      x3 = bdry(3,1,k)
      y3 = bdry(3,2,k)

      x4 = bdry(2,1,k)
      y4 = bdry(2,2,k)


      elseif(ssw .eq. 3 .or. ssw .eq. -3) then
      x1 = poly(ivert-2,1,k)
      y1 = poly(ivert-2,2,k)

      x2 = bdry(1,1,k)
      y2 = bdry(1,2,k)

      x3 = bdry(4,1,k)
      y3 = bdry(4,2,k)

      x6 = bdry(2,1,k)
      y6 = bdry(2,2,k)

      x7 = bdry(3,1,k)
      y7 = bdry(3,2,k)
      else
        print *,"problem in irregfluxgauss"
      endif

         ix0 = ix(k)
         iy0 = iy(k)

         x0   = xcirr_ho(k)
         y0   = ycirr_ho(k)

         dlen = 0.d0
      do 33 nn = 1,nlinequad_ho
         r = (rline_ho(nn)+1.d0)/2.d0
         if(ssw .eq. 2  .or. ssw .eq. -2) then

         bxpt = x2*r*(-1.d0+2.d0*r)
     .         +x3*(1.d0-r)*(1.d0-2.d0*r)
     .         +4.d0*x4*r*(1.d0-r)
         bypt = y2*r*(-1.d0+2.d0*r)
     .         +y3*(1.d0-r)*(1.d0-2.d0*r)
     .         +4.d0*y4*r*(1.d0-r)

         Tx = x2*(-1.d0+2.d0*r)+2.d0*x2*r-x3*(1.d0-2.d0*r)
     .-2.d0*x3*(1.d0-r)+4.d0*x4*(1.d0-r)-4.d0*x4*r
         Ty = y2*(-1.d0+2.d0*r)+2.d0*y2*r-y3*(1.d0-2.d0*r)
     .-2.d0*y3*(1.d0-r)+4.d0*y4*(1.d0-r)-4.d0*y4*r


         dl = dsqrt(Tx**2.d0 + Ty**2.d0)

         dNx =  Ty/dl
         dNy = -Tx/dl
         elseif(ssw .eq. 3 .or. ssw .eq. -3) then
      bxpt = (-dble(9 * r ** 2 * (1 - r)) + dble((-9 * (1 - r) ** 2 + 11
     . - 9 * r) * r) / 0.2D1) * x2 + (-0.9D1 / 0.2D1 * dble(r ** 2) * db
     .le(1 - r) + dble((-18 * (1 - r) ** 2 + 9 - 9 * r) * r) / 0.2D1 + 0
     ..1D1 - dble(r)) * x3 + (0.27D2 / 0.2D1 * dble(r ** 2) * dble(1 - r
     .) - 0.9D1 / 0.2D1 * dble(r) * dble(1 - r)) * x6 + dble((27 * (1 -
     .r) ** 2 - 9 + 9 * r) * r * x7) / 0.2D1
      bypt = (-dble(9 * r ** 2 * (1 - r)) + dble((-9 * (1 - r) ** 2 + 11
     . - 9 * r) * r) / 0.2D1) * y2 + (-0.9D1 / 0.2D1 * dble(r ** 2) * db
     .le(1 - r) + dble((-18 * (1 - r) ** 2 + 9 - 9 * r) * r) / 0.2D1 + 0
     ..1D1 - dble(r)) * y3 + (0.27D2 / 0.2D1 * dble(r ** 2) * dble(1 - r
     .) - 0.9D1 / 0.2D1 * dble(r) * dble(1 - r)) * y6 + dble((27 * (1 -
     .r) ** 2 - 9 + 9 * r) * r * y7) / 0.2D1

      Tx = (-dble(18 * r * (1 - r)) + dble(9 * r ** 2) + dble((9 - 18 *
     .r) * r) / 0.2D1 - 0.9D1 / 0.2D1 * dble((1 - r) ** 2) + 0.11D2 / 0.
     .2D1 - 0.9D1 / 0.2D1 * dble(r)) * x2 + (-dble(9 * r * (1 - r)) + 0.
     .9D1 / 0.2D1 * dble(r ** 2) + dble((27 - 36 * r) * r) / 0.2D1 - dbl
     .e(9 * (1 - r) ** 2) + 0.7D1 / 0.2D1 - 0.9D1 / 0.2D1 * dble(r)) * x
     .3 + (dble(27 * r * (1 - r)) - 0.27D2 / 0.2D1 * dble(r ** 2) - 0.9D
     .1 / 0.2D1 + dble(9 * r)) * x6 + dble((-45 + 54 * r) * r * x7) / 0.
     .2D1 + dble((27 * (1 - r) ** 2 - 9 + 9 * r) * x7) / 0.2D1

      Ty = (-dble(18 * r * (1 - r)) + dble(9 * r ** 2) + dble((9 - 18 *
     .r) * r) / 0.2D1 - 0.9D1 / 0.2D1 * dble((1 - r) ** 2) + 0.11D2 / 0.
     .2D1 - 0.9D1 / 0.2D1 * dble(r)) * y2 + (-dble(9 * r * (1 - r)) + 0.
     .9D1 / 0.2D1 * dble(r ** 2) + dble((27 - 36 * r) * r) / 0.2D1 - dbl
     .e(9 * (1 - r) ** 2) + 0.7D1 / 0.2D1 - 0.9D1 / 0.2D1 * dble(r)) * y
     .3 + (dble(27 * r * (1 - r)) - 0.27D2 / 0.2D1 * dble(r ** 2) - 0.9D
     .1 / 0.2D1 + dble(9 * r)) * y6 + dble((-45 + 54 * r) * r * y7) / 0.
     .2D1 + dble((27 * (1 - r) ** 2 - 9 + 9 * r) * y7) / 0.2D1

         dl = dsqrt(Tx**2.d0 + Ty**2.d0)
         dNx =  Ty/dl
         dNy = -Tx/dl

         else
         print *, "problem in irregflux_gauss"
         endif

!         xdif = bxpt-x0
!         ydif = bypt-y0
!
!        call evalU(Uout,xdif, ydif, dx, dy,
!     .                 q, qx, qy,
!     .                 qxx, qxy, qyy,
!     .                 qxxx, qxxy, qxyy, qyyy,
!     .                 ix0, iy0, k, mitot, mjtot, 4, 1)
!
!
!      rho = Uout(1)
!
!
!        xright = xlow + mitot*dx
!        ytop   = ylow + mjtot*dy
!
!      call f(val, bxpt,bypt, time, iprob)
!      call getdir(dir,bxpt,bypt)
!
!      dot = dNx*dir(1) +dNy*dir(2)
!
!
!         firreg(k,1) = firreg(k,1) + ( 0.5d0*dot * ( val(1) + rho )
!     .            +0.5d0*abs(dot)*(  val(1) -rho  ) )*dl*0.5d0
!         firreg(k,2) = 0.d0
!         firreg(k,3) = 0.d0
!         firreg(k,4) = 0.d0

         firreg(k,:) = 0.d0
  33  continue
 100  continue







      write(6,*) '>>> ERROR in irreg3... ik > mtsize'
      stop
 120  continue


      ikmax = ik-1

!      print *, "L1 boundary flux error ", err
      return
      end subroutine

