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

         xdif = bxpt-x0
         ydif = bypt-y0

        do m = 1,4
        call evalU(Uout,xdif, ydif, dx, dy,
     .                 q, qx, qy,
     .                 qxx, qxy, qyy,
     .                 qxxx, qxxy, qxyy, qyyy,
     .                 ix0, iy0, k, mitot, mjtot, 4, m)
        end do

        if(quad) then
          rhol = Uout(1)
          ul =   Uout(2)
          vl =   Uout(3)
          pl =   Uout(4)
        else
          rhol =    Uout(1)
          rhoul =   Uout(2)
          rhovl =   Uout(3)
          El =      Uout(4)

          ul = rhoul / rhol
          vl = rhovl / rhol
          pl = (gamma - 1.d0)*(El-0.5d0*(rhoul**2.d0+rhovl**2.d0)/rhol)
        endif


!       pressure BC
!         firreg(k,1) = firreg(k,1) + wline_ho(nn)*(  0.d0         )
!         firreg(k,2) = firreg(k,2) + wline_ho(nn)*(  -dl* dNx *pl  )
!         firreg(k,3) = firreg(k,3) + wline_ho(nn)*(  -dl* dNy *pl  )
!         firreg(k,4) = firreg(k,4) + wline_ho(nn)*(  0.d0         )


!        exact BC-1
!         call f(stater, bxpt, bypt, 0, 28)
!         call eulerflux(fllf, stater, dNx, dNy)
!         err(1) = err(1) + dl* wline_ho(nn)*dabs(fllf(1))
!         err(2) = err(2) + dl* wline_ho(nn)*dabs(fllf(2)- dNx *pl)
!         err(3) = err(3) + dl* wline_ho(nn)*dabs(fllf(3)- dNy *pl)
!         err(4) = err(4) + dl* wline_ho(nn)*dabs(fllf(4))

!         firreg(k,1) = firreg(k,1) + wline_ho(nn)*(  0.d0         )
!         firreg(k,2) = firreg(k,2) + wline_ho(nn)*(  -dl* dNx *pr  )
!         firreg(k,3) = firreg(k,3) + wline_ho(nn)*(  -dl* dNy *pr  )
!         firreg(k,4) = firreg(k,4) + wline_ho(nn)*(  0.d0         )

!        exact BC-2
!        statel(1) = rhol
!        statel(2) = rhol * ul
!        statel(3) = rhol * vl
!        statel(4) = pl/(gamma-1.d0) + 0.5d0*rhol*( ul**2.d0+ vl**2.d0 )
!         call f(stater, bxpt, bypt, 0, 28)
!         call riemann(fllf, statel, stater, dNx, dNy)
!         firreg(k,:) = firreg(k,:) - wline_ho(nn)* fllf(:) * dl


!       LLF BC
        statel(1) = rhol
        statel(2) = rhol * ul
        statel(3) = rhol * vl
        statel(4) = pl/(gamma-1.d0) + 0.5d0*rhol*( ul**2.d0+ vl**2.d0 )
        stater(:) = statel(:)
        dot = statel(2) * dNx + statel(3) * dNy;
        stater(2) = statel(2) - 2.d0*dot*dNx;
        stater(3) = statel(3) - 2.d0*dot*dNy;
         call riemann(fllf, statel, stater, dNx, dNy)
         firreg(k,:) = firreg(k,:) - wline_ho(nn)* fllf(:) * dl



  33  continue
 100  continue







      write(6,*) '>>> ERROR in irreg3... ik > mtsize'
      stop
 120  continue


      ikmax = ik-1

!      print *, "L1 boundary flux error ", err
      return
      end subroutine



      subroutine irregflux_gauss_lo(q,firreg,irr,mitot,mjtot,
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
      dimension err(4)
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
c     # compute (alf,beta) = unit normal to boundary pointing out.
         hsx1 = x2
         hsx2 = x1
         hsy1 = y2
         hsy2 = y1
         rlen = dsqrt((hsy1-hsy2)**2 + (hsx1-hsx2)**2)
         alf = -(hsy1-hsy2)/rlen
         beta = -(hsx2-hsx1)/rlen
         xnrml(ik) = alf
         ynrml(ik) = beta
         rleng(ik) = rlen
c
c     # compute boundary data  - q already has primitive variables
c


         x0   = xcirr(k)
         y0   = ycirr(k)

      do 33 nn = 1,nlinequad
            bxpt = 0.5d0*(1+rline(nn))*hsx1+0.5d0*(1-rline(nn))*hsx2
            bypt = 0.5d0*(1+rline(nn))*hsy1+0.5d0*(1-rline(nn))*hsy2


         xdif = (bxpt-x0)
         ydif = (bypt-y0)

        do m = 1,4
        call evalU(Uout,xdif, ydif, dx, dy,
     .                 q, qx, qy,
     .                 qxx, qxy, qyy,
     .                 qxxx, qxxy, qxyy, qyyy,
     .                 ix0, iy0, k, mitot, mjtot, 4, m)
        end do

        if(quad) then
          rhol = Uout(1)
          ul =   Uout(2)
          vl =   Uout(3)
          pl =   Uout(4)
        else
          rhol =    Uout(1)
          rhoul =   Uout(2)
          rhovl =   Uout(3)
          El =      Uout(4)

          ul = rhoul / rhol
          vl = rhovl / rhol
          pl = (gamma - 1.d0)*(El-0.5d0*(rhoul**2.d0+rhovl**2.d0)/rhol)
        endif

        statel(1) = rhol
        statel(2) = rhol * ul
        statel(3) = rhol * vl
        statel(4) = pl/(gamma-1.d0) + 0.5d0*rhol*( ul**2.d0+ vl**2.d0 )

!      ialg = 3
!
!      dot = sqrt(bxpt**2.d0 + bypt**2.d0);
!      dNx = bxpt / dot
!      dNy = bypt / dot
!      if(alf * dNx + beta * dNy < 0 ) then
!         dNx = -dNx
!         dNy = -dNy
!      endif


         firreg(k,1) = firreg(k,1) + wline(nn)*(  0.d0              )
         firreg(k,2) = firreg(k,2) + wline(nn)*(  rlen*(-alf) *pl   )
         firreg(k,3) = firreg(k,3) + wline(nn)*(  rlen*(-beta)*pl   )
         firreg(k,4) = firreg(k,4) + wline(nn)*(  0.d0              )



!        exact BC-1
!         call f(stater, bxpt, bypt, 0, 28)
!         call eulerflux(fllf, stater, alf, beta)
!         err(1) = err(1) + rlen* wline(nn)*dabs(fllf(1))
!         err(2) = err(2) + rlen* wline(nn)*dabs(fllf(2) - alf *pl)
!         err(3) = err(3) + rlen* wline(nn)*dabs(fllf(3) - beta *pl)
!         err(4) = err(4) + rlen* wline(nn)*dabs(fllf(4))
  33  continue
 100  continue


      write(6,*) '>>> ERROR in irreg3... ik > mtsize'
      stop
 120  continue


      ikmax = ik-1

!      print *, "L1 boundary flux error ", err
      return
      end subroutine

      subroutine riemann(fllf,statel, stater, alf,beta)
      implicit double precision (a-h,o-z)
      dimension
     &   fl(4),fr(4),
     &   gl(4),gr(4), fllf(4),
     &   statel(4),stater(4)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold

          rhol =    statel(1)
          rhoul =   statel(2)
          rhovl =   statel(3)
          El =      statel(4)

          ul = rhoul / rhol
          vl = rhovl / rhol
          pl = (gamma - 1.d0)*(El-0.5d0*(rhoul**2.d0+rhovl**2.d0)/rhol)


         rhor = stater(1)
         ur = stater(2)/stater(1)
         vr = stater(3)/stater(1)
         pr = (gamma- 1.d0) * (stater(4)
     .       - (stater(2)**2.d0 + stater(3)**2.d0) / 2.d0 / stater(1))

         cl = dsqrt(gamma * pl / rhol)
         cr = dsqrt(gamma * pr / rhor)

         vell = alf * ul  + beta * vl;
         velr = alf * ur  + beta * vr;

         sl = max(dabs(vell + cl), dabs(vell - cl) )
         sr = max(dabs(velr + cr), dabs(velr - cr) )

         sp = max(sl,sr)

         fl(1) = statel(2)
         fl(2) = statel(2)**2 / statel(1)+pl
         fl(3) = statel(2)*statel(3)/statel(1)
         fl(4) = (statel(2)/statel(1))*(statel(4) + pl)

         fr(1) = stater(2)
         fr(2) = stater(2)**2 / stater(1)+pr
         fr(3) = stater(2)*stater(3)/stater(1)
         fr(4) = (stater(2)/stater(1))*(stater(4) + pr)


         gl(1) = statel(3)
         gl(2) = statel(3)*statel(2)/statel(1)
         gl(3) = statel(3)**2 / statel(1)+pl
         gl(4) = (statel(3)/statel(1))*(statel(4) + pl)

         gr(1) = stater(3)
         gr(2) = stater(3)*stater(2)/stater(1)
         gr(3) = stater(3)**2 / stater(1)+pr
         gr(4) = (stater(3)/stater(1))*(stater(4) + pr)

         fllf(:) = 0.5d0*(alf*(fl(:) + fr(:))
     .                      +beta*(gl(:) + gr(:))
     .                 + sp*(statel(:)-stater(:))  )


      end subroutine

      subroutine eulerflux(fllf,statel, alf,beta)
      implicit double precision (a-h,o-z)
      dimension
     &   fl(4),fr(4),
     &   gl(4),gr(4), fllf(4),
     &   statel(4),stater(4)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold

          rhol =    statel(1)
          rhoul =   statel(2)
          rhovl =   statel(3)
          El =      statel(4)

          ul = rhoul / rhol
          vl = rhovl / rhol
          pl = (gamma - 1.d0)*(El-0.5d0*(rhoul**2.d0+rhovl**2.d0)/rhol)

         fl(1) = statel(2)
         fl(2) = statel(2)**2 / statel(1)+pl
         fl(3) = statel(2)*statel(3)/statel(1)
         fl(4) = (statel(2)/statel(1))*(statel(4) + pl)


         gl(1) = statel(3)
         gl(2) = statel(3)*statel(2)/statel(1)
         gl(3) = statel(3)**2 / statel(1)+pl
         gl(4) = (statel(3)/statel(1))*(statel(4) + pl)

         fllf(:) = alf*fl(:)+beta*gl(:)


      end subroutine
!      if(ssw .eq. 1 .or. ssw .eq. -1) then
!         firreg(k,1) = firreg(k,1) + wline(nn)*(  0.d0              )
!         firreg(k,2) = firreg(k,2) + wline(nn)*(  rlen*(-alf) *pl   )
!         firreg(k,3) = firreg(k,3) + wline(nn)*(  rlen*(-beta)*pl   )
!         firreg(k,4) = firreg(k,4) + wline(nn)*(  0.d0              )
!      elseif(ialg .eq. 1) then ! Algorithm 1
!         rhog = rhol
!         ug =  (ul * dNy - vl * dNx) * dNy
!         vg = -(ul * dNy - vl * dNx) * dNx
!         pg =  pl
!         Eg =  pg / (gamma - 1.d0) + 0.5d0*rhog*(ug**2.d0+ vg**2.d0)
!
!       vbdotn = ug * alf + vg * beta
!       firreg(k,1) = firreg(k,1) - wline(nn)*(rhog      * vbdotn )*rlen
!       firreg(k,2) = firreg(k,2) - wline(nn)*(rhog * ug * vbdotn +
!     . alf*pg)*rlen
!       firreg(k,3) = firreg(k,3) - wline(nn)*(rhog * vg * vbdotn  +
!     . beta*pg)*rlen
!       firreg(k,4) = firreg(k,4) - wline(nn)*((Eg + pg) * vbdotn )*rlen
!
!      elseif(ialg .eq. 2) then ! Algorithm 2
!         rhog = rhol
!         ug = ul * (dNy**2.d0 - dNx**2.d0) - 2.d0*dNx*dNy*vl
!         vg = vl * (dNx**2.d0 - dNy**2.d0) - 2.d0*dNx*dNy*ul
!         pg = pl
!
!         stater(1) =  statel(1)
!         stater(2) =  rhog * ug
!         stater(3) =  rhog * vg
!         stater(4) =  pg/(gamma-1.d0)+0.5d0*rhol*(ug**2.d0+vg**2.d0)
!
!!         ddd = alf * (ul+ug)/2.d0 + beta * (vl+vg)/2.d0
!
!         call riemann(fllf, statel, stater, alf, beta)
!         firreg(k,:) = firreg(k,:) - wline(nn)* fllf(:) * rlen
!
!      elseif(ialg .eq. 3) then ! Algorithm 3
!         x = alf*dNx+beta*dNy
!!         angle = dacos(alf*dNx+beta*dNy)
!         vdotn = ul * alf + vl * beta
!
!         dsign = 1.d0
!         if(vdotn < 0.d0) dsign = -1.d0
!
!         v_T = ul * dNy - vl * dNx
!!         v_N = vdotn - 2.d0*(vdotn-dsign*dabs(v_T)*dtan(angle))
!         v_N = vdotn-2.d0*(vdotn-dsign*dabs(v_T)*dsqrt(1.d0-x**2.d0)/x)
!
!         ug = v_N*alf  + v_T*beta
!         vg = v_N*beta - v_T*alf
!
!!         ddd = dNx * (ul+ug)/2.d0 + dNy * (vl + vg)/2.d0
!
!         stater(1) = rhol
!         stater(2) = rhol * ug
!         stater(3) = rhol * vg
!         stater(4) = pl/(gamma-1.d0) + 0.5d0*rhol*( ug**2.d0+ vg**2.d0 )
!
!         call riemann(fllf, statel, stater, alf, beta)
!         firreg(k,:) = firreg(k,:) - wline(nn)* fllf(:) * rlen
!
!
!      else ! exact BC
!!     EXACT BC
!         call f(stater, bxpt, bypt, 0, 28)
!         call riemann(fllf, statel, stater, alf, beta)
!         firreg(k,:) = firreg(k,:) - wline(nn)* fllf(:) * rlen
!      endif
