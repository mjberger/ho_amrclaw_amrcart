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
     &   qlnrml(mtsize,4),qrnrml(mtsize,4), fnrml(mtsize,4)
     
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
          rho = Uout(1)
          u =   Uout(2)
          v =   Uout(3)
          pr =  Uout(4)
        else
          rho =    Uout(1)
          rhou =   Uout(2)
          rhov =   Uout(3)
          E =      Uout(4)

          u = rhou / rho
          v = rhov / rho
          pr = (gamma - 1.d0) * (E-0.5d0*(rhou**2 + rhov**2)/ rho)
        endif

c     
c        # compute rotated data inside cell for b.r.p. and save
         qrnrml(ik,1) = rho
         qrnrml(ik,2) = alf*u + beta*v
         qrnrml(ik,3) = -beta*u + alf*v
         qrnrml(ik,4) = pr

         qlnrml(ik,1) = qrnrml(ik,1)
         qlnrml(ik,2) = -qrnrml(ik,2)
         qlnrml(ik,3) = qrnrml(ik,3)
         qlnrml(ik,4) = qrnrml(ik,4)

c     NEW WAY: NO RP, JUST PRESSURE AT BNDRY
c     
         firreg(k,1) = firreg(k,1) + wline(nn)*(  0.d0          )
         firreg(k,2) = firreg(k,2) + wline(nn)*(  rlen*alf*pr   )
         firreg(k,3) = firreg(k,3) + wline(nn)*(  rlen*beta*pr  )
         firreg(k,4) = firreg(k,4) + wline(nn)*(  0.d0          )
33    continue


 100  continue
      write(6,*) '>>> ERROR in irreg3... ik > mtsize'
      stop
 120  continue


      ikmax = ik-1


      return
      end


c     OLD WAY: SOLVE RP AT BNDRY
!--c     # solve Riemann problems at all boundary faces:
!--      call vrm(qrnrml,qlnrml,fnrml,1,ikmax,xrp,mtsize)
!--
!--c     # rotate resulting fluxes back to x-y and save in firreg:
!--c     # Note that only pressure is nonzero.
!--c
!--      k = lstgrd
!--      do 150 ik=1,ikmax
!--         k = iabs(nxtirr(k))
!--         firreg(k,1) = 0.d0
!--         firreg(k,2) = rleng(ik)*xnrml(ik)*fnrml(ik,2)
!--         firreg(k,3) = rleng(ik)*ynrml(ik)*fnrml(ik,2)
!--         firreg(k,4) = 0.d0
!-- 150  continue
c
!         rho_c       = q(ix0,iy0,1)
!         u_c         = q(ix0,iy0,2)
!         v_c         = q(ix0,iy0,3)
!         pr_c        = q(ix0,iy0,4)
!
!         rho = rho_c + xdif*qx(ix0,iy0,1) + ydif*qy(ix0,iy0,1)
!     .   + 0.5d0*qxx(ix0,iy0,1)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,1)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,1)*( ydif**2 -(dy**2)*poly(10,1,k) )
!
!         u   = u_c +   xdif*qx(ix0,iy0,2) + ydif*qy(ix0,iy0,2)
!     .   + 0.5d0*qxx(ix0,iy0,2)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,2)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,2)*( ydif**2 -(dy**2)*poly(10,1,k) )
!
!         v   = v_c +   xdif*qx(ix0,iy0,3) + ydif*qy(ix0,iy0,3)
!     .   + 0.5d0*qxx(ix0,iy0,3)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,3)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,3)*( ydif**2 -(dy**2)*poly(10,1,k) )
!
!         pr  = pr_c +  xdif*qx(ix0,iy0,4) + ydif*qy(ix0,iy0,4)
!     .   + 0.5d0*qxx(ix0,iy0,4)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,4)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,4)*( ydif**2 -(dy**2)*poly(10,1,k) )

      !         rho_c       = q(ix0,iy0,1)
!         rhou_c         = q(ix0,iy0,2)
!         rhov_c         = q(ix0,iy0,3)
!         E_c        = q(ix0,iy0,4)
!
!         rho = rho_c + xdif*qx(ix0,iy0,1) + ydif*qy(ix0,iy0,1)
!     .   + 0.5d0*qxx(ix0,iy0,1)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,1)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,1)*( ydif**2 -(dy**2)*poly(10,1,k) )
!
!         rhou   = rhou_c +   xdif*qx(ix0,iy0,2) + ydif*qy(ix0,iy0,2)
!     .   + 0.5d0*qxx(ix0,iy0,2)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,2)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,2)*( ydif**2 -(dy**2)*poly(10,1,k) )
!
!         rhov   = rhov_c +   xdif*qx(ix0,iy0,3) + ydif*qy(ix0,iy0,3)
!     .   + 0.5d0*qxx(ix0,iy0,3)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,3)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,3)*( ydif**2 -(dy**2)*poly(10,1,k) )
!
!         E  = E_c +  xdif*qx(ix0,iy0,4) + ydif*qy(ix0,iy0,4)
!     .   + 0.5d0*qxx(ix0,iy0,4)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,4)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,4)*( ydif**2 -(dy**2)*poly(10,1,k) )






!         rho_c       = q(ix0,iy0,1)
!         u_c         = q(ix0,iy0,2)
!         v_c         = q(ix0,iy0,3)
!         pr_c        = q(ix0,iy0,4)
!
!         rho = rho_c + xdif*qx(ix0,iy0,1) + ydif*qy(ix0,iy0,1)
!     .   + 0.5d0*qxx(ix0,iy0,1)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,1)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,1)*( ydif**2 -(dy**2)*poly(10,1,k) )
!
!         u   = u_c +   xdif*qx(ix0,iy0,2) + ydif*qy(ix0,iy0,2)
!     .   + 0.5d0*qxx(ix0,iy0,2)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,2)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,2)*( ydif**2 -(dy**2)*poly(10,1,k) )
!
!         v   = v_c +   xdif*qx(ix0,iy0,3) + ydif*qy(ix0,iy0,3)
!     .   + 0.5d0*qxx(ix0,iy0,3)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,3)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,3)*( ydif**2 -(dy**2)*poly(10,1,k) )
!
!         pr  = pr_c +  xdif*qx(ix0,iy0,4) + ydif*qy(ix0,iy0,4)
!     .   + 0.5d0*qxx(ix0,iy0,4)*( xdif**2-(dx**2)*poly(8,1,k)   )
!     .   +       qxy(ix0,iy0,4)*( xdif*ydif-dx*dy*poly(9,1,k) )
!     .   + 0.5d0*qyy(ix0,iy0,4)*( ydif**2 -(dy**2)*poly(10,1,k) )
