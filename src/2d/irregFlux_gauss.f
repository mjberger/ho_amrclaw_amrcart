c
c---------------------------------------------------------------------
c
c
      subroutine irregFlux_gauss(q,firreg,irr,mitot,mjtot,
     &                 dx,dy,lstgrd,xlow,ylow,mptr,qx,qy,qxx,qxy,qyy,
     &                 lwidth, nvar, time)

      use amr_module

      implicit double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot), irr(mitot,mjtot)
      dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
      dimension qxx(nvar,mitot,mjtot),qyy(nvar,mitot,mjtot)
      dimension qxy(nvar,mitot,mjtot)
      dimension firreg(nvar,-1:irrsize)
c
c     # Compute irregular flux in cut cells. Solve Riemann problems
c     # assume there are less than mtsize irregular cell sides on 
c     # any particular grid.
c
c ###
      logical  debug,quad, nolimiter
      integer  xrp, yrp
      integer  icell,jcele

      common /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold,pwconst
      common  /order2/  ssw, quad, nolimiter

      data       xrp/1/, yrp/0/
      data  debug/.false./
c
c       ========================================================
c
c     # loop around boundary: compute pressure at boundary
c     # no RP in this version
c
c     # later will add curved boundary options
c
c       =======================================================
c
      write(*,*)"Resetting pr for debugging in irregFlux_gauss"

      k = iabs(nxtirr(lstgrd))
      do while (k .ne. 0) 
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
         ix0 = ixg(k)
         iy0 = iyg(k)
c     # compute (alf,beta) = unit normal to boundary pointing in.
         hsx1 = x2
         hsx2 = x1
         hsy1 = y2
         hsy2 = y1
         rlen = dsqrt((hsy1-hsy2)**2 + (hsx1-hsx2)**2)
         alf = (hsy1-hsy2)/rlen
         beta = (hsx2-hsx1)/rlen
c     
c     # compute boundary data  - q already has primitive variables 
c     
         w1 = 0.5d0*(1-dsqrt(3.d0)/3.d0)
         w2 = 0.5d0*(1+dsqrt(3.d0)/3.d0)
         x0   = xcirr(k)
         y0   = ycirr(k)

         firreg(:,k) = 0.d0   ! initialize since will accumulate at gauss points
         do 33 nn = 1,2

           if(nn .eq. 1) then
             bxpt = w1*hsx1+(1.d0-w1)*hsx2
             bypt = w1*hsy1+(1.d0-w1)*hsy2
           else
             bxpt = w2*hsx1+(1.d0-w2)*hsx2
             bypt = w2*hsy1+(1.d0-w2)*hsy2
           endif


           xdif = (bxpt-x0)
           ydif = (bypt-y0)

           rho_c    = q(1,ix0,iy0)
           u_c      = q(2,ix0,iy0)
           v_c      = q(3,ix0,iy0)
           pr_c     = q(4,ix0,iy0)

c          comment out other fields only need density in this approach
c          if want to go back to solving RP and rotating see older codes

c          rho = rho_c + xdif*qx(1,ix0,iy0) + ydif*qy(1,ix0,iy0)
c      .   + 0.5d0*qxx(1,ix0,iy0)*( xdif**2-(dx**2)*poly(8,1,k)   )
c      .   +       qxy(1,ix0,iy0)*( xdif*ydif-dx*dy*poly(9,1,k) )
c      .   + 0.5d0*qyy(1,ix0,iy0)*( ydif**2 -(dy**2)*poly(10,1,k) )

c          u   = u_c +   xdif*qx(2,ix0,iy0,2) + ydif*qy(2,ix0,iy0,2)
c      .   + 0.5d0*qxx(2,ix0,iy0)*( xdif**2-(dx**2)*poly(8,1,k)   )
c      .   +       qxy(2,ix0,iy0)*( xdif*ydif-dx*dy*poly(9,1,k) )
c      .   + 0.5d0*qyy(2,ix0,iy0)*( ydif**2 -(dy**2)*poly(10,1,k) )

c          v   = v_c +   xdif*qx(3,ix0,iy0) + ydif*qy(3,ix0,iy0)
c      .   + 0.5d0*qxx(3,ix0,iy0)*( xdif**2-(dx**2)*poly(8,1,k)   )
c      .   +       qxy(3,ix0,iy0)*( xdif*ydif-dx*dy*poly(9,1,k) )
c      .   + 0.5d0*qyy(3,ix0,iy0)*( ydif**2 -(dy**2)*poly(10,1,k) )

           pr  = pr_c +  xdif*qx(4,ix0,iy0) + ydif*qy(4,ix0,iy0)
     &      + 0.5d0*qxx(4,ix0,iy0)*( xdif**2-(dx**2)*poly(8,1,k)   )
     &      +       qxy(4,ix0,iy0)*( xdif*ydif-dx*dy*poly(9,1,k) )
     &      + 0.5d0*qyy(4,ix0,iy0)*( ydif**2 -(dy**2)*poly(10,1,k) )

c     
c          ##  NEW WAY: NO RP, JUST PRESSURE AT BNDRY
c     
           pr = 1.d0

           firreg(1,k) = firreg(1,k) + 0.5d0*(  0.d0          )
           firreg(2,k) = firreg(2,k) + 0.5d0*(  rlen*alf*pr   )
           firreg(3,k) = firreg(3,k) + 0.5d0*(  rlen*beta*pr  )
           firreg(4,k) = firreg(4,k) + 0.5d0*(  0.d0          )
 33     continue
        k = iabs(nxtirr(k))
      end do


 120  continue

      return
      end      
