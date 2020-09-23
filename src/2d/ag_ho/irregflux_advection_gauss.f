c
c---------------------------------------------------------------------
c
c
      subroutine irregflux_gauss(q,firreg,irr,mitot,mjtot,
     &                 dx,dy,lstgrd,xlow,ylow,mptr,qx,qy,qxx,qxy,qyy,
     &                 lwidth, msize, time)
      implicit double precision (a-h,o-z)
      parameter ( mtsize=9033)

      include "cirr.i"

      dimension q(mitot,mjtot,4), irr(mitot,mjtot)
      dimension qx(mitot,mjtot,4),qy(mitot,mjtot,4)
      dimension qxx(mitot,mjtot,4),qxy(mitot,mjtot,4),qyy(mitot,mjtot,4)
      dimension firreg(-1:irrsize,4), val(4)
      dimension dir(2)
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
         rho_c       = q(ix0,iy0,1)
         u_c         = q(ix0,iy0,2)
         v_c         = q(ix0,iy0,3)
         pr_c        = q(ix0,iy0,4)

         bxpt = .5d0*(hsx1+hsx2)
         bypt = .5d0*(hsy1+hsy2)
         x0   = xcirr(k)
         y0   = ycirr(k)
         w1 = 0.5d0*(1-dsqrt(3.d0)/3.d0)
         w2 = 0.5d0*(1+dsqrt(3.d0)/3.d0)

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


         rho = rho_c + xdif*qx(ix0,iy0,1) + ydif*qy(ix0,iy0,1)
     .   + 0.5d0*qxx(ix0,iy0,1)*( xdif**2-(dx**2)*poly(8,1,k)   )
     .   +       qxy(ix0,iy0,1)*( xdif*ydif-dx*dy*poly(9,1,k) )
     .   + 0.5d0*qyy(ix0,iy0,1)*( ydif**2 -(dy**2)*poly(10,1,k) )

         u   = u_c +   xdif*qx(ix0,iy0,2) + ydif*qy(ix0,iy0,2)
     .   + 0.5d0*qxx(ix0,iy0,2)*( xdif**2-(dx**2)*poly(8,1,k)   )
     .   +       qxy(ix0,iy0,2)*( xdif*ydif-dx*dy*poly(9,1,k) )
     .   + 0.5d0*qyy(ix0,iy0,2)*( ydif**2 -(dy**2)*poly(10,1,k) )

         v   = v_c +   xdif*qx(ix0,iy0,3) + ydif*qy(ix0,iy0,3)
     .   + 0.5d0*qxx(ix0,iy0,3)*( xdif**2-(dx**2)*poly(8,1,k)   )
     .   +       qxy(ix0,iy0,3)*( xdif*ydif-dx*dy*poly(9,1,k) )
     .   + 0.5d0*qyy(ix0,iy0,3)*( ydif**2 -(dy**2)*poly(10,1,k) )

         pr  = pr_c +  xdif*qx(ix0,iy0,4) + ydif*qy(ix0,iy0,4)
     .   + 0.5d0*qxx(ix0,iy0,4)*( xdif**2-(dx**2)*poly(8,1,k)   )
     .   +       qxy(ix0,iy0,4)*( xdif*ydif-dx*dy*poly(9,1,k) )
     .   + 0.5d0*qyy(ix0,iy0,4)*( ydif**2 -(dy**2)*poly(10,1,k) )





        xright = xlow + mitot*dx
        ytop   = ylow + mjtot*dy

      call f(val, bxpt,bypt, time, iprob)
      call getdir(dir)

      dot = alf*dir(1) +beta*dir(2)


         firreg(k,1) = firreg(k,1) + ( 0.5d0*dot * ( val(1) + rho )
     .            +0.5d0*abs(dot)*(  val(1) -rho  ) )*rlen*0.5d0
         firreg(k,2) = 0.d0
         firreg(k,3) = 0.d0
         firreg(k,4) = 0.d0


 33   continue




 100  continue
      write(6,*) '>>> ERROR in irreg3... ik > mtsize'
      stop
 120  continue


      ikmax = ik-1


      return
      end      
