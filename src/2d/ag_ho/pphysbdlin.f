c
c ------------------------------------------------------------------
c
      subroutine pphysbdlin(xleft,xright,ybot,ytop,level,nrow,ncol,nvar,
     1                  val,time,hx,hy,qx,qy,irr,lstgrd)
 
c     This routine takes an (enlarged) grid (or grid patch)
c     with mesh widths hx,hy, and sets the values of any piece of
c     of the patch which extends outside the physical domain using the
c     values given by the boundary conditions. 
c
c  # use linear extrap.  variables in primitive already
c
      implicit double precision (a-h,o-z)
      include "../cirr.i"

      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
      dimension val(nrow,ncol,nvar)
      data  pi/3.1415926535d0/
      dimension qx(nrow,ncol,nvar), qy(nrow,ncol,nvar)
      dimension irr(nrow,ncol)

 
      pr = 1.d0
      u = .1d0
      v = .01d0
c
      hxmarg = hx*.01
      hymarg = hy*.01

c     left boundary
 
      if (xleft .lt. -hxmarg) then
 
        nxl = (hxmarg-xleft)/hx
 
           do 400 i = 1,nxl
           do 400 j = 1,ncol
            kuse = irr(i,j)
            if (kuse .ne. -1 .and. kuse .ne. lstgrd) then
                xcen = xcirr(kuse)
                ycen = ycirr(kuse)
            else
              ycen = ybot  + (dfloat(j)-.5d0)* hy
              xcen = xleft + (dfloat(i)-.5d0)* hx
            endif
            rho =  ycen - .1d0*xcen + .5d0
c            rho = 1.d0
            val(i,j,1) = rho
            val(i,j,2) = u
            val(i,j,3) = v
            val(i,j,4) = pr
400        continue
c
      endif
 
 
c     top boundary.  
 
      if (ytop .gt. yprob+hymarg) then
 
        nyt = (ytop - yprob + hymarg)/hy
        jbeg = max0(ncol-nyt+1, 1)
 
           do 100 j= jbeg,ncol
           do 100 i    = 1, nrow
           do 100 ivar=1, nvar
             val(i,j,ivar) = val(i,jbeg-1,ivar)

100        continue
c
      endif
 
c     right boundary. 
      if (xright .gt. xprob+hxmarg) then
 
        nxr = (xright - xprob + hxmarg)/hx
        nxr = nxr - 1
 
        ibeg = max0(nrow-nxr, 1)

c start extrap at bottom of grid, not including ghost cells
           if (ybot .lt. -hymarg) then
              jbeg = (hymarg-ybot)/hy + 1
           else
              jbeg = 1
           endif

           do 300 j = jbeg, ncol
           do 300 i = ibeg, nrow
            kuse = irr(i,j)
            if (kuse .ne. -1 .and. kuse .ne. lstgrd) then
                xcen = xcirr(kuse)
                ycen = ycirr(kuse)
            else
              ycen = ybot  + (dfloat(j)-.5d0)* hy
              xcen = xleft + (dfloat(i)-.5d0)* hx
            endif
            rho =  ycen - .1d0*xcen + .5d0
c            rho = 1.d0
            val(i,j,1) = rho
            val(i,j,2) = u
            val(i,j,3) = v
            val(i,j,4) = pr
300        continue
c
      endif
 
 
c     bottom boundary. 
 
      if (ybot .lt. -hymarg) then
        nyb = (hymarg-ybot)/hy
           do 200 j = 1,nyb
           do 200 i = nxl+1,nrow   ! the first cells set in left bc above
c$$$               x =  xleft + (dfloat(i)-.5d0)* hx
c$$$               rho = val(i,2*nyb+1-j,1)
c$$$               u   = val(i,2*nyb+1-j,2)
c$$$               v   = val(i,2*nyb+1-j,3)
c$$$               vel2 = u*u + v*v
c$$$               pr  = val(i,2*nyb+1-j,4)
c$$$               val(i,j,1) =  rho
c$$$               if (x .gt. xprob/3.) then   ! start of viscous wall, no slip
c$$$                  val(i,j,2) = -u
c$$$               else
c$$$                  val(i,j,2) =  u   ! otherwise inviscis, extrap tang. flow
c$$$               endif
c$$$           val(i,j,3) = -v   ! negate for no normal flow
c$$$           val(i,j,4) = val(i,2*nyb+1-j,4)    ! pw const pressure

             val(i,j,1) = rho
             val(i,j,2) = u
             val(i,j,3) = v
             val(i,j,4) = pr

200        continue
c
      endif
c
 99   return
      end
