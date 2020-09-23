c
c ------------------------------------------------------------------
c
      subroutine physbd(xleft,xright,ybot,ytop,level,nrow,ncol,nvar,
     1                  val,time,hx,hy)
 
c     This routine takes an (enlarged) grid (or grid patch)
c     with mesh widths hx,hy, and sets the values of any piece of
c     of the patch which extends outside the physical domain using the
c     values given by the boundary conditions. 
c
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      dimension val(nrow,ncol,nvar), poly(10,2)
      data  pi/3.1415926535d0/
      logical viscous/.false./
 
      rmach = .2d0
      alf = alpha*pi/180.
      rho = 1.
c      rho =  1.+cos(alf)*xcen + sin(alf)*ycen
      pr = 1./gamma   ! to match cart3d normalization
      u = rmach *  dcos(alf)
      v = rmach *  dsin(alf)
c
      hxmarg = hx*.01
      hymarg = hy*.01

c     left boundary
 
      if (xleft .lt. -hxmarg) then
 
        nxl = (hxmarg-xleft)/hx
 
           do 400 i = 1,nxl
           do 400 j = 1,ncol
              

c            ycen = ybot  + (dfloat(j)-.5d0)* hy
c            xcen = xleft + (dfloat(i)-.5d0)* hx
c            rho =  4.+.1*(cos(alf)*ycen - sin(alf)*xcen)
               val(i,j,1) = rho
 	       val(i,j,2) = rho*u
 	       val(i,j,3) = rho*v
 	       vel2 = u*u + v*v
 	       val(i,j,4) = pr/gamma1 + .5d0*rho*vel2
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

c  set pressure, compute rho by ratio, extrap rest equiv to dT/dn=0 at outflow
c  vars are in conserved variables
           rhoin  = val(ibeg-1,j,1)
           uin    = val(ibeg-1,j,2)/rhoin
           vin    = val(ibeg-1,j,3)/rhoin
           prin = gamma1*(val(ibeg-1,j,4)-.5d0*rhoin*(uin*uin+vin*vin))
           do 300 i = ibeg, nrow
             val(i,j,1) = pr*rhoin/prin
             rhoout = val(i,j,1)
             val(i,j,2) = rhoout * uin  ! extrap velocities
             val(i,j,3) = rhoout * vin
             val(i,j,4) = pr/gamma1+.5d0*rhoout*(uin*uin+vin*vin)
300        continue
c
      endif
 
 
c     bottom boundary. 
 
      if (ybot .lt. -hymarg) then
        nyb = (hymarg-ybot)/hy
           do 200 j = 1,nyb
           do 200 i = 1,nrow   ! the first lwidth set in left bc above
c$$$               x =  xleft + (dfloat(i)-.5d0)* hx
c$$$               rho = val(i,2*nyb+1-j,1)
c$$$               u   = val(i,2*nyb+1-j,2)/rho
c$$$               v   = val(i,2*nyb+1-j,3)/rho
c$$$               vel2 = u*u + v*v
c$$$               pr  = gamma1*(val(i,2*nyb+1-j,4) - .5d0*rho*vel2)
c$$$               val(i,j,1) =  rho
c$$$               if (x .gt. xprob/3.) then   ! start of viscous wall, no slip
c$$$ 	          val(i,j,2) = -rho*u
c$$$               else
c$$$                  val(i,j,2) =  rho*u   ! otherwise inviscis, extrap tang. flow
c$$$               endif
c$$$ 	       val(i,j,3) = -rho*v   ! negate for no normal flow
c$$$ 	       val(i,j,4) = val(i,2*nyb+1-j,4)    ! pw const pressure

               val(i,j,1) = rho
 	       val(i,j,2) = rho*u
 	       val(i,j,3) = rho*v
 	       vel2 = u*u + v*v
 	       val(i,j,4) = pr/gamma1 + .5d0*rho*vel2
200        continue
c
      endif
c
 99   return
      end
