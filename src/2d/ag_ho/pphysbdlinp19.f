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
c 
c
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      dimension val(nrow,ncol,nvar)
      data  pi/3.1415926535d0/
      logical viscous/.false./
      dimension qx(nrow,ncol,nvar), qy(nrow,ncol,nvar)
      dimension irr(nrow,ncol)
      include "cirr.i"

 

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
              call makep(poly(1,1,kuse),i,j,xlow,ylow,hx,hy)
            endif
            call p19tru(xcen,ycen,rhot,ut,vt,pt,poly(1,1,kuse),kuse)
c           call p19old(xcen,ycen,rhot,ut,vt,pt)

               val(i,j,1) = rhot
 	       val(i,j,2) = ut
 	       val(i,j,3) = vt
 	       val(i,j,4) = pt
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
            rho =  ycen - .1*xcen + .5
               val(i,j,1) = rho
 	       val(i,j,2) = u
 	       val(i,j,3) = v
 	       val(i,j,4) = pr
300        continue
c
      endif
 
 
c     bottom boundary. 
c
c for supersonic vortex, are currently setting exact solution
c might want to extrapolate outflow instead
 
      if (ybot .lt. -hymarg) then
        nyb = (hymarg-ybot)/hy
           do 200 j = 1,nyb
           do 200 i = 1,nrow   
            kuse = irr(i,j)
            if (kuse .ne. -1 .and. kuse .ne. lstgrd) then
                xcen = xcirr(kuse)
                ycen = ycirr(kuse)
            else
              ycen = ybot  + (dfloat(j)-.5d0)* hy
              xcen = xleft + (dfloat(i)-.5d0)* hx
              call makep(poly(1,1,kuse),i,j,xlow,ylow,hx,hy)
            endif
            call p19tru(xcen,ycen,rhot,ut,vt,pt,poly(1,1,kuse),kuse)
c           call p19old(xcen,ycen,rhot,ut,vt,pt)

               val(i,j,1) = rhot
 	       val(i,j,2) = ut
 	       val(i,j,3) = vt
 	       val(i,j,4) = pt

200        continue
c
      endif
c
 99   return
      end
