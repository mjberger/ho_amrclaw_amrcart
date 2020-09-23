c
c --------------------------------------------------------------
c
       subroutine linfcn(q,qx,qy,irr,mitot,mjtot,mptr,ix0,iy0,x,y,
     &            qatxy,slnorm,meqn)
c
      implicit double precision (a-h,o-z)

      include "cirr.i"

      dimension q(mitot,mjtot,meqn),qx(mitot,mjtot,meqn),
     &           qy(mitot,mjtot,meqn),irr(mitot,mjtot)
      dimension qatxy(meqn)
      logical slnorm
      logical    graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder, mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
c  
c  #  evaluate the solution at the point (x,y) in cell (ix0,iy0) using
c  #  a linear reconstruction (if slnorm true, otherwise pw const)  
c  #  of the function, with slopes in qx,qy
c  #  from the cell centroid. 
c
c  # NOTE: for now, q contains primitive variables, slopes are in
c          primitive variables
c
       lstgrd = node(17,mptr)
       level  = node(4,mptr)
       hx     = hxposs(level)
       hy     = hyposs(level)
       k      = irr(ix0,iy0)
c
       if (k .eq. lstgrd) then
c         # regular cell - compute centroid
          xc = rnode(1,mptr) + (ix0-lwidth-.5d0)*hx
          yc = rnode(2,mptr) + (iy0-lwidth-.5d0)*hy
       else
c         # centroid precomputed in irregular data structure
          xc = xcirr(k)
          yc = ycirr(k)
       endif
c
c if euler, no need to convert to primitive variables
c
       do m = 1, meqn
          qatxy(m) =  q(ix0,iy0,m)
       end do
c
c  reconstruct linear function: all slopes (both regular and
c  irregular) should be in qx,qy; evaluate at (x,y).
c
       if (slnorm) then
         do m = 1, meqn
           qatxy(m) = qatxy(m) + (x-xc)*qx(ix0,iy0,m) + 
     &                           (y-yc)*qy(ix0,iy0,m)
         end do
       endif
c

       return
       end
