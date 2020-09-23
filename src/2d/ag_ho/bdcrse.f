c
c --------------------------------------------------------------
c
      subroutine bdcrse( xleft,xright,ybot,ytop,hx,hy,time,intime,
     1                   level,nvar,lw2,valdub,midub,mjdub,valbig,
     2                   mi2tot,mj2tot,irr,mitot,mjtot)
c 
c  bdcrse: first set double the stencil size worth of boundary values,
c          then coarsen them for the giant step integration.
c
      implicit double precision (a-h,o-z)
      include "cirr.i"
      dimension valdub(midub,mjdub,nvar), valbig(mi2tot,mj2tot,nvar),
     1          irr(mitot,mjtot)
 
c     left boundary
 
      xl = xleft - lw2*hx
      yb = ybot - lw2*hy
      xr = xleft
      yt = ytop
 
      call filpc(xl,xr,yb,yt,level,nvar,valdub,time,intime,midub,
     1           mjdub,1,lw2,1,mjdub-lw2,irr,mitot,mjtot)
 
c     right boundary
 
      xl = xright
      yb = ybot
      xr = xright + lw2*hx
      yt = ytop + lw2*hy
 
      call filpc(xl,xr,yb,yt,level,nvar,valdub,time,intime,midub,
     1           mjdub,midub-lw2+1,midub,lw2+1,mjdub,irr,mitot,mjtot)
 
c     bottom boundary
 
      xl = xleft
      yb = ybot - lw2*hy
      xr = xright + lw2*hx
      yt = ybot
 
      call filpc(xl,xr,yb,yt,level,nvar,valdub,time,intime,midub,
     1           mjdub,lw2+1,midub,1,lw2,irr,mitot,mjtot)
 
c     top boundary
 
      xl = xleft - lw2*hx
      yb = ytop
      xr = xright
      yt = ytop + lw2*hy
 
      call filpc(xl,xr,yb,yt,level,nvar,valdub,time,intime,midub,
     1           mjdub,1,midub-lw2,mjdub-lw2+1,mjdub,irr,mitot,mjtot)
 
c
c now coarsen by 2 in every direction - conservatively, so that
c irregular and exterior cells count appropriately.
c irr array only exists for lwidth, not 2*lwidth sized grid.
c boundary cells not done correctly here.
c
      lwidth = lw2 / 2
      ar(-1) = 0.d0
      do 100 ivar = 1,nvar
      do 100 j = 1, mj2tot
         jfine = 2*(j - 1) + 1
         do 100 i = 1, mi2tot
            ifine = 2*(i - 1) + 1
	    tarea = 0.d0
	    if ((ifine .le. lwidth) .or. (jfine .le. lwidth) .or.
     .          (ifine+1 .gt. midub-lwidth) .or. 
     .          (jfine+1 .gt. mjdub-lwidth)) then
	       a11 = 1.d0
	       a21 = 1.d0
	       a12 = 1.d0
	       a22 = 1.d0
	    else
	       ieff = ifine - lwidth
	       jeff = jfine - lwidth
	       a11   = ar(irr(ieff,jeff))
	       a21   = ar(irr(ieff+1,jeff))
	       a12   = ar(irr(ieff,jeff+1))
	       a22   = ar(irr(ieff+1,jeff+1))
	    endif
c	
	    tarea = a11+a12+a21+a22
c
            if (tarea .eq. 0.d0) then
	      valbig(i,j,ivar) = valdub(ifine,jfine,ivar)
	    else
              valbig(i,j,ivar)=(valdub(ifine,jfine,ivar)*a11 +
     &                    valdub(ifine+1,jfine,ivar)*a21+
     &                    valdub(ifine,jfine+1,ivar)*a12 +
     &                    valdub(ifine+1,jfine+1,ivar)*a22)/tarea
	    endif
100   continue
c
c  external boundary conditions
c 
       lwidth = lw2 / 2
       hx2    = 2.d0*hx
       hy2    = 2.d0*hy
       xl = xleft  - lwidth*hx2
       yb = ybot   - lwidth*hy2
       xr = xright + lwidth*hx2
       yt = ytop   + lwidth*hy2
       call physbd(xl,xr,yb,yt,level,mi2tot,mj2tot,nvar,valbig,
     1             time,hx2,hy2)
c 
      return
      end
