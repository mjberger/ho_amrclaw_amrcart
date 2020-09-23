c
c -------------------------------------------------------------
c
      subroutine bound(xleft,xright,ybot,ytop,hx,hy,time,intime,
     1                 level,nvar,lwidth,valbig,mitot,mjtot,irr,lstgrd)
 
      implicit double precision (a-h,o-z)
      dimension valbig(mitot,mjtot,nvar), irr(mitot,mjtot)
 
c     This routine sets the boundary values  for the grid mptr
c     at level level.
c     We are setting the values for a strip lwidth zones wide all
c     the way around the border, in 4 rectangular strips.
c
c     Inputs to this routine:
c     xleft, xright, ybot, ytop = the location in physical space of
c     corners of the grid.
c     hx, hy = the spatial increments of the grid
c
c     Outputs from this routine:
c     The values around the border of the grid are inserted 
c     directly into the enlarged valbig array.
c
c     This routine calls the routine filpatch
c     which for any block of mesh points on a given level,
c     intersects that block with all grids on that level and with
c     the physical boundaries, copies the values into the
c     appropriate intersecting regions, and interpolates the remaining
c     cells from coarser grids as required.
 
c     left boundary
 
      xl = xleft - lwidth*hx
      yb = ybot - lwidth*hy
      xr = xleft
      yt = ytop
 
      call filpatch(xl,xr,yb,yt,level,nvar,valbig,time,intime,
     1              mitot,mjtot,1,lwidth,1,mjtot-lwidth,irr)
 
c     right boundary
 
      xl = xright
      yb = ybot
      xr = xright + lwidth*hx
      yt = ytop + lwidth*hy
 
      call filpatch(xl,xr,yb,yt,level,nvar,valbig,time,intime,
     1              mitot,mjtot,mitot-lwidth+1,mitot,lwidth+1,mjtot,irr) 
 
c     bottom boundary
 
      xl = xleft
      yb = ybot - lwidth*hy
      xr = xright + lwidth*hx
      yt = ybot
 
      call filpatch(xl,xr,yb,yt,level,nvar,valbig,time,intime,
     1              mitot,mjtot,lwidth+1,mitot,1,lwidth,irr)
 
c     top boundary
 
      xl = xleft - lwidth*hx
      yb = ytop
      xr = xright
      yt = ytop + lwidth*hy
 
      call filpatch(xl,xr,yb,yt,level,nvar,valbig,time,intime,
     1              mitot,mjtot,1,mitot-lwidth,mjtot-lwidth+1,mjtot,irr) 
c
c external boundary conditions
c
      xl = xleft - lwidth*hx
      yb = ybot - lwidth*hy
      xr = xright + lwidth*hx
      yt = ytop + lwidth*hy
      call physbd(xl,xr,yb,yt,level,mitot,
     1            mjtot,nvar,valbig,time,hx,hy)
c
      return
      end
