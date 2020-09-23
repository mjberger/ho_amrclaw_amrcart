c
c---------------------------------------------------------------------
c
      subroutine makebx(ix0,iy0,ixadj,iyadj,ivert,x1,y1,x2,y2,ddx,ddy,
     &                  box,areabx,ix1,iy1,ix2,iy2)
      implicit double precision (a-h,o-z)
      dimension box(10,2)
c
c     # create box, a list of vertices corresponding to a box starting at
c     # (x1,y1), (x2,y2) and going out (ddx,ddy).  
c     # On input, 
c     #    (ix0,iy0) is index of current cell, 
c     #               with (x1,y1),(x2,y2) being one side of this cell
c     #    (ixadj,iyadj) is the adjacent cell across this side
c     #    ivert = 1 if this is a vertical side
c     #            0 if this is a horizontal side
c     # On output,
c     #    box is the box
c     #    areabx is the area of the box
c     #    (ix1,iy1) is the cell that the box first intersects, going
c     #             away from the initial side
c     #    (ix2,iy2) is the other cell this box intersects, if any.
c
      ihoriz = 1-ivert
      if (ddx*(y2-y1)-ddy*(x2-x1) .gt. 0.d0) then
	  is = 1
	  ix1 = ix0
	  iy1 = iy0
	else
	  is = -1
	  ix1 = ixadj
	  iy1 = iyadj
	endif
c
      ix2 = ix1 + ihoriz*isig(ddx)
      iy2 = iy1 + ivert*isig(ddy)
c     ix2 = ix1 + dsign(dfloat(ihoriz),(ddx))
c     iy2 = iy1 + dsign(dfloat(ivert),(ddy))
      areabx = ivert*dabs((y2-y1)*ddx) + ihoriz*dabs((x2-x1)*ddy)
c
      box(1,1) = x1
      box(1,2) = y1
      box(3-is,1) = x2
      box(3-is,2) = y2
      box(3,1) = x2 + ddx
      box(3,2) = y2 + ddy
      box(3+is,1) = x1 + ddx
      box(3+is,2) = y1 + ddy
      box(5,1) = x1
      box(5,2) = y1
      box(6,1) = -11
c
      return
      end
