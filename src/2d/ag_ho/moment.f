c
c ----------------------------------------------------------
c
      subroutine moment (rect,badpts,npt,usage,buffx,buffy,hx,hy)
c
      implicit double precision (a-h,o-z)
      dimension     rect(12),badpts(2,npt)
c
c input parameters:
c     badpts      = x,y coords of flagged badpts grouped into clusters
c                   are in the first two rows
c     npt         = num. of badpts. in the cluster.
c     buffx,y     = expand new grids with buffer zone of size buffer.
c                   this also prevents degenerate point rectangles.
c output parameters:
c     usage       = ratio of flagged to unflagged badpts. in new grid
c                   measures goodness of fit and clustering
c    rect( )      = stores some info. for grid created herein.
c                   sometimes rect = rnode, sometimes = temp. array.
c                   depending on calling prog. (grdfit or expand)
c
c
c  moment = compute rectangle around flagged points.
c  save some info., even tho. usage might be low and rect. scrapped.
c
      rn = dfloat(npt)
c
c compute length of enclosing rectangles to include all flagged badpts.
c expand on all sides by 'buffer' zone.  this also prevents
c case of degenerate rectangle (line).
c
      emx1 = badpts(1,1)
      emn1 = emx1
      emx2 = badpts(2,1)
      emn2 = emx2
      do 80 ipt = 1, npt
          if (badpts(1,ipt) .gt. emx1) emx1 = badpts(1,ipt)
          if (badpts(1,ipt) .lt. emn1) emn1 = badpts(1,ipt)
          if (badpts(2,ipt) .gt. emx2) emx2 = badpts(2,ipt)
          if (badpts(2,ipt) .lt. emn2) emn2 = badpts(2,ipt)
 80   continue
      emx1  = emx1 + buffx
      emx2  = emx2 + buffy
      emn1  = emn1 - buffx
      emn2  = emn2 - buffy
c
c from length of the sides, determine the 4 rect. corners
c number the corners clockwise
c
      rect(1) = emn1
      rect(2) = emn2
      rect(3) = emn1
      rect(4) = emx2
      rect(5) = emx1
      rect(6) = emx2
      rect(7) = emx1
      rect(8) = emn2
c
c compute volume cutoff ratio
c
      iside1 = (hx+emx1-emn1-2.0*buffx+.00001) / hx
      iside2 = (hy+emx2-emn2-2.0*buffy+.00001) / hy
      gpall  = iside1 * iside2
      usage  = rn / gpall
c
      return
      end
