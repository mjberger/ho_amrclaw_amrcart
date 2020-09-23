c
c -----------------------------------------------------------
c
       subroutine getCellCentroid(i,j,xcent,ycent,irr,mitot,mjtot,
     .                            xlow,ylow,dx,dy,lstgrd)

       implicit double precision (a-h,o-z)
       include "cirr.i"
       dimension irr(mitot,mjtot)

       k = irr(i,j)
       if (k .eq. lstgrd) then   ! regular cell
            xcent = xlow + (i-.5d0)*dx 
            ycent = ylow + (j-.5d0)*dy
       else if (k .eq. -1) then   ! solid cell. set to 0.
            xcent = 0.d0
            ycent = 0.d0
       else  ! this has been precomputed
          xcent = xcirr(k)
          ycent = ycirr(k)
       endif

       return
       end
