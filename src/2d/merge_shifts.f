      subroutine merge_shifts(irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,
     . lstgrd, numHoods) 

! !   INPUTS:
! !   REQUIRED GRID DATA
! !
! !   irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,lstgrd
! !
! !   REQUIRED MERGING DATA
! !   volMerge, xcmerge, ycmerge, numHoods,
! !   mi, mj, ncount
! !
! !   OUTPUT
! !   qmshifts
! !
! !   NOTE qmshifts(1,k) = shift for xx
! !        qmshifts(2,k) = shift for xy
! !        qmshifts(3,k) = shift for yy
! !   Called dmergeshift in AG code


      use amr_module
      implicit double precision (a-h, o-z)

      dimension irr(mitot,mjtot), numHoods(mitot,mjtot)

      logical IS_OUTSIDE

      include "quadrature.i"

      IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                   y .lt. ylower .or. y .gt. yupper)



      ! recompute the shifts
      do 719 j = 1, mjtot    ! iterate over each cell in the grid
      do 719 i = 1, mitot    ! iterate over each cell in the grid
         k = irr(i,j)


         if (k .eq. -1 .or. k .eq. lstgrd) cycle ! skip solid and whole cells
         call getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)
         if (IS_OUTSIDE(xc,yc)) cycle

          shiftmxx = 0.d0
          shiftmxy = 0.d0
          shiftmyy = 0.d0

          x0 = xcentMerge(k)
          y0 = ycentMerge(k)
          arr = volMerge(k)

          do 340 ic = 0, ncount(k) ! my count 1 less than AG
!           AG data structure)
!           ioff = mi(k,ic)
!           joff = mj(k,ic)

!           my data structure
            if (ic .eq. 0) then
               ioff = i  ! first do cell itself
               joff = j
            else ! now pointing to nbor
               ioff = iidx(ic,k)
               joff = jidx(ic,k)
            endif

            koff = irr(ioff, joff)

            if(koff .eq. lstgrd) then ! full cell
               itri = 2
               call makep(poly(1,1,lstgrd),ioff,joff,xlow,ylow,hx,hy)
            else                      ! cut cells
              ivert = 1
              do 239 while (poly(ivert+1,1,koff) .ne. -11.)
                ivert = ivert + 1
 239          continue
              itri = ivert - 3
            endif
            idx1 = 1
            do 219 it = 1, itri ! for each  triangle
             idx2 = it + 1
             idx3 = it + 2

             x1 = poly(idx1,1,koff)
             y1 = poly(idx1,2,koff)

             x2 = poly(idx2,1,koff)
             y2 = poly(idx2,2,koff)

             x3 = poly(idx3,1,koff)
             y3 = poly(idx3,2,koff)

             ! YOU MUST WEIGHT THE TRIANGLE AREA BY NUMHOODS
             artri = triangle_area(x1, x2, x3,  y1, y2, y3) / 
     .               numHoods(ioff, joff)

             do 229 itq = 1,ntriquad

                 xval = x1 * rtri(itq) + x2 * stri(itq)
     .               +  x3 * (1.d0-rtri(itq)-stri(itq))
                 yval = y1 * rtri(itq) + y2 * stri(itq)
     .               +  y3 * (1.d0-rtri(itq)-stri(itq))

                 shiftmxx = shiftmxx + (artri/arr)*wtri(itq)*
     .                       (xval-x0)**2 /(hx**2)
                 shiftmxy = shiftmxy + (artri/arr)*wtri(itq)*
     .                       (xval-x0)*(yval-y0)/(hx * hy)
                 shiftmyy = shiftmyy + (artri/arr)*wtri(itq)*
     .                       (yval-y0)**2 /(hy**2)

 229         continue ! for each quadrature point on each triangle
 219        continue ! for each triangle

 340   continue  ! iterate over each cell in the neighborhood

        ! save for use in merge nhood reconstruction
        qmshifts(1,k)  = shiftmxx
        qmshifts(2,k)  = shiftmxy
        qmshifts(3,k)  = shiftmyy

 719  continue ! iterate over each cell in the grid

      return
      end
