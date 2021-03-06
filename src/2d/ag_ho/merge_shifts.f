      subroutine merge_shifts(irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,
     . lstgrd,
     . vmerge, xcmerge, ycmerge, nhoods,
     . mi, mj, nc,
     . qmshifts)
  !   INPUTS:
  !   REQUIRED GRID DATA
  !
  !   irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,lstgrd
  !
  !   REQUIRED MERGING DATA
  !   vmerge, xcmerge, ycmerge, nhoods,
  !   mi, mj, nc
  !
  !   OUTPUT
  !   qmshifts
  !
  !   NOTE qmshifts(1,kirr) = shift for xx
  !        qmshifts(2,kirr) = shift for xy
  !        qmshifts(3,kirr) = shift for yy


      implicit double precision (a-h, o-z)
      include "cirr.i"
      dimension irr(mitot,mjtot), nhoods(5200,5200)
      dimension vmerge(irrsize)
      dimension xcmerge(irrsize), ycmerge(irrsize)
      dimension mi(irrsize,10),mj(irrsize,10),nc(irrsize)
      dimension qmshifts(3,irrsize)
      dimension rtri(3), stri(3), wtri(3)
      dimension rquad(4), squad(4), wquad(4)

      qmshifts = 0.d0

      ! QUADRATURE RULE ON TRIANGLES
      ! This quadrature rule integrates quadratics exactly
      rtri = (/ 1.d0/6.d0, 2.d0/3.d0,1.d0/6.d0 /)
      stri = (/ 1.d0/6.d0, 1.d0/6.d0,2.d0/3.d0 /)
      wtri = (/ 1.d0/3.d0, 1.d0/3.d0,1.d0/3.d0 /)
      ntriquad = 3

      ! QUADRATURE RULE ON CARTESIAN CELLS
      ! This quadrature rule integrates cubics exactly
      rquad = (/-dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0),
     .           dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0)/)
      squad = (/ -dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0),
     .            dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0) /)
      wquad = (/ 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0 /)
      nquadquad = 4


      ! recompute the shifts
      do 719 j = lwidth+1, mjtot-lwidth ! iterate over each cell in the grid
      do 719 i = lwidth+1, mitot-lwidth ! iterate over each cell in the grid
       kirr = irr(i,j)




       if (kirr .eq. -1 .or. kirr .eq. lstgrd) cycle ! skip solid and whole cells

           shiftmxx = 0.d0
           shiftmxy = 0.d0
           shiftmyy = 0.d0

           x0 = xcmerge(kirr)
           y0 = ycmerge(kirr)

           arr = vmerge(kirr)


      do 340 ic = 1, nc(kirr) ! iterate over each cell in the neighborhood
            ioff = mi(kirr,ic)
            joff = mj(kirr,ic)

            koff = irr(ioff, joff)

      if(koff .eq. lstgrd) then ! full cell
            x1 = xlow + (dfloat(ioff)-1.d0)*hx
            y1 = ylow + (dfloat(joff)-1.d0)*hy

            x3 = xlow + (dfloat(ioff)-0.d0)*hx
            y3 = ylow + (dfloat(joff)-0.d0)*hy

            ! YOU MUST WEIGHT THE QUADRILATERAL AREA BY NUMHOODS
            arcell = hx * hy / nhoods(ioff, joff)

      do 11 iq = 1, nquadquad
          xval = x1 * (1.d0 - rquad(iq))/2. + x3 * (1.d0 + rquad(iq))/2.
          yval = y1 * (1.d0 - squad(iq))/2. + y3 * (1.d0 + squad(iq))/2.

      shiftmxx = shiftmxx + (arcell/arr)*wquad(iq)*((xval-x0)/hx)**2

      shiftmxy = shiftmxy
     .       + (arcell/arr)*wquad(iq)*((xval-x0)/hx)*((yval-y0)/hy)

      shiftmyy = shiftmyy + (arcell/arr)*wquad(iq)*((yval-y0)/hy)**2
 11         continue ! for each quadrature point on full cell




      else                      ! cut cells




          ivert = 1
          do 239 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
 239     continue


          itri = ivert - 3
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
            artri = triangle_area(x1, x2, x3,
     .                            y1, y2, y3) / nhoods(ioff, joff)

            do 229 itq = 1,ntriquad

                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

      shiftmxx = shiftmxx + (artri/arr)*wtri(itq)*(xval-x0)**2 /(hx**2)
      shiftmxy = shiftmxy
     .          + (artri/arr)*wtri(itq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftmyy = shiftmyy + (artri/arr)*wtri(itq)*(yval-y0)**2 /(hy**2)

 229      continue ! for each quadrature point on each triangle
 219      continue ! for each triangle

      endif
 340     continue  ! iterate over each cell in the neighborhood


           qmshifts(1,kirr)  = shiftmxx
           qmshifts(2,kirr)  = shiftmxy
           qmshifts(3,kirr)  = shiftmyy

 719      continue ! iterate over each cell in the grid

      end


 !    computes the area of the triangle with the vertices
 !    (x1,y1), (x2,y2), (x3,y3)
      double precision function triangle_area(x1, x2, x3, y1, y2, y3)
      implicit double precision (a-h,o-z)

      triangle_area = abs(0.5d0 * ( (x2 - x1) * (y3 - y1)
     .                    - (x3 - x1) * (y2 - y1) ))

      return
      end
