      subroutine merge_shifts_ho(irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,
     . lstgrd, ssw)
  !   INPUTS:
  !   REQUIRED GRID DATA
  !
  !   irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,lstgrd
  !
  !
  !   OUTPUT
  !   volmerge_ho, xcentmerge_ho, ycentmerge_ho, qmshifts_ho

  !   NOTE qmshifts_ho(1,kirr) = shift for xx
  !        qmshifts_ho(2,kirr) = shift for xy
  !        qmshifts_ho(3,kirr) = shift for yy
  !        qmshifts_ho(4,kirr) = shift for xxx
  !        qmshifts_ho(5,kirr) = shift for xxy
  !        qmshifts_ho(6,kirr) = shift for xyy
  !        qmshifts_ho(7,kirr) = shift for yyy


      implicit double precision (a-h, o-z)
      include "cirr.i"
      include "quadrature.i"
      dimension irr(mitot,mjtot)



          qmshifts_ho = 0.d0
          volMerge_ho = 0.d0
          xcentmerge_ho = 0.d0
          ycentmerge_ho = 0.d0


      ar_ho(lstgrd) = hx*hy
      ! compute xcmerge, ycmerge, volmerge
      do 30 j = lwidth+1, mjtot-lwidth
      do 30 i = lwidth+1, mitot-lwidth
         k = irr(i,j)

         if(k .eq. -1 .or. k .eq. lstgrd) goto 30

!         tv = 0.d0
!         xct = 0.d0
!         yct = 0.d0

         do 20 ic = 1, ncount(k)
            icurr = iidx(k,ic)
            jcurr = jidx(k,ic)
            kcurr = irr(icurr,jcurr)

            if (kcurr .eq. lstgrd) then
               xc = xlow + (icurr-0.5d0)*hx
               yc = ylow + (jcurr-0.5d0)*hy
            else
               xc = xcirr_ho(kcurr)
               yc = ycirr_ho(kcurr)
            endif


            volMerge_ho(k) = volMerge_ho(k)
     .                    +    ar_ho(kcurr) / numHoods(icurr,jcurr)
            xcentmerge_ho(k) = xcentmerge_ho(k)
     .                    + xc*ar_ho(kcurr) / numHoods(icurr,jcurr)
            ycentmerge_ho(k) = ycentmerge_ho(k)
     .                    + yc*ar_ho(kcurr) / numHoods(icurr,jcurr)


!            tv = tv  +    ar(kcurr) / numHoods(icurr,jcurr)
!            xct = xct+ xc*ar(kcurr) / numHoods(icurr,jcurr)
!            yct = yct+ yc*ar(kcurr) / numHoods(icurr,jcurr)

 20      continue

         xcentmerge_ho(k) = xcentmerge_ho(k) / volMerge_ho(k)
         ycentmerge_ho(k) = ycentmerge_ho(k) / volMerge_ho(k)

!          xct = xct / tv
!         yct = yct / tv
 30   continue


      ! recompute the shifts
      do 719 j = lwidth+1, mjtot-lwidth ! iterate over each cell in the grid
      do 719 i = lwidth+1, mitot-lwidth ! iterate over each cell in the grid
       kirr = irr(i,j)




       if (kirr .eq. -1 .or. kirr .eq. lstgrd) cycle ! skip solid and whole cells

           shiftmxx  = 0.d0
           shiftmxy  = 0.d0
           shiftmyy  = 0.d0
           shiftmxxx = 0.d0
           shiftmxxy = 0.d0
           shiftmxyy = 0.d0
           shiftmyyy = 0.d0

           x0 = xcentmerge_ho(kirr)
           y0 = ycentmerge_ho(kirr)

           arr = volmerge_ho(kirr)


      do 340 ic = 1, ncount(kirr) ! iterate over each cell in the neighborhood
            ioff = iidx(kirr,ic)
            joff = jidx(kirr,ic)

            koff = irr(ioff, joff)

      if(koff .eq. lstgrd) then ! full cell
            x1 = xlow + (dfloat(ioff)-1.d0)*hx
            y1 = ylow + (dfloat(joff)-1.d0)*hy

            x3 = xlow + (dfloat(ioff)-0.d0)*hx
            y3 = ylow + (dfloat(joff)-0.d0)*hy

            ! YOU MUST WEIGHT THE QUADRILATERAL AREA BY NUMHOODS
            arcell = hx * hy / numHoods(ioff, joff)

      do 11 iq = 1, nquadquad
          xval = x1 * (1.d0 - rquad(iq))/2. + x3 * (1.d0 + rquad(iq))/2.
          yval = y1 * (1.d0 - squad(iq))/2. + y3 * (1.d0 + squad(iq))/2.

      shiftmxx = shiftmxx + (arcell/arr)*wquad(iq)*((xval-x0)/hx)**2

      shiftmxy = shiftmxy
     .       + (arcell/arr)*wquad(iq)*((xval-x0)/hx)*((yval-y0)/hy)

      shiftmyy = shiftmyy + (arcell/arr)*wquad(iq)*((yval-y0)/hy)**2


      shiftmxxx = shiftmxxx
     .       + (arcell/arr)*wquad(iq)*(xval-x0)**3 /(hx**3)
      shiftmxxy = shiftmxxy
     .       + (arcell/arr)*wquad(iq)*(yval-y0)*(xval-x0)**2/(hx*hx*hy)
      shiftmxyy = shiftmxyy
     .       + (arcell/arr)*wquad(iq)*(xval-x0)*(yval-y0)**2/(hx*hy*hy)
      shiftmyyy = shiftmyyy
     .       + (arcell/arr)*wquad(iq)*(yval-y0)**3/(hy**3)
 11         continue ! for each quadrature point on full cell




      else                      ! cut cells




          ivert = 1
          do 239 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
 239     continue


          itri = ivert - 3
          idx1 = 1

          do 219 it = 1, itri-1 ! for each  triangle except the last one
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            ! YOU MUST WEIGHT THE TRIANGLE AREA BY NUMHOODS
            artri = tri_area(x1, x2, x3,
     .                            y1, y2, y3) / numHoods(ioff, joff)

            do 229 itq = 1,ntriquad

                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

      shiftmxx = shiftmxx + (artri/arr)*wtri(itq)*(xval-x0)**2 /(hx**2)
      shiftmxy = shiftmxy
     .          + (artri/arr)*wtri(itq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftmyy = shiftmyy + (artri/arr)*wtri(itq)*(yval-y0)**2 /(hy**2)


      shiftmxxx = shiftmxxx
     .                   + (artri/arr)*wtri(itq)*(xval-x0)**3 /(hx**3)
      shiftmxxy = shiftmxxy
     .        + (artri/arr)*wtri(itq)*(yval-y0)*(xval-x0)**2/(hy*hx**2)
      shiftmxyy = shiftmxyy
     .        + (artri/arr)*wtri(itq)*(xval-x0)*(yval-y0)**2 /(hx*hy**2)
      shiftmyyy = shiftmyyy
     .        + (artri/arr)*wtri(itq)*(yval-y0)**3 /(hy**3)

 229      continue ! for each quadrature point on each triangle
 219      continue ! for each triangle
      if(ssw .eq. 2  .or. ssw .eq. -2) then
      x1 = poly(ivert-2,1,koff)
      y1 = poly(ivert-2,2,koff)

      x2 = bdry(1,1,koff)
      y2 = bdry(1,2,koff)

      x3 = bdry(3,1,koff)
      y3 = bdry(3,2,koff)

      x4 = bdry(2,1,koff)
      y4 = bdry(2,2,koff)

      do nq = 1,ntriquad_ho
      artri = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)/ numHoods(ioff, joff)

      xval = rs2xy_2(x1,x2,x3,x4,nq)
      yval = rs2xy_2(y1,y2,y3,y4,nq)

      shiftmxx = shiftmxx
     .          + (artri/arr)*wtri_ho(nq)*(xval-x0)**2 /(hx**2)
      shiftmxy = shiftmxy
     .          + (artri/arr)*wtri_ho(nq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftmyy = shiftmyy
     .          + (artri/arr)*wtri_ho(nq)*(yval-y0)**2 /(hy**2)


      shiftmxxx = shiftmxxx +
     .     (artri/arr)*wtri_ho(nq)*(xval-x0)**3 /(hx**3)
      shiftmxxy = shiftmxxy +
     .     (artri/arr)*wtri_ho(nq)*(yval-y0)*(xval-x0)**2 /(hy*hx**2)
      shiftmxyy = shiftmxyy +
     .     (artri/arr)*wtri_ho(nq)*(xval-x0)*(yval-y0)**2 /(hx*hy**2)
      shiftmyyy =shiftmyyy
     .   + (artri/arr)*wtri_ho(nq)*(yval-y0)**3 /(hy**3)

      enddo
      elseif(ssw .eq. 3 .or. ssw .eq. -3) then
      x1 = poly(ivert-2,1,koff)
      y1 = poly(ivert-2,2,koff)

      x2 = bdry(1,1,koff)
      y2 = bdry(1,2,koff)

      x3 = bdry(4,1,koff)
      y3 = bdry(4,2,koff)

      x4 = bdry(2,1,koff)
      y4 = bdry(2,2,koff)

      x5 = bdry(3,1,koff)
      y5 = bdry(3,2,koff)

      do nq = 1,ntriquad_ho
      artri = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)/ numHoods(ioff,joff)

      xval = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yval = rs2xy_3(y1,y2,y3,y4,y5,nq)

      shiftmxx = shiftmxx
     .       + (artri/arr)*wtri_ho(nq)*(xval-x0)**2 /(hx**2)
      shiftmxy = shiftmxy
     .       + (artri/arr)*wtri_ho(nq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftmyy = shiftmyy
     .       + (artri/arr)*wtri_ho(nq)*(yval-y0)**2 /(hy**2)


      shiftmxxx = shiftmxxx
     .     + (artri/arr)*wtri_ho(nq)*(xval-x0)**3 /(hx**3)
      shiftmxxy = shiftmxxy
     .     + (artri/arr)*wtri_ho(nq)*(yval-y0)*(xval-x0)**2 /(hy*hx**2)
      shiftmxyy = shiftmxyy
     .     + (artri/arr)*wtri_ho(nq)*(xval-x0)*(yval-y0)**2 /(hx*hy**2)
      shiftmyyy =shiftmyyy
     .     + (artri/arr)*wtri_ho(nq)*(yval-y0)**3 /(hy**3)

      end do

      endif


      endif
 340     continue  ! iterate over each cell in the neighborhood


           qmshifts_ho(1,kirr)  = shiftmxx
           qmshifts_ho(2,kirr)  = shiftmxy
           qmshifts_ho(3,kirr)  = shiftmyy
           qmshifts_ho(4,kirr)  = shiftmxxx
           qmshifts_ho(5,kirr)  = shiftmxxy
           qmshifts_ho(6,kirr)  = shiftmxyy
           qmshifts_ho(7,kirr)  = shiftmyyy

 719      continue ! iterate over each cell in the grid

      end


 !    computes the area of the triangle with the vertices
 !    (x1,y1), (x2,y2), (x3,y3)
      double precision function tri_area(x1, x2, x3, y1, y2, y3)
      implicit double precision (a-h,o-z)

      tri_area = dabs(0.5d0 * ( (x2 - x1) * (y3 - y1)
     .                    - (x3 - x1) * (y2 - y1) ))

      return
      end
