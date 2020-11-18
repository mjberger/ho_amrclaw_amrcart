c
c ------------------------------------------------------------
c
       subroutine computeCubicShifts(shiftxxx,shiftxxy,
     &                   shiftxyy,shiftyyy,hx,hy,kirr)

       use amr_module

       implicit none
       include "quadrature.i"

       double precision shiftxxx,shiftxxy
       double precision shiftxyy,shiftyyy
       double precision x0,y0,arr,artri
       double precision x1,y1,x2,y2,x3,y3
       double precision xval,yval,hx,hy
       double precision triangle_area
       integer idx1,idx2,idx3,ivert,it,itq,itri
       integer kirr

c      compute integrals of cubics over irreular cells
c      (if kirr is regular then integrals are zero and this 
c       routine not called).

       shiftxxx = 0.d0
       shiftxxy = 0.d0
       shiftxyy = 0.d0
       shiftyyy = 0.d0

       x0 = xcirr(kirr)
       y0 = xcirr(kirr)
       arr = ar(kirr)

       ivert = 1
       do while (poly(ivert+1,1,kirr) .ne. -11)
          ivert = ivert + 1
       end do

       itri = ivert - 3
       idx1 = 1
       x1 = poly(idx1,1,kirr)
       y1 = poly(idx1,2,kirr)

       do it = 1, itri
          idx2 = it + 1
          idx3 = it + 3
          x2 = poly(idx2,1,kirr)
          y2 = poly(idx2,2,kirr)
          x3 = poly(idx3,1,kirr)
          y3 = poly(idx3,2,kirr)
          artri = triangle_area(x1,x2,x3,y1,y2,y3)

          do  itq = 1, ntriquad
            xval = x1 * rtri(itq) + x2 * stri(itq)
     &          +  x3 * (1.d0-rtri(itq)-stri(itq))
            yval = y1 * rtri(itq) + y2 * stri(itq)
     &          +  y3 * (1.d0-rtri(itq)-stri(itq))

            shiftxxx = shiftxxx + (artri/arr)*wtri(itq) *
     &                            (xval-x0)**3 / (hx**3)
             shiftxxy = shiftxxy + (artri/arr)*wtri(itq) *
     &                            (yval-y0)*(xval-x0)**2 / (hy*hx**2)
             shiftxyy = shiftxyy + (artri/arr)*wtri(itq) *
     &                            (xval-x0)*(yval-y0)**2 / (hx*hy**2)
             shiftyyy =shiftyyy + (artri/arr)*wtri(itq) *
     &                           (yval-y0)**3 / (hy**3)
          end do  ! end loop over quadrature points
       end do ! end loop over each triangle in polygon

       return
       end
