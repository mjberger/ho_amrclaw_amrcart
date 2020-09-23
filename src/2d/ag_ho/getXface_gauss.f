c
c -----------------------------------------------------------
c
       subroutine getXface_gauss(irow,j,xface,yface,irr,mitot,mjtot,
     .                     xlow,ylow,dx,dy,lstgrd,num, missing)

       implicit double precision (a-h,o-z)
       include "cirr.i"
       include "./quadrature.i"
       dimension irr(mitot,mjtot)
       logical IS_CUT, missing

       IS_CUT(k) = ((k .ne. lstgrd) .and. (k .ne. -1))

c need midpoint of x face between 2 cells. 
c    if bottom is irregular need top face
c    (or if top is irregular need bottom face)
c    if at least one is regular, easy case
c    if one is solid, return face isnt used.
c
      missing = .false.

      w1 = 0.5d0*(1-dsqrt(3.d0)/3.d0)
      w2 = 0.5d0*(1+dsqrt(3.d0)/3.d0)
      if (irr(irow,j) .eq. lstgrd .or. irr(irow+1,j) .eq. lstgrd) then  ! face must be regular
             xface = xlow + dfloat(irow)*dx
c             yface = ylow + (dfloat(j)-.5d0)*dy
             y1 = ylow + (dfloat(j)+0.d0)*dy
             y2 = ylow + (dfloat(j)-1.d0)*dy

             yface = 0.5d0*(1+rline(num)) * y1
     .             + 0.5d0*(1-rline(num)) * y2

!             if(num .eq. 1) then
!              yface = w1*y1 + (1.d0-w1)*y2
!             else
!              yface = w2*y1 + (1.d0-w2)*y2     ! midpoint of face
!             endif

      else if (IS_CUT(irr(irow,j)) .and. IS_CUT(irr(irow+1,j))) then !use bottom cell to find top face
              k = irr(irow,j)
              xreg =  xlow + irow*dx
              do 10 kside = 1, 6
                if (poly(kside+2,1,k) .eq. -11) then
!                   write(*,*)"missing side in getXface for ",irow,j
!                   stop
                    missing = .true.
                    return
                endif
                x1 = poly(kside,1,k)
                y1 = poly(kside,2,k)
                x2 = poly(kside+1,1,k)
                y2 = poly(kside+1,2,k)
                if (y1 .eq. y2) go to 10  ! not an x face
                if (x1 .ne. x2)  then
                  write(*,*)"screwy side in getXface for cell ",irow,j
                  stop
                endif

                if (dabs(x1-xreg) .lt. dx/1000.)  then
                   xface = x1
                   yface = 0.5d0*(1+rline(num)) * y1
     .                   + 0.5d0*(1-rline(num)) * y2
!                   if(num .eq. 1) then
!                    yface = w1*y1 + (1.d0-w1)*y2
!                   else
!                    yface = w2*y1 + (1.d0-w2)*y2     ! midpoint of face
!                   endif
                   go to 99
                endif
 10           continue

      else  ! one or both must be solid, no need for real computation

          xface = 0.d0
          yface = 0.d0
          missing = .true.
      endif

 99   return
      end
