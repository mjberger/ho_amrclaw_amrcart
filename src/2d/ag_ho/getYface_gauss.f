c
c -----------------------------------------------------------
c
       subroutine getYface_gauss(i,jcol,xface,yface,irr,mitot,mjtot,
     .                     xlow,ylow,dx,dy,lstgrd,num, missing)

       implicit double precision (a-h,o-z)
       include "cirr.i"
       include "./quadrature.i"
       dimension irr(mitot,mjtot)
       logical IS_CUT, missing

       IS_CUT(k) = ((k .ne. lstgrd) .and. (k .ne. -1))

c need midpoint of y face between 2 cells. 
c    if bottom is irregular need top face
c    (or if top is irregular need bottom face)
c    if at least one is regular, easy case
c    if one is solid, return face isnt used.
c
      missing = .false.

      w1 = 0.5d0*(1-dsqrt(3.d0)/3.d0)
      w2 = 0.5d0*(1+dsqrt(3.d0)/3.d0)
      if (irr(i,jcol) .eq. lstgrd .or. irr(i,jcol+1) .eq. lstgrd) then  ! face must be regular
             yface = ylow + dfloat(jcol)*dy

c             xface = xlow + (dfloat(i)-.5d0)*dx
             x1 = xlow + (dfloat(i)+0.d0)*dx
             x2 = xlow + (dfloat(i)-1.d0)*dx
             xface = 0.5d0*(1+rline(num)) * x1
     .             + 0.5d0*(1-rline(num)) * x2
!             if(num .eq. 1) then
!              xface = w1*x1 + (1.d0-w1)*x2
!             else
!              xface = w2*x1 + (1.d0-w2)*x2
!             endif
      else if (IS_CUT(irr(i,jcol)) .and. IS_CUT(irr(i,jcol+1))) then !use bottom cell to find top face
              k = irr(i,jcol)
              yreg =  ylow + jcol*dy
              do 10 kside = 1, 6
                if (poly(kside+2,1,k) .eq. -11) then
!                   write(*,*)"missing side in getYface for ",i,jcol
!                   stop
                    missing = .true.
                    return
                endif
                x1 = poly(kside,1,k)
                y1 = poly(kside,2,k)
                x2 = poly(kside+1,1,k)
                y2 = poly(kside+1,2,k)
                if (x1 .eq. x2) go to 10  ! not a y face
                if (y1 .ne. y2)  then
                  write(*,*)"screwy side in getYface for cell ",i,jcol
                  stop
                endif

                if (dabs(y1-yreg) .lt. dy/1000.)  then
                   yface = y1
                   xface = 0.5d0*(1+rline(num)) * x1
     .                   + 0.5d0*(1-rline(num)) * x2
!                   if(num .eq. 1) then
!                    xface = w1*x1 + (1.d0-w1)*x2
!                   else
!                    xface = w2*x1 + (1.d0-w2)*x2     ! midpoint of face
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
