!       todo: separate the p = 1 and p = 2 reconstructions ...

c
c------------------------------------------------------------
c
       subroutine slopes(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
       common /order2/ ssw, quad, nolimiter
       logical IS_GHOST
       dimension ata(5,5), rhs(5,nvar), a(5)
       dimension b(5), w(5,nvar)
       dimension inuf(mitot,mjtot)
       include "cirr.i"
       include "./quadrature.i"



      qx = 0.d0
      qy = 0.d0
      qxx = 0.d0
      qxy = 0.d0
      qyy = 0.d0
      inuf = 0

        if(ssw .eq. 0) return

       do 33 i = lwidth, mitot-lwidth+1
       do 33 j = lwidth, mjtot-lwidth+1
            k = irr(i,j)
            if (k .ne. lstgrd) goto 33
            inuf(i,j) = 1

            do 333 ioff = -1,1
            do 333 joff = -1,1
                if(irr(i+ioff, j+joff) .ne. lstgrd) then
                    inuf(i,j) = 0
                endif
 333        continue
 33    continue


       do 34 i = lwidth, mitot-lwidth+1
       do 34 j = lwidth, mjtot-lwidth+1
            if(inuf(i,j) .eq. 0) cycle

            if(ssw .eq. 2) then
            qyy(i,j,:) =  (qp(i,j+1,:)-2*qp(i,j,:)+qp(i,j-1,:))/hy**2
            qxy(i,j,:) =  (qp(i+1,j+1,:)-qp(i+1,j-1,:)-qp(i-1,j+1,:)
     .                     +qp(i-1,j-1,:))/(4*hx*hy)
            qxx(i,j,:) =  (qp(i+1,j,:)-2*qp(i,j,:)+qp(i-1,j,:))/hx**2
            endif

            qy(i,j,:)  =  (qp(i,j+1,:) - qp(i,j-1,:))/(2*hy)
            qx(i,j,:)  =  (qp(i+1,j,:) - qp(i-1,j,:))/(2*hx)
 34    continue

       do 35 i = lwidth, mitot-lwidth+1
       do 35 j = lwidth, mjtot-lwidth+1
            kirr = irr(i,j)
            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
            if(kirr .eq. lstgrd) then
                x0 = xlow + (i-.5d0)*hx
                y0 = ylow + (j-.5d0)*hy
            else
                x0 = xcirr(kirr)
                y0 = ycirr(kirr)
            endif

            ata = 0.d0
            rhs = 0.d0



            do 394 moff = 1,3

            do 334 ioff = -moff,moff
            do 334 joff = -moff,moff
                idxmax = max(abs(ioff),abs(joff))
                if(idxmax .ne. moff) goto 334


                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 334





                if(koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*hx
                    yoff = ylow + (j+joff-.5d0)*hy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

                deltax = xoff - x0
                deltay = yoff - y0

                a(1) = deltax/hx
                a(2) = deltay/hy

                ata(1,1) = ata(1,1) + a(1) * a(1)
                ata(1,2) = ata(1,2) + a(1) * a(2)
                ata(2,2) = ata(2,2) + a(2) * a(2)
                rhs(1,:) = rhs(1,:) + a(1) *
     .                     (qp(i+ioff,j+joff,:) - qp(i,j,:))
                rhs(2,:) = rhs(2,:) + a(2) *
     .                     (qp(i+ioff,j+joff,:) - qp(i,j,:))




                if( ssw .eq. 1) then
                    cycle
                endif



         ! load the offsets
         a(3) = -poly(8 ,1,kirr)
         a(4) = -poly(9 ,1,kirr)
         a(5) = -poly(10,1,kirr)



        if(koff .ne. lstgrd) then
          arr = ar(koff)
          ivert = 1
          do 20 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
  20      continue
           itri = ivert - 3
        else
          arr = hx*hy
          itri = 2
          poly(1,1,koff) = xlow + (dfloat(i+ioff)-1.d0)*hx
          poly(1,2,koff) = ylow + (dfloat(j+joff)-1.d0)*hy

          poly(2,1,koff) = xlow + (dfloat(i+ioff)-0.d0)*hx
          poly(2,2,koff) = ylow + (dfloat(j+joff)-1.d0)*hy

          poly(3,1,koff) = xlow + (dfloat(i+ioff)-0.d0)*hx
          poly(3,2,koff) = ylow + (dfloat(j+joff)-0.d0)*hy

          poly(4,1,koff) = xlow + (dfloat(i+ioff)-1.d0)*hx
          poly(4,2,koff) = ylow + (dfloat(j+joff)-0.d0)*hy
        endif







          idx1 = 1
          do 21 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 22 itq = 1,ntriquad
                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

                a(3) = a(3) + (artri/arr) * wtri(itq) *( (xval-x0)**2)
     .                                                / (hx**2)

                a(4) = a(4) + (artri/arr) * wtri(itq) *( (xval-x0)/hx )
     .                                                *( (yval-y0)/hy )

                a(5) = a(5) + (artri/arr) * wtri(itq) *( (yval-y0)**2 )
     .                                                / (hy**2)
  22        continue ! for each quadrature point on each triangle
  21      continue ! for each triangle


            ata(1,3) = ata(1,3) + a(1)*a(3)
            ata(1,4) = ata(1,4) + a(1)*a(4)
            ata(1,5) = ata(1,5) + a(1)*a(5)
            ata(2,3) = ata(2,3) + a(2)*a(3)
            ata(2,4) = ata(2,4) + a(2)*a(4)
            ata(2,5) = ata(2,5) + a(2)*a(5)
            ata(3,3) = ata(3,3) + a(3)*a(3)
            ata(3,4) = ata(3,4) + a(3)*a(4)
            ata(3,5) = ata(3,5) + a(3)*a(5)
            ata(4,4) = ata(4,4) + a(4)*a(4)
            ata(4,5) = ata(4,5) + a(4)*a(5)
            ata(5,5) = ata(5,5) + a(5)*a(5)

            rhs(3,:) = rhs(3,:) + a(3) *
     .                     (qp(i+ioff,j+joff,:) - qp(i,j,:))
            rhs(4,:) = rhs(4,:) + a(4) *
     .                     (qp(i+ioff,j+joff,:) - qp(i,j,:))
            rhs(5,:) = rhs(5,:) + a(5) *
     .                     (qp(i+ioff,j+joff,:) - qp(i,j,:))

 334         continue


                 if( ssw .eq. 1) then
                    ! check for singularity of A
                     if( ata(1,1) < 1d-8 ) cycle
                     c11 = dsqrt(ata(1,1))
                     c12 = ata(1,2)/c11


                     if( ata(2,2)-c12**2 < 1d-8 ) cycle
                     c22 = dsqrt(ata(2,2)-c12**2)

                     ! quit the double for loop if the matrix is nonsingular.
                     goto 435
                endif






            ! check for singularity of A
             if( ata(1,1) < 1d-8 ) cycle
             c11 = dsqrt(ata(1,1))
             c12 = ata(1,2)/c11
             c13 = ata(1,3)/c11
             c14 = ata(1,4)/c11
             c15 = ata(1,5)/c11

             if( ata(2,2)-c12**2 < 1d-8 ) cycle
             c22 = dsqrt(ata(2,2)-c12**2)
             c23 = (ata(2,3)-c12*c13)/c22
             c24 = (ata(2,4)-c12*c14)/c22
             c25 = (ata(2,5)-c12*c15)/c22

             if( ata(3,3)-c13**2 - c23**2 < 1d-8 ) cycle
             c33 = dsqrt(ata(3,3)-c13**2 - c23**2)
             c34 = (ata(3,4)-c13*c14-c23*c24)/c33
             c35 = (ata(3,5)-c13*c15-c23*c25)/c33

             if( ata(4,4)-c14**2-c24**2-c34**2 < 1d-8 ) cycle
             c44 = dsqrt(ata(4,4)-c14**2-c24**2-c34**2)
             c45 = (ata(4,5)-c14*c15-c24*c25-c34*c35)
     &            /c44

            if( ata(5,5)-c15**2-c25**2-c35**2-c45**2 < 1d-8 ) cycle
             c55 = dsqrt(ata(5,5)-c15**2-
     &                      c25**2-c35**2-c45**2)

             ! quit the double for loop if the matrix is nonsingular.
             goto 335

 394         continue





      ! SOLVE THE NORMAL EQUATIONS.

 435         continue ! quit the double for loop if the matrix is nonsingular.
             c11 = sqrt(ata(1,1))
             c12 = ata(1,2) / c11
             c22 = sqrt(ata(2,2) - c12**2)

             do mm = 1, nvar
               b(1) = rhs(1,mm) / c11
               b(2) = (rhs(2,mm) - c12*b(1))/c22
               qy(i,j,mm) = b(2) / c22
               qx(i,j,mm) = (b(1) - c12*qy(i,j,mm))/c11
             end do

             qx(i,j,:) = qx(i,j,:) / hx
             qy(i,j,:) = qy(i,j,:) / hy
!                  if(i .eq. 30 .and. j .eq. 5) then
!          print *, "here"
!          endif
             cycle

 335         continue ! quit the double for loop if the matrix is nonsingular.

             c11 = dsqrt(ata(1,1))
             c12 = ata(1,2)/c11
             c13 = ata(1,3)/c11
             c14 = ata(1,4)/c11
             c15 = ata(1,5)/c11

             c22 = dsqrt(ata(2,2)-c12**2)
             c23 = (ata(2,3)-c12*c13)/c22
             c24 = (ata(2,4)-c12*c14)/c22
             c25 = (ata(2,5)-c12*c15)/c22

             c33 = dsqrt(ata(3,3)-c13**2 - c23**2)
             c34 = (ata(3,4)-c13*c14-c23*c24)/c33
             c35 = (ata(3,5)-c13*c15-c23*c25)/c33

             c44 = dsqrt(ata(4,4)-c14**2-c24**2-c34**2)
             c45 = (ata(4,5)-c14*c15-c24*c25-c34*c35)
     &            /c44
             c55 = dsqrt(ata(5,5)-c15**2-
     &                      c25**2-c35**2-c45**2)


             do mm = 1, nvar
                b(1) = rhs(1,mm) / c11
                b(2) = (rhs(2,mm) - c12*b(1)) / c22
                b(3) = (rhs(3,mm) - c13*b(1)-c23*b(2))/c33
                b(4) = (rhs(4,mm) - c14*b(1)-c24*b(2)
     &                    -c34*b(3))/c44
                b(5) = (rhs(5,mm) - c15*b(1)-c25*b(2)
     &                    -c35*b(3)-c45*b(4))/c55

                w(5,mm) = b(5) / c55
                w(4,mm) = (b(4)-c45*w(5,mm))/c44
                w(3,mm) = (b(3)-c35*w(5,mm)-c34*w(4,mm))/c33
                w(2,mm) = (b(2)-c25*w(5,mm)-c24*w(4,mm)
     &               -c23*w(3,mm))/c22
                w(1,mm) = (b(1)-c15*w(5,mm)-c14*w(4,mm)
     &               -c13*w(3,mm)-c12*w(2,mm))/c11




                qyy(i,j,mm) =  2.d0*w(5,mm)/(hy**2)
                qxy(i,j,mm) =  w(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*w(3,mm)/(hx**2)
                qy(i,j,mm)  =  w(2,mm)/hy
                qx(i,j,mm)  =  w(1,mm)/hx

            end do


 35    continue

       return
       end






             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
             ! will have to add robustness checks
!             c11 = sqrt(a(1,1))
!             c12 = a(1,2) / c11
!             c22 = sqrt(a(2,2) - c12**2)
!
!
!             do mm = 1, nvar
!               b(1) = rhs(1,mm) / c11
!               b(2) = (rhs(2,mm) - c12*b(1))/c22
!               qy(i,j,mm) = b(2) / c22
!               qx(i,j,mm) = (b(1) - c12*qy(i,j,mm))/c11
!             end do
!
!
!
!             do mm = 1, nvar
!               qy(i,j,mm) = qy(i,j,mm) / hy
!               qx(i,j,mm) = qx(i,j,mm) / hx
!             end do

!               a(1,1) = a(1,1) + deltax*deltax / hx / hx
!               a(1,2) = a(1,2) + deltax*deltay / hx / hy
!               a(2,2) = a(2,2) + deltay*deltay / hy / hy
!               rhs(1,:) = rhs(1,:) + (deltax / hx) *
!     .                    (qp(i+ioff,j+joff,:) - qp(i,j,:))
!               rhs(2,:) = rhs(2,:) + (deltay / hy) *
!     .                    (qp(i+ioff,j+joff,:) - qp(i,j,:))
!
!
!               a(1,3) = a(1,3) + (deltax / hx) ** 3
!               a(1,4) = a(1,4) + (deltay / hy) * (deltax / hx) ** 2
!               a(1,5) = a(1,5) + (deltax / hx) * (deltay / hy) ** 2
!               a(2,3) = a(2,3) + (deltay / hy) * (deltax / hx) ** 2
!               a(2,4) = a(2,4) + (deltax / hx) * (deltay / hy) ** 2
!               a(2,5) = a(2,5) + (deltay / hy) ** 3
!               a(3,3) = a(3,3) + (deltax / hx) ** 4
!               a(3,4) = a(3,4) + (deltay / hy) * (deltax / hx) ** 3
!               a(3,5) = a(3,5) + ( (deltax / hx) ** 2 )
!     .                                  * (  (deltay / hy) ** 2  )
!               a(4,4) = a(4,4) + ( (deltax / hx) ** 2 )
!     .                                  * (  (deltay / hy) ** 2  )
!               a(4,5) = a(4,5) + (deltax / hx) * (deltay / hy) ** 3
!               a(5,5) = a(5,5) + (deltay / hy) ** 4
!
!               rhs(3,:) = rhs(3,:) + (qp(i+ioff,j+joff,:) - qp(i,j,:))
!     .                               *(deltax / hx) ** 2
!               rhs(4,:) = rhs(4,:) + (qp(i+ioff,j+joff,:) - qp(i,j,:))
!     .                               *(deltax / hx) * (deltay / hy)
!               rhs(5,:) = rhs(5,:) + (qp(i+ioff,j+joff,:) - qp(i,j,:))
!     .                               *(deltay / hy) ** 2

!
!             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
!             ! will have to add robustness checks
!             c11 = sqrt(a(1,1))
!             c12 = a(1,2) / c11
!             c22 = sqrt(a(2,2) - c12**2)
!
!             ! now back solve (C^t C = rhs of A^tdu ) to get x and y gradient for all variables
!             do mm = 1, nvar
!               b(1) = rhs(1,mm) / c11
!               b(2) = (rhs(2,mm) - c12*b(1))/c22
!               qy(i,j,mm) = b(2) / c22
!               qx(i,j,mm) = (b(1) - c12*qy(i,j,mm))/c11
!             end do
!


!        if(i .eq. 33 .and. j .eq. 10) then
!            print *,0,0, kirr
!!            print *, x0,y0, xoff, yoff
!!            print *, xoff, qp(i+ioff,j+joff,1)
!             write(*, 212,advance="no") '{'
! 212          format(A)
!            do ii = 1,3
!!             write(*, 90,advance="no") '{', poly(ii,1,koff),',',
!!     .                                      poly(ii,2,koff)
!! 90          format(A,e8.3,A,e8.3)
!             write(*, 31,advance="no") '{', poly(ii,1,kirr),',',
!     .                                      poly(ii,2,kirr),'}'
! 31          format(A,e25.16,A,e25.16,A)
!
!            end do
!             write(*, 32,advance="no") '};'
! 32          format(A)
!            print *, ""
!        endif


!            if(i .eq. 47 .and. j .eq. 34) then
!            print *, ata(1,:)
!            print *, ata(2,:)
!            print *, ata(3,:)
!            print *, ata(4,:)
!            print *, ata(5,:)
!            print *," "
!            print *,rhs(:,1)
!            print *," "
!            print *, qx(i,j,1), qy(i,j,1), qxx(i,j,1), qxy(i,j,1),
!     .               qyy(i,j,1)
!            print *," "
!            endif
!        if(i .eq. 33 .and. j .eq. 10) then
!            print *,ioff,joff, koff
!!            print *, x0,y0, xoff, yoff
!!            print *, xoff, qp(i+ioff,j+joff,1)
!             write(*, 91,advance="no") '{'
! 91          format(A)
!            do ii = 1,itri+2
!!             write(*, 90,advance="no") '{', poly(ii,1,koff),',',
!!     .                                      poly(ii,2,koff)
!! 90          format(A,e8.3,A,e8.3)
!             write(*, 90,advance="no") '{', poly(ii,1,koff),',',
!     .                                      poly(ii,2,koff),'}'
! 90          format(A,e25.16,A,e25.16,A)
!
!            end do
!             write(*, 92,advance="no") '};'
! 92          format(A)
!            print *, ""
!        endif
!          if(i .eq. 33 .and. j .eq. 10) then
!            print *,ioff,joff
!            print *,a
!            print *,qp(i+ioff,j+joff,1) - qp(i,j,1)
!          endif
