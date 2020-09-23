
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
       dimension qxxx(mitot,mjtot,4), qxxy(mitot,mjtot,4),
     .           qxyy(mitot,mjtot,4), qyyy(mitot,mjtot,4)
       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
       common /order2/ ssw, quad, nolimiter
       logical IS_GHOST
       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension b(9), w(9,nvar)
       dimension inuf(mitot,mjtot)
       include "cirr.i"
       include "./quadrature.i"



      qx = 0.d0
      qy = 0.d0

      qxx = 0.d0
      qxy = 0.d0
      qyy = 0.d0

      qxxx = 0.d0
      qxxy = 0.d0
      qxyy = 0.d0
      qyyy = 0.d0


      inuf = 0
      inuforder = 0

      if(ssw .eq. 1) inuforder = 1
      if(ssw .eq. 2) inuforder = 1
      if(ssw .eq. 3) inuforder = 2


       do 33 i = lwidth, mitot-lwidth+1
       do 33 j = lwidth, mjtot-lwidth+1
            k = irr(i,j)
            if (k .ne. lstgrd) goto 33
            inuf(i,j) = 1

            do 333 ioff = -inuforder,inuforder
            do 333 joff = -inuforder,inuforder
                if(irr(i+ioff, j+joff) .ne. lstgrd) then
                    inuf(i,j) = 0
                endif
 333        continue
 33    continue

      if(ssw .eq. 0) then

      ! do nothing

      return
      elseif(ssw .eq. 1) then
       do 34 i = lwidth, mitot-lwidth+1
       do 34 j = lwidth, mjtot-lwidth+1
            if(inuf(i,j) .eq. 0) cycle
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



 334         continue


      ! check for singularity of A
       if( ata(1,1) < 1d-8 ) cycle
       c11 = dsqrt(ata(1,1))
       c12 = ata(1,2)/c11


       if( ata(2,2)-c12**2 < 1d-8 ) cycle
       c22 = dsqrt(ata(2,2)-c12**2)

       ! quit the double for loop if the matrix is nonsingular.
       goto 435




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

 35   continue
      return




      elseif(ssw .eq. 2) then





       do 234 i = lwidth, mitot-lwidth+1
       do 234 j = lwidth, mjtot-lwidth+1
            if(inuf(i,j) .eq. 0) cycle

            qyy(i,j,:) =  (qp(i,j+1,:)-2*qp(i,j,:)+qp(i,j-1,:))/hy**2
            qxy(i,j,:) =  (qp(i+1,j+1,:)-qp(i+1,j-1,:)-qp(i-1,j+1,:)
     .                     +qp(i-1,j-1,:))/(4*hx*hy)
            qxx(i,j,:) =  (qp(i+1,j,:)-2*qp(i,j,:)+qp(i-1,j,:))/hx**2
            qy(i,j,:)  =  (qp(i,j+1,:) - qp(i,j-1,:))/(2*hy)
            qx(i,j,:)  =  (qp(i+1,j,:) - qp(i-1,j,:))/(2*hx)
 234    continue

       do 235 i = lwidth, mitot-lwidth+1
       do 235 j = lwidth, mjtot-lwidth+1
            kirr = irr(i,j)
            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
            if(kirr .eq. lstgrd) then
                x0 = xlow + (i-.5d0)*hx
                y0 = ylow + (j-.5d0)*hy
            else
                x0 = xcirr(kirr)
                y0 = ycirr(kirr)
            endif

            a = 0.d0
            ata = 0.d0
            rhs = 0.d0



            do 2394 moff = 1,3

            do 2334 ioff = -moff,moff
            do 2334 joff = -moff,moff
                idxmax = max(abs(ioff),abs(joff))
                if(idxmax .ne. moff) goto 2334


                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 2334





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


         ! load the offsets
         a(3) = -poly(8 ,1,kirr)
         a(4) = -poly(9 ,1,kirr)
         a(5) = -poly(10,1,kirr)



        if(koff .ne. lstgrd) then
          arr = ar(koff)
          ivert = 1
          do 220 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
  220      continue
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
          do 221 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 222 itq = 1,ntriquad
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
  222        continue ! for each quadrature point on each triangle
  221      continue ! for each triangle


            do ii = 1,5
            do jj = 1,5
                ata(ii,jj) = ata(ii,jj) + a(ii)*a(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .                      + a(ii) * (qp(i+ioff,j+joff,:) - qp(i,j,:))
            enddo

 2334         continue


            call cholesky(9, 5, ata, G, ipd)
            if(ipd .eq. 1) then
             ! quit the double for loop if the matrix is nonsingular.
                goto 2335
            endif

 2394         continue

      ! SOLVE THE NORMAL EQUATIONS.

 2335         continue ! quit the double for loop if the matrix is nonsingular.


             do mm = 1, nvar
                call solve_d(9,5,G,rhs(:,mm))

                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2)
                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2)
                qy(i,j,mm)  =  rhs(2,mm)/hy
                qx(i,j,mm)  =  rhs(1,mm)/hx

            end do


 235    continue
      return

      elseif(ssw .eq. 3) then


       do 434 i = lwidth, mitot-lwidth+1
       do 434 j = lwidth, mjtot-lwidth+1
            if(inuf(i,j) .eq. 0) cycle
            qxxx(i,j,:) = (qp(i+2, j,:) - 2*qp(i+1, j,:) + 2*qp(i-1,j,:)
     .                     -qp(i-2,j,:)) * 6.d0 / (2 * hx**3)
            qxxy(i,j,:) =2.d0*(qp(i+1,j+1,:)-2*qp(i,j+1,:)+qp(i-1,j+1,:)
     .    -qp(i+1,j-1,:)+2*qp(i,j-1,:)-qp(i-1,j-1,:))/(2*hy*hx**2)
            qxyy(i,j,:) = 2.d0*(qp(i,j+1,:)-2*qp(i,j,:)+qp(i,j-1,:)
     .            -qp(i,j+1,:)+2*qp(i,j,:)-qp(i,j-1,:))/(2 * hx * hy**2)
            qyyy(i,j,:) = (qp(i, j+2,:)-2 * qp(i, j+1,:) + 2*qp(i,j-1,:)
     .                     -qp(i,j-2,:))*6.d0 / (2 * hx**3)
            qyy(i,j,:) =  (qp(i,j+1,:)-2*qp(i,j,:)+qp(i,j-1,:))/hy**2
            qxy(i,j,:) =  (qp(i+1,j+1,:)-qp(i+1,j-1,:)-qp(i-1,j+1,:)
     .                     +qp(i-1,j-1,:))/(4*hx*hy)
            qxx(i,j,:) =  (qp(i+1,j,:)-2*qp(i,j,:)+qp(i-1,j,:))/hx**2
            qy(i,j,:)  =  (qp(i,j+1,:) - qp(i,j-1,:))/(2*hy)
            qx(i,j,:)  =  (qp(i+1,j,:) - qp(i-1,j,:))/(2*hx)
 434    continue








       do 735 i = lwidth, mitot-lwidth+1
       do 735 j = lwidth, mjtot-lwidth+1
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



            do 7394 moff = 1,4

            do 7334 ioff = -moff,moff
            do 7334 joff = -moff,moff
                idxmax = max(abs(ioff),abs(joff))
                if(idxmax .ne. moff) goto 7334


                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 7334





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

         ! load the offsets
         a(3) = -poly(8 ,1,kirr)
         a(4) = -poly(9 ,1,kirr)
         a(5) = -poly(10,1,kirr)
         a(6) = -dcubicshifts(1,kirr)
         a(7) = -dcubicshifts(2,kirr)
         a(8) = -dcubicshifts(3,kirr)
         a(9) = -dcubicshifts(4,kirr)

        if(koff .ne. lstgrd) then
          arr = ar(koff)
          ivert = 1
          do 720 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
  720      continue
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
          do 721 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 722 itq = 1,ntriquad
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

      a(6) = a(6) + (artri/arr) * wtri(itq) *( (xval-x0)**3 )/ (hx**3)
      a(7) = a(7) + (artri/arr) * wtri(itq) *( (yval-y0)*(xval-x0)**2 )
     .                                                   / (hy * hx**2)
      a(8) = a(8) + (artri/arr) * wtri(itq) *((xval-x0)*(yval-y0)**2)
     .                                        / (hx*hy**2)
      a(9) = a(9) + (artri/arr) * wtri(itq) *( (yval-y0)**3 )/ (hy**3)

  722        continue ! for each quadrature point on each triangle
  721      continue ! for each triangle

            do ii = 1,9
            do jj = 1,9
                ata(ii,jj) = ata(ii,jj) + a(ii)*a(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .                      + a(ii) * (qp(i+ioff,j+joff,:) - qp(i,j,:))
            enddo

 7334         continue

            call cholesky(9, 9, ata, G, ipd)
            if(ipd .eq. 1) then
             ! quit the double for loop if the matrix is nonsingular.
                goto 7335
            endif

 7394         continue





      ! SOLVE THE NORMAL EQUATIONS.

 7335         continue ! quit the double for loop if the matrix is nonsingular.
        do mm = 1, nvar
            call solve_d(9,9, G, rhs(:,mm))
                qyyy(i,j,mm) =   rhs(9,mm)/(hy**3)
                qxyy(i,j,mm) =   rhs(8,mm)/(hx*hy**2)
                qxxy(i,j,mm) =   rhs(7,mm)/(hy*hx**2)
                qxxx(i,j,mm) =   rhs(6,mm)/(hx**3)
                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2)
                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2)
                qy(i,j,mm)  =  rhs(2,mm)/hy
                qx(i,j,mm)  =  rhs(1,mm)/hx
        enddo


 735    continue





















      return
      endif

      end






