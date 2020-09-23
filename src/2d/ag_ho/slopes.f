
c
c------------------------------------------------------------
c
       subroutine slopes(qp,qx,qy,qxx,qxy,qyy,
     .                   qxxx, qxxy, qxyy, qyyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
       dimension qxxx(mitot,mjtot,4), qxxy(mitot,mjtot,4),
     .           qxyy(mitot,mjtot,4), qyyy(mitot,mjtot,4)
       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold, ilts,ibc
       common /order2/ ssw, quad, nolimiter
       dimension sigma(4)

      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)


       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"

      qx = 0.d0
      qy = 0.d0

      qxx = 0.d0
      qxy = 0.d0
      qyy = 0.d0

      qxxx = 0.d0
      qxxy = 0.d0
      qxyy = 0.d0
      qyyy = 0.d0





      if(ssw .eq. 0) then
      ! do nothing
      return
      elseif(ssw .eq. 1) then
        call slopes1(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
      return

      elseif(ssw .eq. -1) then
        call slopes1q(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)

      return
      elseif(ssw .eq. -10) then
        call slopes1qp(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)

      return
      elseif(ssw .eq. 2 ) then

      select case(ihob)
        case(0)
        call slopes2(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
        case(1)
        call slopes2_ho(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)

        call slopes2_ho_qr(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
      end select

      return

      elseif(ssw .eq. -2) then
      select case(ihob)
        case(0)
        call slopes3(qp,qx,qy,qxx,qxy,qyy,
     .                   qxxx, qxxy, qxyy, qyyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
        case(1)
        call slopes2q_ho(qp,qx,qy,qxx,qxy,qyy,
     .                   qxxx, qxxy, qxyy, qyyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
      end select

      qxxx = 0.d0
      qxxy = 0.d0
      qxyy = 0.d0
      qyyy = 0.d0

      return


      elseif(ssw .eq. 3 .or. ssw .eq. -3) then
      select case(ihob)
        case(0)
        call slopes3(qp,qx,qy,qxx,qxy,qyy,
     .                   qxxx, qxxy, qxyy, qyyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
        case(1)
        call slopes3_ho(qp,qx,qy,qxx,qxy,qyy,
     .                   qxxx, qxxy, qxyy, qyyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
      end select

      return

      endif

      end









      subroutine slopes2q_ho(qp,qx,qy,qxx,qxy,qyy,
     .                   qxxx, qxxy, qxyy, qyyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
       dimension qxxx(mitot,mjtot,4), qxxy(mitot,mjtot,4),
     .           qxyy(mitot,mjtot,4), qyyy(mitot,mjtot,4)
      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)
       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"





       do 434 i = 1, mitot
       do 434 j = 1, mjtot
            if(inuf(i,j) .eq. 0) cycle




      qyyy(i,j,:) = (2*qp(i,j-1,:) - qp(i,j-2,:) - 2*qp(i,j+1,:)
     .+ qp(i,j+2,:) + 2*qp(i-1,j-1,:) - qp(i-1,j-2,:) - 2*qp(i-1,j+1,:)
     .+ qp(i-1,j+2,:) + 2*qp(i-2,j-1,:) - qp(i-2,j-2,:)
     . - 2*qp(i-2,j+1,:) + qp(i-2,j+2,:) +  2*qp(i+1,j-1,:)
     . - qp(i+1,j-2,:) - 2*qp(i+1,j+1,:) + qp(i+1,j+2,:)
     .+  2*qp(i+2,j-1,:) - qp(i+2,j-2,:) - 2*qp(i+2,j+1,:)
     . + qp(i+2,j+2,:))/(10.*hy**3)

      qxyy(i,j,:) = (2*qp(i-1,j,:) + qp(i-1,j-1,:) - 2*qp(i-1,j-2,:)
     . + qp(i-1,j+1,:) - 2*qp(i-1,j+2,:) + 4*qp(i-2,j,:)
     .+ 2*qp(i-2,j-1,:) - 4*qp(i-2,j-2,:) + 2*qp(i-2,j+1,:)
     .- 4*qp(i-2,j+2,:) - 2*qp(i+1,j,:) - qp(i+1,j-1,:) +2*qp(i+1,j-2,:)
     .- qp(i+1,j+1,:) + 2*qp(i+1,j+2,:) - 4*qp(i+2,j,:)
     .-2*qp(i+2,j-1,:) + 4*qp(i+2,j-2,:) - 2*qp(i+2,j+1,:)
     .+ 4*qp(i+2,j+2,:))/(70.*hx*hy**2)


      qxxy(i,j,:) = (2*qp(i,j-1,:) + 4*qp(i,j-2,:) - 2*qp(i,j+1,:)
     .- 4*qp(i,j+2,:) +qp(i-1,j-1,:) + 2*qp(i-1,j-2,:) - qp(i-1,j+1,:)
     . - 2*qp(i-1,j+2,:) - 2*qp(i-2,j-1,:) - 4*qp(i-2,j-2,:)
     . + 2*qp(i-2,j+1,:) + 4*qp(i-2,j+2,:)+ qp(i+1,j-1,:)
     . + 2*qp(i+1,j-2,:) - qp(i+1,j+1,:)
     .- 2*qp(i+1,j+2,:) -2*qp(i+2,j-1,:) - 4*qp(i+2,j-2,:)
     . + 2*qp(i+2,j+1,:) + 4*qp(i+2,j+2,:))/(70.*hx**2*hy)

!            qyyy(i,j,:) = (qp(i,j+2,:)-2.d0*qp(i,j+1,:)
!     .                    +2.d0*qp(i,j-1,:)-qp(i,j-2,:))/(2*hy**3)

      qxxx(i,j,:) = (2*qp(i-1,j,:) + 2*qp(i-1,j-1,:) + 2*qp(i-1,j-2,:)
     .+ 2*qp(i-1,j+1,:) + 2*qp(i-1,j+2,:) - qp(i-2,j,:) - qp(i-2,j-1,:)
     . - qp(i-2,j-2,:) - qp(i-2,j+1,:) - qp(i-2,j+2,:) - 2*qp(i+1,j,:)
     . - 2*qp(i+1,j-1,:) - 2*qp(i+1,j-2,:) - 2*qp(i+1,j+1,:)
     . - 2*qp(i+1,j+2,:) + qp(i+2,j,:) + qp(i+2,j-1,:) + qp(i+2,j-2,:)
     . + qp(i+2,j+1,:) +  qp(i+2,j+2,:))/(10.*hx**3)


      qyy(i,j,:) = (-350*qp(i,j,:) + 17*qp(i,j-1,:) + 68*qp(i,j-2,:)
     . + 17*qp(i,j+1,:) + 68*qp(i,j+2,:) - 10*qp(i-1,j,:)
     . + 7*qp(i-1,j-1,:) + 58*qp(i-1,j-2,:)  + 7*qp(i-1,j+1,:)
     .+ 58*qp(i-1,j+2,:) - 40*qp(i-2,j,:) -  23*qp(i-2,j-1,:)
     . + 28*qp(i-2,j-2,:) - 23*qp(i-2,j+1,:) + 28*qp(i-2,j+2,:)
     . - 10*qp(i+1,j,:) + 7*qp(i+1,j-1,:) + 58*qp(i+1,j-2,:)
     . + 7*qp(i+1,j+1,:) + 58*qp(i+1,j+2,:) -  40*qp(i+2,j,:)
     . - 23*qp(i+2,j-1,:) + 28*qp(i+2,j-2,:) - 23*qp(i+2,j+1,:)
     . + 28*qp(i+2,j+2,:))/(945.*hy**2)

      qxy(i,j,:) = (qp(i-1,j-1,:) + 2*qp(i-1,j-2,:) - qp(i-1,j+1,:)
     .- 2*qp(i-1,j+2,:) + 2*qp(i-2,j-1,:) + 4*qp(i-2,j-2,:)
     .- 2*qp(i-2,j+1,:) - 4*qp(i-2,j+2,:)- qp(i+1,j-1,:)
     . - 2*qp(i+1,j-2,:) + qp(i+1,j+1,:) + 2*qp(i+1,j+2,:)
     .- 2*qp(i+2,j-1,:) - 4*qp(i+2,j-2,:) + 2*qp(i+2,j+1,:)
     .+ 4*qp(i+2,j+2,:))/(100.*hx*hy)

      qxx(i,j,:) = (-350*qp(i,j,:) - 10*qp(i,j-1,:) - 40*qp(i,j-2,:)
     . - 10*qp(i,j+1,:) - 40*qp(i,j+2,:) + 17*qp(i-1,j,:)
     . + 7*qp(i-1,j-1,:) - 23*qp(i-1,j-2,:) + 7*qp(i-1,j+1,:)
     . - 23*qp(i-1,j+2,:) + 68*qp(i-2,j,:) + 58*qp(i-2,j-1,:)
     . + 28*qp(i-2,j-2,:) + 58*qp(i-2,j+1,:) + 28*qp(i-2,j+2,:)
     .+ 17*qp(i+1,j,:) + 7*qp(i+1,j-1,:) - 23*qp(i+1,j-2,:)
     . + 7*qp(i+1,j+1,:) - 23*qp(i+1,j+2,:) + 68*qp(i+2,j,:)
     . + 58*qp(i+2,j-1,:) + 28*qp(i+2,j-2,:) + 58*qp(i+2,j+1,:)
     . + 28*qp(i+2,j+2,:))/(945.*hx**2)

         qy(i,j,:)  = (-288*qp(i,j-1,:) - 65*qp(i,j-2,:)
     .+ 288*qp(i,j+1,:) + 65*qp(i,j+2,:) - 263*qp(i-1,j-1,:)
     . - 15*qp(i-1,j-2,:) + 263*qp(i-1,j+1,:) + 15*qp(i-1,j+2,:)
     .- 188*qp(i-2,j-1,:) + 135*qp(i-2,j-2,:) + 188*qp(i-2,j+1,:)
     .- 135*qp(i-2,j+2,:) - 263*qp(i+1,j-1,:) - 15*qp(i+1,j-2,:)
     .+ 263*qp(i+1,j+1,:) + 15*qp(i+1,j+2,:) - 188*qp(i+2,j-1,:)
     . + 135*qp(i+2,j-2,:) + 188*qp(i+2,j+1,:)
     . - 135*qp(i+2,j+2,:))/(1680.*hy)
         qx(i,j,:)  = (-288*qp(i-1,j,:) - 263*qp(i-1,j-1,:)
     .- 188*qp(i-1,j-2,:) - 263*qp(i-1,j+1,:) - 188*qp(i-1,j+2,:)
     .- 65*qp(i-2,j,:) - 15*qp(i-2,j-1,:) + 135*qp(i-2,j-2,:)
     .- 15*qp(i-2,j+1,:) + 135*qp(i-2,j+2,:) + 288*qp(i+1,j,:)
     .+ 263*qp(i+1,j-1,:) + 188*qp(i+1,j-2,:) + 263*qp(i+1,j+1,:)
     .+ 188*qp(i+1,j+2,:) + 65*qp(i+2,j,:) + 15*qp(i+2,j-1,:)
     .- 135*qp(i+2,j-2,:) + 15*qp(i+2,j+1,:)
     .- 135*qp(i+2,j+2,:))/(1680.*hx)

!            if(i .eq. 10 .and. j .eq. 10) then
!            print *, "here"
!            endif
 434    continue








       do 735 i = 1, mitot
       do 735 j = 1, mjtot
            kirr = irr(i,j)
            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
            if(kirr .eq. lstgrd) then
                x0 = xlow + (i-.5d0)*hx
                y0 = ylow + (j-.5d0)*hy
            else
                x0 = xcirr_ho(kirr)
                y0 = ycirr_ho(kirr)
            endif


      a = 0.d0
      ata = 0.d0
      rhs = 0.d0

            do 7334 ioff = -iir(i,j),iir(i,j)
            do 7334 joff = -jjr(i,j),jjr(i,j)
                if (.not. IS_REAL(i+ioff,j+joff)) go to 7334
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 7334


                if(koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*hx
                    yoff = ylow + (j+joff-.5d0)*hy
                else
                    xoff = xcirr_ho(koff)
                    yoff = ycirr_ho(koff)
                endif

                deltax = xoff - x0
                deltay = yoff - y0

                a(1) = deltax/hx
                a(2) = deltay/hy

                 ! load the offsets
                 a(3) = -poly_ho(8 ,1,kirr)
                 a(4) = -poly_ho(9 ,1,kirr)
                 a(5) = -poly_ho(10,1,kirr)
                 a(6) = -dcubicshifts_ho(1,kirr)
                 a(7) = -dcubicshifts_ho(2,kirr)
                 a(8) = -dcubicshifts_ho(3,kirr)
                 a(9) = -dcubicshifts_ho(4,kirr)

                 if(koff .ne. lstgrd) then
                   arr = ar_ho(koff)
                   ivert = 1
                   do 720 while (poly(ivert+1,1,koff) .ne. -11.)
                     ivert = ivert + 1
  720               continue
                    itri = ivert - 3
                 else
                   arr = hx*hy
                   itri = 2+1 ! added an extra one artificially
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
          do 721 it = 1, itri-1 ! for each  triangle
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



      if(koff .ne. lstgrd) then
      x1 = poly(ivert-2,1,koff)
      y1 = poly(ivert-2,2,koff)

      x2 = bdry(1,1,koff)
      y2 = bdry(1,2,koff)

      x3 = bdry(3,1,koff)
      y3 = bdry(3,2,koff)

      x4 = bdry(2,1,koff)
      y4 = bdry(2,2,koff)

      do nq = 1,ntriquad_ho
      artri = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)

      xval = rs2xy_2(x1,x2,x3,x4,nq)
      yval = rs2xy_2(y1,y2,y3,y4,nq)

                a(3) = a(3) + (artri/arr) * wtri_ho(nq) *( (xval-x0)**2)
     .                                                / (hx**2)

                a(4) = a(4) + (artri/arr) * wtri_ho(nq) *( (xval-x0)/hx)
     .                                                *( (yval-y0)/hy )

                a(5) = a(5) + (artri/arr) * wtri_ho(nq) *( (yval-y0)**2)
     .                                                / (hy**2)

      a(6) = a(6) + (artri/arr) * wtri_ho(nq) *( (xval-x0)**3 )/ (hx**3)
      a(7) = a(7) + (artri/arr) * wtri_ho(nq) *( (yval-y0)*(xval-x0)**2)
     .                                                   / (hy * hx**2)
      a(8) = a(8) + (artri/arr) * wtri_ho(nq) *((xval-x0)*(yval-y0)**2)
     .                                        / (hx*hy**2)
      a(9) = a(9) + (artri/arr) * wtri_ho(nq) *( (yval-y0)**3 )/ (hy**3)


      enddo
      endif







            do ii = 1,9
            do jj = 1,9
                ata(ii,jj) = ata(ii,jj) + a(ii)*a(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .                      + a(ii) * (qp(i+ioff,j+joff,:) - qp(i,j,:))
            enddo

 7334         continue


      ! SOLVE THE NORMAL EQUATIONS.
        call cholesky(9, 9, ata, G)
        do mm = 1, nvar
            call solve(9,9, G, rhs(:,mm))
                qyyy(i,j,mm) = 6.d0*rhs(9,mm)/(hy**3)
                qxyy(i,j,mm) = 2.d0*rhs(8,mm)/(hx*hy**2)
                qxxy(i,j,mm) = 2.d0*rhs(7,mm)/(hy*hx**2)
                qxxx(i,j,mm) = 6.d0*rhs(6,mm)/(hx**3)
                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2)
                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2)
                qy(i,j,mm)  =  rhs(2,mm)/hy
                qx(i,j,mm)  =  rhs(1,mm)/hx

!                if(mm .eq. 1) then
!                print *, "here"
!                endif
        enddo


 735    continue











      end subroutine






























      subroutine slopes3_ho(qp,qx,qy,qxx,qxy,qyy,
     .                   qxxx, qxxy, qxyy, qyyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
       dimension qxxx(mitot,mjtot,4), qxxy(mitot,mjtot,4),
     .           qxyy(mitot,mjtot,4), qyyy(mitot,mjtot,4)
      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)
       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"





       do 434 i = 1, mitot
       do 434 j = 1, mjtot
            if(inuf(i,j) .eq. 0) cycle




      qyyy(i,j,:) = (2*qp(i,j-1,:) - qp(i,j-2,:) - 2*qp(i,j+1,:)
     .+ qp(i,j+2,:) + 2*qp(i-1,j-1,:) - qp(i-1,j-2,:) - 2*qp(i-1,j+1,:)
     .+ qp(i-1,j+2,:) + 2*qp(i-2,j-1,:) - qp(i-2,j-2,:)
     . - 2*qp(i-2,j+1,:) + qp(i-2,j+2,:) +  2*qp(i+1,j-1,:)
     . - qp(i+1,j-2,:) - 2*qp(i+1,j+1,:) + qp(i+1,j+2,:)
     .+  2*qp(i+2,j-1,:) - qp(i+2,j-2,:) - 2*qp(i+2,j+1,:)
     . + qp(i+2,j+2,:))/(10.*hy**3)

      qxyy(i,j,:) = (2*qp(i-1,j,:) + qp(i-1,j-1,:) - 2*qp(i-1,j-2,:)
     . + qp(i-1,j+1,:) - 2*qp(i-1,j+2,:) + 4*qp(i-2,j,:)
     .+ 2*qp(i-2,j-1,:) - 4*qp(i-2,j-2,:) + 2*qp(i-2,j+1,:)
     .- 4*qp(i-2,j+2,:) - 2*qp(i+1,j,:) - qp(i+1,j-1,:) +2*qp(i+1,j-2,:)
     .- qp(i+1,j+1,:) + 2*qp(i+1,j+2,:) - 4*qp(i+2,j,:)
     .-2*qp(i+2,j-1,:) + 4*qp(i+2,j-2,:) - 2*qp(i+2,j+1,:)
     .+ 4*qp(i+2,j+2,:))/(70.*hx*hy**2)


      qxxy(i,j,:) = (2*qp(i,j-1,:) + 4*qp(i,j-2,:) - 2*qp(i,j+1,:)
     .- 4*qp(i,j+2,:) +qp(i-1,j-1,:) + 2*qp(i-1,j-2,:) - qp(i-1,j+1,:)
     . - 2*qp(i-1,j+2,:) - 2*qp(i-2,j-1,:) - 4*qp(i-2,j-2,:)
     . + 2*qp(i-2,j+1,:) + 4*qp(i-2,j+2,:)+ qp(i+1,j-1,:)
     . + 2*qp(i+1,j-2,:) - qp(i+1,j+1,:)
     .- 2*qp(i+1,j+2,:) -2*qp(i+2,j-1,:) - 4*qp(i+2,j-2,:)
     . + 2*qp(i+2,j+1,:) + 4*qp(i+2,j+2,:))/(70.*hx**2*hy)

!            qyyy(i,j,:) = (qp(i,j+2,:)-2.d0*qp(i,j+1,:)
!     .                    +2.d0*qp(i,j-1,:)-qp(i,j-2,:))/(2*hy**3)

      qxxx(i,j,:) = (2*qp(i-1,j,:) + 2*qp(i-1,j-1,:) + 2*qp(i-1,j-2,:)
     .+ 2*qp(i-1,j+1,:) + 2*qp(i-1,j+2,:) - qp(i-2,j,:) - qp(i-2,j-1,:)
     . - qp(i-2,j-2,:) - qp(i-2,j+1,:) - qp(i-2,j+2,:) - 2*qp(i+1,j,:)
     . - 2*qp(i+1,j-1,:) - 2*qp(i+1,j-2,:) - 2*qp(i+1,j+1,:)
     . - 2*qp(i+1,j+2,:) + qp(i+2,j,:) + qp(i+2,j-1,:) + qp(i+2,j-2,:)
     . + qp(i+2,j+1,:) +  qp(i+2,j+2,:))/(10.*hx**3)


      qyy(i,j,:) = (-350*qp(i,j,:) + 17*qp(i,j-1,:) + 68*qp(i,j-2,:)
     . + 17*qp(i,j+1,:) + 68*qp(i,j+2,:) - 10*qp(i-1,j,:)
     . + 7*qp(i-1,j-1,:) + 58*qp(i-1,j-2,:)  + 7*qp(i-1,j+1,:)
     .+ 58*qp(i-1,j+2,:) - 40*qp(i-2,j,:) -  23*qp(i-2,j-1,:)
     . + 28*qp(i-2,j-2,:) - 23*qp(i-2,j+1,:) + 28*qp(i-2,j+2,:)
     . - 10*qp(i+1,j,:) + 7*qp(i+1,j-1,:) + 58*qp(i+1,j-2,:)
     . + 7*qp(i+1,j+1,:) + 58*qp(i+1,j+2,:) -  40*qp(i+2,j,:)
     . - 23*qp(i+2,j-1,:) + 28*qp(i+2,j-2,:) - 23*qp(i+2,j+1,:)
     . + 28*qp(i+2,j+2,:))/(945.*hy**2)

      qxy(i,j,:) = (qp(i-1,j-1,:) + 2*qp(i-1,j-2,:) - qp(i-1,j+1,:)
     .- 2*qp(i-1,j+2,:) + 2*qp(i-2,j-1,:) + 4*qp(i-2,j-2,:)
     .- 2*qp(i-2,j+1,:) - 4*qp(i-2,j+2,:)- qp(i+1,j-1,:)
     . - 2*qp(i+1,j-2,:) + qp(i+1,j+1,:) + 2*qp(i+1,j+2,:)
     .- 2*qp(i+2,j-1,:) - 4*qp(i+2,j-2,:) + 2*qp(i+2,j+1,:)
     .+ 4*qp(i+2,j+2,:))/(100.*hx*hy)

      qxx(i,j,:) = (-350*qp(i,j,:) - 10*qp(i,j-1,:) - 40*qp(i,j-2,:)
     . - 10*qp(i,j+1,:) - 40*qp(i,j+2,:) + 17*qp(i-1,j,:)
     . + 7*qp(i-1,j-1,:) - 23*qp(i-1,j-2,:) + 7*qp(i-1,j+1,:)
     . - 23*qp(i-1,j+2,:) + 68*qp(i-2,j,:) + 58*qp(i-2,j-1,:)
     . + 28*qp(i-2,j-2,:) + 58*qp(i-2,j+1,:) + 28*qp(i-2,j+2,:)
     .+ 17*qp(i+1,j,:) + 7*qp(i+1,j-1,:) - 23*qp(i+1,j-2,:)
     . + 7*qp(i+1,j+1,:) - 23*qp(i+1,j+2,:) + 68*qp(i+2,j,:)
     . + 58*qp(i+2,j-1,:) + 28*qp(i+2,j-2,:) + 58*qp(i+2,j+1,:)
     . + 28*qp(i+2,j+2,:))/(945.*hx**2)

         qy(i,j,:)  = (-288*qp(i,j-1,:) - 65*qp(i,j-2,:)
     .+ 288*qp(i,j+1,:) + 65*qp(i,j+2,:) - 263*qp(i-1,j-1,:)
     . - 15*qp(i-1,j-2,:) + 263*qp(i-1,j+1,:) + 15*qp(i-1,j+2,:)
     .- 188*qp(i-2,j-1,:) + 135*qp(i-2,j-2,:) + 188*qp(i-2,j+1,:)
     .- 135*qp(i-2,j+2,:) - 263*qp(i+1,j-1,:) - 15*qp(i+1,j-2,:)
     .+ 263*qp(i+1,j+1,:) + 15*qp(i+1,j+2,:) - 188*qp(i+2,j-1,:)
     . + 135*qp(i+2,j-2,:) + 188*qp(i+2,j+1,:)
     . - 135*qp(i+2,j+2,:))/(1680.*hy)
         qx(i,j,:)  = (-288*qp(i-1,j,:) - 263*qp(i-1,j-1,:)
     .- 188*qp(i-1,j-2,:) - 263*qp(i-1,j+1,:) - 188*qp(i-1,j+2,:)
     .- 65*qp(i-2,j,:) - 15*qp(i-2,j-1,:) + 135*qp(i-2,j-2,:)
     .- 15*qp(i-2,j+1,:) + 135*qp(i-2,j+2,:) + 288*qp(i+1,j,:)
     .+ 263*qp(i+1,j-1,:) + 188*qp(i+1,j-2,:) + 263*qp(i+1,j+1,:)
     .+ 188*qp(i+1,j+2,:) + 65*qp(i+2,j,:) + 15*qp(i+2,j-1,:)
     .- 135*qp(i+2,j-2,:) + 15*qp(i+2,j+1,:)
     .- 135*qp(i+2,j+2,:))/(1680.*hx)

!            if(i .eq. 10 .and. j .eq. 10) then
!            print *, "here"
!            endif
 434    continue








       do 735 i = 1, mitot
       do 735 j = 1, mjtot
            kirr = irr(i,j)
            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
            if(kirr .eq. lstgrd) then
                x0 = xlow + (i-.5d0)*hx
                y0 = ylow + (j-.5d0)*hy
            else
                x0 = xcirr_ho(kirr)
                y0 = ycirr_ho(kirr)
            endif


      a = 0.d0
      ata = 0.d0
      rhs = 0.d0

            do 7334 ioff = -iir(i,j),iir(i,j)
            do 7334 joff = -jjr(i,j),jjr(i,j)
                if (.not. IS_REAL(i+ioff,j+joff)) go to 7334
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 7334


                if(koff .eq. lstgrd) then
                    xoff = xlow + (dfloat(i+ioff)-.5d0)*hx
                    yoff = ylow + (dfloat(j+joff)-.5d0)*hy
                else
                    xoff = xcirr_ho(koff)
                    yoff = ycirr_ho(koff)
                endif

                deltax = xoff - x0
                deltay = yoff - y0

                a(1) = deltax/hx
                a(2) = deltay/hy

                 ! load the offsets
                 a(3) = -poly_ho(8 ,1,kirr)
                 a(4) = -poly_ho(9 ,1,kirr)
                 a(5) = -poly_ho(10,1,kirr)
                 a(6) = -dcubicshifts_ho(1,kirr)
                 a(7) = -dcubicshifts_ho(2,kirr)
                 a(8) = -dcubicshifts_ho(3,kirr)
                 a(9) = -dcubicshifts_ho(4,kirr)

                 if(koff .ne. lstgrd) then
                   arr = ar_ho(koff)
                   ivert = 1
                   do 720 while (poly(ivert+1,1,koff) .ne. -11.)
                     ivert = ivert + 1
  720               continue
                    itri = ivert - 3
                 else
                   arr = hx*hy
                   itri = 2+1 ! added an extra one artificially
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
          do 721 it = 1, itri-1 ! for each  triangle
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



      if(koff .ne. lstgrd) then
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
      artri = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)

      xval = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yval = rs2xy_3(y1,y2,y3,y4,y5,nq)

                a(3) = a(3) + (artri/arr) * wtri_ho(nq) *( (xval-x0)**2)
     .                                                / (hx**2)

                a(4) = a(4) + (artri/arr) * wtri_ho(nq) *( (xval-x0)/hx)
     .                                                *( (yval-y0)/hy )

                a(5) = a(5) + (artri/arr) * wtri_ho(nq) *( (yval-y0)**2)
     .                                                / (hy**2)

      a(6) = a(6) + (artri/arr) * wtri_ho(nq) *( (xval-x0)**3 )/ (hx**3)
      a(7) = a(7) + (artri/arr) * wtri_ho(nq) *( (yval-y0)*(xval-x0)**2)
     .                                                   / (hy * hx**2)
      a(8) = a(8) + (artri/arr) * wtri_ho(nq) *((xval-x0)*(yval-y0)**2)
     .                                        / (hx*hy**2)
      a(9) = a(9) + (artri/arr) * wtri_ho(nq) *( (yval-y0)**3 )/ (hy**3)




!
!      a(3) = a(3) + (artri/arr) * wtri_ho(nq) *( (xval-x0)**2)
!     .                                                / (hx**2)
!
!      a(4) = a(4) + (artri/arr) * wtri_ho(nq) *( (xval-x0)/hx )
!     .                                                *( (yval-y0)/hy )
!
!      a(5) = a(5) + (artri/arr) * wtri_ho(nq) *( (yval-y0)**2 )
!     .                                                / (hy**2)
      enddo
      endif







            do ii = 1,9
            do jj = 1,9
                ata(ii,jj) = ata(ii,jj) + a(ii)*a(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .                      + a(ii) * (qp(i+ioff,j+joff,:) - qp(i,j,:))
            enddo

 7334         continue


      ! SOLVE THE NORMAL EQUATIONS.
        call cholesky(9, 9, ata, G)
!        if(i .eq. 9 .and. j .eq. 24) then
!            print *, ata
!            print *, G
!            print *, rhs(:,1)
!        endif
        do mm = 1, nvar
            call solve(9,9, G, rhs(:,mm))
                qyyy(i,j,mm) = 6.d0*rhs(9,mm)/(hy**3)
                qxyy(i,j,mm) = 2.d0*rhs(8,mm)/(hx*hy**2)
                qxxy(i,j,mm) = 2.d0*rhs(7,mm)/(hy*hx**2)
                qxxx(i,j,mm) = 6.d0*rhs(6,mm)/(hx**3)
                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2)
                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2)
                qy(i,j,mm)  =  rhs(2,mm)/hy
                qx(i,j,mm)  =  rhs(1,mm)/hx

!                if(mm .eq. 1) then
!                print *, "here"
!                endif
        enddo


 735    continue











      end subroutine








      subroutine slopes3(qp,qx,qy,qxx,qxy,qyy,
     .                   qxxx, qxxy, qxyy, qyyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
       dimension qxxx(mitot,mjtot,4), qxxy(mitot,mjtot,4),
     .           qxyy(mitot,mjtot,4), qyyy(mitot,mjtot,4)
      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)
       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"

       do 434 i = 1, mitot
       do 434 j = 1, mjtot
            if(inuf(i,j) .eq. 0) cycle




      qyyy(i,j,:) = (2*qp(i,j-1,:) - qp(i,j-2,:) - 2*qp(i,j+1,:)
     .+ qp(i,j+2,:) + 2*qp(i-1,j-1,:) - qp(i-1,j-2,:) - 2*qp(i-1,j+1,:)
     .+ qp(i-1,j+2,:) + 2*qp(i-2,j-1,:) - qp(i-2,j-2,:)
     . - 2*qp(i-2,j+1,:) + qp(i-2,j+2,:) +  2*qp(i+1,j-1,:)
     . - qp(i+1,j-2,:) - 2*qp(i+1,j+1,:) + qp(i+1,j+2,:)
     .+  2*qp(i+2,j-1,:) - qp(i+2,j-2,:) - 2*qp(i+2,j+1,:)
     . + qp(i+2,j+2,:))/(10.*hy**3)

      qxyy(i,j,:) = (2*qp(i-1,j,:) + qp(i-1,j-1,:) - 2*qp(i-1,j-2,:)
     . + qp(i-1,j+1,:) - 2*qp(i-1,j+2,:) + 4*qp(i-2,j,:)
     .+ 2*qp(i-2,j-1,:) - 4*qp(i-2,j-2,:) + 2*qp(i-2,j+1,:)
     .- 4*qp(i-2,j+2,:) - 2*qp(i+1,j,:) - qp(i+1,j-1,:) +2*qp(i+1,j-2,:)
     .- qp(i+1,j+1,:) + 2*qp(i+1,j+2,:) - 4*qp(i+2,j,:)
     .-2*qp(i+2,j-1,:) + 4*qp(i+2,j-2,:) - 2*qp(i+2,j+1,:)
     .+ 4*qp(i+2,j+2,:))/(70.*hx*hy**2)


      qxxy(i,j,:) = (2*qp(i,j-1,:) + 4*qp(i,j-2,:) - 2*qp(i,j+1,:)
     .- 4*qp(i,j+2,:) +qp(i-1,j-1,:) + 2*qp(i-1,j-2,:) - qp(i-1,j+1,:)
     . - 2*qp(i-1,j+2,:) - 2*qp(i-2,j-1,:) - 4*qp(i-2,j-2,:)
     . + 2*qp(i-2,j+1,:) + 4*qp(i-2,j+2,:)+ qp(i+1,j-1,:)
     . + 2*qp(i+1,j-2,:) - qp(i+1,j+1,:)
     .- 2*qp(i+1,j+2,:) -2*qp(i+2,j-1,:) - 4*qp(i+2,j-2,:)
     . + 2*qp(i+2,j+1,:) + 4*qp(i+2,j+2,:))/(70.*hx**2*hy)

!            qyyy(i,j,:) = (qp(i,j+2,:)-2.d0*qp(i,j+1,:)
!     .                    +2.d0*qp(i,j-1,:)-qp(i,j-2,:))/(2*hy**3)

      qxxx(i,j,:) = (2*qp(i-1,j,:) + 2*qp(i-1,j-1,:) + 2*qp(i-1,j-2,:)
     .+ 2*qp(i-1,j+1,:) + 2*qp(i-1,j+2,:) - qp(i-2,j,:) - qp(i-2,j-1,:)
     . - qp(i-2,j-2,:) - qp(i-2,j+1,:) - qp(i-2,j+2,:) - 2*qp(i+1,j,:)
     . - 2*qp(i+1,j-1,:) - 2*qp(i+1,j-2,:) - 2*qp(i+1,j+1,:)
     . - 2*qp(i+1,j+2,:) + qp(i+2,j,:) + qp(i+2,j-1,:) + qp(i+2,j-2,:)
     . + qp(i+2,j+1,:) +  qp(i+2,j+2,:))/(10.*hx**3)


      qyy(i,j,:) = (-350*qp(i,j,:) + 17*qp(i,j-1,:) + 68*qp(i,j-2,:)
     . + 17*qp(i,j+1,:) + 68*qp(i,j+2,:) - 10*qp(i-1,j,:)
     . + 7*qp(i-1,j-1,:) + 58*qp(i-1,j-2,:)  + 7*qp(i-1,j+1,:)
     .+ 58*qp(i-1,j+2,:) - 40*qp(i-2,j,:) -  23*qp(i-2,j-1,:)
     . + 28*qp(i-2,j-2,:) - 23*qp(i-2,j+1,:) + 28*qp(i-2,j+2,:)
     . - 10*qp(i+1,j,:) + 7*qp(i+1,j-1,:) + 58*qp(i+1,j-2,:)
     . + 7*qp(i+1,j+1,:) + 58*qp(i+1,j+2,:) -  40*qp(i+2,j,:)
     . - 23*qp(i+2,j-1,:) + 28*qp(i+2,j-2,:) - 23*qp(i+2,j+1,:)
     . + 28*qp(i+2,j+2,:))/(945.*hy**2)

      qxy(i,j,:) = (qp(i-1,j-1,:) + 2*qp(i-1,j-2,:) - qp(i-1,j+1,:)
     .- 2*qp(i-1,j+2,:) + 2*qp(i-2,j-1,:) + 4*qp(i-2,j-2,:)
     .- 2*qp(i-2,j+1,:) - 4*qp(i-2,j+2,:)- qp(i+1,j-1,:)
     . - 2*qp(i+1,j-2,:) + qp(i+1,j+1,:) + 2*qp(i+1,j+2,:)
     .- 2*qp(i+2,j-1,:) - 4*qp(i+2,j-2,:) + 2*qp(i+2,j+1,:)
     .+ 4*qp(i+2,j+2,:))/(100.*hx*hy)

      qxx(i,j,:) = (-350*qp(i,j,:) - 10*qp(i,j-1,:) - 40*qp(i,j-2,:)
     . - 10*qp(i,j+1,:) - 40*qp(i,j+2,:) + 17*qp(i-1,j,:)
     . + 7*qp(i-1,j-1,:) - 23*qp(i-1,j-2,:) + 7*qp(i-1,j+1,:)
     . - 23*qp(i-1,j+2,:) + 68*qp(i-2,j,:) + 58*qp(i-2,j-1,:)
     . + 28*qp(i-2,j-2,:) + 58*qp(i-2,j+1,:) + 28*qp(i-2,j+2,:)
     .+ 17*qp(i+1,j,:) + 7*qp(i+1,j-1,:) - 23*qp(i+1,j-2,:)
     . + 7*qp(i+1,j+1,:) - 23*qp(i+1,j+2,:) + 68*qp(i+2,j,:)
     . + 58*qp(i+2,j-1,:) + 28*qp(i+2,j-2,:) + 58*qp(i+2,j+1,:)
     . + 28*qp(i+2,j+2,:))/(945.*hx**2)

         qy(i,j,:)  = (-288*qp(i,j-1,:) - 65*qp(i,j-2,:)
     .+ 288*qp(i,j+1,:) + 65*qp(i,j+2,:) - 263*qp(i-1,j-1,:)
     . - 15*qp(i-1,j-2,:) + 263*qp(i-1,j+1,:) + 15*qp(i-1,j+2,:)
     .- 188*qp(i-2,j-1,:) + 135*qp(i-2,j-2,:) + 188*qp(i-2,j+1,:)
     .- 135*qp(i-2,j+2,:) - 263*qp(i+1,j-1,:) - 15*qp(i+1,j-2,:)
     .+ 263*qp(i+1,j+1,:) + 15*qp(i+1,j+2,:) - 188*qp(i+2,j-1,:)
     . + 135*qp(i+2,j-2,:) + 188*qp(i+2,j+1,:)
     . - 135*qp(i+2,j+2,:))/(1680.*hy)
         qx(i,j,:)  = (-288*qp(i-1,j,:) - 263*qp(i-1,j-1,:)
     .- 188*qp(i-1,j-2,:) - 263*qp(i-1,j+1,:) - 188*qp(i-1,j+2,:)
     .- 65*qp(i-2,j,:) - 15*qp(i-2,j-1,:) + 135*qp(i-2,j-2,:)
     .- 15*qp(i-2,j+1,:) + 135*qp(i-2,j+2,:) + 288*qp(i+1,j,:)
     .+ 263*qp(i+1,j-1,:) + 188*qp(i+1,j-2,:) + 263*qp(i+1,j+1,:)
     .+ 188*qp(i+1,j+2,:) + 65*qp(i+2,j,:) + 15*qp(i+2,j-1,:)
     .- 135*qp(i+2,j-2,:) + 15*qp(i+2,j+1,:)
     .- 135*qp(i+2,j+2,:))/(1680.*hx)

!            if(i .eq. 10 .and. j .eq. 10) then
!            print *, "here"
!            endif
 434    continue








       do 735 i = 1, mitot
       do 735 j = 1, mjtot
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

            do 7334 ioff = -iir(i,j),iir(i,j)
            do 7334 joff = -jjr(i,j),jjr(i,j)
                if (.not. IS_REAL(i+ioff,j+joff)) go to 7334
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 7334


                if(koff .eq. lstgrd) then
                    xoff = xlow + (dfloat(i+ioff)-.5d0)*hx
                    yoff = ylow + (dfloat(j+joff)-.5d0)*hy
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
  720               continue
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


      ! SOLVE THE NORMAL EQUATIONS.
        call cholesky(9, 9, ata, G)
!        if(i .eq. 9 .and. j .eq. 24) then
!            print *, ata
!            print *, G
!            print *, rhs(:,1)
!        endif
        do mm = 1, nvar
            call solve(9,9, G, rhs(:,mm))
                qyyy(i,j,mm) = 6.d0*rhs(9,mm)/(hy**3)
                qxyy(i,j,mm) = 2.d0*rhs(8,mm)/(hx*hy**2)
                qxxy(i,j,mm) = 2.d0*rhs(7,mm)/(hy*hx**2)
                qxxx(i,j,mm) = 6.d0*rhs(6,mm)/(hx**3)
                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2)
                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2)
                qy(i,j,mm)  =  rhs(2,mm)/hy
                qx(i,j,mm)  =  rhs(1,mm)/hx

!                if(mm .eq. 1) then
!                print *, "here"
!                endif
        enddo


 735    continue
      end subroutine





      subroutine slopes2(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)

      logical IS_GHOST
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)
       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)
       dimension temprhs(5), tempata(5,5)
       common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold, ilts,ibc


       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"
       do 234 i = 1, mitot
       do 234 j = 1, mjtot
            if(inuf(i,j) .eq. 0) cycle

            qyy(i,j,:) =  (qp(i,j+1,:)-2*qp(i,j,:)+qp(i,j-1,:))/hy**2
            qxy(i,j,:) =  (qp(i+1,j+1,:)-qp(i+1,j-1,:)-qp(i-1,j+1,:)
     .                     +qp(i-1,j-1,:))/(4*hx*hy)
            qxx(i,j,:) =  (qp(i+1,j,:)-2*qp(i,j,:)+qp(i-1,j,:))/hx**2
            qy(i,j,:)  =  (qp(i,j+1,:) - qp(i,j-1,:))/(2*hy)
            qx(i,j,:)  =  (qp(i+1,j,:) - qp(i-1,j,:))/(2*hx)
 234    continue

!      qyy = 0.d0
!      qxy = 0.d0
!      qxx = 0.d0


!       do 235 i = lwidth, mitot-lwidth+1
!       do 235 j = lwidth, mjtot-lwidth+1
       do 235 i = 1, mitot
       do 235 j = 1, mjtot
            if (ibc .eq. 1 .and. IS_GHOST(i,j)) go to 235

            kirr = irr(i,j)
            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
            if(kirr .eq. lstgrd) then
                x0 = xlow + (dfloat(i)-.5d0)*hx
                y0 = ylow + (dfloat(j)-.5d0)*hy
            else
                x0 = xcirr(kirr)
                y0 = ycirr(kirr)
            endif


      a = 0.d0
      ata = 0.d0
      rhs = 0.d0

!            do 2394 mmoff = 1,3

            do 2334 ioff = -iir(i,j),iir(i,j)
            do 2334 joff = -jjr(i,j),jjr(i,j)
!                idxmax = max(abs(ioff),abs(joff))
!                if(idxmax .ne. mmoff) goto 2334
            if (.not. IS_REAL(i+ioff,j+joff)) go to 2334

                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 2334

               ! skip ghost cells when using better BCs
               if(ibc .eq. 1 .and. IS_GHOST(i+ioff,j+joff))  go to 2334




                if(koff .eq. lstgrd) then
                    xoff = xlow + (dfloat(i+ioff)-.5d0)*hx
                    yoff = ylow + (dfloat(j+joff)-.5d0)*hy
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

               a(3) = a(3)
     .+ (artri/arr) * wtri(itq) *( ((xval-x0)/hx)**2.d0 )

                a(4) = a(4) + (artri/arr) * wtri(itq) *( (xval-x0)/hx )
     .                                                *( (yval-y0)/hy )

         a(5) = a(5) + (artri/arr) * wtri(itq) *( ((yval-y0)/hy)**2.d0 )

  222        continue ! for each quadrature point on each triangle
  221      continue ! for each triangle


            do ii = 1,5
            do jj = 1,5
                ata(ii,jj) = ata(ii,jj) + a(ii)*a(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .                      + a(ii) * (qp(i+ioff,j+joff,:) - qp(i,j,:))
            enddo

!            if( i .eq. 4 .and. j .eq. 23) then
!              print *,"4,23", qp(i+ioff,j+joff,1) - qp(i,j,1),qp(i,j,1)
!            endif
!            if( i .eq. 5 .and. j .eq. 23) then
!              print *,"5,23", qp(i+ioff,j+joff,1) - qp(i,j,1),qp(i,j,1)
!            endif


 2334         continue

!      if( i .eq. 4 .and. j .eq. 23) then
!      print *,""
!      temprhs(1:5) = rhs(1:5,1)
!      tempata(1:5, 1:5) = ata(1:5,1:5)
!      endif
!      if( i .eq. 5 .and. j .eq. 23) then
!      print *,"here"
!      temprhs(1) = rhs(1,1) + temprhs(1)
!      temprhs(2) = rhs(2,1) - temprhs(2)
!      temprhs(3) = rhs(3,1) - temprhs(3)
!      temprhs(4) = rhs(4,1) + temprhs(4)
!      temprhs(5) = rhs(5,1) - temprhs(5)
!
!      tempata = tempata - ata(1:5,1:5)
!      endif

      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 5, ata, G)
             do mm = 1, nvar
                call solve(9,5,G,rhs(:,mm))

                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2.d0)
                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2.d0)
                qy(i,j,mm)  =  rhs(2,mm)/hy
                qx(i,j,mm)  =  rhs(1,mm)/hx


!                if(abs(qy(i,j,mm)) > 1e-10) then
!                print *,""
!                endif
            end do


 235    continue



      end subroutine



      subroutine slopes2_ho(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
      logical IS_REAL, IS_GHOST
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)

      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)

       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)
       dimension qp0(4), qm0(4), q0p(4), q0m(4)
       dimension qpp(4), qmp(4), qpm(4), qmm(4)
       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"
       common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold, ilts,ibc

       do 454 i = 1, mitot
       do 454 j = 1, mjtot
            if(inuf(i,j) .eq. 0) cycle

            qyy(i,j,:) =  (qp(i,j+1,:)-2*qp(i,j,:)+qp(i,j-1,:))/hy**2
            qxy(i,j,:) =  (qp(i+1,j+1,:)-qp(i+1,j-1,:)-qp(i-1,j+1,:)
     .                     +qp(i-1,j-1,:))/(4*hx*hy)
            qxx(i,j,:) =  (qp(i+1,j,:)-2*qp(i,j,:)+qp(i-1,j,:))/hx**2
            qy(i,j,:)  =  (qp(i,j+1,:) - qp(i,j-1,:))/(2*hy)
            qx(i,j,:)  =  (qp(i+1,j,:) - qp(i-1,j,:))/(2*hx)
!      qp0 = qp(i+1,j  ,:)-qp(i,j,:)
!      qm0 = qp(i-1,j  ,:)-qp(i,j,:)
!      q0p = qp(i  ,j+1,:)-qp(i,j,:)
!      q0m = qp(i  ,j-1,:)-qp(i,j,:)
!      qpp = qp(i+1,j+1,:)-qp(i,j,:)
!      qmp = qp(i-1,j+1,:)-qp(i,j,:)
!      qpm = qp(i+1,j-1,:)-qp(i,j,:)
!      qmm = qp(i-1,j-1,:)-qp(i,j,:)

!      qx(i,j,:) = (1.d0/6.d0)*qp0/hx-(1.d0/6.d0)*qm0/hx
!     .+(1.d0/6.d0)*qpp/hx+(1.d0/6.d0)*qpm/hx
!     .-(1.d0/6.d0)*qmm/hx-(1.d0/6.d0)*qmp/hx
!
!      qxx(i,j,:) =(3.d0/5.d0)*qp0/hx**2.d0+(3.d0/5.d0)*qm0/hx**2.d0
!     .-(2.d0/5.d0)*q0p/hx**2.d0-(2.d0/5.d0)*q0m/hx**2.d0
!     .+(1.d0/5.d0)*qpp/hx**2.d0
!     .+(1.d0/5.d0)*qpm/hx**2.d0+(1.d0/5.d0)*qmm/hx**2.d0
!     .+(1.d0/5.d0)*qmp/hx**2.d0
!
!
!      qxy(i,j,:) = (1.d0/4.d0)*qpp/(hx*hy)-(1.d0/4.d0)*qpm/(hx*hy)
!     .+(1.d0/4.d0)*qmm/(hx*hy)-(1.d0/4.d0)*qmp/(hx*hy)
!
!      qy(i,j,:) = (1.d0/6.d0)*q0p/hy-(1.d0/6.d0)*q0m/hy
!     .+(1.d0/6.d0)*qpp/hy-(1.d0/6.d0)*qpm/hy
!     .-(1.d0/6.d0)*qmm/hy+(1.d0/6.d0)*qmp/hy
!
!      qyy(i,j,:) = -(2.d0/5.d0)*qp0/hy**2.d0-(2.d0/5.d0)*qm0/hy**2.d0
!     .+(3.d0/5.d0)*q0p/hy**2.d0
!     .+(3.d0/5.d0)*q0m/hy**2.d0+(1.d0/5.d0)*qpp/hy**2.d0
!     .+(1.d0/5.d0)*qpm/hy**2.d0+(1.d0/5.d0)*qmm/hy**2.d0
!     .+(1.d0/5.d0)*qmp/hy**2.d0





 454    continue

       do 435 i = 1, mitot
       do 435 j = 1, mjtot
            if (ibc .eq. 1 .and. IS_GHOST(i,j)) go to 435

            kirr = irr(i,j)
            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
            if(kirr .eq. lstgrd) then
                x0 = xlow + (i-.5d0)*hx
                y0 = ylow + (j-.5d0)*hy
            else
                x0 = xcirr_ho(kirr)
                y0 = ycirr_ho(kirr)
            endif


      a = 0.d0
      ata = 0.d0
      rhs = 0.d0

            do 4334 ioff = -iir(i,j),iir(i,j)
            do 4334 joff = -jjr(i,j),jjr(i,j)
!                idxmax = max(abs(ioff),abs(joff))
!                if(idxmax .ne. mmoff) goto 2334
            if (.not. IS_REAL(i+ioff,j+joff)) go to 4334

                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 4334


            if(ibc .eq. 1 .and. IS_GHOST(i+ioff,j+joff))  go to 4334


                if(koff .eq. lstgrd) then
                    xoff = xlow + (dfloat(i+ioff)-.5d0)*hx
                    yoff = ylow + (dfloat(j+joff)-.5d0)*hy
                else
                    xoff = xcirr_ho(koff)
                    yoff = ycirr_ho(koff)
                endif

                deltax = xoff - x0
                deltay = yoff - y0

                a(1) = deltax/hx
                a(2) = deltay/hy


         ! load the offsets
         a(3) = -poly_ho(8 ,1,kirr)
         a(4) = -poly_ho(9 ,1,kirr)
         a(5) = -poly_ho(10,1,kirr)



        if(koff .ne. lstgrd) then
          arr = ar_ho(koff)
          ivert = 1
          do 420 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
  420      continue
           itri = ivert - 3
        else
          arr = hx*hy
          itri = 2+1 ! added an extra one artificially
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
          do 421 it = 1, itri-1 ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 422 itq = 1,ntriquad
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
  422        continue ! for each quadrature point on each triangle
  421      continue ! for each triangle

      if(koff .ne. lstgrd) then
      x1 = poly(ivert-2,1,koff)
      y1 = poly(ivert-2,2,koff)

      x2 = bdry(1,1,koff)
      y2 = bdry(1,2,koff)

      x3 = bdry(3,1,koff)
      y3 = bdry(3,2,koff)

      x4 = bdry(2,1,koff)
      y4 = bdry(2,2,koff)

      do nq = 1,ntriquad_ho
      artri = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)

      xval = rs2xy_2(x1,x2,x3,x4,nq)
      yval = rs2xy_2(y1,y2,y3,y4,nq)

      a(3) = a(3) + (artri/arr) * wtri_ho(nq) *( (xval-x0)**2)
     .                                                / (hx**2)

      a(4) = a(4) + (artri/arr) * wtri_ho(nq) *( (xval-x0)/hx )
     .                                                *( (yval-y0)/hy )

      a(5) = a(5) + (artri/arr) * wtri_ho(nq) *( (yval-y0)**2 )
     .                                                / (hy**2)
      enddo
      endif





      do ii = 1,5
      do jj = 1,5
          ata(ii,jj) = ata(ii,jj) + a(ii)*a(jj)
      enddo

          rhs(ii,:) = rhs(ii,:)
     .                + a(ii) * (qp(i+ioff,j+joff,:) - qp(i,j,:))
      enddo


 4334         continue




      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 5, ata, G)
             do mm = 1, nvar
                call solve(9,5,G,rhs(:,mm))

                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2)
                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2)
                qy(i,j,mm)  =  rhs(2,mm)/hy
                qx(i,j,mm)  =  rhs(1,mm)/hx


!                if(abs(qy(i,j,mm)) > 1e-10) then
!                print *,""
!                endif
            end do


 435    continue
      end subroutine






























      subroutine slopes2_ho_qr(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
      logical IS_REAL, IS_GHOST
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)

      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)

       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension dmat(100,5), vec(100,nvar), sol(5)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)
       dimension qp0(4), qm0(4), q0p(4), q0m(4)
       dimension qpp(4), qmp(4), qpm(4), qmm(4)
       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"
       common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold, ilts,ibc

       do 454 i = 1, mitot
       do 454 j = 1, mjtot
            if(inuf(i,j) .eq. 0) cycle

            qyy(i,j,:) =  (qp(i,j+1,:)-2*qp(i,j,:)+qp(i,j-1,:))/hy**2
            qxy(i,j,:) =  (qp(i+1,j+1,:)-qp(i+1,j-1,:)-qp(i-1,j+1,:)
     .                     +qp(i-1,j-1,:))/(4*hx*hy)
            qxx(i,j,:) =  (qp(i+1,j,:)-2*qp(i,j,:)+qp(i-1,j,:))/hx**2
            qy(i,j,:)  =  (qp(i,j+1,:) - qp(i,j-1,:))/(2*hy)
            qx(i,j,:)  =  (qp(i+1,j,:) - qp(i-1,j,:))/(2*hx)


 454    continue

       do 435 i = 1, mitot
       do 435 j = 1, mjtot
            if (ibc .eq. 1 .and. IS_GHOST(i,j)) go to 435

            kirr = irr(i,j)
            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
            if(kirr .eq. lstgrd) then
                x0 = xlow + (i-.5d0)*hx
                y0 = ylow + (j-.5d0)*hy
            else
                x0 = xcirr_ho(kirr)
                y0 = ycirr_ho(kirr)
            endif


      a = 0.d0
      ata = 0.d0
      rhs = 0.d0

      dmat = 0.d0
      vec = 0.d0

      neighnum = 0

            do 4334 ioff = -iir(i,j),iir(i,j)
            do 4334 joff = -jjr(i,j),jjr(i,j)
!                idxmax = max(abs(ioff),abs(joff))
!                if(idxmax .ne. mmoff) goto 2334
            if (.not. IS_REAL(i+ioff,j+joff)) go to 4334

                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 4334


            if(ibc .eq. 1 .and. IS_GHOST(i+ioff,j+joff))  go to 4334

            neighnum = neighnum + 1

                if(koff .eq. lstgrd) then
                    xoff = xlow + (dfloat(i+ioff)-.5d0)*hx
                    yoff = ylow + (dfloat(j+joff)-.5d0)*hy
                else
                    xoff = xcirr_ho(koff)
                    yoff = ycirr_ho(koff)
                endif

                deltax = xoff - x0
                deltay = yoff - y0


                dmat(neighnum,1) = deltax/hx
                dmat(neighnum,2) = deltay/hy


         ! load the offsets
         dmat(neighnum,3) =  -poly_ho(8 ,1,kirr)
         dmat(neighnum,4) =  -poly_ho(9 ,1,kirr)
         dmat(neighnum,5) =  -poly_ho(10,1,kirr)



        if(koff .ne. lstgrd) then
          arr = ar_ho(koff)
          ivert = 1
          do 420 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
  420      continue
           itri = ivert - 3
        else
          arr = hx*hy
          itri = 2+1 ! added an extra one artificially
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
          do 421 it = 1, itri-1 ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 422 itq = 1,ntriquad
                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))


        dmat(neighnum,3) = dmat(neighnum,3)
     .+ (artri/arr) * wtri(itq) *( (xval-x0)**2)/ (hx**2)
        dmat(neighnum,4) = dmat(neighnum,4)
     .+ (artri/arr) * wtri(itq) *( (xval-x0)/hx )*( (yval-y0)/hy )
        dmat(neighnum,5) = dmat(neighnum,5)
     .+ (artri/arr) * wtri(itq) *( (yval-y0)**2 )/ (hy**2)



  422        continue ! for each quadrature point on each triangle
  421      continue ! for each triangle

      if(koff .ne. lstgrd) then
      x1 = poly(ivert-2,1,koff)
      y1 = poly(ivert-2,2,koff)

      x2 = bdry(1,1,koff)
      y2 = bdry(1,2,koff)

      x3 = bdry(3,1,koff)
      y3 = bdry(3,2,koff)

      x4 = bdry(2,1,koff)
      y4 = bdry(2,2,koff)

      do nq = 1,ntriquad_ho
      artri = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)

      xval = rs2xy_2(x1,x2,x3,x4,nq)
      yval = rs2xy_2(y1,y2,y3,y4,nq)

      dmat(neighnum,3) = dmat(neighnum,3)
     .+ (artri/arr) * wtri_ho(nq) *( (xval-x0)**2)/ (hx**2)

      dmat(neighnum,4) = dmat(neighnum,4)
     .+ (artri/arr) * wtri_ho(nq) *( (xval-x0)/hx )*( (yval-y0)/hy )

      dmat(neighnum,5) = dmat(neighnum,5)
     . + (artri/arr) * wtri_ho(nq) *( (yval-y0)**2 )/ (hy**2)
      enddo
      endif



      vec(neighnum,:) = qp(i+ioff,j+joff,:) - qp(i,j,:)




 4334         continue


      sol = 0.d0
      do mm = 1, nvar
      call svd_solve ( neighnum, 5,
     . dmat(1:neighnum,1:5), vec(:,mm), sol )
      qyy(i,j,mm) =  2.d0*sol(5)/(hy**2)
      qxy(i,j,mm) =  sol(4)/(hx*hy)
      qxx(i,j,mm) =  2.d0*sol(3)/(hx**2)
      qy(i,j,mm)  =  sol(2)/hy
      qx(i,j,mm)  =  sol(1)/hx
      end do


      ! SOLVE THE NORMAL EQUATIONS.
!             call cholesky(9, 5, ata, G)
!             do mm = 1, nvar
!                call solve(9,5,G,rhs(:,mm))
!
!                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2)
!                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
!                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2)
!                qy(i,j,mm)  =  rhs(2,mm)/hy
!                qx(i,j,mm)  =  rhs(1,mm)/hx
!
!            end do


 435    continue
      end subroutine







      subroutine slopes1(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)
       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"









       do 34 i = 1, mitot
       do 34 j = 1, mjtot
            if(inuf(i,j) .eq. 0) cycle
            qy(i,j,:)  =  (qp(i,j+1,:) - qp(i,j-1,:))/(2*hy)
            qx(i,j,:)  =  (qp(i+1,j,:) - qp(i-1,j,:))/(2*hx)
 34    continue

!     standard approach

!       do 35 i = 1, mitot
!       do 35 j = 1, mjtot
!            kirr = irr(i,j)
!            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
!            if(kirr .eq. lstgrd) then
!                x0 = xlow + (i-.5d0)*hx
!                y0 = ylow + (j-.5d0)*hy
!            else
!                x0 = xcirr(kirr)
!                y0 = ycirr(kirr)
!            endif
!
!
!      a = 0.d0
!      ata = 0.d0
!      rhs = 0.d0
!
!!     nothing special, first order slopes
!        do 334 ioff = -iir(i,j),iir(i,j)
!        do 334 joff = -jjr(i,j),jjr(i,j)
!                if (.not. IS_REAL(i+ioff,j+joff)) go to 334
!                koff = irr(i+ioff,j+joff)
!                if (koff .eq. -1) go to 334
!
!                if(koff .eq. lstgrd) then
!                    xoff = xlow + (i+ioff-.5d0)*hx
!                    yoff = ylow + (j+joff-.5d0)*hy
!                else
!                    xoff = xcirr(koff)
!                    yoff = ycirr(koff)
!                endif
!
!                deltax = xoff - x0
!                deltay = yoff - y0
!
!                a(1) = deltax/hx
!                a(2) = deltay/hy
!
!                ata(1,1) = ata(1,1) + a(1) * a(1)
!                ata(1,2) = ata(1,2) + a(1) * a(2)
!                ata(2,2) = ata(2,2) + a(2) * a(2)
!                rhs(1,:) = rhs(1,:) + a(1) *
!     .                     (qp(i+ioff,j+joff,:) - qp(i,j,:))
!                rhs(2,:) = rhs(2,:) + a(2) *
!     .                     (qp(i+ioff,j+joff,:) - qp(i,j,:))
!
!
!!        print *, "top", deltax, deltay, qp(i+ioff,j+joff,1) - qp(i,j,1)
!
!
!
! 334         continue
!!      print *,""
!
!      ! SOLVE THE NORMAL EQUATIONS.
!             c11 = sqrt(ata(1,1))
!             c12 = ata(1,2) / c11
!             c22 = sqrt(ata(2,2) - c12**2)
!
!             do mm = 1, nvar
!               b(1) = rhs(1,mm) / c11
!               b(2) = (rhs(2,mm) - c12*b(1))/c22
!               qy(i,j,mm) = b(2) / c22
!               qx(i,j,mm) = (b(1) - c12*qy(i,j,mm))/c11
!             end do
!
!             qx(i,j,:) = qx(i,j,:) / hx
!             qy(i,j,:) = qy(i,j,:) / hy
! 35   continue


!     fancy approach, on cuts and neighbors of cuts

       do 35 i = 1, mitot
       do 35 j = 1, mjtot
            kirr = irr(i,j)
            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
            if(kirr .eq. lstgrd) then
                x0 = xlow + (i-.5d0)*hx
                y0 = ylow + (j-.5d0)*hy
            else
                x0 = xcirr(kirr)
                y0 = ycirr(kirr)
            endif














      AA = 0.d0
      a = 0.d0
      ata = 0.d0
      rhs = 0.d0
      rhsA = 0.d0
      ! found normal of edge that is closest to me
      ! now reconstruct
        do 335 ioff = -iir(i,j),iir(i,j)
        do 335 joff = -jjr(i,j),jjr(i,j)
                if (.not. IS_REAL(i+ioff,j+joff)) go to 335
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 335


                if(koff .eq. lstgrd) then
                    xoff = xlow + (dfloat(i+ioff)-.5d0)*hx
                    yoff = ylow + (dfloat(j+joff)-.5d0)*hy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif


                ! standard LSQ for rhou and rhov
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


 335  continue

!      ! SOLVE THE NORMAL EQUATIONS.
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






      end subroutine


      subroutine slopes1qp(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)
       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"
       do 634 i = 1, mitot
       do 634 j = 1, mjtot
            if(inuf(i,j) .eq. 0) cycle
            qy(i,j,:)  =  (qp(i,j+1,:) - qp(i,j-1,:))/(2*hy)
            qx(i,j,:)  =  (qp(i+1,j,:) - qp(i-1,j,:))/(2*hx)
 634    continue

       do 635 i = 1, mitot
       do 635 j = 1, mjtot
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


            do 6334 ioff = -iir(i,j),iir(i,j)
            do 6334 joff = -jjr(i,j),jjr(i,j)
                if (.not. IS_REAL(i+ioff,j+joff)) go to 6334
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 6334


                if(koff .eq. lstgrd) then
                    xoff = xlow + (dfloat(i+ioff)-.5d0)*hx
                    yoff = ylow + (dfloat(j+joff)-.5d0)*hy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

                deltax = xoff - x0
                deltay = yoff - y0

            a(1) =  deltax/hx
            a(2) =  deltay/hy
            a(3) = (deltax/hx)**2.d0
            a(4) = (deltax/hx)*(deltay/hy)
            a(5) = (deltay/hy)**2.d0

            do ii = 1,5
            do jj = 1,5
                ata(ii,jj) = ata(ii,jj) + a(ii)*a(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .                      + a(ii) * (qp(i+ioff,j+joff,:) - qp(i,j,:))
            enddo

 6334         continue




      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 5, ata, G)
             do mm = 1, nvar
                call solve(9,5,G,rhs(:,mm))

                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2)
                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2)
                qy(i,j,mm)  =  rhs(2,mm)/hy
                qx(i,j,mm)  =  rhs(1,mm)/hx


            end do


 635    continue

                qyy = 0.d0
                qxy = 0.d0
                qxx = 0.d0
      end subroutine




      subroutine slopes1q(qp,qx,qy,qxx,qxy,qyy,
     .                   mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension qpr(4)
       dimension qxx(mitot,mjtot,4), qxy(mitot,mjtot,4),
     .           qyy(mitot,mjtot,4)
      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
       dimension ata(9,9), rhs(9,nvar), a(9), g(9,9)
       dimension AA(9,9), rhsA(9,nvar)
       dimension b(9), w(9,nvar)
       include "cirr.i"
       include "./quadrature.i"
       include "./reconinfo.i"
       do 634 i = 1, mitot
       do 634 j = 1, mjtot
            if(inuf(i,j) .eq. 0) cycle
            qy(i,j,:)  =  (qp(i,j+1,:) - qp(i,j-1,:))/(2*hy)
            qx(i,j,:)  =  (qp(i+1,j,:) - qp(i-1,j,:))/(2*hx)
 634    continue

       do 635 i = 1, mitot
       do 635 j = 1, mjtot
            kirr = irr(i,j)
            if(inuf(i,j) .eq. 1 .or. kirr .eq. -1) cycle
            if(kirr .eq. lstgrd) then
                x0 = xlow + (dfloat(i)-.5d0)*hx
                y0 = ylow + (dfloat(j)-.5d0)*hy
            else
                x0 = xcirr(kirr)
                y0 = ycirr(kirr)
            endif


      a = 0.d0
      ata = 0.d0
      rhs = 0.d0


            do 6334 ioff = -iir(i,j),iir(i,j)
            do 6334 joff = -jjr(i,j),jjr(i,j)

            if (.not. IS_REAL(i+ioff,j+joff)) go to 6334

                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 6334





                if(koff .eq. lstgrd) then
                    xoff = xlow + (dfloat(i+ioff)-.5d0)*hx
                    yoff = ylow + (dfloat(j+joff)-.5d0)*hy
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
          do 620 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
  620      continue
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
          do 621 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 622 itq = 1,ntriquad
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
  622        continue ! for each quadrature point on each triangle
  621      continue ! for each triangle


            do ii = 1,5
            do jj = 1,5
                ata(ii,jj) = ata(ii,jj) + a(ii)*a(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .                      + a(ii) * (qp(i+ioff,j+joff,:) - qp(i,j,:))
            enddo

 6334         continue




      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 5, ata, G)
             do mm = 1, nvar
                call solve(9,5,G,rhs(:,mm))

                qyy(i,j,mm) =  2.d0*rhs(5,mm)/(hy**2)
                qxy(i,j,mm) =  rhs(4,mm)/(hx*hy)
                qxx(i,j,mm) =  2.d0*rhs(3,mm)/(hx**2)
                qy(i,j,mm)  =  rhs(2,mm)/hy
                qx(i,j,mm)  =  rhs(1,mm)/hx


            end do


 635    continue

                qyy = 0.d0
                qxy = 0.d0
                qxx = 0.d0
      end subroutine
