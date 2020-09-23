c
c ----------------------------------------------------------
c

       subroutine vctoprm(mitot,mjtot,meqn,
     . hx, hy, lstgrd, xlow, ylow, irr,
     . q,qp,qx,qy,qxx, qxy, qyy, qxxx, qxxy, qxyy, qyyy)
c
c      convert to primitie variables. note for future
c      that may be overwriting conserved variables
c
       implicit double precision(a-h,o-z)
       include "./quadrature.i"
       include "./cirr.i"
      common  /order2/  ssw, quad, nolimiter
      logical quad,nolimiter
      dimension q(mitot,mjtot,4),qp(mitot,mjtot,4)
      dimension qx(mitot,mjtot,meqn),qy(mitot,mjtot,meqn)
      dimension qxx(mitot,mjtot,meqn),qyy(mitot,mjtot,meqn)
      dimension qxy(mitot,mjtot,meqn)
      dimension qxxx(mitot,mjtot,meqn),qxyy(mitot,mjtot,meqn)
      dimension qxxy(mitot,mjtot,meqn),qyyy(mitot,mjtot,meqn)
       dimension Uout(4)
       dimension irr(mitot,mjtot)
       common   /userdt/  gamma,dgamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold


        select case(ihob)
        case(0)
      if(ssw .eq. -1 .or. ssw .eq. -10 .or. ssw .eq. 1) then
           call vctoprm_centroids(q,qp,mitot,mjtot, irr)
      else
           call vctoprm_lo(mitot,mjtot,meqn,
     . hx, hy, lstgrd, xlow, ylow, irr,
     . q,qp,qx,qy,qxx, qxy, qyy, qxxx, qxxy, qxyy, qyyy)
      endif
        case(1)
           call vctoprm_ho(mitot,mjtot,meqn,
     . hx, hy, lstgrd, xlow, ylow, irr,
     . q,qp,qx,qy,qxx, qxy, qyy, qxxx, qxxy, qxyy, qyyy)

        end select

       return
       end





       subroutine vctoprm_lo(mitot,mjtot,meqn,
     . hx, hy, lstgrd, xlow, ylow, irr,
     . q,qp,qx,qy,qxx, qxy, qyy, qxxx, qxxy, qxyy, qyyy)
c
c      convert to primitie variables. note for future
c      that may be overwriting conserved variables
c
       implicit double precision(a-h,o-z)
       include "./quadrature.i"
       include "./cirr.i"
       dimension q(mitot,mjtot,4),qp(mitot,mjtot,4)
      dimension qx(mitot,mjtot,meqn),qy(mitot,mjtot,meqn)
      dimension qxx(mitot,mjtot,meqn),qyy(mitot,mjtot,meqn)
      dimension qxy(mitot,mjtot,meqn)
      dimension qxxx(mitot,mjtot,meqn),qxyy(mitot,mjtot,meqn)
      dimension qxxy(mitot,mjtot,meqn),qyyy(mitot,mjtot,meqn)
       dimension Uout(4)
       dimension irr(mitot,mjtot)
       common   /userdt/  gamma,dgamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c


       gamma1 = .4d0
       do 10 j = 1, mjtot
       do 10 i = 1, mitot
          kirr = irr(i,j)
          if(kirr .eq. -1) goto 10
          if(kirr .ne. lstgrd) then
            arr = ar(kirr)
            ivert = 1
            do 720 while (poly(ivert+1,1,kirr) .ne. -11.)
              ivert = ivert + 1
  720        continue
             itri = ivert - 3

            x0 = xcirr(kirr)
            y0 = ycirr(kirr)
          else
            arr = hx*hy
            itri = 2
            poly(1,1,kirr) = xlow + (dfloat(i)-1.d0)*hx
            poly(1,2,kirr) = ylow + (dfloat(j)-1.d0)*hy

            poly(2,1,kirr) = xlow + (dfloat(i)-0.d0)*hx
            poly(2,2,kirr) = ylow + (dfloat(j)-1.d0)*hy

            poly(3,1,kirr) = xlow + (dfloat(i)-0.d0)*hx
            poly(3,2,kirr) = ylow + (dfloat(j)-0.d0)*hy

            poly(4,1,kirr) = xlow + (dfloat(i)-1.d0)*hx
            poly(4,2,kirr) = ylow + (dfloat(j)-0.d0)*hy


             x0 = xlow + (dfloat(i)-.5d0)*hx
             y0 = ylow + (dfloat(j)-.5d0)*hy


          endif


          rho = 0.d0
          u = 0.d0
          v = 0.d0
          pr = 0.d0

          idx1 = 1
          do 721 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,kirr)
            y1 = poly(idx1,2,kirr)

            x2 = poly(idx2,1,kirr)
            y2 = poly(idx2,2,kirr)

            x3 = poly(idx3,1,kirr)
            y3 = poly(idx3,2,kirr)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 722 itq = 1,ntriquad
                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

                xdif = xval -x0
                ydif = yval -y0
        do m = 1,4
        call evalU(Uout,xdif, ydif, hx, hy,
     .                 q, qx, qy,
     .                 qxx, qxy, qyy,
     .                 qxxx, qxxy, qxyy, qyyy,
     .                 i, j, kirr, mitot, mjtot, 4, m)
        end do

!      rho = rho + (artri/arr) * wtri(itq) * Uout(1)
      u   = u   + (artri/arr) * wtri(itq) * Uout(2)/Uout(1)
      v   = v   + (artri/arr) * wtri(itq) * Uout(3)/Uout(1)
      pr  = pr
     . + (artri/arr) * wtri(itq) * gamma1*(Uout(4)-.5d0*(Uout(2)**2+
     .                                 Uout(3)**2)/Uout(1))

  722        continue ! for each quadrature point on each triangle
  721      continue ! for each triangle







          qp(i,j,1) = q(i,j,1)
          qp(i,j,2) = u
          qp(i,j,3) = v
          qp(i,j,4) = pr
 10    continue
c
       return
       end


       subroutine vctoprm_ho(mitot,mjtot,meqn,
     . hx, hy, lstgrd, xlow, ylow, irr,
     . q,qp,qx,qy,qxx, qxy, qyy, qxxx, qxxy, qxyy, qyyy)
c
c      convert to primitie variables. note for future
c      that may be overwriting conserved variables
c
       implicit double precision(a-h,o-z)
       include "./quadrature.i"
       include "./cirr.i"
       dimension q(mitot,mjtot,4),qp(mitot,mjtot,4)
      dimension qx(mitot,mjtot,meqn),qy(mitot,mjtot,meqn)
      dimension qxx(mitot,mjtot,meqn),qyy(mitot,mjtot,meqn)
      dimension qxy(mitot,mjtot,meqn)
      dimension qxxx(mitot,mjtot,meqn),qxyy(mitot,mjtot,meqn)
      dimension qxxy(mitot,mjtot,meqn),qyyy(mitot,mjtot,meqn)
       dimension Uout(4)
       dimension irr(mitot,mjtot)
       common   /userdt/  gamma,dgamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c

      common  /order2/  ssw, quad, nolimiter
      logical quad, nolimiter

       gamma1 = .4d0
       do 10 j = 1, mjtot
       do 10 i = 1, mitot
          kirr = irr(i,j)
          if(kirr .eq. -1) goto 10
          if(kirr .ne. lstgrd) then
            arr = ar_ho(kirr)
            ivert = 1
            do 720 while (poly(ivert+1,1,kirr) .ne. -11.)
              ivert = ivert + 1
  720        continue

             ! artificially reduce by 1
             itri = ivert - 3 - 1


            x0 = xcirr_ho(kirr)
            y0 = ycirr_ho(kirr)
          else
            arr = hx*hy
            itri = 2
            poly(1,1,kirr) = xlow + (dfloat(i)-1.d0)*hx
            poly(1,2,kirr) = ylow + (dfloat(j)-1.d0)*hy

            poly(2,1,kirr) = xlow + (dfloat(i)-0.d0)*hx
            poly(2,2,kirr) = ylow + (dfloat(j)-1.d0)*hy

            poly(3,1,kirr) = xlow + (dfloat(i)-0.d0)*hx
            poly(3,2,kirr) = ylow + (dfloat(j)-0.d0)*hy

            poly(4,1,kirr) = xlow + (dfloat(i)-1.d0)*hx
            poly(4,2,kirr) = ylow + (dfloat(j)-0.d0)*hy


             x0 = xlow + (i-.5d0)*hx
             y0 = ylow + (j-.5d0)*hy


          endif


          rho = 0.d0
          u = 0.d0
          v = 0.d0
          pr = 0.d0

          idx1 = 1
          do 721 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,kirr)
            y1 = poly(idx1,2,kirr)

            x2 = poly(idx2,1,kirr)
            y2 = poly(idx2,2,kirr)

            x3 = poly(idx3,1,kirr)
            y3 = poly(idx3,2,kirr)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 722 itq = 1,ntriquad
                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

                xdif = xval -x0
                ydif = yval -y0
        do m = 1,4
        call evalU(Uout,xdif, ydif, hx, hy,
     .                 q, qx, qy,
     .                 qxx, qxy, qyy,
     .                 qxxx, qxxy, qxyy, qyyy,
     .                 i, j, kirr, mitot, mjtot, 4, m)
        end do

      rho = rho + (artri/arr) * wtri(itq) * Uout(1)
      u   = u   + (artri/arr) * wtri(itq) * Uout(2)/Uout(1)
      v   = v   + (artri/arr) * wtri(itq) * Uout(3)/Uout(1)
      pr  = pr
     . + (artri/arr) * wtri(itq) * gamma1*(Uout(4)-.5d0*(Uout(2)**2+
     .                                 Uout(3)**2)/Uout(1))

  722        continue ! for each quadrature point on each triangle
  721      continue ! for each triangle


      if((ssw .eq. 2  .or. ssw .eq. -2) .and. kirr .ne. lstgrd) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)
      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)
      x3 = bdry(3,1,kirr)
      y3 = bdry(3,2,kirr)
      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)

      do nq = 1,ntriquad_ho
      artri = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)
      xval = rs2xy_2(x1,x2,x3,x4,nq)
      yval = rs2xy_2(y1,y2,y3,y4,nq)



        xdif = xval -x0
        ydif = yval -y0
        do m = 1,4
        call evalU(Uout,xdif, ydif, hx, hy,
     .                 q, qx, qy,
     .                 qxx, qxy, qyy,
     .                 qxxx, qxxy, qxyy, qyyy,
     .                 i, j, kirr, mitot, mjtot, 4, m)
        end do

      rho = rho + (artri/arr) * wtri_ho(nq) * Uout(1)
      u   = u   + (artri/arr) * wtri_ho(nq) * Uout(2)/Uout(1)
      v   = v   + (artri/arr) * wtri_ho(nq) * Uout(3)/Uout(1)
      pr  = pr
     . + (artri/arr) * wtri_ho(nq) * gamma1*(Uout(4)-.5d0*(Uout(2)**2+
     .                                 Uout(3)**2)/Uout(1))
      enddo
!      print *,du


      elseif((ssw .eq. 3 .or. ssw .eq. -3) .and. kirr .ne. lstgrd) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)
      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)
      x3 = bdry(4,1,kirr)
      y3 = bdry(4,2,kirr)
      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)
      x5 = bdry(3,1,kirr)
      y5 = bdry(3,2,kirr)

      do nq = 1,ntriquad_ho
      artri = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)
      xval = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yval = rs2xy_3(y1,y2,y3,y4,y5,nq)

        xdif = xval -x0
        ydif = yval -y0
        do m = 1,4
        call evalU(Uout,xdif, ydif, hx, hy,
     .                 q, qx, qy,
     .                 qxx, qxy, qyy,
     .                 qxxx, qxxy, qxyy, qyyy,
     .                 i, j, kirr, mitot, mjtot, 4, m)
        end do

      rho = rho + (artri/arr) * wtri_ho(nq) * Uout(1)
      u   = u   + (artri/arr) * wtri_ho(nq) * Uout(2)/Uout(1)
      v   = v   + (artri/arr) * wtri_ho(nq) * Uout(3)/Uout(1)
      pr  = pr
     . + (artri/arr) * wtri_ho(nq) * gamma1*(Uout(4)-.5d0*(Uout(2)**2+
     .                                 Uout(3)**2)/Uout(1))


      enddo
      endif




          qp(i,j,1) = q(i,j,1)
          qp(i,j,2) = u
          qp(i,j,3) = v
          qp(i,j,4) = pr
 10    continue
c
       return
       end


c
c ----------------------------------------------------------
c
       subroutine vctoprm_centroids(q,qp,mitot,mjtot, irr)
c
c      convert to primitie variables. note for future
c      that may be overwriting conserved variables
c
       implicit double precision(a-h,o-z)
       dimension q(mitot,mjtot,4),qp(mitot,mjtot,4)
       dimension irr(mitot,mjtot)
       common   /userdt/  gamma,dgamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold

       gamma1 = .4d0
       do 10 j = 1, mjtot
       do 10 i = 1, mitot
          kirr = irr(i,j)
          if(kirr .eq. -1) goto 10

          rho       = q(i,j,1)
          u         = q(i,j,2)/rho
          v         = q(i,j,3)/rho
c         ei        = q(i,j,4)/rho - .5d0*(u**2+v**2)
c         pr        = gamma1*rho*ei
          pr        = gamma1*(q(i,j,4)-.5d0*(q(i,j,2)**2+
     .                q(i,j,3)**2)/rho)
c         qp(i,j,1) = rho
          qp(i,j,2) = u
          qp(i,j,3) = v
          qp(i,j,4) = pr
 10    continue
c
       return
       end
