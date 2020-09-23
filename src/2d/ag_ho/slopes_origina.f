c
c------------------------------------------------------------
c
       subroutine slopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension dumax(nvar),dumin(nvar),phimin(nvar), graddot(nvar)
       dimension dalpha(nvar), recon(nvar), b(2)
       logical  regular, quad, nolimiter
       logical MC/.true./
       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
       common /order2/ ssw, quad, nolimiter
       logical IS_GHOST
       dimension a(2,2), rhs(2,nvar)
       include "cirr.i"
c      regular(i,j) = ((i.gt. lwidth).and.(i.le.mitot-lwidth).and.
c    &                 (j.gt. lwidth).and.(j.le.mjtot-lwidth))

c
c      # ssw = slope switch (1. for slopes, 0 for donor cell 0 slopes)
c      # now set in amrcart
c
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   compute slopes using muscl limiter
c   qp contains the primitive variables
c
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)
c
c  initialize all slopes to 0 first
c
      qx = 0.d0
      qy = 0.d0

      isloperecon = 1
      ! 0 is no slopes
      ! 1 is linear
      ! 2 is quadratic

      inumneigh = 4
      ! 4 neighbors
      ! 8 neighbors

      if(isloperecon .eq. 0) then ! no slopes
       do 5 i = 1, mitot
       do 5 j = 1, mjtot
       do 5 m = 1, nvar 
          qx(i,j,m) = 0.d0
          qy(i,j,m) = 0.d0
 5     continue


      elseif(isloperecon .eq. 2) then
      ! LSQ 4 neighbors

       do 34 i = 2, mitot-1
       do 34 j = 2, mjtot-1
            if (irr(i,j) .ne. lstgrd) goto 34

            rhs = 0.d0 ! initialize for accumulation
            a = 0.d0
            k = irr(i,j)
            if(k .eq. lstgrd) then
                x0 = xlow + (i-.5d0)*hx
                y0 = ylow + (j-.5d0)*hy
            else
                x0 = xcirr(k)
                y0 = ycirr(k)
            endif


            do 22 joff = -1, 1
            do 22 ioff = -1, 1
                if (ioff .eq. 0 .and. joff .eq. 0) goto 22 ! no eqn to solve

                if(inumneigh .eq. 4) then ! only 4 neighbors
                    if(abs(joff) .eq. 1 .and. abs(ioff) .eq. 1) goto 22
                endif

                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) goto 22

                if(koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*hx
                    yoff = ylow + (j+joff-.5d0)*hy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

               deltax =  xoff - x0
               deltay =  yoff - y0
               a(1,1) = a(1,1) + deltax*deltax
               a(1,2) = a(1,2) + deltax*deltay
               a(2,2) = a(2,2) + deltay*deltay
               rhs(1,:) = rhs(1,:) + deltax *
     .                    (qp(i+ioff,j+joff,:) - qp(i,j,:))
               rhs(2,:) = rhs(2,:) + deltay *
     .                    (qp(i+ioff,j+joff,:) - qp(i,j,:))
 22          continue

             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
             ! will have to add robustness checks
             c11 = sqrt(a(1,1))
             c12 = a(1,2) / c11
             c22 = sqrt(a(2,2) - c12**2)

             ! now back solve (C^t C = rhs of A^tdu ) to get x and y gradient for all variables
             do mm = 1, nvar
               b(1) = rhs(1,mm) / c11
               b(2) = (rhs(2,mm) - c12*b(1))/c22
               qy(i,j,mm) = b(2) / c22
               qx(i,j,mm) = (b(1) - c12*qy(i,j,mm))/c11
             end do

 34    continue


      endif





 99    continue
       return
       end
