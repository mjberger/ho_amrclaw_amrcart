c
c ---------------------------------------------------------------------
c
       subroutine SRD_cellMerge_U(q,dlimit,nvar,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
!
       implicit double precision (a-h, o-z)
       include "cirr.i"
       include "./quadrature.i"
       common /order2/ ssw, quad, nolimiter
       logical    quad, nolimiter
       dimension q(mitot,mjtot,nvar)
       dimension irr(mitot,mjtot)
       dimension qx(mitot,mjtot,nvar), qy(mitot,mjtot,nvar)
       dimension gradmx(irrsize,nvar), gradmy(irrsize,nvar)
       dimension gradmxx(irrsize,nvar), gradmxy(irrsize,nvar)
       dimension gradmyy(irrsize,nvar)
       dimension gradmxxx(irrsize,nvar), gradmxxy(irrsize,nvar)
       dimension gradmxyy(irrsize,nvar), gradmyyy(irrsize,nvar)
       dimension qMerge(mitot,mjtot, nvar), qm(nvar)

       dimension qmx(irrsize,nvar),qmy(irrsize,nvar)
       dimension a(9,9), rhs(9,nvar), b(9), f(9), w(9,nvar)
       dimension G(9,9)
       dimension qpr(4)
       dimension valnew(mitot,mjtot,nvar)
       dimension dumax(nvar),dumin(nvar),phimin(nvar)
       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                  ismp,gradThreshold

       logical IS_GHOST, verbose/.true./
       IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)


!       nco = 1 ! reconstruct slope using 3x3 neighborhood.

!
!
            call outgeom(lstgrd, mitot, mjtot, irr,xlow, ylow, dx, dy)
!
!
!
      TOL = 1.d-10

       if (verbose) then
          totmass =  bigconck(q,irr,mitot,mjtot,lwidth,nvar)
          write(*,911) totmass
 911      format(/,"         SRDmass before redistribution is ",e30.20)
       endif






      qMerge(:,:,:) = 0.d0
      qMerge(1:lwidth,:,:) =  q(1:lwidth,:,:)
      qMerge(mitot-lwidth+1:mitot,:,:) =  q(mitot-lwidth+1:mitot,:,:)
      qMerge(:,1:lwidth,:) =  q(:,1:lwidth,:)
      qMerge(:,mjtot-lwidth+1:mjtot,:) =  q(:,mjtot-lwidth+1:mjtot,:)


c       form qMerge vals, these are the same regardless of the order of accuracy!
      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth
            k = irr(i,j)

            if (k .eq. -1 .or. k .eq. lstgrd) then
              qMerge(i,j,:) = q(i,j,:)
              go to 10
            endif


             do 27 ic = 1, ncount(k)
                icurr = iidx(k,ic)
                jcurr = jidx(k,ic)
                koff = irr(icurr,jcurr)
                qMerge(i,j,:) = qMerge(i,j,:) + ar(koff)*
     .                        q(icurr,jcurr,:)/numHoods(icurr,jcurr)

 27          continue
             qMerge(i,j,:) = qMerge(i,j,:) / volMerge(k)
 10     continue





         ! gradient of merge neighborhoods, initialized to 0. set using neighboring merged tiels
         gradmx = 0.d0
         gradmy = 0.d0
         gradmxx = 0.d0
         gradmxy = 0.d0
         gradmyy = 0.d0

        if( ssw .eq. 1) then
        do 20 j = lwidth+1, mjtot-lwidth
        do 20 i = lwidth+1, mitot-lwidth
             k = irr(i,j)
             if (k .eq. -1 .or. k .eq. lstgrd) go to 20
             ! solid/regular cell doesnt participate in anything
             ! can reconstruct high order poly here, but for now do linear reconstruction
             ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save


             rhs = 0.d0 ! initialize for accumulation
             a = 0.d0
             x0 = xcentMerge(k)
             y0 = ycentMerge(k)
!             do 223 nco = 1,3
             do 22 joff = -mjoff(i,j), mjoff(i,j)
             do 22 ioff = -mioff(i,j), mioff(i,j)
!                 if (ioff .eq. 0 .and. joff .eq. 0) go to 22 ! no eqn to solve
                 if (IS_GHOST(i+ioff,j+joff)) go to 22
                 koff = irr(i+ioff,j+joff)
                 if (koff .eq. -1) go to 22

                 if(koff .eq. lstgrd) then
!                     call getCellCentroid(lstgrd,i+ioff,j+joff,
!     .                               xcoff,ycoff,xlow,ylow,dx,dy,koff)
                     xcoff = xlow + (i+ioff-0.5d0)*dx
                     ycoff = ylow + (j+joff-0.5d0)*dy
                     deltax = xcoff - x0
                     deltay = ycoff - y0
                 else
                     deltax = xcentMerge(koff) - x0
                     deltay = ycentMerge(koff) - y0
                 endif

                 a(1,1) = a(1,1) + deltax*deltax
                 a(1,2) = a(1,2) + deltax*deltay
                 a(2,2) = a(2,2) + deltay*deltay
                 rhs(1,:) = rhs(1,:) + deltax *
     .                      (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
                 rhs(2,:) = rhs(2,:) + deltay *
     .                      (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))



!           if(ireflect .eq. 1 .and. ioff .eq. 0 .and. joff .eq. 0
!     .         .and. irr(i,j) .ne. lstgrd) then
!            ! find irregular boundary of (i,j)
!                 k = irr(i,j)
!                 do 221 kside=1,6
!                    if (poly(kside+2,1,k).eq.-11.) then
!                       x1 = poly(kside,1,  k)
!                       y1 = poly(kside,2,  k)
!                       x2 = poly(kside+1,1,k)
!                       y2 = poly(kside+1,2,k)
!                       go to 251
!                    endif
!221      continue
!251      continue
!               rlen = dsqrt((y1-y2)**2 + (x1-x2)**2)
!               dnx = (y1-y2)/rlen
!               dny = (x2-x1)/rlen
!
!!                 dm = (y2-y1)/(x2-x1)
!                 if(koff .eq. lstgrd) then
!                     xoff = xlow + (i+ioff-0.5d0)*dx
!                     yoff = ylow + (j+joff-0.5d0)*dy
!                 else
!                     xoff = xcentMerge(koff)
!                     yoff = ycentMerge(koff)
!                 endif
!
!!                xoff0 = xoff - x1
!!                yoff0 = yoff - y1
!!
!!          xoff = ((1-dm**2)*xoff0 + 2.d0*dm *yoff0)/(1+dm**2) + x1
!!          yoff = (2.d0*dm * xoff0 + (dm**2-1.d0) *yoff0)/(1+dm**2) + y1
!
!!               rlen = dsqrt((y1-y2)**2 + (x1-x2)**2)
!!               dnx = (y1-y2)/rlen
!!               dny = (x2-x1)/rlen
!
!
!
!               distx = x0-x1
!               disty = y0-y1
!               dot = distx * dnx + disty * dny
!
!               xref = x0 - 2.d0 * dot * dnx
!               yref = y0 - 2.d0 * dot * dny
!
!
!               qpr = qMerge(i+ioff,j+joff,:)
!               dot = qMerge(i+ioff,j+joff,2) * dnx
!     .             + qMerge(i+ioff,j+joff,3) * dny
!               qpr(2) = qpr(2) - 2.d0 * dot * dnx
!               qpr(3) = qpr(3) - 2.d0 * dot * dny
!
!                deltax = xref - x0
!                deltay = yref - y0
!
!                 a(1,1) = a(1,1) + deltax*deltax
!                 a(1,2) = a(1,2) + deltax*deltay
!                 a(2,2) = a(2,2) + deltay*deltay
!                 rhs(1,:) = rhs(1,:) + deltax*(qpr - qMerge(i,j,:))
!                 rhs(2,:) = rhs(2,:) + deltay*(qpr - qMerge(i,j,:))
!
!!                ata(1,1) = ata(1,1) + a(1) * a(1)
!!                ata(1,2) = ata(1,2) + a(1) * a(2)
!!                ata(2,2) = ata(2,2) + a(2) * a(2)
!!                rhs(1,:) = rhs(1,:) + a(1) * (qpr - qp(i,j,:))
!!                rhs(2,:) = rhs(2,:) + a(2) * (qpr - qp(i,j,:))
!!        print *, "bot", deltax, deltay, qpr(1) - qp(i,j,1)
!            endif ! end of reflection






 22           continue

!              if(a(1,1) > TOL .and. (a(2,2) - c12**2) > TOL) then
                ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
                ! will have to add robustness checks
!                c11 = sqrt(a(1,1))
!                c12 = a(1,2) / c11
!                c22 = sqrt(a(2,2) - c12**2)

!                goto 125 ! backsolve
!              endif


! 223        continue

! 125        continue

                c11 = sqrt(a(1,1))
                c12 = a(1,2) / c11
                c22 = sqrt(a(2,2) - c12**2)
              ! now back solve (C^t C = rhs of A^tdu ) to get x and y gradient for all variables
              do m = 1, nvar
                b(1) = rhs(1,m) / c11
                b(2) = (rhs(2,m) - c12*b(1))/c22
                gradmy(k,m) = b(2) / c22
                gradmx(k,m) = (b(1) - c12*gradmy(k,m))/c11
              end do
!            print *,nco
  20     continue






        elseif (ssw .eq. 2) then






        do 120 j = lwidth+1, mjtot-lwidth
        do 120 i = lwidth+1, mitot-lwidth
             k = irr(i,j)
             if (k .eq. -1 .or. k .eq. lstgrd) go to 120
             ! solid/regular cell doesnt participate in anything
             ! can reconstruct high order poly here, but for now do linear reconstruction
             ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save

             rhs = 0.d0 ! initialize for accumulation
             a = 0.d0
             x0 = xcentMerge(k)
             y0 = ycentMerge(k)


            do 334 ioff = -mioff(i,j), mioff(i,j)
            do 334 joff = -mjoff(i,j), mjoff(i,j)
                 if (IS_GHOST(i+ioff,j+joff)) go to 334
                 koff = irr(i+ioff,j+joff)
                 if (koff .eq. -1) goto 334

                 if(koff .eq. lstgrd) then
                     xcoff = xlow + (i+ioff-0.5d0)*dx
                     ycoff = ylow + (j+joff-0.5d0)*dy
                     deltax = xcoff - x0
                     deltay = ycoff - y0
                 else
                     deltax = xcentMerge(koff) - x0
                     deltay = ycentMerge(koff) - y0
                 endif

                 f(1) = deltax/dx
                 f(2) = deltay/dy

                 f(3) = -dmergeshifts(1,k)
                 f(4) = -dmergeshifts(2,k)
                 f(5) = -dmergeshifts(3,k)

                 if(koff .eq. lstgrd) then ! set up data on whole cell
                    ncount(koff)  = 1
                    iidx(koff, 1) = i+ioff
                    jidx(koff, 1) = j+joff
                    volmerge(koff) = dx * dy / numhoods(i+ioff, j+joff)
                 endif


                 cvm = volmerge(koff) ! off vol merge
      do 297 ic = 1, ncount(koff) ! compute weighted inner product of monomials on this neighborhood
                    icurr = iidx(koff,ic)
                    jcurr = jidx(koff,ic)
                    kcurr = irr(icurr,jcurr)
                    nhc = numhoods(icurr, jcurr) ! num hoods current
                    if(kcurr .eq. lstgrd) then
                      itri = 2
                      poly(1,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(1,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(2,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(2,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(3,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(3,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy

                      poly(4,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(4,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy
                    else
                     ivert = 1
                      do 290 while (poly(ivert+1,1,kcurr) .ne. -11.)
                        ivert = ivert + 1
  290                   continue
                       itri = ivert - 3
                    endif

      ! computing the inner product on each triangle of neighborhood member
                  idx1 = 1
                  do 291 it = 1, itri ! for each  triangle
                    idx2 = it + 1
                    idx3 = it + 2

                    x1 = poly(idx1,1,kcurr)
                    y1 = poly(idx1,2,kcurr)

                    x2 = poly(idx2,1,kcurr)
                    y2 = poly(idx2,2,kcurr)

                    x3 = poly(idx3,1,kcurr)
                    y3 = poly(idx3,2,kcurr)

                    artri = area(x1, x2, x3, y1, y2, y3)

                    do 292 itq = 1,ntriquad
                        xval = x1 * rtri(itq) + x2 * stri(itq)
     .                     +  x3 * (1.d0-rtri(itq)-stri(itq))
                        yval = y1 * rtri(itq) + y2 * stri(itq)
     .                     +  y3 * (1.d0-rtri(itq)-stri(itq))

            f(3) = f(3) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)**2 )
     .                                                / (dx**2)

            f(4) = f(4) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)/dx )
     .                                                *( (yval-y0)/dy )

            f(5) = f(5) + (artri/nhc/cvm) * wtri(itq) *( (yval-y0)**2 )
     .                                                / (dy**2)
 292        continue ! for each quadrature point on each triangle
 291        continue ! for each triangle
 297        continue ! for each subelement of the current neighborhood
        ! integration is complete, let's accumulate

            do ii = 1,5
            do jj = 1,5
                a(ii,jj) = a(ii,jj) + f(ii)*f(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .               + f(ii) * (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
            enddo

 2334         continue


 334           continue



      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 5, a, G)
             do mm = 1, nvar
                call solve(9,5,G,rhs(:,mm))

                gradmyy(k,mm) =  2.d0*rhs(5,mm)/(dy**2)
                gradmxy(k,mm) =  rhs(4,mm)/(dx*dy)
                gradmxx(k,mm) =  2.d0*rhs(3,mm)/(dx**2)
                gradmy(k,mm)  =  rhs(2,mm)/dy
                gradmx(k,mm)  =  rhs(1,mm)/dx
            end do




 120    continue ! iterate over merging tiles on entire grid
!                gradmyy =  0.d0
!                gradmxy =  0.d0
!                gradmxx =  0.d0

        elseif (ssw .eq. -1) then






        do 820 j = lwidth+1, mjtot-lwidth
        do 820 i = lwidth+1, mitot-lwidth
             k = irr(i,j)
             if (k .eq. -1 .or. k .eq. lstgrd) go to 820
             ! solid/regular cell doesnt participate in anything
             ! can reconstruct high order poly here, but for now do linear reconstruction
             ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save

             rhs = 0.d0 ! initialize for accumulation
             a = 0.d0
             x0 = xcentMerge(k)
             y0 = ycentMerge(k)


            do 834 ioff = -mioff(i,j), mioff(i,j)
            do 834 joff = -mjoff(i,j), mjoff(i,j)
                 if (IS_GHOST(i+ioff,j+joff)) go to 834
                 koff = irr(i+ioff,j+joff)
                 if (koff .eq. -1) goto 834

                 if(koff .eq. lstgrd) then
                     xcoff = xlow + (i+ioff-0.5d0)*dx
                     ycoff = ylow + (j+joff-0.5d0)*dy
                     deltax = xcoff - x0
                     deltay = ycoff - y0
                 else
                     deltax = xcentMerge(koff) - x0
                     deltay = ycentMerge(koff) - y0
                 endif

                 f(1) = deltax/dx
                 f(2) = deltay/dy

                 f(3) = -dmergeshifts(1,k)
                 f(4) = -dmergeshifts(2,k)
                 f(5) = -dmergeshifts(3,k)

                 if(koff .eq. lstgrd) then ! set up data on whole cell
                    ncount(koff)  = 1
                    iidx(koff, 1) = i+ioff
                    jidx(koff, 1) = j+joff
                    volmerge(koff) = dx * dy / numhoods(i+ioff, j+joff)
                 endif


                 cvm = volmerge(koff) ! off vol merge
      do 897 ic = 1, ncount(koff) ! compute weighted inner product of monomials on this neighborhood
                    icurr = iidx(koff,ic)
                    jcurr = jidx(koff,ic)
                    kcurr = irr(icurr,jcurr)
                    nhc = numhoods(icurr, jcurr) ! num hoods current
                    if(kcurr .eq. lstgrd) then
                      itri = 2
                      poly(1,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(1,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(2,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(2,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(3,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(3,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy

                      poly(4,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(4,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy
                    else
                     ivert = 1
                      do 890 while (poly(ivert+1,1,kcurr) .ne. -11.)
                        ivert = ivert + 1
  890                   continue
                       itri = ivert - 3
                    endif

      ! computing the inner product on each triangle of neighborhood member
                  idx1 = 1
                  do 891 it = 1, itri ! for each  triangle
                    idx2 = it + 1
                    idx3 = it + 2

                    x1 = poly(idx1,1,kcurr)
                    y1 = poly(idx1,2,kcurr)

                    x2 = poly(idx2,1,kcurr)
                    y2 = poly(idx2,2,kcurr)

                    x3 = poly(idx3,1,kcurr)
                    y3 = poly(idx3,2,kcurr)

                    artri = area(x1, x2, x3, y1, y2, y3)

                    do 892 itq = 1,ntriquad
                        xval = x1 * rtri(itq) + x2 * stri(itq)
     .                     +  x3 * (1.d0-rtri(itq)-stri(itq))
                        yval = y1 * rtri(itq) + y2 * stri(itq)
     .                     +  y3 * (1.d0-rtri(itq)-stri(itq))

            f(3) = f(3) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)**2 )
     .                                                / (dx**2)

            f(4) = f(4) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)/dx )
     .                                                *( (yval-y0)/dy )

            f(5) = f(5) + (artri/nhc/cvm) * wtri(itq) *( (yval-y0)**2 )
     .                                                / (dy**2)
 892        continue ! for each quadrature point on each triangle
 891        continue ! for each triangle
 897        continue ! for each subelement of the current neighborhood
        ! integration is complete, let's accumulate

            do ii = 1,5
            do jj = 1,5
                a(ii,jj) = a(ii,jj) + f(ii)*f(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .               + f(ii) * (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
            enddo

 8334         continue


 834           continue



      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 5, a, G)
             do mm = 1, nvar
                call solve(9,5,G,rhs(:,mm))

                gradmyy(k,mm) =  2.d0*rhs(5,mm)/(dy**2)
                gradmxy(k,mm) =  rhs(4,mm)/(dx*dy)
                gradmxx(k,mm) =  2.d0*rhs(3,mm)/(dx**2)
                gradmy(k,mm)  =  rhs(2,mm)/dy
                gradmx(k,mm)  =  rhs(1,mm)/dx
            end do




 820    continue ! iterate over merging tiles on entire grid
        gradmyy =  0.d0
        gradmxy =  0.d0
        gradmxx =  0.d0

!       irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,
!     . lstgrd,
!     . vmerge, xcmerge, ycmerge, nhoods,numcount,
!     . mi, mj, nc,
!     . mreconi, mreconj,
!     . qMerge, qmshifts,
!     . nvar,
!     . qmx,qmy

!      call qmslopes(irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,
!     . lstgrd,
!     . volmerge, xcentmerge, ycentmerge, numhoods,
!     . iidx, jidx, ncount,
!     . mioff, mjoff,
!     . qMerge, qmshift,
!     . nvar,
!     . qmx,qmy)
!      print *, "here!\n"










        elseif (ssw .eq. -10) then






        do 320 j = lwidth+1, mjtot-lwidth
        do 320 i = lwidth+1, mitot-lwidth
             k = irr(i,j)
             if (k .eq. -1 .or. k .eq. lstgrd) go to 320
             ! solid/regular cell doesnt participate in anything
             ! can reconstruct high order poly here, but for now do linear reconstruction
             ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save

             rhs = 0.d0 ! initialize for accumulation
             a = 0.d0
             x0 = xcentMerge(k)
             y0 = ycentMerge(k)


            do 434 ioff = -mioff(i,j), mioff(i,j)
            do 434 joff = -mjoff(i,j), mjoff(i,j)
                 if (IS_GHOST(i+ioff,j+joff)) go to 434
                 koff = irr(i+ioff,j+joff)
                 if (koff .eq. -1) goto 434

                 if(koff .eq. lstgrd) then
                     xcoff = xlow + (i+ioff-0.5d0)*dx
                     ycoff = ylow + (j+joff-0.5d0)*dy
                     deltax = xcoff - x0
                     deltay = ycoff - y0
                 else
                     deltax = xcentMerge(koff) - x0
                     deltay = ycentMerge(koff) - y0
                 endif

                 f(1) =  deltax/dx
                 f(2) =  deltay/dy
                 f(3) = (deltax/dx)**2
                 f(4) = (deltax/dx)*(deltay/dy)
                 f(5) = (deltay/dy)**2

            do ii = 1,5
            do jj = 1,5
                a(ii,jj) = a(ii,jj) + f(ii)*f(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .               + f(ii) * (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
            enddo

 434           continue



      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 5, a, G)
             do mm = 1, nvar
                call solve(9,5,G,rhs(:,mm))

                gradmyy(k,mm) =  2.d0*rhs(5,mm)/(dy**2)
                gradmxy(k,mm) =  rhs(4,mm)/(dx*dy)
                gradmxx(k,mm) =  2.d0*rhs(3,mm)/(dx**2)
                gradmy(k,mm)  =  rhs(2,mm)/dy
                gradmx(k,mm)  =  rhs(1,mm)/dx
            end do




 320    continue ! iterate over merging tiles on entire grid
        gradmyy =  0.d0
        gradmxy =  0.d0
        gradmxx =  0.d0

        elseif(ssw .eq. 3) then









        do 620 j = lwidth+1, mjtot-lwidth
        do 620 i = lwidth+1, mitot-lwidth
             k = irr(i,j)
             if (k .eq. -1 .or. k .eq. lstgrd) go to 620
             ! solid/regular cell doesnt participate in anything
             ! can reconstruct high order poly here, but for now do linear reconstruction
             ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save

             rhs = 0.d0 ! initialize for accumulation
             a = 0.d0
             x0 = xcentMerge(k)
             y0 = ycentMerge(k)


            do 634 ioff = -mioff(i,j), mioff(i,j)
            do 634 joff = -mjoff(i,j), mjoff(i,j)
                 if (IS_GHOST(i+ioff,j+joff)) go to 634
                 koff = irr(i+ioff,j+joff)
                 if (koff .eq. -1) goto 634

                 if(koff .eq. lstgrd) then
                     xcoff = xlow + (i+ioff-0.5d0)*dx
                     ycoff = ylow + (j+joff-0.5d0)*dy
                     deltax = xcoff - x0
                     deltay = ycoff - y0
                 else
                     deltax = xcentMerge(koff) - x0
                     deltay = ycentMerge(koff) - y0
                 endif

                 f(1) = deltax/dx
                 f(2) = deltay/dy
                 f(3) = -dmergeshifts(1,k)
                 f(4) = -dmergeshifts(2,k)
                 f(5) = -dmergeshifts(3,k)
                 f(6) = -dmergeshifts(4,k)
                 f(7) = -dmergeshifts(5,k)
                 f(8) = -dmergeshifts(6,k)
                 f(9) = -dmergeshifts(7,k)

                 if(koff .eq. lstgrd) then ! set up data on whole cell
                    ncount(koff)  = 1
                    iidx(koff, 1) = i+ioff
                    jidx(koff, 1) = j+joff
                    volmerge(koff) = dx * dy / numhoods(i+ioff, j+joff)
                 endif


                 cvm = volmerge(koff) ! off vol merge
      do 697 ic = 1, ncount(koff) ! compute weighted inner product of monomials on this neighborhood
                    icurr = iidx(koff,ic)
                    jcurr = jidx(koff,ic)
                    kcurr = irr(icurr,jcurr)
                    nhc = numhoods(icurr, jcurr) ! num hoods current
                    if(kcurr .eq. lstgrd) then
                      itri = 2
                      poly(1,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(1,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(2,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(2,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(3,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(3,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy

                      poly(4,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(4,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy
                    else
                     ivert = 1
                      do 690 while (poly(ivert+1,1,kcurr) .ne. -11.)
                        ivert = ivert + 1
  690                   continue
                       itri = ivert - 3
                    endif

      ! computing the inner product on each triangle of neighborhood member
                  idx1 = 1
                  do 691 it = 1, itri ! for each  triangle
                    idx2 = it + 1
                    idx3 = it + 2

                    x1 = poly(idx1,1,kcurr)
                    y1 = poly(idx1,2,kcurr)

                    x2 = poly(idx2,1,kcurr)
                    y2 = poly(idx2,2,kcurr)

                    x3 = poly(idx3,1,kcurr)
                    y3 = poly(idx3,2,kcurr)

                    artri = area(x1, x2, x3, y1, y2, y3)

                    do 692 itq = 1,ntriquad
                        xval = x1 * rtri(itq) + x2 * stri(itq)
     .                     +  x3 * (1.d0-rtri(itq)-stri(itq))
                        yval = y1 * rtri(itq) + y2 * stri(itq)
     .                     +  y3 * (1.d0-rtri(itq)-stri(itq))

            f(3) = f(3) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)**2 )
     .                                                / (dx**2)

            f(4) = f(4) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)/dx )
     .                                                *( (yval-y0)/dy )

            f(5) = f(5) + (artri/nhc/cvm) * wtri(itq) *( (yval-y0)**2 )
     .                                                / (dy**2)





            f(6) = f(6) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)**3 )
     .                                                / (dx**3)

            f(7) = f(7) + (artri/nhc/cvm) * wtri(itq)
     .                   *( (xval-x0)**2 / dx**2)*( (yval-y0)/dy )

            f(8) = f(8) + (artri/nhc/cvm) * wtri(itq)
     .                   *( (xval-x0) / dx)*( (yval-y0)**2/dy**2 )

            f(9) = f(9) + (artri/nhc/cvm) * wtri(itq) *( (yval-y0)**3 )
     .                                                / (dy**3)




 692        continue ! for each quadrature point on each triangle
 691        continue ! for each triangle
 697        continue ! for each subelement of the current neighborhood
        ! integration is complete, let's accumulate

            do ii = 1,9
            do jj = 1,9
                a(ii,jj) = a(ii,jj) + f(ii)*f(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .               + f(ii) * (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
            enddo

 6334         continue


 634           continue



      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 9, a, G)
             do mm = 1, nvar
                call solve(9,9,G,rhs(:,mm))
                gradmyyy(k,mm) = 6.d0*rhs(9,mm)/(dy**3)
                gradmxyy(k,mm) = 2.d0*rhs(8,mm)/(dx*dy**2)
                gradmxxy(k,mm) = 2.d0*rhs(7,mm)/(dy*dx**2)
                gradmxxx(k,mm) = 6.d0*rhs(6,mm)/(dx**3)
                gradmyy(k,mm) =  2.d0*rhs(5,mm)/(dy**2)
                gradmxy(k,mm) =  rhs(4,mm)/(dx*dy)
                gradmxx(k,mm) =  2.d0*rhs(3,mm)/(dx**2)
                gradmy(k,mm)  =  rhs(2,mm)/dy
                gradmx(k,mm)  =  rhs(1,mm)/dx
            end do


 620    continue ! iterate over merging tiles on entire grid












        endif


















      if(nolimiter .eqv. .false.) then
        call limiter_srd(qMerge,gradmx,gradmy,mitot,mjtot,irr, nvar,
     .                   dx, dy,lwidth, lstgrd, xlow, ylow)
      endif






















      valnew = 0.d0  !  all cells initialized to 0

      if(ssw .eq. 0) then
      do 55 j = lwidth+1, mjtot-lwidth
      do 55 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1) then
             valnew(i,j,:) = q(i,j,:)
             go to 55  ! does not contribute
          endif

          if(k .eq. lstgrd) then
           valnew(i,j,:) = valnew(i,j,:) + qMerge(i,j,:) / numHoods(i,j)
          else
          do 44 ic = 1, ncount(k)
                 ioff = iidx(k,ic)
                 joff = jidx(k,ic)
                 koff = irr(ioff,joff)
             call getCellCentroid(lstgrd,ioff,joff,xc,yc,xlow,
     .                            ylow,dx,dy,koff)
             qm(:) = qMerge(i,j,:)
             valnew(ioff,joff,:) = valnew(ioff,joff,:) +
     .              qm(:)/numHoods(ioff,joff)

 44       continue
        endif
 55    continue


      elseif(ssw .eq. 1) then


      do 50 j = lwidth+1, mjtot-lwidth
      do 50 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1) then
             valnew(i,j,:) = q(i,j,:)
             go to 50  ! does not contribute
          endif

          if(k .eq. lstgrd) then
           valnew(i,j,:) = valnew(i,j,:) + qMerge(i,j,:) / numHoods(i,j)
          else
          do 40 ic = 1, ncount(k)
                 ioff = iidx(k,ic)
                 joff = jidx(k,ic)
                 koff = irr(ioff,joff)
             call getCellCentroid(lstgrd,ioff,joff,xc,yc,xlow,
     .                            ylow,dx,dy,koff)
             qm(:) = qMerge(i,j,:) + (xc-xcentMerge(k))*gradmx(k,:)
     .                             + (yc-ycentMerge(k))*gradmy(k,:)
             valnew(ioff,joff,:) = valnew(ioff,joff,:) +
     .              qm(:)/numHoods(ioff,joff)

 40       continue
        endif
 50    continue

      elseif(ssw .eq. -10) then


      do 150 j = lwidth+1, mjtot-lwidth
      do 150 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1) then
             valnew(i,j,:) = q(i,j,:)
             go to 150  ! does not contribute
          endif

          if(k .eq. lstgrd) then
           valnew(i,j,:) = valnew(i,j,:) + qMerge(i,j,:) / numHoods(i,j)
          else
          do 140 ic = 1, ncount(k)
                 ioff = iidx(k,ic)
                 joff = jidx(k,ic)
                 koff = irr(ioff,joff)
             call getCellCentroid(lstgrd,ioff,joff,xc,yc,xlow,
     .                            ylow,dx,dy,koff)
             qm(:) = qMerge(i,j,:) + (xc-xcentMerge(k))*gradmx(k,:)
     .                             + (yc-ycentMerge(k))*gradmy(k,:)
             valnew(ioff,joff,:) = valnew(ioff,joff,:) +
     .              qm(:)/numHoods(ioff,joff)

 140       continue
        endif
 150    continue



      elseif(ssw .eq. -1) then
      do 650 j = lwidth+1, mjtot-lwidth
      do 650 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1) then
             valnew(i,j,:) = q(i,j,:)
             go to 650  ! does not contribute
          endif

          if(k .eq. lstgrd) then
           valnew(i,j,:) = valnew(i,j,:) + qMerge(i,j,:) / numHoods(i,j)
          else
          do 640 ic = 1, ncount(k)
                 ioff = iidx(k,ic)
                 joff = jidx(k,ic)
                 koff = irr(ioff,joff)
             call getCellCentroid(lstgrd,ioff,joff,xc,yc,xlow,
     .                            ylow,dx,dy,koff)
             qm(:) = qMerge(i,j,:) + (xc-xcentMerge(k))*gradmx(k,:)
     .                             + (yc-ycentMerge(k))*gradmy(k,:)
             valnew(ioff,joff,:) = valnew(ioff,joff,:) +
     .              qm(:)/numHoods(ioff,joff)

 640       continue
        endif
 650    continue
!      do 516 j = lwidth+1, mjtot-lwidth
!      do 516 i = lwidth+1, mitot-lwidth
!          k = irr(i,j)
!          if (k .eq. -1) then
!             valnew(i,j,:) = q(i,j,:)
!             go to 516  ! does not contribute
!          endif
!          if(k .eq. lstgrd) then
!           valnew(i,j,:) = valnew(i,j,:) + qMerge(i,j,:) / numHoods(i,j)
!           cycle
!          endif
!
!      ! if here, I am on a cut cell
!
!      do 416 ic = 1, ncount(k)
!                 ioff = iidx(k,ic)
!                 joff = jidx(k,ic)
!                 koff = irr(ioff,joff)
!                 call getCellCentroid(lstgrd,ioff,joff,xc,yc,xlow,
!     .                            ylow,dx,dy,koff)
!      qm = 0.d0
!      if(koff .eq. lstgrd) then
!          arr = dx*dy
!          itri = 2
!          poly(1,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
!          poly(1,2,koff) = ylow + (dfloat(joff)-1.d0)*dy
!
!          poly(2,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
!          poly(2,2,koff) = ylow + (dfloat(joff)-1.d0)*dy
!
!          poly(3,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
!          poly(3,2,koff) = ylow + (dfloat(joff)-0.d0)*dy
!
!          poly(4,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
!          poly(4,2,koff) = ylow + (dfloat(joff)-0.d0)*dy
!      else
!          arr = ar(koff)
!          ivert = 1
!          do 246 while (poly(ivert+1,1,koff) .ne. -11.)
!            ivert = ivert + 1
!  246      continue
!           itri = ivert - 3
!      endif
!
!      idx1 = 1
!      do 26 it = 1, itri ! for each  triangle
!            idx2 = it + 1
!            idx3 = it + 2
!
!            x1 = poly(idx1,1,koff)
!            y1 = poly(idx1,2,koff)
!
!            x2 = poly(idx2,1,koff)
!            y2 = poly(idx2,2,koff)
!
!            x3 = poly(idx3,1,koff)
!            y3 = poly(idx3,2,koff)
!
!            artri = area(x1, x2, x3, y1, y2, y3)
!
!            do 226 itq = 1,ntriquad
!                xval = x1 * rtri(itq) + x2 * stri(itq)
!     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
!                yval = y1 * rtri(itq) + y2 * stri(itq)
!     .              +  y3 * (1.d0-rtri(itq)-stri(itq))
!
!                diffx = xval-xcentMerge(k)
!                diffy = yval-ycentMerge(k)
!
!            qm(:) = qm(:) + (artri/arr) * wtri(itq) *
!     .  (qMerge(i,j,:)+ gradmx(k,:) * diffx + gradmy(k,:) * diffy
!     . + 0.5d0*gradmxx(k,:)*(diffx*diffx - (dx**2)*dmergeshifts(1,k))
!     . +       gradmxy(k,:)*(diffx*diffy - dx*dy*dmergeshifts(2,k))
!     . + 0.5d0*gradmyy(k,:)*(diffy*diffy - (dy**2)*dmergeshifts(3,k))  )
!  226   continue ! for each quadrature point on each triangle
!  26   continue ! for each triangle
!
!        valnew(ioff,joff,:) = valnew(ioff,joff,:) +
!     .              qm(:)/numHoods(ioff,joff)
!
! 416       continue ! cycle through all the subcells
!
!
!
! 516    continue ! onto the next cell in the grid






      elseif(ssw .eq. 2) then


      do 510 j = lwidth+1, mjtot-lwidth
      do 510 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1) then
             valnew(i,j,:) = q(i,j,:)
             go to 510  ! does not contribute
          endif
          if(k .eq. lstgrd) then
           valnew(i,j,:) = valnew(i,j,:) + qMerge(i,j,:) / numHoods(i,j)
           cycle
          endif

      ! if here, I am on a cut cell

      do 410 ic = 1, ncount(k)
                 ioff = iidx(k,ic)
                 joff = jidx(k,ic)
                 koff = irr(ioff,joff)
                 call getCellCentroid(lstgrd,ioff,joff,xc,yc,xlow,
     .                            ylow,dx,dy,koff)
      qm = 0.d0
      if(koff .eq. lstgrd) then
          arr = dx*dy
          itri = 2
          poly(1,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
          poly(1,2,koff) = ylow + (dfloat(joff)-1.d0)*dy

          poly(2,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
          poly(2,2,koff) = ylow + (dfloat(joff)-1.d0)*dy

          poly(3,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
          poly(3,2,koff) = ylow + (dfloat(joff)-0.d0)*dy

          poly(4,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
          poly(4,2,koff) = ylow + (dfloat(joff)-0.d0)*dy
      else
          arr = ar(koff)
          ivert = 1
          do 240 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
  240      continue
           itri = ivert - 3
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

            do 225 itq = 1,ntriquad
                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

                diffx = xval-xcentMerge(k)
                diffy = yval-ycentMerge(k)

            qm(:) = qm(:) + (artri/arr) * wtri(itq) *
     .  (qMerge(i,j,:)+ gradmx(k,:) * diffx + gradmy(k,:) * diffy
     . + 0.5d0*gradmxx(k,:)*(diffx*diffx - (dx**2)*dmergeshifts(1,k))
     . +       gradmxy(k,:)*(diffx*diffy - dx*dy*dmergeshifts(2,k))
     . + 0.5d0*gradmyy(k,:)*(diffy*diffy - (dy**2)*dmergeshifts(3,k))  )
  225   continue ! for each quadrature point on each triangle
  21   continue ! for each triangle

        valnew(ioff,joff,:) = valnew(ioff,joff,:) +
     .              qm(:)/numHoods(ioff,joff)

 410       continue ! cycle through all the subcells



 510    continue ! onto the next cell in the grid










      elseif(ssw .eq. 3) then


      do 1510 j = lwidth+1, mjtot-lwidth
      do 1510 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1) then
             valnew(i,j,:) = q(i,j,:)
             go to 1510  ! does not contribute
          endif
          if(k .eq. lstgrd) then
           valnew(i,j,:) = valnew(i,j,:) + qMerge(i,j,:) / numHoods(i,j)
           cycle
          endif

      ! if here, I am on a cut cell

      do 1410 ic = 1, ncount(k)
                 ioff = iidx(k,ic)
                 joff = jidx(k,ic)
                 koff = irr(ioff,joff)
                 call getCellCentroid(lstgrd,ioff,joff,xc,yc,xlow,
     .                            ylow,dx,dy,koff)
      qm = 0.d0
      if(koff .eq. lstgrd) then
          arr = dx*dy
          itri = 2
          poly(1,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
          poly(1,2,koff) = ylow + (dfloat(joff)-1.d0)*dy

          poly(2,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
          poly(2,2,koff) = ylow + (dfloat(joff)-1.d0)*dy

          poly(3,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
          poly(3,2,koff) = ylow + (dfloat(joff)-0.d0)*dy

          poly(4,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
          poly(4,2,koff) = ylow + (dfloat(joff)-0.d0)*dy
      else
          arr = ar(koff)
          ivert = 1
          do 1240 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
1240      continue
           itri = ivert - 3
      endif

      idx1 = 1
      do 121 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 1225 itq = 1,ntriquad
                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

                diffx = xval-xcentMerge(k)
                diffy = yval-ycentMerge(k)

            qm(:) = qm(:) + (artri/arr) * wtri(itq) *
     .  (qMerge(i,j,:)+ gradmx(k,:) * diffx + gradmy(k,:) * diffy
     . + 0.5d0*gradmxx(k,:)*(diffx*diffx - (dx**2)*dmergeshifts(1,k))
     . +       gradmxy(k,:)*(diffx*diffy - dx*dy*dmergeshifts(2,k))
     . + 0.5d0*gradmyy(k,:)*(diffy*diffy - (dy**2)*dmergeshifts(3,k))
     .+(1.d0/6.d0)*gradmxxx(k,:)*(diffx**3 - (dx**3)*dmergeshifts(4,k))
     .+ 0.5d0*gradmxxy(k,:)*(diffy*diffx**2-dx**2*dy*dmergeshifts(5,k))
     .+ 0.5d0*gradmxyy(k,:)*(diffy**2*diffx-dx*dy**2*dmergeshifts(6,k))
     .+(1.d0/6.d0)*gradmyyy(k,:)*(diffy**3 - (dy**3)*dmergeshifts(7,k))
     .)
 1225   continue ! for each quadrature point on each triangle
 121   continue ! for each triangle

        valnew(ioff,joff,:) = valnew(ioff,joff,:) +
     .              qm(:)/numHoods(ioff,joff)

!        recon = qm(:)
!        velsq = recon(2)**2+recon(3)**2
!        press = .4d0*(recon(4)-0.5d0*velsq/recon(1))
!        if( press .lt. 0.d0) then
!        print *, "uh oh"
!        endif

 1410       continue ! cycle through all the subcells



 1510    continue ! onto the next cell in the grid






      endif

      valnew(1:lwidth,:,:) =  q(1:lwidth,:,:)
      valnew(mitot-lwidth+1:mitot,:,:) =  q(mitot-lwidth+1:mitot,:,:)
      valnew(:,1:lwidth,:) =  q(:,1:lwidth,:)
      valnew(:,mjtot-lwidth+1:mjtot,:) =  q(:,mjtot-lwidth+1:mjtot,:)

      q = valnew



       if (verbose) then
          totmassafter =  bigconck(q,irr,mitot,mjtot,lwidth,nvar)
          write(*,921) totmassafter
 921      format(/,"         SRDmass after redistribution is  ",e30.20)
          print *, "diff is ",totmassafter-totmass
       endif







      end
!
!
!
!
!
!c
!c -------------------------------------------------------------------------
!c
      subroutine getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)

      implicit double precision (a-h,o-z)
      include "cirr.i"

      if (k .eq. lstgrd) then
         xc = xlow + (i-0.5d0)*dx
         yc = ylow + (j-0.5d0)*dy
      else if (k .eq. -1) then
         xc = 0.d0
         yc = 0.d0
      else
         if(ihob .eq. 0) then
             xc = xcirr(k)
             yc = ycirr(k)
         else
             xc = xcirr_ho(k)
             yc = ycirr_ho(k)
         endif
      endif

      return
      end
c
c ----------------------------------------------------------------------------
c
      function bigconck(q,irr,mitot,mjtot,lwidth,nvar)
c
      implicit double precision (a-h,o-z)
      include "cirr.i"

      dimension q(mitot,mjtot,nvar), irr(mitot,mjtot)

      totmass = 0.d0
      ar(-1)  = 0.d0


      do 10 iv = 1, nvar
      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth

         k = irr(i,j)
         totmass = totmass + q(i,j,iv)*ar(k)

 10   continue

      bigconck = totmass

      return
      end




















































c     MINIMUM DENSITY
c ----------------------------------------------------------------------------
c
c      subroutine determineDirection(q,areaTol,irr,mitot,mjtot,lwidth,
c     .nvar,i,j,idir)
c      implicit double precision (a-h, o-z)
c      include "cirr.i"
c      dimension irr(mitot,mjtot), q(mitot,mjtot,nvar)
c      dimension area(4), inum(4), avgrho(4)
c
c        k = irr(i,j)
c        ! compute volume in the four directions
c        ! ignore volumes that are < tolerance
c        ! choose direction that has the minimum density
c
c        do 30 id = 1,4
c            inum(id) = 0
c            icurr = i
c            jcurr = j
c            kcurr = irr(icurr,jcurr)
c            area(id) = ar(kcurr)
c            do while (area(id) < areaTol)
c                if(id .eq. 1) then
c                    icurr = icurr - 1
c                elseif(id .eq. 2) then
c                    icurr = icurr + 1
c                elseif(id .eq. 3) then
c                    jcurr = jcurr - 1
c                elseif(id .eq. 4) then
c                    jcurr = jcurr + 1
c                endif
c                kcurr = irr(icurr,jcurr)
c
c                if(kcurr .eq. -1) then
c                    exit
c                endif
c
c                area(id) = area(id) + ar(kcurr)
c                inum(id) = inum(id) + 1
c            end do
c  30    continue
c
c
c
c        ! find the direction that has the minimum density such that
c        ! area(id) > areaTOL.
c
c        ifirstimethru = 1
c        minrho = -1.d0
c        minidx = -1
c        do 50 id = 1,4
c            if(area(id) > areaTOL) then
c            icurr = i
c            jcurr = j
c            kcurr = irr(icurr,jcurr)
c            avgrho(id) = q(icurr,jcurr,1)*ar(kcurr)
c            do 60 ic = 1,inum(id)
c                if(id .eq. 1) then
c                    icurr = icurr - 1
c                elseif(id .eq. 2) then
c                    icurr = icurr + 1
c                elseif(id .eq. 3) then
c                    jcurr = jcurr - 1
c                elseif(id .eq. 4) then
c                    jcurr = jcurr + 1
c                endif
c                kcurr = irr(icurr,jcurr)
c                avgrho(id) = avgrho(id) + q(icurr,jcurr,1)*ar(kcurr)
c  60        continue
c            avgrho(id) = avgrho(id) / area(id)
c
c            if(ifirstimethru .eq. 1) then
c                minrho = avgrho(id)
c                minidx = id
c                ifirstimethru = 0
c            else
c                if(minrho > avgrho(id)) then
c                    minrho = avgrho(id)
c                    minidx = id
c                endif
c            endif
c
c            endif
c  50    continue
c
c        if(minidx .eq. -1) print *, "PROBLEM with direction."
c        idir = minidx
c
c
c      return
c      end




c       prioritize minimum shear
c
c ----------------------------------------------------------------------------
c
c     subroutine determineDirection(q,areaTol,irr,mitot,mjtot,lwidth,
c    .nvar,i,j,idir)
c     implicit double precision (a-h, o-z)
c     include "cirr.i"
c     dimension irr(mitot,mjtot), q(mitot,mjtot,nvar)
c     dimension area(4), inum(4), avgq(4,nvar), shear(4)
c     logical IS_GHOST
c     IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
c    .                 j .le. lwidth .or. j .gt. mjtot-lwidth)
c
c       k = irr(i,j)
c       ! compute volume in the four directions
c       ! ignore volumes that are < tolerance
c       ! choose direction that has the minimum density
c
c       do 30 id = 1,4
c           inum(id) = 0
c           icurr = i
c           jcurr = j
c           kcurr = irr(icurr,jcurr)
c           area(id) = ar(kcurr)
c           do while (area(id) < areaTol)
c               if(id .eq. 1) then
c                   icurr = icurr - 1
c               elseif(id .eq. 2) then
c                   icurr = icurr + 1
c               elseif(id .eq. 3) then
c                   jcurr = jcurr - 1
c               elseif(id .eq. 4) then
c                   jcurr = jcurr + 1
c               endif
c               kcurr = irr(icurr,jcurr)
c
c               if(kcurr .eq. -1) then
c                   exit
c               endif
c
c               area(id) = area(id) + ar(kcurr)
c               inum(id) = inum(id) + 1
c           end do
c 30    continue
c
c
c
c       ! find the direction that has the minimum density such that
c       ! area(id) > areaTOL.
c
c       ifirstimethru = 1
c       minshear = -1.d0
c       minidx = -1
c       shear(:) = -1.d0
c       do 50 id = 1,4
c           if(area(id) > areaTOL) then
c           icurr = i
c           jcurr = j
c           kcurr = irr(icurr,jcurr)
c           avgq(id,:) = q(icurr,jcurr,:)*ar(kcurr)
c           do 60 ic = 1,inum(id)
c               if(id .eq. 1) then
c                   icurr = icurr - 1
c               elseif(id .eq. 2) then
c                   icurr = icurr + 1
c               elseif(id .eq. 3) then
c                   jcurr = jcurr - 1
c               elseif(id .eq. 4) then
c                   jcurr = jcurr + 1
c               endif
c               kcurr = irr(icurr,jcurr)
c               avgq(id,:) = avgq(id,:) + q(icurr,jcurr,:)*ar(kcurr)
c 60        continue
c           avgq(id,:) = avgq(id,:) / area(id)
c
c
c
c           icurr  = i
c           icurr1 = i
c           icurr2 = i
c           icurr3 = i
c           jcurr  = j
c           jcurr1 = j
c           jcurr2 = j
c           jcurr3 = j
c
c
c
c       if(id .eq. 1) then
c           jcurr1 = jcurr-1
c           jcurr2 = jcurr+1
c
c           icurr3 = icurr+1
c       elseif(id .eq. 2) then
c           jcurr1 = jcurr-1
c           jcurr2 = jcurr+1
c
c           icurr3 = icurr-1
c       elseif(id .eq. 3) then
c           icurr1 = icurr-1
c           icurr2 = icurr+1
c
c           jcurr3 = jcurr+1
c       elseif(id .eq. 4) then
c           icurr1 = icurr-1
c           icurr2 = icurr+1
c
c           jcurr3 = jcurr-1
c       endif
c
c           kcurr  = irr(icurr,jcurr)
c           kcurr1 = irr(icurr1,jcurr1)
c           kcurr2 = irr(icurr2,jcurr2)
c           kcurr3 = irr(icurr3,jcurr3)
c
c           shear1 = -1.d0
c           shear2 = -1.d0
c           shear3 = -1.d0
c
c           if(id .eq. 1 .or. id .eq. 2) then
c               speed  = avgq(id,2)/avgq(id,1)
c               if(kcurr1 > -1 .and. ar(kcurr1) > areaTOL) then
c                   speed1 = q(icurr1,jcurr1,2)/q(icurr1,jcurr1,1)
c                   shear1 = abs(speed1-speed)
c               endif
c               if(kcurr2 > -1 .and. ar(kcurr2) > areaTOL) then
c                   speed2 = q(icurr2,jcurr2,2)/q(icurr2,jcurr2,1)
c                   shear2 = abs(speed2-speed)
c               endif
c           elseif(id .eq. 3 .or. id .eq. 4) then
c               speed  = avgq(id,3)/avgq(id,1)
c               if(kcurr1 > -1 .and. ar(kcurr1) > areaTOL) then
c                   speed1 = q(icurr1,jcurr1,3)/q(icurr1,jcurr1,1)
c                   shear1 = abs(speed1-speed)
c               endif
c               if(kcurr2 > -1 .and. ar(kcurr2) > areaTOL) then
c                   speed2 = q(icurr2,jcurr2,3)/q(icurr2,jcurr2,1)
c                   shear2 = abs(speed2-speed)
c               endif
c           endif
c           if(id .eq. 1 .or. id .eq. 2) then
c               speed  = avgq(id,3)/avgq(id,1)
c               if(kcurr3 > -1 .and. ar(kcurr3) > areaTOL) then
c                   speed3 = q(kcurr3,kcurr3,3)/q(kcurr3,jcurr3,1)
c                   shear3 = abs(speed3-speed)
c               endif
c           elseif(id .eq. 3 .or. id .eq. 4) then
c               speed  = avgq(id,2)/avgq(id,1)
c               if(kcurr3 > -1 .and. ar(kcurr3) > areaTOL) then
c                   speed3 = q(kcurr3,kcurr3,2)/q(kcurr3,jcurr3,1)
c                   shear3 = abs(speed3-speed)
c               endif
c           endif
c
c
c
c           shear(id) = max(shear1, shear2, shear3)
c
c
c           icurr  = i
c           icurr1 = i
c           icurr2 = i
c           icurr3 = i
c           jcurr  = j
c           jcurr1 = j
c           jcurr2 = j
c           jcurr3 = j
c
c
c
c       ! compute max shear using this direction's stencil
c           do 70 ic = 1,inum(id)
c               icurrprev = icurr
c               jcurrprev = jcurr
c               if(id .eq. 1) then
c                   icurr  = icurrprev-1
c                   icurr1 = icurrprev-1
c                   icurr2 = icurrprev-1
c                   icurr3 = icurrprev-2
c
c                   jcurr1 = jcurrprev-1
c                   jcurr2 = jcurrprev+1
c               elseif(id .eq. 2) then
c                   icurr  = icurrprev+1
c                   icurr1 = icurrprev+1
c                   icurr2 = icurrprev+1
c                   icurr3 = icurrprev+2
c
c                   jcurr1 = jcurrprev-1
c                   jcurr2 = jcurrprev+1
c               elseif(id .eq. 3) then
c                   icurr1 = icurrprev+1
c                   icurr2 = icurrprev-1
c
c                   jcurr1 = jcurrprev-1
c                   jcurr2 = jcurrprev-1
c                   jcurr3 = jcurrprev-2
c                   jcurr =  jcurrprev-1
c               elseif(id .eq. 4) then
c                   icurr1 = icurrprev+1
c                   icurr2 = icurrprev-1
c
c                   jcurr1 = jcurrprev+1
c                   jcurr2 = jcurrprev+1
c                   jcurr3 = jcurrprev+2
c                   jcurr  = jcurrprev+1
c               endif
c               kcurr  = irr(icurr,jcurr)
c               kcurr1 = irr(icurr1,jcurr1)
c               kcurr2 = irr(icurr2,jcurr2)
c               kcurr3 = irr(icurr3,jcurr3)
c
c               shear1 = -1.d0
c               shear2 = -1.d0
c               shear3 = -1.d0
c
c               if(id .eq. 1 .or. id .eq. 2) then
c                   speed  = avgq(id,2)/avgq(id,1)
c                   if(kcurr1 > -1 .and. ar(kcurr1) > areaTOL) then
c                       speed1 = q(icurr1,jcurr1,2)/q(icurr1,jcurr1,1)
c                       shear1 = abs(speed1-speed)
c                   endif
c                   if(kcurr2 > -1 .and. ar(kcurr2) > areaTOL) then
c                       speed2 = q(icurr2,jcurr2,2)/q(icurr2,jcurr2,1)
c                       shear2 = abs(speed2-speed)
c                   endif
c               elseif(id .eq. 3 .or. id .eq. 4) then
c                   speed  = avgq(id,3)/avgq(id,1)
c                   if(kcurr1 > -1 .and. ar(kcurr1) > areaTOL) then
c                       speed1 = q(icurr1,jcurr1,3)/q(icurr1,jcurr1,1)
c                       shear1 = abs(speed1-speed)
c                   endif
c                   if(kcurr2 > -1 .and. ar(kcurr2) > areaTOL) then
c                       speed2 = q(icurr2,jcurr2,3)/q(icurr2,jcurr2,1)
c                       shear2 = abs(speed2-speed)
c                   endif
c               endif
c
c           if(ic .eq. inum(id)) then
c               if(id .eq. 1 .or. id .eq. 2) then
c                   speed  = avgq(id,3)/avgq(id,1)
c                   if(kcurr3 > -1 .and. ar(kcurr3) > areaTOL) then
c                       speed3 = q(icurr3,jcurr3,3)/q(icurr3,jcurr3,1)
c                       shear3 = abs(speed3-speed)
c                   endif
c               elseif(id .eq. 3 .or. id .eq. 4) then
c                   speed  = avgq(id,2)/avgq(id,1)
c                   if(kcurr3 > -1 .and. ar(kcurr3) > areaTOL) then
c                       speed3 = q(icurr3,jcurr3,2)/q(icurr3,jcurr3,1)
c                       shear3 = abs(speed3-speed)
c                   endif
c               endif
c           endif
c
c               shear(id) = max(shear(id),shear1, shear2, shear3)
c70        continue
c
c          if(ifirstimethru .eq. 1) then
c              minshear = shear(id)
c              minidx = id
c              ifirstimethru = 0
c          else
c              if(minshear > shear(id)) then
c                  minshear = shear(id)
c                  minidx = id
c              endif
c          endif
c
c       endif
c 50    continue
c
c       if(minidx .eq. -1) print *, "PROBLEM with direction."
c       idir = minidx
c
c
c     return
c     end
c
c        do 70 ioff = -nco, nco
c        do 70 joff = -nco, nco
c              if (ioff .eq. 0 .and. joff .eq. 0) go to 31
c              if (IS_GHOST(i+ioff,j+joff)) go to 31
c              koff = irr(i+ioff,j+joff)
c              if (koff .eq. -1) go to 31
c              if (ar(koff) < areaTOL) go to 31
c
c
c
c 70     continue

