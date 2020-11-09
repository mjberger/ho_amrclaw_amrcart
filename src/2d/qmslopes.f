      subroutine qmslopes(irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,
     &                     lstgrd, numHoods,  mioff, mjoff,
     &                     qMerge, nvar, qmx,qmy,qmxx,qmxy,qmyy)

! !   INPUTS:
! !   REQUIRED GRID DATA
! !
! !   irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,lstgrd
! !
! !   REQUIRED MERGING DATA
! !   volMerge, xcentMerge, ycentMerge, numHoods, ncount
! !
! !   mioff, mjoff: RECONSTRUCTION NEIGHBORHOODS ON MERGING NEIGHS
! !
! !   qmshifts, qMerge
! !
! !   nvar: number of equations
! !
! !   OUTPUT
! !   qmx,qmy: second order merging slopes


      use amr_module
      implicit double precision (a-h, o-z)
      include "quadrature.i"
      include "cuserdt.i"

      dimension irr(mitot,mjtot), numHoods(mitot,mjtot)

      dimension qMerge(nvar,mitot,mjtot)
      dimension qmx(nvar,irrsize), qmy(nvar,irrsize)
      dimension qmxx(nvar,irrsize), qmyy(nvar,irrsize)
      dimension qmxy(nvar,irrsize)

      dimension a(9,9), rhs(9,nvar), b(9), f(9), w(9,nvar),G(9,9)

      dimension mioff(mitot,mjtot), mjoff(mitot,mjtot)

      logical OUT_OF_RANGE, NOT_TRUSTED
      logical  debug


       OUT_OF_RANGE(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     .                      j .lt. 1 .or. j .gt. mjtot)

       NOT_TRUSTED(i,j) = (i .eq. 1 .or. i .eq. mitot-1 .or.
     .                     j .eq. 1 .or. j .eq. mjtot-1)

      areaMin = areaFrac*dx*dy
      debug = .true.

      ist = 2
      jst = 2
      iend = mitot-1
      jend = mjtot-1

        do 821 j = jst, jend
        do 820 i = ist, iend
             k = irr(i,j)
             if (k .eq. -1 .or. k .eq. lstgrd) cycle
             if (ar(k) .gt. areaMin) cycle

             rhs = 0.d0 ! initialize for accumulation
             a = 0.d0
             x0 = xcentMerge(k)
             y0 = ycentMerge(k)


            do 834 joff = -mjoff(i,j), mjoff(i,j)
            do 834 ioff = -mioff(i,j), mioff(i,j)
                 if (ioff .eq. 0 .and. joff .eq. 0) go to 834
                 if (OUT_OF_RANGE(i+ioff,j+joff)) go to 834
                 ! cant use 1st and last cell in nhood since no
                 ! good update from method
                 if (NOT_TRUSTED(i+ioff,j+joff)) go to 834
                 koff = irr(i+ioff,j+joff)
                 if (koff .eq. -1) goto 834

                 if (koff .eq. lstgrd) then
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

                 f(3) = -qmshifts(1,k) !AG calls them dmergeshifts
                 f(4) = -qmshifts(2,k)
                 f(5) = -qmshifts(3,k)

                 if (koff .eq. lstgrd) then ! set up data on whole cell
                    cvm = dx * dy / numHoods(i+ioff, j+joff)
                 else
                    cvm = volMerge(koff) ! off vol merge
                 endif

                 if (koff.eq. lstgrd .and. ncount(koff) .ne. 0) then
                   write(*,*)"found problem"
                   stop
                 endif
                 ! compute weighted inner product of monomials on this neighborhood
                 do 897 ic = 0, ncount(koff) 
                    if (ic .eq. 0) then !nhood is itself
                       icurr = i+ioff
                       jcurr = j+joff
                       kcurr = irr(icurr,jcurr)
                    else ! nbor had its own nhood
                       icurr = iidx(ic,koff)
                       jcurr = jidx(ic,koff)
                       kcurr = irr(icurr,jcurr)
                    endif

                    nhc = numHoods(icurr, jcurr) ! num hoods current
                    if (kcurr .eq. lstgrd) then
                      itri = 2
                      call makep(poly(1,1,lstgrd),icurr,jcurr,xlow,ylow,
     &                           dx,dy)
                    else
                      ivert = 1
                      do 890 while (poly(ivert+1,1,kcurr) .ne. -11.)
                        ivert = ivert + 1
  890                   continue
                      itri = ivert - 3
                    endif

      ! computing the inner product on each triangle of neighborhood member
                  indx1 = 1
                  do 891 it = 1, itri ! for each  triangle
                    indx2 = it + 1
                    indx3 = it + 2

                    x1 = poly(indx1,1,kcurr)
                    y1 = poly(indx1,2,kcurr)

                    x2 = poly(indx2,1,kcurr)
                    y2 = poly(indx2,2,kcurr)

                    x3 = poly(indx3,1,kcurr)
                    y3 = poly(indx3,2,kcurr)

                    artri = triangle_area(x1, x2, x3, y1, y2, y3)

                    do 892 itq = 1,ntriquad
                        xval = x1 * rtri(itq) + x2 * stri(itq)
     .                     +  x3 * (1.d0-rtri(itq)-stri(itq))
                        yval = y1 * rtri(itq) + y2 * stri(itq)
     .                     +  y3 * (1.d0-rtri(itq)-stri(itq))

                        f(3) = f(3) + (artri/nhc/cvm) * wtri(itq) *
     .                          ( (xval-x0)**2 ) / (dx**2)

                        f(4) = f(4) + (artri/nhc/cvm) * wtri(itq) *
     .                          ( (xval-x0)/dx )*( (yval-y0)/dy )

                        f(5) = f(5) + (artri/nhc/cvm) * wtri(itq) *
     .                          ( (yval-y0)**2 ) / (dy**2)
 892        continue ! for each quadrature point on each triangle
 891        continue ! for each triangle
 897        continue ! for each subelement of the current neighborhood
        ! integration is complete, let's accumulate

            do ii = 1, 5
              do jj = 1, 5
                a(ii,jj) = a(ii,jj) + f(ii)*f(jj)
              enddo

              rhs(ii,:) = rhs(ii,:)
     .               + f(ii) * (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
            enddo

 8334       continue


 834        continue



      ! SOLVE THE NORMAL EQUATIONS USING THE CHOLESKY FACTORIZATION.
            call cholesky(9, 5, a, G)
            do mm = 1, nvar
               call solve(9,5,G,rhs(:,mm))
               qmyy(mm,k)  =  2.d0*rhs(5,mm)/dy**2  ! not sure about 0.5
               qmxy(mm,k)  =  rhs(4,mm)/(dx*dy)
               qmxx(mm,k)  =  2.d0*rhs(3,mm)/dx**2  ! not sure about 0.5
               qmy(mm,k)  =  rhs(2,mm)/dy
               qmx(mm,k)  =  rhs(1,mm)/dx
            end do

 820    continue 
 821    continue ! iterate over merging tiles on entire grid

      return
      end
