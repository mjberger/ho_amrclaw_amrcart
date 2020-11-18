c
c ---------------------------------------------------------------------
c
       subroutine qslopes3(qp,qx,qy,qxx,qxy,qyy,mitot,mjtot,
     &                  irr,lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar,
     &                  iir,jjr,istage) 

      use amr_module
      implicit double precision(a-h,o-z)
      include "quadrature.i"

      dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &          qy(nvar,mitot,mjtot),irr(mitot,mjtot)
      dimension qxx(nvar,mitot,mjtot),qyy(nvar,mitot,mjtot),
     &          qxy(nvar,mitot,mjtot)
      dimension qxxx(nvar,mitot,mjtot),qxxy(nvar,mitot,mjtot),
     &          qxyy(nvar,mitot,mjtot),qyyy(nvar,mitot,mjtot)
      dimension iir(mitot,mjtot), jjr(mitot,mjtot)
      logical nolimiter
      common /order2/ ssw, quad, nolimiter
      include "cuserdt.i"


      dimension ata(9,9)
      dimension ag(9), rhs(nvar,9), G(9,9)
      logical   prflag, quad, enufNbor
      logical all_nbors_exist
      logical IS_REAL, IS_GHOST
      data      prflag/.false./

      IS_REAL(i,j) = (i .gt. 0 .and. i .le. mitot .and.
     &                j .gt. 0 .and. j .le. mjtot)

      IS_GHOST(i,j) = (i .le. lwidth .or. i .le. mitot-lwidth .or.
     &                 j .le. lwidth .or. j .le. mjtot-lwidth)


c   ##########
c   #  compute slopes for cut cells using least squares approach
c   #  code for igradChoice  = 5 - fit higher order accurate quadratic using iir,jjr  
c   ##########

c
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   qp contains primitive variables, regular cell slopes already set 
c   all slopes initialized to zero in regular slopes subr.
c
c
c     first save everything into unlim (even though limited already
c     on flow cells. will add unlimited cut gradients

      nTerms = 5 ! for quadratic through cell avg
      quad = .false.

      do 110 iy0 = 1, mjtot
      do 110 ix0 = 1, mitot
         k = irr(ix0,iy0)
         if (k .eq. -1) then  ! set to 0 for easier debugging
            qx(:,ix0,iy0) = 0.d0
            qy(:,ix0,iy0) = 0.d0
            qxx(:,ix0,iy0) = 0.d0
            qxy(:,ix0,iy0) = 0.d0
            qyy(:,ix0,iy0) = 0.d0
            go to 110
         endif

         if (k .eq. lstgrd) then
            nco = 2  ! testing new AG formulation of regular cells that uses 2
         else
            nco = 2
         endif
         if (all_nbors_exist(ix0,iy0,nco,irr,lstgrd,mitot,mjtot)) cycle
c     
         if (ar(k)/ar(lstgrd) .lt. gradThreshold) then
            qx(:,ix0,iy0) = 0.d0
            qy(:,ix0,iy0) = 0.d0
            qxx(:,ix0,iy0) = 0.d0
            qxy(:,ix0,iy0) = 0.d0
            qyy(:,ix0,iy0) = 0.d0
            go to 110    ! leave 0 gradient:  more stable for teeny cells w/o slopes
         endif
c     
c      # this cell needs derivatives
         if (k .ne. lstgrd) then 
            x0 = xcirr(k)
            y0 = ycirr(k)
         else
            x0 = xlow + (ix0-.5d0)*hx
            y0 = ylow + (iy0-.5d0)*hy
         endif

         ata = 0.d0
         rhs = 0.d0
         do 22 joff = -jjr(ix0,iy0), jjr(ix0,iy0)
         do 22 ioff = -iir(ix0,iy0), iir(ix0,iy0)
            if (ioff .eq. 0 .and. joff .eq. 0) cycle
            iyn = iy0 + joff
            ixn = ix0 + ioff
            if (.not. IS_REAL(ixn,iyn)) cycle
            kn = irr(ixn,iyn)
            if (kn .eq. -1) cycle

            if (kn .ne. lstgrd) then
               xn = xcirr(kn)
               yn = ycirr(kn)
            else
               xn = xlow + (ixn-.5d0)*hx
               yn = ylow + (iyn-.5d0)*hy
            endif

            deltax = xn - x0
            deltay = yn - y0
            ag(1) = deltax/hx
            ag(2) = deltay/hy
            ag(3) = -poly(8,1,k)
            ag(4) = -poly(9,1,k)
            ag(5) = -poly(10,1,k)
            ag(6) = -dcubicshifts(1,k)
            ag(7) = -dcubicshifts(2,k)
            ag(8) = -dcubicshifts(3,k)
            ag(9) = -dcubicshifts(4,k)

              ! handle cut cells first
              if (kn .ne. lstgrd) then
                arr = ar(kn)
                ivert = 1
                do while (poly(ivert+1,1,kn) .ne. -11)
                   ivert = ivert+1
                end do
                itri = ivert - 3
              else
                arr = hx*hy
                call makep(poly(1,1,lstgrd),ixn,iyn,xlow,ylow,hx,hy)
                itri = 2
              endif

              idx1 = 1
              x1 = poly(idx1,1,kn)
              y1 = poly(idx1,2,kn)
              do it = 1, itri
                idx2 = it + 1
                idx3 = it + 2
                x2 = poly(idx2,1,kn)
                y2 = poly(idx2,2,kn)
                x3 = poly(idx3,1,kn)
                y3 = poly(idx3,2,kn)
                artri = triangle_area(x1,x2,x3,y1,y2,y3)

                do itq = 1, ntriquad
                   xval = x1 * rtri(itq) + x2 * stri(itq)
     &                 +  x3 * (1.d0-rtri(itq)-stri(itq))
                   yval = y1 * rtri(itq) + y2 * stri(itq)
     &                 +  y3 * (1.d0-rtri(itq)-stri(itq))

                   ag(3) = ag(3) + (artri/arr) * wtri(itq) *
     &                               ((xval-x0)**2) / (hx**2)

                   ag(4) = ag(4) + (artri/arr) * wtri(itq) *
     &                               ((xval-x0)/hx)*((yval-y0)/hy)

                   ag(5) = ag(5) + (artri/arr) * wtri(itq) *
     &                              ((yval-y0)**2) / (hy**2)

                   ag(6) = ag(6) + (artri/arr) * wtri(itq) *
     &                             ((xval-x0)**3) / (hx**3)
                   ag(7) = ag(7) + (artri/arr) * wtri(itq) *
     &                           ((yval-y0)*(xval-x0)**2) / (hy*hx**2)
                   ag(8) = ag(8) + (artri/arr) * wtri(itq) *
     &                           ((xval-x0)*(yval-y0)**2) / (hx*hy**2)
                   ag(9) = ag(9) + (artri/arr) * wtri(itq) *
     &                           ((yval-y0)**3) / (hy**3)
                end do ! end loop over quadrature points
              end do ! end loop over triangles

c
c        ## form normal equations and solve via cholesky
c
          do 31 ii = 1, 9
          do 30 jj = 1, 9
             ata(jj,ii) = ata(jj,ii) + ag(ii)*ag(jj)
 30       continue
             rhs(:,ii) = rhs(:,ii) + ag(ii)*
     &                  (qp(:,ixn,iyn) - qp(:,ix0,iy0))
 31       continue

 22       continue

          call cholesky(9,9,ata,G)

             do 60 mm = 1, nvar
                call solve(9,9,G,rhs(mm,:))

                qyyy(mm,ix0,iy0) = 6.d0*rhs(mm,9)/(hy**3)
                qxyy(mm,ix0,iy0) = 2.d0*rhs(mm,8)/(hx*hy**2)
                qxxy(mm,ix0,iy0) = 2.d0*rhs(mm,7)/(hy*hx**2)
                qxxx(mm,ix0,iy0) = 6.d0*rhs(mm,6)/(hx**3)
                qyy(mm,ix0,iy0) =  2.d0*rhs(mm,5)/(hy**2)
                qxy(mm,ix0,iy0) =  rhs(mm,4)/(hx*hy)
                qxx(mm,ix0,iy0) =  2.d0*rhs(mm,3)/(hx**2)
                qy(mm,ix0,iy0)  =  rhs(mm,2)/hy
                qx(mm,ix0,iy0)  =  rhs(mm,1)/hx
c     
 60       continue
c  
 110  continue
c
 120  if (prflag) then
         write(21,*)' qx '
         call prDeriv(qx,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qy '
         call prDeriv(qy,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qxx'
         call prDeriv(qxx,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qxy'
         call prDeriv(qxy,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qyy'
         call prDeriv(qyy,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qxxx'
         call prDeriv(qxxx,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qxxy'
         call prDeriv(qxxy,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qxyy'
         call prDeriv(qxyy,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qyyy'
         call prDeriv(qyyy,irr,mitot,mjtot,nvar,lstgrd)

      endif
c
 99   return
      end

