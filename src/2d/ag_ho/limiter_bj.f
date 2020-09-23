      subroutine limiter(qp, qx, qy, mitot, mjtot, irr, nvar, hx, hy,
     . lwidth, lstgrd, xlow, ylow)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)

       common /order2/ ssw, primitive, nolimiter
       dimension nlist(25,2)
       logical primitive
       logical quadrecon
       include "./reconinfo.i"
       include "cirr.i"

      logical  LP_success, output_debug, tricell
      double precision neighb_x(8),neighb_y(8), neighb_u(8)
      double precision Phi_x, Phi_y, LP_Dx,LP_Dy, LS_Dx, LS_Dy
      integer row_A
      logical do_pos_constr ! are we really a cutcell?
      integer my_iters, tot_iters, max_iters   ! maximum no of iterations needed in Simplex
      double precision avg_iters    ! avg number of iterations needed in Simplex
      integer count_cells
      logical   prflag, enufNbor, outside, almost_outside
      logical dump/.false./


       dimension dumax(4),dumin(4),phimin(4)
       dimension graddot(4),alpha(4),recon(4)
       logical IS_GHOST, verbose/.true./
       IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)




      outside(i,j) = ((i .lt. 1) .or. (j.lt. 1) .or.
     .                (i .gt. mitot) .or. (j.gt. mjtot))


       if(primitive .eqv. .false.) then
        print *, "PROBLEM: THE CURRENT LIMITER IMPLEMENTATION ONLY WORKS
     .            WHEN RECONSTRUCTING IN PRIMITIVE VARIABLES."
       endif

      ! quadratic limiting does not work right now.
      quadrecon = .false.

      ! limit using the MC limiter everywhere that you can
      do 50 i = 1, mitot
      do 50 j = 1, mjtot
        k = irr(i,j)
        if(k .eq. lstgrd .and. .not. outside(i+1,j)) then
        if(irr(i+1,j) .eq. lstgrd) then

            ! simple MC limiter
            do m = 1,nvar
                dfx = 2.d0*( qp(i+1,j,m) - qp(i,j,m) ) / hx
                qx(i,j,m) = dminmod(qx(i,j,m) ,dfx)
            enddo
        endif
        endif

        if(k .eq. lstgrd .and. .not. outside(i-1,j) ) then
        if( irr(i-1,j) .eq. lstgrd) then
            ! simple MC limiter
            do m = 1,nvar
                dbx = 2.d0*( qp(i,j,m) - qp(i-1,j,m) ) / hx
                qx(i,j,m) = dminmod(dbx, qx(i,j,m) )
            enddo
        endif
        endif

        if(k .eq. lstgrd .and. .not. outside(i,j+1)) then
        if( irr(i,j+1) .eq. lstgrd ) then

            ! simple MC limiter
            do m = 1,nvar
                dfy = 2.d0*( qp(i,j+1,m) - qp(i,j,m) ) / hy
                qy(i,j,m) = dminmod(qy(i,j,m) ,dfy)
            enddo
        endif
        endif

        if(k .eq. lstgrd .and. .not. outside(i,j-1) ) then
        if( irr(i,j-1) .eq. lstgrd ) then

            ! simple MC limiter
            do m = 1,nvar
                dby = 2.d0*( qp(i,j,m) - qp(i,j-1,m) ) / hy
                qy(i,j,m) = dminmod(dby, qy(i,j,m) )
            enddo
        endif
        endif



50    continue


      ! I still must
      ! limit the cut cells
      ! AND
      ! whole cells adjacent the cut cells

c
c  **************************************
c
c   LIMIT CUT CELLS USING LINEARITY PRESERVING BARTH-JESPERSEN OVER RECONSTRUCTION NEIGHBORHOOD
c

      do 30 j = 1, mjtot
      do 30 i = 1, mitot
          k = irr(i,j)
          if (k .eq. -1 .or. k .eq. lstgrd) goto 30 ! solid/full cells have no need for a gradient

           call getCellCentroid(lstgrd,i,j,xcent,ycent,xlow,
     .                            ylow,hx,hy,k)


          ! find max and min needed for BJ limiting
          dumax = 0.d0
          dumin = 0.d0
          do 31 joff = -jjr(i,j),jjr(i,j)
          do 31 ioff = -iir(i,j),iir(i,j)
            if (ioff .eq. 0 .and. joff .eq. 0) go to 31
            if (IS_GHOST(i+ioff,j+joff)) go to 31

            koff = irr(i+ioff,j+joff)
            if (koff .eq. -1) go to 31
            dumax = max(dumax,qp(i+ioff,j+joff,:)-qp(i,j,:))
            dumin = min(dumin,qp(i+ioff,j+joff,:)-qp(i,j,:))
  31      continue

          phimin = 1.d0
          do 32 joff = -jjr(i,j), jjr(i,j)
          do 32 ioff = -iir(i,j), iir(i,j)
             if (ioff .eq. 0 .and. joff .eq. 0) go to 32 ! no eqn to solve
             if (IS_GHOST(i+ioff,j+joff)) go to 32
             koff = irr(i+ioff,j+joff)
             if (koff .eq. -1) go to 32
             call getCellCentroid(lstgrd,i+ioff,j+joff,xc,yc,xlow,
     .                            ylow,hx,hy,koff)
              diffx = xc-xcent
              diffy = yc-ycent
              graddot  = qx(i,j,:)*diffx + qy(i,j,:)*diffy
              recon = qp(i,j,:) + graddot
                do m = 1,4
                   if (graddot(m) > 1.d-10) then
                      alpha(m) = min(1.d0, dumax(m)/graddot(m))
                   else if (graddot(m) < -1.d-10) then
                      alpha(m) = min(1.d0, dumin(m)/graddot(m))
                   else
                      alpha(m) = 1.d0
                   endif
                end do
                if (recon(1) .le. 0.d0) alpha(1) = 0.d0
                if (recon(4) .le. 0.d0) alpha(4) = 0.d0
                phimin = min(phimin, alpha)
  32      continue
          qx(i,j,:) = qx(i,j,:)*phimin(:)
          qy(i,j,:) = qy(i,j,:)*phimin(:)


      ! positivity at quadrature points is required, midpoints because p is always = 1.
         do 510 kside=1,6


               x1 = poly(kside,1,k)
               y1 = poly(kside,2,k)
               x2 = poly(kside+1,1,k)
               y2 = poly(kside+1,2,k)

               d_bxp = (x1+x2)/2.d0
               d_byp = (y1+y2)/2.d0

               diffx = d_bxp-xcent
               diffy = d_byp-ycent
               graddot  = qx(i,j,:)*diffx + qy(i,j,:)*diffy
               recon = qp(i,j,:) + graddot


              alpha = 1.d0
              if (recon(1) .le. 0.d0) alpha(1) = 0.d0
              if (recon(4) .le. 0.d0) alpha(4) = 0.d0

!                if (alpha(1) < 0.5d0) then
!                print *, "here"
!                endif

              qx(i,j,1) = qx(i,j,1)*alpha(1)
              qy(i,j,4) = qy(i,j,4)*alpha(4)
              if (poly(kside+2,1,k).eq.-11.) exit
 510      continue

!      alpha = 1.d0
!      if (recon(1) .le. 0.d0 .and. graddot(1) > 0) then
!      alpha(1) = (d_eps-qp(i,j,1)) / graddot(1)
!      alpha(1) = min(alpha(1), 1.d0)
!      endif
!
!      if (recon(4) .le. 0.d0 .and. graddot(4) > 0) then
!      alpha(4) = (d_eps-qp(i,j,4)) / graddot(4)
!      alpha(4) = min(alpha(4), 1.d0)
!      endif
!
!                if (alpha(1) < 1.d0) then
!                print *, "here"
!                endif
!
!              qx(i,j,1) = qx(i,j,1)*alpha(1)
!              qy(i,j,4) = qy(i,j,4)*alpha(4)





  30  continue

c
c  **************************************
c
c   LIMIT FLOW CELLS ADJACENT TO CUTS USING RECENTERING
c
      do 210 ix0 = 1, mitot
      do 210 iy0 = 1, mjtot
        k = irr(ix0,iy0)
        if (k .eq. -1) go to 210
        if (k .ne. lstgrd) go to 210


            x0 = xlow + (ix0-.5d0)*hx
            y0 = ylow + (iy0-.5d0)*hy

            do m = 1, nvar
c
c  for each irregular side, use recentering. adjacent cant be solid since this cell regular
       if(.not. outside(ix0+1,iy0)) then
       if (irr(ix0+1,iy0) .ne. lstgrd) then
            kp1 = irr(ix0+1,iy0)
            xcut = xcirr(kp1)
            ycut = ycirr(kp1)
            recentM =   qp(ix0+1,iy0,m) + (y0-ycut)*qy(ix0+1,iy0,m)
            cutmin= (recentM - qp(ix0,iy0,m))/(xcut-x0)

             othermin = qx(ix0,iy0,m)
             qx(ix0,iy0,m) = .5*min(abs(cutmin),abs(othermin)) *
     &                         (sign(1.d0,othermin)+sign(1.d0,cutmin))
!            endif
       endif
       endif
       if(.not. outside(ix0-1,iy0)) then
       if (irr(ix0-1,iy0) .ne. lstgrd) then
            km1 = irr(ix0-1,iy0)
            xcut = xcirr(km1)
            ycut = ycirr(km1)
            recentM =   qp(ix0-1,iy0,m) + (y0-ycut)*qy(ix0-1,iy0,m)
            cutmin= (recentM - qp(ix0,iy0,m))/(xcut-x0)

            othermin = qx(ix0,iy0,m)
            qx(ix0,iy0,m) = .5*min(abs(cutmin),abs(othermin))*
     &                     (sign(1.d0,othermin)+sign(1.d0,cutmin))
       endif
       endif
c  for each irregular side, use recentering.
       if(.not. outside(ix0,iy0+1)) then
       if (irr(ix0,iy0+1) .ne. lstgrd) then
            kp1 = irr(ix0,iy0+1)
            xcut = xcirr(kp1)
            ycut = ycirr(kp1)
            recentM =   qp(ix0,iy0+1,m) + (x0-xcut)*qx(ix0,iy0+1,m)
            cutmin= (recentM - qp(ix0,iy0,m))/(ycut-y0)

            othermin = qy(ix0,iy0,m)
            qy(ix0,iy0,m) = .5*min(abs(cutmin),abs(othermin))*
     &                          (sign(1.d0,othermin)+sign(1.d0,cutmin))
!            endif
       endif
       endif
       if(.not. outside(ix0,iy0-1)) then
       if (irr(ix0,iy0-1) .ne. lstgrd) then
            km1 = irr(ix0,iy0-1)
            xcut = xcirr(km1)
            ycut = ycirr(km1)
            recentM =   qp(ix0,iy0-1,m) + (x0-xcut)*qx(ix0,iy0-1,m)
            cutmin= (recentM - qp(ix0,iy0,m))/(ycut-y0)

            othermin = qy(ix0,iy0,m)
            qy(ix0,iy0,m) = .5* min(abs(cutmin),abs(othermin))*
     &                      (sign(1.d0,othermin)+sign(1.d0,cutmin))

       endif
       endif

       end do
 210  continue







      return
      end

      function dminmod(a,b)
      implicit double precision(a-h,o-z)
      if( a > 0 .and. b > 0 ) then
        dminmod = min(a,b)
      elseif( a < 0 .and. b < 0 ) then
        dminmod = max(a,b)
      else
        dminmod = 0.d0
      endif
      end





        ! apparently recentering + Sandra's limiter does not work

      ! limit the cut cells first
!      do 110 ix0 = 1, mitot
!      do 110 iy0 = 1, mjtot
!         k = irr(ix0,iy0)
!         if (k .eq. -1) go to 110
!         if (k .eq. lstgrd) go to 110
!
!
!         x0 = xcirr(k)
!         y0 = ycirr(k)
!
!
!         do 37 kside = 1, 6
!            if (poly(kside+2,1,k).eq.-11)then
!               sidex2 = poly(kside,1,k)
!               sidey2 = poly(kside,2,k)
!               sidex1 = poly(kside+1,1,k)
!               sidey1 = poly(kside+1,2,k)
!               go to 39
!            endif
! 37      continue
! 39      bxpt = .5d0*(sidex1+sidex2)
!         bypt = .5d0*(sidey1+sidey2)
!
!
!
!
!
!
!
!         nlist(1,1) = ix0
!         nlist(1,2) = iy0
!         nst        = 1
!         nend       = 1
!         call addneighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
!     .        lstgrd,quadrecon,xlow,ylow,hx,hy)
!
!c
!c convert from my routine to Sandras cut cell lp routine
!c
!          center_x = x0
!          center_y = y0
!
!          num_neighb = newend-1  ! not counting cell itself
!          output_debug = .false.
!
!
!          do m = 1, 4   ! loop over the 4 variables
!            row_A     = 2*num_neighb + 4
!            do_pos_constr = .false.
!            if (irr(ix0,iy0) .ne. lstgrd) then
!               if (m .eq. 1 .OR. m .eq. 4) then
!                  do_pos_constr = .true.
!                  row_A = row_A + 1 ! additional pos constraint for rho, p
!               endif
!            endif
!
!            do inbor = 1, num_neighb
!              ixn = nlist(inbor+1,1)
!              iyn = nlist(inbor+1,2)
!              kn =  irr(ixn,iyn)
!              if (kn .ne. lstgrd) then
!                 xn = xcirr(kn)
!                 yn = ycirr(kn)
!              else
!                 xn = xlow + (ixn-.5d0)*hx
!                 yn = ylow + (iyn-.5d0)*hy
!              endif
!
!              neighb_x(inbor) = xn ! neighboring centroid, not diffs
!              neighb_y(inbor) = yn ! neighboring centroid, not diffs
!
!              neighb_u(inbor) = qp(ixn,iyn, m)  ! neighboring u, need to loop over m
!            end do  ! do inbor = 1, num_neighb
!              center_u = qp(ix0,iy0,m)
!              LS_Dx = qx(ix0,iy0,m)
!              LS_Dy = qy(ix0,iy0,m)
!
!            LP_success = .false.
!
!
!
!             call vector_lim_allineq(num_neighb,neighb_x,neighb_y,
!     &               neighb_u,center_x,center_y,center_u,output_debug,
!     &               row_A,Phi_x,Phi_y, LP_Dx,LP_Dy,
!     &         LP_success,LS_Dx,LS_Dy,bxpt,bypt,do_pos_constr,my_iters)
!
!c           Increase the count of number of iterations
!            tot_iters = tot_iters + my_iters
!            if (my_iters .GT. max_iters) then
!               max_iters = my_iters
!            endif
!            count_cells = count_cells + 1
!
!
! 402         format("sandras lsq ",2e15.7,10x," limited ",2e15.7)
!
!             if (LP_success) then
!                  qx(ix0,iy0,m) = LP_Dx
!                  qy(ix0,iy0,m) = LP_Dy
!
! 400              format("limiter values for m = ",i3," are ",2e15.7)
!             else
!                  write(*,*)" LP failure for cell ",ix0,iy0
!c                  stop
!             endif
!
!          end do   ! do m=1,4
!
!
!
!
!
!110   continue
!
!
!c
!c  **************************************
!c
!c   LIMIT FLOW CELLS ADJACENT TO CUTS USING RECENTERING
!c
!      do 210 ix0 = 1, mitot
!      do 210 iy0 = 1, mjtot
!        k = irr(ix0,iy0)
!        if (k .eq. -1) go to 210
!        if (k .ne. lstgrd) go to 210
!
!
!            x0 = xlow + (ix0-.5d0)*hx
!            y0 = ylow + (iy0-.5d0)*hy
!
!            do m = 1, nvar
!c
!c  for each irregular side, use recentering. adjacent cant be solid since this cell regular
!       if(.not. outside(ix0+1,iy0)) then
!       if (irr(ix0+1,iy0) .ne. lstgrd) then
!            kp1 = irr(ix0+1,iy0)
!            xcut = xcirr(kp1)
!            ycut = ycirr(kp1)
!            recentM =   qp(ix0+1,iy0,m) + (y0-ycut)*qy(ix0+1,iy0,m)
!            cutmin= (recentM - qp(ix0,iy0,m))/(xcut-x0)
!
!             othermin = qx(ix0,iy0,m)
!             qx(ix0,iy0,m) = .5*min(abs(cutmin),abs(othermin)) *
!     &                         (sign(1.d0,othermin)+sign(1.d0,cutmin))
!!            endif
!       endif
!       endif
!       if(.not. outside(ix0-1,iy0)) then
!       if (irr(ix0-1,iy0) .ne. lstgrd) then
!            km1 = irr(ix0-1,iy0)
!            xcut = xcirr(km1)
!            ycut = ycirr(km1)
!            recentM =   qp(ix0-1,iy0,m) + (y0-ycut)*qy(ix0-1,iy0,m)
!            cutmin= (recentM - qp(ix0,iy0,m))/(xcut-x0)
!
!            othermin = qx(ix0,iy0,m)
!            qx(ix0,iy0,m) = .5*min(abs(cutmin),abs(othermin))*
!     &                     (sign(1.d0,othermin)+sign(1.d0,cutmin))
!       endif
!       endif
!c  for each irregular side, use recentering.
!       if(.not. outside(ix0,iy0+1)) then
!       if (irr(ix0,iy0+1) .ne. lstgrd) then
!            kp1 = irr(ix0,iy0+1)
!            xcut = xcirr(kp1)
!            ycut = ycirr(kp1)
!            recentM =   qp(ix0,iy0+1,m) + (x0-xcut)*qx(ix0,iy0+1,m)
!            cutmin= (recentM - qp(ix0,iy0,m))/(ycut-y0)
!
!            othermin = qy(ix0,iy0,m)
!            qy(ix0,iy0,m) = .5*min(abs(cutmin),abs(othermin))*
!     &                          (sign(1.d0,othermin)+sign(1.d0,cutmin))
!!            endif
!       endif
!       endif
!       if(.not. outside(ix0,iy0-1)) then
!       if (irr(ix0,iy0-1) .ne. lstgrd) then
!            km1 = irr(ix0,iy0-1)
!            xcut = xcirr(km1)
!            ycut = ycirr(km1)
!            recentM =   qp(ix0,iy0-1,m) + (x0-xcut)*qx(ix0,iy0-1,m)
!            cutmin= (recentM - qp(ix0,iy0,m))/(ycut-y0)
!
!            othermin = qy(ix0,iy0,m)
!            qy(ix0,iy0,m) = .5* min(abs(cutmin),abs(othermin))*
!     &                      (sign(1.d0,othermin)+sign(1.d0,cutmin))
!
!       endif
!       endif
!
!       end do
! 210  continue

