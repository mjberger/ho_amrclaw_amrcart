c
c -------------------------------------------------------------------
c
      subroutine makeFVReconHood(irr,mitot,mjtot,lwidth,lstgrd,dx,dy,
     &                           xlow,ylow,iir,jjr,igradChoice)

      use amr_module

      implicit double precision (a-h, o-z)
      dimension irr(mitot,mjtot)
      dimension iir(mitot,mjtot),jjr(mitot,mjtot)
      logical IS_REAL, IS_BORDER
      logical is_reg

      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)

      IS_BORDER(i,j) = (i .eq. 1 .or. i .eq. mitot .or.
     .                  j .eq. 1 .or. j .eq.  mjtot)

c  igradChoie = 4 fits ordinary quadratic
c  igradChoie = 5 fits high-order accurate  quadratic
c                 by fitting cubic but then ignoring cubic terms

      if (igradChoice .eq. 4) then
        reconTOLx = 1.5d0*dx  ! for cut cells
        reconTOLy = 1.5d0*dy
        iir = 1
        jjr = 1
      else if (igradChoice .eq. 5) then
        reconTOLx = 2.5d0*dx  ! for cut cells
        reconTOLy = 2.5d0*dy
        iir = 2
        jjr = 2
      endif


!      set up for quadratic slopes on 5 by 5 neighborhood for cut cells
!      and 3 by 3 neighborhood for neighbor of cut cells.  Check for
!      stability and enlarge nhood as needed.

      do j = 1, mjtot   ! does first/last cell need a gradient?
      do i = 1, mitot
         k = irr(i,j)
         if (k .eq. -1) cycle ! solid so do nothing
         if ((k .eq. lstgrd) .and. 
     &       is_reg(irr,mitot,mjtot,i,j,lstgrd,3)) cycle 
         if (k .ne. lstgrd) then
            iir(i,j) = iir(i,j) + 1
            jjr(i,j) = jjr(i,j) + 1
         endif

         if( k .eq. lstgrd) then
             xi = xlow + (i-.5d0)*dx
             yi = ylow + (j-.5d0)*dy
         else
             xi = xcirr(k)
             yi = ycirr(k)
         endif

         icont = 1
         do while(icont .eq. 1)
            icont = 0

            diffx = -1.d0
            diffy = -1.d0

            do 31 joff = -jjr(i,j), jjr(i,j)
            do 30 ioff = -iir(i,j), iir(i,j)

                if (ioff .eq. 0 .and. joff .eq. 0) go to 30
                if (.not. IS_REAL(i+ioff,j+joff)) go to 30
                koff = irr(i+ioff, j+joff)
                if (koff .eq. -1) goto 30

                if( koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*dx
                    yoff = ylow + (j+joff-.5d0)*dy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

                diffx = max(diffx, dabs(xoff- xi) )
                diffy = max(diffy, dabs(yoff- yi) )
30          continue
31          continue

              if(diffx < reconTOLx) then
                 iir(i,j) = iir(i,j) + 1
                 icont = 1
              endif
              if(diffy < reconTOLy) then
                 jjr(i,j) = jjr(i,j) + 1
                 icont = 1
              endif

         end do  ! end while loop

      end do 
      end do

      return
      end subroutine
c
c -------------------------------------------------------------------
c
      logical function is_reg(irr,mitot,mjtot,i,j,lstgrd,nhoodSize)

      integer irr(mitot,mjtot)
      logical IS_REAL

      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
      
      ! nhoodSize is total size, so divide by 2 (and truncate) to get halfwidth
      nhsize = nhoodSize/2
      is_reg = .true.

      do joff = -nhsize, nhsize
      do ioff = -nhsize, nhsize
         if (.not. IS_REAL(i+ioff,j+joff)) then
            is_reg = .false.
            return
         endif
         if (irr(i+ioff,j+joff) .ne. lstgrd) then
            is_reg = .false.
            return
         endif
      end do
      end do

      return
      end
