c
c ----------------------------------------------------------------------
c
        subroutine  makeMergeGradHood(irr,lwidth,mitot,mjtot,lstgrd,
     &                                dx,dy,xlow,ylow,mptr,mioff,mjoff)

        use amr_module
        implicit double precision (a-h,o-z)
        include "cuserdt.i"
        dimension irr(mitot,mjtot)
        dimension mioff(mitot,mjtot),mjoff(mitot,mjtot)
        logical NOT_OK_GHOST, OKstencil, IS_OUTSIDE, OUT_OF_RANGE
        logical debug/.false./

        IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                     y .lt. ylower .or. y .gt. yupper)

        NOT_OK_GHOST(i,j) = (i .lt. 3 .or. 
     .                       i .gt. mitot-2 .or.
     .                       j .lt. 3 .or. 
     .                       j .gt. mjtot-2)


        OUT_OF_RANGE(i,j) =  (i .lt. 1 .or. i .gt. mitot .or.
     .                        j .lt. 1 .or. j .gt. mjtot)


        ! to get same as previous behavior uncomment next line
        ! if (numMergeTerms .eq. 2) return  

        numMergeTerms = 5
        if (numMergeTerms .eq. 2) then
          initVal = 1   ! initial neighborhood is 1 to each side
          tolx = 0.5d0 * dx
          toly = 0.5d0 * dy
        else ! numMergeTerms is 5, want 2nd gradient
          tolx = 1.5d0*dx 
          toly = 1.5d0*dy 
          initVal = 2   ! to fit quadratic need more cells 
        endif
        mioff = initVal
        mjoff = initVal

        !do j = lwidth+1, mjtot-lwidth
        !do i = lwidth+1, mitot-lwidth
        do j = 1, mjtot
        do i = 1, mitot
           k = irr(i,j)
           if (k .eq. -1 .or. k .eq. lstgrd) cycle !default mi/mjoff works

           OKstencil = .false.
           do while (.not. OKstencil) 
             diffx = -1.d0
             diffy = -1.d0
             ! form nhood, see if stable or need to expand
             OKstencil = .true.
             do joff = -mjoff(i,j), mjoff(i,j)
             do ioff = -mioff(i,j), mioff(i,j)
                if (OUT_OF_RANGE(i+ioff,j+joff)) cycle
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) cycle
                if (koff .eq. lstgrd) then
                   xm = xlow + (i+ioff-0.5d0)*dx
                   ym = ylow + (j+joff-0.5d0)*dy
                else
                   xm = xcentMerge(koff)
                   ym = ycentMerge(koff)
                endif
                if (IS_OUTSIDE(xm,ym) .or.
     &              NOT_OK_GHOST(i+ioff,j+joff)) cycle
                diffx = max(diffx, dabs(xm-xcentMerge(k)))
                diffy = max(diffy, dabs(ym-ycentMerge(k)))
             end do
             end do

             if (diffx .lt. tolx) then
                mioff(i,j) = mioff(i,j) + 1
                if (debug) write(*,900) i, j
 900            format("increasing mioff for cell ",2i5)
                OKstencil = .false.
             endif
             if (diffy .lt. toly) then
                mjoff(i,j) = mjoff(i,j) + 1
                if (debug) write(*,901) i, j
 901            format("increasing mjoff for cell ",2i5)
                OKstencil = .false.
             endif
             end do   ! end while loop checking for good stencil
        end do
        end do

        return
        end
