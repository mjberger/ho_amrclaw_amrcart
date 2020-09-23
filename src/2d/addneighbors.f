c
c ---------------------------------------------------------------
c
       subroutine addneighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .            lstgrd,quad,xlow,ylow,hx,hy)

       use amr_module
       implicit double precision (a-h, o-z)

       include "cuserdt.i"
       dimension nlist(25,2), irr(mitot,mjtot)
       logical OUT_OF_BOUNDS, quad
       integer cutoffDist

       OUT_OF_BOUNDS(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     &                       j .lt. 1 .or. j .gt. mjtot)

       L1_LEN(ixn,iyn) = iabs(ixn-nlist(1,1))+iabs(iyn-nlist(1,2))


         newend = nend

c        ::: add neighbors of this point to nlist
         ix0 = nlist(1,1)
         iy0 = nlist(1,2)
         k0  = irr(ix0,iy0)

         if (k0 .eq. lstgrd) then
c          # use 3 by 3 neighborhood
            nco = 1 
         else
            nco = 2 ! try larger stencil for cut cells since half missing roughly
         endif

 8     continue
c
c        # for each cell on the list, find if its suitable - not solid, not
c        #  and not already on the list
c        # should really preprocess and save

         do  15 joff = -nco,nco
         do  14 ioff = -nco,nco
           if (ioff .eq. 0 .and. joff .eq. 0) cycle  !dont count cell itself
c           ::: prerequisite for edge sharing if one of the offsets = 0
c           ## if next line commented out, then diagonal nbors are allwoed
           ixn = ix0 + ioff
           iyn = iy0 + joff
           if (OUT_OF_BOUNDS(ixn,iyn)) cycle   

           kn  = irr(ixn,iyn)
           if (kn .eq. -1) cycle    

c          ## commented out code to check for thin bodies
c          if (kn .eq. lstgrd  .or. k0 .eq. lstgrd) go to 25
c           ::: both cells cut. must check that really share
c           ::: a common edge
c          go to 25  ! skip next check for larger stencil
c          do 23 kside = 1,6
c             if (poly(kside+2,1,k0) .eq. -11) cycle
c             x1 = poly(kside,1,k0)
c             y1 = poly(kside,2,k0)
c             x2 = poly(kside+1,1,k0)
c             y2 = poly(kside+1,2,k0)
c             if (x1 .ne. x2) then
c                ihoriz = 1
c                ivert  = 0
c             else
c                ihoriz = 0
c                ivert  = 1
c             endif
c     compute indices of adjacent cell
c             ixr = ix0 + ivert*isig(y1-y2)
c             iyr = iy0 + ihoriz*isig(x2-x1)
c             if ((ixr .eq. ixn) .and. (iyr .eq. iyn)) go to 25
c23        continue
c          ::: should never fall through
c     
c          ::: make sure new cell not already on list
 25        do 30 ncheck = 1, newend
              if (ixn .eq. nlist(ncheck,1) .and.
     .            iyn .eq. nlist(ncheck,2)) cycle  
 30        continue
           if (kn .ne. lstgrd) then
 1            xn = xcirr(kn)
              yn = ycirr(kn)
           else
              xn = xlow + (ixn-.5d0)*hx
              yn = ylow + (iyn-.5d0)*hy
           endif

c          need to figure out using ghost cells - they need slopes too though
c          ::: dont use new cell if outside domain.  cant test if 
c          ::: ghost cell since might have multiple patches
c          if (xn .gt. xprob .or. xn .lt. 0. .or.  
c    &         yn .gt. yprob .or. yn .lt. 0.) cycle

c          ## use this cell
           newend = newend + 1
           nlist(newend,1) = ixn 
           nlist(newend,2) = iyn 
 14     continue
 15     continue

      if (newend .lt. 7) then  ! repeat with larger offset to find more nbors
         newend = nend  ! reset for another go round with larger nhood
         nco = nco + 1
         !! note that right now using square neighborhood
         !! want to figure out bad dir and only increase that
         go to 8 
      endif

       return
       end
