c
c ---------------------------------------------------------------
c
       subroutine addMoreNeighbors(irr,nlist,nst,nend,newend,mitot,
     .            mjtot,lstgrd,quad,xlow,ylow,hx,hy,lwidth,nhood,istage)

       implicit double precision (a-h, o-z)
       include "cirr.i"
       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                  ismp,gradThreshold
       parameter(nsize=25)
       dimension nlist(nsize,2), irr(mitot,mjtot)

       logical outside, quad, IS_GHOST
       logical verbose/.false./

       outside(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     &                 j .lt. 1 .or. j .gt. mjtot)
       l1_len(ixn,iyn) = iabs(ixn-nlist(1,1))+iabs(iyn-nlist(1,2))

       IS_GHOST(i,j) = (i .lt. ist .or. i .gt. iend .or.
     .                  j .lt. jst .or. j .gt. jend)

c      IS_FAR_GHOST(i,j,ioff,joff) =((i.le.lwidth .and. ioff.lt.0) .or.
c    .                               (j.le.lwidth .and. joff.lt.0) .or.
c    .                     (i .gt. mitot-lwidth .and. ioff .gt. 0) .or.
c    .                     (j .gt. mjtot-lwidth .and. joff .gt. 0))

c this block sets usable cell area, depends on stage 
c
       if (istage .eq. 1) then
           ist = lwidth-1
           jst = lwidth-1
           iend = mitot-lwidth+2
           jend = mjtot-lwidth+2
        else if (istage .eq. 2) then
           ist = lwidth+1
           jst = lwidth+1
           iend = mitot-lwidth
           jend = mjtot-lwidth
        endif
c
c nhood is size of search area around cell to check. 
       newend = nend

c search for neighbors of celsl already on nlist. adapted from qsloeps routine. actually only
c searches nbors of the "one" cell on the list.
c
       do 10 n = nst, nend
c        ::: add neighbors of this point to nlist
         ix0 = nlist(n,1)
         iy0 = nlist(n,2)
         k0  = irr(ix0,iy0)

         do 20 ioff = -nhood, nhood
         do 20 joff = -nhood, nhood
           if (ioff .eq. 0 .and. joff .eq. 0) go to 20  ! this point already on list
c           ::: prerequisite for edge sharing. if commented out then diag nbors allowed.
c           if (ioff.ne.0 .and. joff.ne.0) go to 20
           ixn = ix0 + ioff
           iyn = iy0 + joff
           kn  = irr(ixn,iyn)
           if (kn .eq. -1) go to 20
           if (outside(ixn,iyn)) go to 20
           if (IS_GHOST(ix0+ioff,iy0+joff)) go to 20
c           if (kn .ne. lstgrd) go to 20 ! see if this improves stability
           if (ar(kn).lt. .5d0*ar(lstgrd)) go to 20  ! dont use small cut cells in stencil
c           ::: dont use points 2 away in each coordinate direction
c           if (l1_len(ixn,iyn) .ge. 3) go to 20
           if (verbose .and. l1_len(ixn,iyn).gt.4) then
              write(*,*)" not using ",ixn,iyn," in stencil for ",ix0,iy0
              go to 20 ! changed since using bigger  nhood and no cuts
           endif
           if (kn .eq. lstgrd  .or. k0 .eq. lstgrd) go to 25 ! thin bodies check
c
c  this code did not handle cells diagonally away that you want in the stencil
c  so it is commented out. this means that it cannot handle thin bodies
c  without fixing up this check.
c
!--c           ::: both cells cut. must check that really share
!--c           ::: a common edge
!--           do 23 kside = 1,6
!--              if (poly(kside+2,1,k0) .eq. -11) go to 20
!--              x1 = poly(kside,1,k0)
!--              y1 = poly(kside,2,k0)
!--              x2 = poly(kside+1,1,k0)
!--              y2 = poly(kside+1,2,k0)
!--              if (x1 .ne. x2) then
!--                 ihoriz = 1
!--                 ivert  = 0
!--              else
!--                 ihoriz = 0
!--                 ivert  = 1
!--              endif
!--c     compute indices of adjacent cell
!--              ixr = ix0 + ivert*isig(y1-y2)
!--              iyr = iy0 + ihoriz*isig(x2-x1)
!--              if ((ixr .eq. ixn) .and. (iyr .eq. iyn)) go to 25
!-- 23        continue
!--c     ::: should never fall through
c     
c     ::: make sure new cell not already on list (might have been added 1st time)
 25        do 30 ncheck = 1, newend
              if (ixn .eq. nlist(ncheck,1) .and.
     .            iyn .eq. nlist(ncheck,2)) go to 20
 30        continue
           if (kn .ne. lstgrd) then
 1            xn = xcirr(kn)
              yn = ycirr(kn)
           else
              xn = xlow + (ixn-.5d0)*hx
              yn = ylow + (iyn-.5d0)*hy
           endif
c     ::: dont use new cell if ghost cell, want 1-sided diffs like cart3d
c     if (xn .gt. xprob .or. xn .lt. 0.) go to 20
c     ### dont use new cell if area way too small (same test for no gradients in this cell)
c     in irreg3hbox
           if (ar(kn)/ar(lstgrd) .lt. gradThreshold) then
              write(*,*)"left out ",ixn,iyn," volfrac",ar(kn)/ar(lstgrd)
              write(*,*)"    in neighbor stencil for cell",ix0,iy0
              go to 20
           endif
           newend = newend + 1
           if (newend .gt. nsize) then
              write(*,*) " need to increase nlist size from ",nsize,
     .                   " in addMoreNeighbors.f"
              stop
           endif
           nlist(newend,1) = ixn 
           nlist(newend,2) = iyn 
 20     continue

 10   continue

       return
       end
