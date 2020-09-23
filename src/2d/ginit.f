c
!> Initializes solution on all grids at level of grid **msave**, 
!! by calling qinit().
!! If **first** = true, (first call to init), then allocate the
!! soln storage area too, else was already allocated.
!!
!! \param msave First grid in the grid list that stores all grids at
!! this level
!! \param first Is this the first call to init?
!! \param nvar number of equations for the system
!! \param naux  number of auxiliary variables
!! \param start_time start time of current simulation
!!
c
c -------------------------------------------------------------
c
      subroutine ginit(msave, first, nvar, naux, start_time)
c
      use amr_module
      implicit double precision (a-h,o-z)
      logical first
      logical quad, nolimiter
      common /order2/ ssw, quad, nolimiter
      include "cuserdt.i"

      

c ::::::::::::::::::::::::::::: GINIT ::::::::::::::::::::::::
c
c  initializes soln on all grids at 'level'  by calling qinit
c  if first = true, (first call to init), then allocate the
c  soln storage area too, else was already allocated.
c
c :::::::::::::::::::::::::::::::::::::::;::::::::::::::::::::

      if (msave .eq. 0) go to 99

      level = node(nestlevel,msave)
      hx    = hxposs(level)
      hy    = hyposs(level)
      mptr  = msave
 
 10       nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot   = nx + 2*nghost
          mjtot   = ny + 2*nghost
          corn1   = rnode(cornxlo,mptr)
          corn2   = rnode(cornylo,mptr)
          if(.not. (first)) go to 20
              loc                 = igetsp(mitot*mjtot*nvar)
              node(store1,mptr)   = loc
              if (naux .gt. 0) then
                locaux              = igetsp(mitot*mjtot*naux)
                do k = 1, mitot*mjtot*naux,naux  ! set first component of aux to signal that it
                   alloc(locaux+k-1) = NEEDS_TO_BE_SET ! needs val, wasnt copied from other grids
                end do
                
               call setaux(nghost,nx,ny,corn1,corn2,hx,hy,
     &                    naux,alloc(locaux))

              else 
                locaux = 1
              endif
              node(storeaux,mptr) = locaux
              ! get space for irr, and for ncount and numhoods too
              ! adding space for recon status 
              node(permstore,mptr) = igetsp(5*mitot*mjtot) 
              locirr = node(permstore,mptr)
              lociir = locirr + 3*mitot*mjtot
              locjjr = lociir + mitot*mjtot
              call setirr(mitot,mjtot,mptr,quad,
     &                    gradThreshold,alloc(locirr),
     &                    alloc(lociir),alloc(locjjr))
              lstgrd = node(lstptr,mptr)
              ! do it all the time now, needed for RK
              !if (level .lt. mxnest) then
                loc2              = igetsp(mitot*mjtot*nvar)
                node(store2,mptr) = loc2
              !endif
              rnode(timemult, mptr) = start_time
              go to 30
 20       continue
c
c  if 2nd time through, put initial values in store2 so finer grids
c  can be advanced with interpolation of their boundary values.
c  new time soln should still be in location store1.
c
          loc    = node(store2,mptr)
          locaux = node(storeaux,mptr)
          locirr = node(permstore,mptr)
          lstgrd = node(lstptr,mptr)
c
   30     continue
          call qinit(nvar,nghost,nx,ny,corn1,corn2,hx,hy,
     &               alloc(loc),naux,alloc(locaux),
     &               lstgrd,alloc(locirr))

c
          mptr  = node(levelptr, mptr)
      if (mptr .ne. 0) go to 10
c
c
 99   continue
      return
      end
