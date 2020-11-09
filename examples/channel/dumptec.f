c
c -----------------------------------------------------
c
      subroutine dumptec (lst,lend,nvar,naux,nplot,time)
c
      use amr_module
      implicit double precision (a-h,o-z)
      include "cuserdt.i"
      common /order2/ ssw, quad, nolimiter
      logical         flag,quad,nolimiter
      character*23  filename 
      character*11   filebndry
c
c dumptec = make tecplot file for finest level soln values over entire domain
c (but for now assume only 1 level  -- output finest level)
c
c nplot == plotnum (0 for initial conditions, 1 (or more) for final conditions
c
      filename = 'disjointErrPlanesxx.dat'
      filename(18:18) = '0'
      filename(19:19) = char(ichar('0')+nplot)
      open(13,file=filename,status='unknown',form='formatted')
      write(13,100) 

      filename = 'disjointCutPlanesxx.dat'
      filename(18:18) = '0'
      filename(19:19) = char(ichar('0')+nplot)

      open(14,file=filename,status='unknown',form='formatted')
      write(14,100) 
 100  format('TITLE = "Extracted cutting Planes through mesh"' )

      ! also dump boundary plots at same time
      filebndry = 'bndry0x.dat'
      filebndry(7:7) = char(ichar('0')+nplot)
      ibunit = 24
      open(ibunit,file=filebndry,status='unknown',form='formatted')
      write(ibunit,*)"# writing bndry file at time ",time
      write(ibunit,*)"#  dist     x_bndry  y_bndry  rho    u     v  ",
     &               "   p       i       j      mptr"
      

      nplot = nplot+1

 10   level = lfine

c     initialize for error computation across all grids at this level
      volDenErrorL1   = 0.d0
      volExactDenL1   = 0.d0
      exactVol        = 0.d0
      volDenErrorMax  = 0.d0

      bndryDenErrorL1  = 0.d0
      bndryDenExactL1 = 0.d0
      bndryDenErrorMax = 0.d0
      imax = -1
      jmax = -1


      
c     fill ghost cells if need them for output or to compute gradients
      mptr = lstart(level)
 20       if (mptr .eq. 0) go to 50
              nx = node(ndihi,mptr)-node(ndilo,mptr)+1
              ny = node(ndjhi,mptr)-node(ndjlo,mptr)+1
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              mitot = nx + 2*nghost
              mjtot = ny + 2*nghost
              locnew = node(store1,mptr)
              locaux = node(storeaux,mptr)
              locirr = node(permstore,mptr)
              lstgrd = node(lstptr,mptr)
              xl   = rnode(cornxlo, mptr)
              yb   = rnode(cornylo, mptr)
              xr   = rnode(cornxhi, mptr)
              yt   = rnode(cornyhi, mptr)
              hx   = hxposs(level)
              hy   = hyposs(level)

             !!if (ssw .ne. 0 .and. .not. pwconst)
               ! 1 means to call external boundary conditions
               ! 0 leaves values in ghost cells, easier to debug
c              istage = 1  
c              call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,
c    &                    mptr,alloc(locaux),naux,istage,
c    &                    alloc(locirr),lstgrd)

              ! remember got 3 times size of irr to include other arrays
              locirr = node(permstore,mptr)
              locnumHoods = locirr + mitot*mjtot
              locreconx = locnumHoods + mitot*mjtot
              locrecony = locreconx + mitot*mjtot

              call outtec(alloc(locnew),nvar,mptr,
     1                    alloc(locirr),mitot,mjtot,
     2                    lstgrd,hx,hy,xlow,ylow,time,
     3                    alloc(locnumHoods),
     4                    ibunit,alloc(locreconx),alloc(locrecony),
     5                    volDenErrorL1,volExactDenL1,exactVol,
     6                    volDenErrorMax,bndryDenErorL1,bndryDenExactL1,
     7                    bndryDenErrorMax,imax,jmax)
c
              mptr = node(levelptr,mptr)
          go to 20
50        continue
 
c
      close(14)
      close(24)
      write(*,*)"done writing ",filename," at time",time

      
c
c  output errors
      write(outunit,600) volDenErrorL1,volExactDenL1,
     .     volDenErrorL1/volExactDenL1,exactVol
      
 600  format("L1 density volume error ", e15.7,/,
     .       "L1 density exact soln   ", e15.7,/,
     .       "L1 Relative density error",e15.7,/,
     .       "Computed volume          ",e15.7,//)

      write(outunit,601) bndryDenErrorL1,bndryDenExactL1,
     &                   bndryDenErrorMax,imax,jmax
 601  format("L1 bndry density error ",e15.7,/,
     &       "L1 bndry density exact ",e15.7,/,
     &     "Max bndry density error ",e14.7,' at ',2i5)
      
      
 99   return
      end
