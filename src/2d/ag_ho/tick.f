cc
cc  -------------------------------------------------------------
cc
c      subroutine tick(nvar,iout,nstart,nstop,cut,cdist,
c     1                vtime,work,time,iousr,steady,lfix,quad,nplot)
cc
c      implicit double precision (a-h,o-z)
c      logical   tprint, vtime, gprint, steady, quad
c      logical            graf,lastout
c      parameter  (maxgr = 192, maxlv=12)
c      dimension work(nvar), intinc(maxlv)
c      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
c     *  newstl(maxlv),
c     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
c     *  lfine,iorder,mxnest,kcheck,lwidth,
c     *  maxcl, graf, lhead
c      include "calloc.i"
c      data      tprint /.true./,  gprint /.true./
cc
cc  parameters:
cc     nstop   = # of coarse grid time steps to be taken
cc     iout    = output interval every 'iout' coarse time steps
cc     vtime   = true for variable timestep, calculated each coarse step
cc
cc  main driver routine,  controls:
cc        integration  of all grids.
cc        error estimation / regridding
cc        output counting
cc        updating of fine to coarse grids
cc
cc  in this version, the integration strategy is geared to steady state
cc  problems - each level goes once, then the next finest goes, etc.
cc
cc  intcnt: counts the number of times a level is integrated.
cc  used to determine who should be integrated next. the finest
cc  level counts as 1 each step. the coarser level counts
cc  as intrat(level)**(mxnest-level). this is so we can do an integer
cc  compare to see who goes next.
cc
cc  icheck: counts the number of steps (incrementing
cc  by 1 each step) to keep track of
cc  when that level should have its error estimated
cc  and finer levels should be regridded.
cc
c      nstepc = nstart
c      stoptm = nstop*possk(1)
c      intinc(mxnest) = 1
c      do 5 i = 1, mxnest-1
c      ii = mxnest - i
c 5    intinc(ii) = intinc(ii+1)*intrat(ii)
cc
cc  ------ start of coarse grid integration loop. ------------------
cc
c 20   if (nstepc .ge. nstop) goto 999
c      if (nstepc .eq. nstop-1) then
c          lastout = .true.
c      else
c          lastout = .false.
c      endif
c      if ((iousr .ne. 0) .and. (nstepc .gt. nstart)) then
c        if ((nstepc/iousr)*iousr.eq.nstepc)call check(nstepc,time,lfix)
c      endif
cc
c          level  = 1
c          delt0  = possk(1)
c          dtnew  = 1.e32
c          do 10 i = 1, maxlv
c 10          intcnt(i) = 0
cc
cc     ------------- regridding  time?  ---------
cc
cc check if either
cc   (i)  this level should have its error estimated before being advanced
cc   (ii) this level needs to provide boundary values for either of
cc        next 2 finer levels to have their error estimated.
cc        this only affects two grid levels higher, occurs because
cc        previous time step needs boundary vals for giant step.
cc  no error estimation on finest possible grid level
cc
c 60       continue
c          if (icheck(level) .ge. kcheck) then
c               lbase = level
c          else if (level+1 .ge. mxnest) then
c               go to 90
c          else if (icheck(level+1) .ge. kcheck) then
c               lbase = level+1
c          else if (level+2 .ge. mxnest) then
c               go to 90
c          else if (icheck(level+2) .ge. kcheck) then
c               lbase = level+2
c          else
c               go to 90
c          endif
c          if (lbase .lt. lfix) then
c              lbasef = lfix
c           else
c              lbasef = lbase
c           endif
c          if (lbasef .eq. mxnest .or. lbasef .gt. lfine) go to 70
cc
cc regrid level 'lbasef+1' up to finest level.
cc level 'lbasef' stays fixed.
cc
c              if (gprint) write(6,101) lbasef
c101           format(8h  level ,i5,32h  stays fixed during regridding )
c              call regrid(nvar,lbasef,cut,cdist,work,steady,quad)
c             call valout(0,1,lfine,time,nvar)
cc
cc  maybe finest level in existence has changed. reset counters.
cc
c              if (gprint .and. lbase .lt. lfine)
c     1           call outtre(lstart(lbase+1),.false.,nvar)
c 70           intnow    = intcnt(lbase)
c              do 80  i  = lbase, lfine
c              intcnt(i) = intnow
c 80           icheck(i) = 0
cc
cc  ------- done regridding --------------------
cc
cc integrate all grids at level 'level'.
cc
c 90       continue
c          if (tprint) write(6,100) level
c100       format(29h  integrating grids at level ,i5)
c          call advanc(level,nvar,dtlev,vtime,steady,lastout)
c          dtnew = dmin1(dtnew,dtlev*float(intinc(1)/intinc(level)))
cc
cc done with a level of integration. update counts, decide who next.
cc
c          intcnt(level)  =  intcnt(level) +intinc(level)
c          icheck(level)  =  icheck(level) + 1
cc
c          if (level .lt. lfine) then
c             level = level + 1
c             go to 60
c          endif
cc
c 105      if (level .eq. 1) go to 110
cc#           if (intcnt(level) .lt. intcnt(level-1)) then
cc               same level goes again
cc#              go to 60
cc#           else
c                level = level - 1
c                call update(level,nvar,work)
cc#           endif
c          go to 105
cc
cc  --------------one complete coarse grid integration cycle done. -----
cc
cc      time for output?  done with the whole thing?
cc
c 110      continue
c          time    = time   + delt0
c          nstepc  = nstepc + 1
c          call conck(1,nvar,time)
c          if ((nstepc/iout)*iout.eq.nstepc .and. nstepc.lt.nstop) then
cc            call valout(0,1,lfine,time,nvar)
c             call dumptec(1, lfine,nvar,nplot,time)
c             if (.not. graf) call outtre(mstart,.true.,nvar)
c          endif
c          if (.not. vtime) go to 20
cc
cc  dtnew holds new delta t. passed back from integration routine.
cc
c         if (.not. steady) then   ! reset timestep if its meaningful
c           delt0    = dtnew
c           do 120 i = 1, mxnest
c 120         possk(i) = delt0 / float(intinc(1)/intinc(i))
c         endif
c         go to 20
cc
c999   continue
cc     if ((nstepc/iout)*iout .ne. nstepc) then
cc        call valout(0,1,lfine,time,nvar)
cc        if (.not. graf) call outtre(mstart, .true., nvar)
cc      endif
cc
c      return
c      end
c

c  -------------------------------------------------------------
c
      subroutine tick(nvar,iout,nstart,nstop,cut,cdist,
     1                vtime,work,time,iousr,steady,lfix,quad,nplot)
c
      implicit double precision (a-h,o-z)
      logical   tprint, vtime, gprint, steady, quad
      logical            graf,lastout
      parameter  (maxgr = 192, maxlv=12)
      dimension work(nvar), intinc(maxlv)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      common /runparam/ ftime,outtime,nstopc,iprint
      common/steadydata/steadydiff, isteadysim
      include "calloc.i"
      data      tprint /.true./,  gprint /.true./
      steadyTOL = 1e-10
c
c  parameters:
c     nstop   = # of coarse grid time steps to be taken
c     iout    = output interval every 'iout' coarse time steps
c     vtime   = true for variable timestep, calculated each coarse step
c
c  main driver routine,  controls:
c        integration  of all grids.
c        error estimation / regridding
c        output counting
c        updating of fine to coarse grids
c
c  in this version, the integration strategy is geared to steady state
c  problems - each level goes once, then the next finest goes, etc.
c
c  intcnt: counts the number of times a level is integrated.
c  used to determine who should be integrated next. the finest
c  level counts as 1 each step. the coarser level counts
c  as intrat(level)**(mxnest-level). this is so we can do an integer
c  compare to see who goes next.
c
c  icheck: counts the number of steps (incrementing
c  by 1 each step) to keep track of
c  when that level should have its error estimated
c  and finer levels should be regridded.
c
      nstepc = nstart
      stoptm = nstop*possk(1)
      intinc(mxnest) = 1
      do 5 i = 1, mxnest-1
      ii = mxnest - i
 5    intinc(ii) = intinc(ii+1)*intrat(ii)
c
c  ------ start of coarse grid integration loop. ------------------
c
c20   if (nstepc .ge. nstop) goto 999
 20   if ((nstop .eq. -1 .and. (time .ge. ftime)) .or.
     .    (nstop  >   -1 .and. (nstepc .ge. nstop)) .or.
     .    (nstop .eq. -2 .and. steadydiff < steadyTOL
     .                   .and. nstepc > 1) )  goto 999
      if (nstepc .eq. nstop-1) then
          lastout = .true.
      else
          lastout = .false.
      endif
      if ((iousr .ne. 0) .and. (nstepc .gt. nstart)) then
        if ((nstepc/iousr)*iousr.eq.nstepc)call check(nstepc,time,lfix)
      endif
c
          level  = 1
          delt0  = possk(1)
          dtnew  = 1.e32
          do 10 i = 1, maxlv
 10          intcnt(i) = 0
c
c     ------------- regridding  time?  ---------
c
c check if either
c   (i)  this level should have its error estimated before being advanced
c   (ii) this level needs to provide boundary values for either of
c        next 2 finer levels to have their error estimated.
c        this only affects two grid levels higher, occurs because
c        previous time step needs boundary vals for giant step.
c  no error estimation on finest possible grid level
c
 60       continue
          if (icheck(level) .ge. kcheck) then
               lbase = level
          else if (level+1 .ge. mxnest) then
               go to 90
          else if (icheck(level+1) .ge. kcheck) then
               lbase = level+1
          else if (level+2 .ge. mxnest) then
               go to 90
          else if (icheck(level+2) .ge. kcheck) then
               lbase = level+2
          else
               go to 90
          endif
          if (lbase .lt. lfix) then
              lbasef = lfix
           else
              lbasef = lbase
           endif
          if (lbasef .eq. mxnest .or. lbasef .gt. lfine) go to 70
c
c regrid level 'lbasef+1' up to finest level.
c level 'lbasef' stays fixed.
c
              if (gprint) write(6,101) lbasef
101           format(8h  level ,i5,32h  stays fixed during regridding )
              call regrid(nvar,lbasef,cut,cdist,work,steady,quad)
             call valout(0,1,lfine,time,nvar)
c
c  maybe finest level in existence has changed. reset counters.
c
              if (gprint .and. lbase .lt. lfine) 
     1           call outtre(lstart(lbase+1),.false.,nvar)
 70           intnow    = intcnt(lbase)
              do 80  i  = lbase, lfine
              intcnt(i) = intnow
 80           icheck(i) = 0
c
c  ------- done regridding --------------------
c
c integrate all grids at level 'level'.
c
 90       continue

          if( ((mod(time,outtime) + possk(1)) .gt. (outtime + 1d-10) )
     .       .and. (nstop .eq. -1) .and. (iprint .ne. nstepc-1)) then
             possk(1) = int((time+possk(1))/outtime)*outtime - time
             delt0 = possk(1)
             iprint = nstepc
          endif
          if (tprint) write(6,100) level
100       format(29h  integrating grids at level ,i5)
          call advanc(level,nvar,dtlev,vtime,steady,lastout)
          dtnew = dmin1(dtnew,dtlev*float(intinc(1)/intinc(level)))
c
c done with a level of integration. update counts, decide who next.
c
          intcnt(level)  =  intcnt(level) +intinc(level)
          icheck(level)  =  icheck(level) + 1
c
          if (level .lt. lfine) then
             level = level + 1
             go to 60
          endif
c
 105      if (level .eq. 1) go to 110
c#           if (intcnt(level) .lt. intcnt(level-1)) then
c               same level goes again
c#              go to 60
c#           else
                level = level - 1
                call update(level,nvar,work)
c#           endif
          go to 105
c
c  --------------one complete coarse grid integration cycle done. -----
c
c      time for output?  done with the whole thing?
c
 110      continue
          time    = time   + delt0
          nstepc  = nstepc + 1
!          call conck(1,nvar,time)
c          if(  (abs(mod(time,0.05d0)) .lt. 1.d-10)    ) then
c         if ((nstepc/iout)*iout.eq.nstepc .and. nstepc.lt.nstop) then
c            call valout(0,1,lfine,time,nvar)

      if(  (nstop .eq. -1 .and. abs(mod(time,outtime)) .lt. 1.d-10)
     . .or.(nstop > -1 .and.
     .     (nstepc/iout)*iout.eq.nstepc .and. nstepc.lt.nstop) ) then
             call dumptec(1, lfine,nvar,nplot,time)
             if (.not. graf) call outtre(mstart,.true.,nvar)
          endif
          if (.not. vtime) go to 20
c
c  dtnew holds new delta t. passed back from integration routine.
c
         if (.not. steady) then   ! reset timestep if its meaningful
           delt0    = dtnew
           do 120 i = 1, mxnest
 120         possk(i) = delt0 / float(intinc(1)/intinc(i))
         endif
         go to 20
c
999   continue
c     if ((nstepc/iout)*iout .ne. nstepc) then
c        call valout(0,1,lfine,time,nvar)
c        if (.not. graf) call outtre(mstart, .true., nvar)
c      endif
c
      call dumptec(1, lfine,nvar,nplot,time)
      return
      end




