      subroutine method(q,f,g,irr,mitot,mjtot,lwidth,dtn,dtnewn,
     &                  lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &                  vtime,steady,qx,qy,level,difmax,lastout,
     &                  meqn,time)

      implicit double precision (a-h,o-z)
      parameter ( msize =5200 )
      dimension q(mitot,mjtot,meqn),qtemp(mitot,mjtot,meqn)
      dimension q1(mitot,mjtot,meqn),q2(mitot,mjtot,meqn)
      dimension q3(mitot,mjtot,meqn), qbefore(mitot,mjtot,meqn)
      dimension sk1(mitot,mjtot,meqn),sk2(mitot,mjtot,meqn)
      dimension sk3(mitot,mjtot,meqn),sk4(mitot,mjtot,meqn)
      dimension f(mitot,mjtot,meqn),g(mitot,mjtot,meqn)
      dimension qx(mitot,mjtot,meqn),qy(mitot,mjtot,meqn)
      dimension qxx(mitot,mjtot,meqn),qyy(mitot,mjtot,meqn)
      dimension qxy(mitot,mjtot,meqn)
      dimension ur(msize,meqn),ul(msize,meqn)
!      dimension res(mitot,mjtot,meqn)
      dimension ff(msize,meqn)
      dimension ffluxlen(msize,msize),gfluxlen(msize,msize)
      dimension recon(meqn)

      integer   irr(mitot,mjtot)
      include "cirr.i"
      common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .        ismp,gradThreshold, ilts,ibc,id,id2,res(ijsize,ijsize,4)


      logical    flag, debug, vtime, steady, quad, nolimiter
      logical    inside
      common /order2/ ssw, quad, nolimiter
      common/steadydata/steadydiff, isteadysim
      common /timedata/ stage_time
      dimension  firreg(-1:irrsize,meqn)
      dimension  dU(mitot,mjtot,meqn),dU2(mitot,mjtot,meqn)
      dimension  dlimit(mitot,mjtot,meqn)
      dimension  diff(mitot,mjtot,meqn)

      logical IS_GHOST

      integer    xrp, yrp
      data       debug/.false./
      data       xrp/1/, yrp/0/
      data       pi/3.14159265357989d0/

      data mstage/3/
      common/eqncount/ieqncount
!      integer imethod/2/
      integer imethod
      ! 1 forward euler
      ! 2 RK2 SSP
      ! 3 RK3-SSP
      ! 4 RK4

      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)

      if(0 .le. ssw .and. ssw .le. 3) then
      imethod = ssw + 1 ! use the same order time stepper as the accuracy..
      elseif(ssw .eq. -1) then ! use RK2-SSP
      imethod = 2
      elseif(ssw .eq. -10) then ! use RK2-SSP
      imethod = 2
      elseif(ssw .eq. -2) then ! use RK3-SSP
      imethod = 3
      elseif(ssw .eq. -3) then ! use RK3-SSP
      imethod = 4
      else
      print *, "UH OH PROBLEM WITH TIME INTEGRATOR CHOICE."
      endif

      call outgeom(1,mitot, mjtot, irr,xlow, ylow, dx, dy)
!      if(isteadysim .eq. 1) imethod = 1
!      imethod = 1

      stage_time = time


      if ( msize .lt. max(mitot,mjtot) ) then
          write(6,*) 'Not enough memory allocated for rowwise flux'
          write(6,*) 'Calculations. Allocated Size ',msize
          write(6,*) 'Space required ',max(mitot,mjtot)
          write(6,*) 'Remember also to allocate the same additional'
          write(6,*) 'space in subroutine vrm '
          stop
       endif



c
      ix1 = lwidth + 1
      ixn = mitot - lwidth
      iy1 = lwidth + 1
      iyn = mjtot - lwidth
      ar(-1) = 1.d0
      ar(lstgrd) = dx*dy
c
      if (.true.) then
         totmass1 = bigconck(q,irr,mitot,mjtot,lwidth,meqn)
         write(*,909) totmass1
 909     format("mass at beginning of RK step: ",e30.16)
      endif

c
c     # initialize fluxes:
c
      firreg = 0.d0
      f = 0.d0
      g = 0.d0
      qx = 0.d0
      qy = 0.d0
      qxx = 0.d0
      qxy = 0.d0
      qyy = 0.d0
      qxxx = 0.d0
      qxxy = 0.d0
      qxyy = 0.d0
      qyyy = 0.d0

c need routine to set face lengths and  midpoints
      call getirrlen(irr,mitot,mjtot,dtn,dx,dy,lstgrd,
     &               mptr,meqn,ffluxlen,gfluxlen,msize)
c
c  could turn on for debugging, or move to regridding section of code
c          call cellsClose(ffluxlen,gfluxlen,mitot,mjtot,irr,lstgrd,
c     .                    lwidth,msize)
c
c
      if(imethod .ne. 4) then

!      dtn = 0.07d0
!      dtnewn = 0.07d0

      istage = 1

c     STAGE 1
      istage = 1
      call eval_dU_gauss(q,dU,f,g,irr,mitot,mjtot,lwidth,dtn,dtnewn,
     &             lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &             vtime,steady,qx,qy,level,difmax,lastout,
     &             meqn,time, ffluxlen, gfluxlen, istage)
      q1 = q+dU;

      if(dtn>0) then
      do i = 1,mitot
      do j = 1,mjtot
      kirr = irr(i,j)
      res(i,j,:) = dU(i,j,:)*ar(kirr)/dtn
      enddo
      enddo
      else
      res = 0.d0
      endif



      if(ismp .eq. 1) then
      if(ihob .eq. 0) then
      call SRD_cellMerge_U(q1,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      else
      call SRD_cellMerge_ho(q1,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      endif
      endif





      if(imethod .eq. 1) then
      istage = mstage

      diff = q1-q
      steadydiff = 0.d0
      do 2 i=lwidth+1,mitot-lwidth
      do 2 j=lwidth+1,mjtot-lwidth
      if(irr(i,j) .eq. -1 .or. IS_GHOST(i,j)) cycle

      do 2 m = 1,ieqncount
      steadydiff = max(steadydiff,diff(i,j,m))
 2    continue

      q = q1

      goto 999
      endif

c     STAGE 2
      istage = 2
      call eval_dU_gauss(q1,dU,f,g,irr,mitot,mjtot,lwidth,dtn,dtnewn,
     &             lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &             vtime,steady,qx,qy,level,difmax,lastout,
     &             meqn,time+dtn, ffluxlen, gfluxlen, istage)
      q2 = q1+dU

      if(ismp .eq. 1) then
      if(ihob .eq. 0) then
      call SRD_cellMerge_U(q2,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      else
      call SRD_cellMerge_ho(q2,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      endif
      endif

      if(imethod .eq. 2) then
      istage = mstage

c     FINAL UPDATE
      diff = 0.5d0*(q+q2) - q

      steadydiff = 0.d0
      do 11 i=lwidth+1,mitot-lwidth
      do 11 j=lwidth+1,mjtot-lwidth
      if(irr(i,j) .eq. -1 .or. IS_GHOST(i,j)) cycle

      do 11 m = 1,ieqncount
      steadydiff = max(steadydiff,diff(i,j,m))
 11   continue

!      steadydiff = maxval( abs( diff(lwidth+1:mitot-lwidth,
!     .                              lwidth+1:mjtot-lwidth,:) ) )


      q = 0.5d0*(q+q2)
      goto 999
      endif

c     STAGE 3
      istage = 3

      qtemp = 0.75d0 * q + 0.25d0 * q2
      call eval_dU_gauss(qtemp,dU,f,g,irr,mitot,mjtot,lwidth,dtn,dtnewn,
     &             lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &             vtime,steady,qx,qy,level,difmax,lastout,
     &             meqn,time+dtn/2.d0, ffluxlen, gfluxlen, istage)
      q3 = qtemp+dU

      if(ismp .eq. 1) then
      if(ihob .eq. 0) then
      call SRD_cellMerge_U(q3,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      else
      call SRD_cellMerge_ho(q3,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      endif
      endif

c     FINAL UPDATE
      diff = (1.d0 / 3.d0) * q + (2.d0 / 3.d0) * q3 - q
      steadydiff = 0.d0
      do 111 i=lwidth+1,mitot-lwidth
      do 111 j=lwidth+1,mjtot-lwidth
      if(irr(i,j) .eq. -1 .or. IS_GHOST(i,j)) cycle

      do 111 m = 1,ieqncount
      steadydiff = max(steadydiff,diff(i,j,m))
 111   continue

      q = (1.d0 / 3.d0) * q + (2.d0 / 3.d0) * q3

      else ! RK4
      qbefore = q

      istage = 1
      call eval_dU_gauss(q,sk1,f,g,irr,mitot,mjtot,lwidth,dtn,dtnewn,
     &             lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &             vtime,steady,qx,qy,level,difmax,lastout,
     &             meqn,time, ffluxlen, gfluxlen, istage)

      q1 = q + 0.5d0 * sk1

      if(dtn>0) then
      do i = 1,mitot
      do j = 1,mjtot
      kirr = irr(i,j)
      res(i,j,:) = sk1(i,j,:)*ar(kirr)/dtn
      enddo
      enddo

      else
      res = 0.d0
      endif

      if(ismp .eq. 1) then
      if(ihob .eq. 0) then
      call SRD_cellMerge_U(q1,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      else
      call SRD_cellMerge_ho(q1,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      endif
      endif


!      istage = mstage
!c     FINAL UPDATE
!      q = q1
!      goto 999


      istage = 2
      call eval_dU_gauss(q1,sk2,f,g,irr,mitot,mjtot,lwidth,dtn,dtnewn,
     &             lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &             vtime,steady,qx,qy,level,difmax,lastout,
     &             meqn,time+dtn/2.d0, ffluxlen, gfluxlen, istage)

      q2 = q + 0.5d0 * sk2

      if(ismp .eq. 1) then
      if(ihob .eq. 0) then
      call SRD_cellMerge_U(q2,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      else
      call SRD_cellMerge_ho(q2,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      endif
      endif

      istage = 3
      call eval_dU_gauss(q2,sk3,f,g,irr,mitot,mjtot,lwidth,dtn,dtnewn,
     &             lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &             vtime,steady,qx,qy,level,difmax,lastout,
     &             meqn,time+dtn/2.d0, ffluxlen, gfluxlen, istage)
      q3 = q + sk3

      if(ismp .eq. 1) then
      if(ihob .eq. 0) then
      call SRD_cellMerge_U(q3,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      else
      call SRD_cellMerge_ho(q3,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      endif
      endif

      istage = 4
      call eval_dU_gauss(q3,sk4,f,g,irr,mitot,mjtot,lwidth,dtn,dtnewn,
     &             lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &             vtime,steady,qx,qy,level,difmax,lastout,
     &             meqn,time+dtn, ffluxlen, gfluxlen, istage)

      q = q + (1.d0/6.d0)*(sk1+2.d0*sk2+2.d0*sk3+sk4)

      if(ismp .eq. 1) then
      if(ihob .eq. 0) then
      call SRD_cellMerge_U(q,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      else
      call SRD_cellMerge_ho(q,dlimit,meqn,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
      endif
      endif

      diff = q - qbefore
      steadydiff = 0.d0
      do 1111 i=lwidth+1,mitot-lwidth
      do 1111 j=lwidth+1,mjtot-lwidth
      if(irr(i,j) .eq. -1 .or. IS_GHOST(i,j)) cycle

      do 1111 m = 1,ieqncount
      steadydiff = max(steadydiff,diff(i,j,m))
1111   continue


      endif




999   continue
c
c
c     :: for multistage method need to adjust for fluxes at each stage
c     :: note that not conservative if using amr unless sum fluxes at each stage around boundary

c
c      :: output adjusted mass for multistage scheme
c
       if (.true.) then
          totmass2 = bigconck(q,irr,mitot,mjtot,lwidth,meqn)  ! final mass on grid
          write(*,919) totmass2
 919      format("mass at end of SSP step: ",e30.16)
          print *, "steady difference is ", steadydiff
c         time step adjustment factor depends on stage.
c         continue to accumulate for 2 stage method and change factor
c         :::: next formula ONLY works for 1 and 2 stage methods
       endif





c     # shift fluxes before returning: (due to mismatch with amr)
c
cdir$ ivdep
         do 920 m = 1, meqn
         do 920 j = iy1, iyn
         do 920 i = ix1, ixn+1
                 f(i-lwidth,j-lwidth,m) = f(i,j,m)
 920    continue

cdir$ ivdep
        do 930 m = 1, meqn
        do 930 j = iy1, iyn+1
        do 930 i = ix1, ixn
                 g(i-lwidth,j-lwidth,m) = g(i,j,m)
 930    continue
c
c     estimate new time step for next round. even if not used will give cfl
      if (vtime) then
         arreg = dx*dy  ! get regular cell info
         rlen = dsqrt(arreg)
         dt3 = 1.d10  ! initialize
         do 140 j = lwidth+1, mjtot-lwidth
         do 140 i = lwidth+1, mitot-lwidth
            k = irr(i,j)
            if (k .ne. -1) then
                sp = speed(q, i, j,mitot,mjtot,meqn,irr,lstgrd,
     .               xlow,ylow,dx,dy)

                call getbothspeeds(q, i, j, spx, spy,
     .               mitot,mjtot,meqn,irr,lstgrd,
     .               xlow,ylow,dx,dy)
                if(imethod .eq. 1 .or.
     .             imethod .eq. 2 .or.
     .             imethod .eq. 3 .or.
     .             imethod .eq. 4) then ! dt *(spx/dx + spy/dy) <= 1
                    delT = 1/((spx / dx) + (spy / dy))
                else
                     delT = dx*dy / (2*dx + 2*dy) / sp
                endif

                dt3 =  min(cfl* delT,dt3)
            endif
 140     continue
         dtnewn = dt3     ! output var, need to return someting


!      dtn    = 0.147114317E-01
!      dtnewn = 0.147114317E-01



         print *, "dt is " , dtnewn
      endif

      return
      end


!                if(iprob .ne. 25) then
!                p = gamma1* (q(i,j,4)- .5d0* (q(i,j,2)*q(i,j,2)/
!     &              q(i,j,1) + q(i,j,3)* q(i,j,3)/q(i,j,1)))
!                c2 = gamma*p/q(i,j,1)
!                if (c2 .le. 0.d0) then
!                    print *, "non physical."
!                    go to 140
!                endif
!                c = dsqrt(c2)
!                u = q(i,j,2)/q(i,j,1)
!                v = q(i,j,3)/q(i,j,1)
!c               ::: dt using muscl cfl limit
!                 speed = dsqrt(u*u+v*v) + c  !use max norm for muscl
!                 delT = dx*dy / (2*dx + 2*dy) / speed
!                 dt3 =  min(cfl* delT,dt3)
!                 else
!                 speed = 1.d0
!                 delT = dx*dy / (2*dx + 2*dy) / speed
!                 dt3 =  min(cfl* delT,dt3)
!                 endif




