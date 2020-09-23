      subroutine eval_dU_llf(qin,dU,f,g,irr,mitot,mjtot,lwidth,
     &               dtn,dtnewn,lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &                  vtime,steady,qx,qy,level,difmax,lastout,
     &                  meqn,time, ffluxlen, gfluxlen, istage)

      implicit double precision (a-h,o-z)
      parameter ( msize =5200 )
      dimension qin(mitot,mjtot,meqn)
      dimension q(mitot,mjtot,meqn)
      dimension f(mitot,mjtot,meqn),g(mitot,mjtot,meqn)
      dimension qx(mitot,mjtot,meqn),qy(mitot,mjtot,meqn)
      dimension qxx(mitot,mjtot,meqn),qyy(mitot,mjtot,meqn)
      dimension qxy(mitot,mjtot,meqn)
      dimension ur(msize,meqn),ul(msize,meqn)
      dimension res(mitot,mjtot,meqn)
      dimension ff(msize,meqn)
      dimension ffluxlen(msize,msize),gfluxlen(msize,msize)


      integer   irr(mitot,mjtot)
      include "cirr.i"
      common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                  ismp,gradThreshold

      logical    flag, debug, vtime, steady, quad, nolimiter
      logical    inside
      common /order2/ ssw, quad, nolimiter
      dimension  firreg(-1:irrsize,meqn)
      dimension  dU(mitot,mjtot,meqn),dU2(mitot,mjtot,meqn)
      dimension  dlimit(mitot,mjtot,meqn)
      dimension fakeState(4)

      integer    xrp, yrp
      data       debug/.false./
      data       xrp/1/, yrp/0/
      data       pi/3.14159265357989d0/

       dimension coeff(2)
       data mstage/2/
       data      coeff/0.5d0,1.d0/

      if ( msize .lt. max(mitot,mjtot) ) then
          write(6,*) 'Not enough memory allocated for rowwise flux'
          write(6,*) 'Calculations. Allocated Size ',msize
          write(6,*) 'Space required ',max(mitot,mjtot)
          write(6,*) 'Remember also to allocate the same additional'
          write(6,*) 'space in subroutine vrm '
          stop
       endif

      fakeState(1) = 1.d0
      fakeState(2) = 0.d0
      fakeState(3) = 0.d0
      fakeState(4) = 2.5d0

c
      ix1 = lwidth + 1
      ixn = mitot - lwidth
      iy1 = lwidth + 1
      iyn = mjtot - lwidth
      ar(-1) = 1.d0


c   :::::   rk with linear reconstruction follows ::::::
c
c  store primitive variables in f for now
c
      q = qin




      call vctoprm(q,q,mitot,mjtot)
c     ### call for exterior bcs at each stage so can use slopes
      xhigh= xlow + mitot*dx
      yhigh = ylow + mjtot*dy
      call pphysbdlin(xlow,xhigh,ylow,yhigh,level,mitot,mjtot,
     &                      meqn,q,time,dx,dy,qx,qy,irr,lstgrd)

c
c irregflux computes the cut cell bndry flux. since no flow
c  through bndry use eval pressure there.
      if (ssw .ne. 0.d0) then   ! recalc slopes at each stage
         call slopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlow,ylow,meqn)
         call qslopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &                 xlow,ylow,qxx,qxy,qyy,mptr,meqn)
      endif

      firreg = 0.d0
      call irregflux(q,firreg,irr,mitot,mjtot,dx,dy,lstgrd,
     .                   xlow,ylow,mptr,qx,qy,lwidth,msize)



      call vprmtoc(q,mitot,mjtot,mitot,mjtot,meqn)

      if (ssw .ne. 0.d0) then   ! recalc slopes at each stage
         call slopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlow,ylow,meqn)
         call qslopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &                 xlow,ylow,qxx,qxy,qyy,mptr,meqn)
      endif

c
c
c  loop through rows of q calculating fluxes one row at time
c  vertical riemann problem first
c
      do 800 jcol = lwidth-2, mjtot-lwidth+3
c
         do 511 i = lwidth-2, mitot-lwidth+3
            call getYface(i,jcol,xface,yface,irr,mitot,mjtot,
     .               xlow,ylow,dx,dy,lstgrd)
            call getCellCentroid(lstgrd,i,jcol+1,xcentp,ycentp,
     .                      xlow,ylow,dx,dy,irr(i,jcol+1))
            call getCellCentroid(lstgrd,i,jcol,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(i,jcol))
            do 512 m = 1, meqn
             if (gfluxlen(i,jcol+1) .ne. 0.d0) then  ! real face
                ur(i,m) = q(i,jcol+1,m) + (xface-xcentp)*qx(i,jcol+1,m)
     .                                  + (yface-ycentp)*qy(i,jcol+1,m)
                ul(i,m) = q(i,jcol,m)   + (xface-xcent)*qx(i,jcol,m)
     .                                  + (yface-ycent)*qy(i,jcol,m)
             else
               ur(i,m) = fakeState(m)
               ul(i,m) = fakeState(m)
             endif
  512    continue

c               check physical
                inp = 1
                if (ur(i,1) .le. 0.d0) inp = 0
                velsq = ur(i,2)**2+ur(i,3)**2
                press = .4d0*(ur(i,4)-0.5d0*velsq/ur(i,1))
                if (press .le. 0.d0) inp = 0

                if(inp .eq. 0) then
                    ur(i,:) = q(i,jcol+1,:)
                    qx(i,jcol+1,:) = 0.d0
                    qy(i,jcol+1,:) = 0.d0
                    print *, "oops nonphysical in LLF"
                endif

                inp = 1
                if (ul(i,1) .le. 0.d0) inp = 0
                velsq = ul(i,2)**2+ul(i,3)**2
                press = .4d0*(ul(i,4)-0.5d0*velsq/ul(i,1))
                if (press .le. 0.d0) inp = 0

                if(inp .eq. 0) then
                    ul(i,:) = q(i,jcol,:)
                    qx(i,jcol,:) = 0.d0
                    qy(i,jcol,:) = 0.d0
                    print *, "oops nonphysical in LLF"
                endif
  511    continue
c
c store fluxes in ff vector, copy into permanent flux array
c
         call vrm_llf(ur,ul,ff,lwidth-2,mitot-lwidth+3,yrp,msize)
c
         do 720 i = lwidth-2, mitot-lwidth+3
         do 720 m = 1, meqn

            g(i,jcol+1,m) = ff(i,m)
  720    continue
c
c

 800  continue
c
c
c    Horizontal riemann problems next
c
      do 900 irow = lwidth-2, mitot-lwidth+3
c
         xright = xlow + (dfloat(irow)-.5d0)*dx
         do 611 j = lwidth-2, mjtot-lwidth+3
            call getXface(irow,j,xface,yface,irr,mitot,mjtot,
     .               xlow,ylow,dx,dy,lstgrd)
            call getCellCentroid(lstgrd,irow+1,j,xcentp,ycentp,
     .                      xlow,ylow,dx,dy,irr(irow+1,j))
            call getCellCentroid(lstgrd,irow,j,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(irow,j))
            do 612 m = 1, meqn
             if (ffluxlen(irow+1,j) .ne. 0.) then ! real face
                ur(j,m) = q(irow+1,j,m) + (xface-xcentp)*qx(irow+1,j,m)
     .                                  + (yface-ycentp)*qy(irow+1,j,m)
                ul(j,m) = q(irow,j,m)   + (xface-xcent)*qx(irow,j,m)
     .                                  + (yface-ycent)*qy(irow,j,m)
             else
c              ur(j,m) = q(irow+1,j,m)
c              ul(j,m) = q(irow,j,m)
               ur(j,m) = fakeState(m)
               ul(j,m) = fakeState(m)
             endif
  612        continue
c            check physical
             inp = 1
             if (ur(j,1) .le. 0.d0) inp = 0
             velsq = ur(j,2)**2+ur(j,3)**2
             press = .4d0*(ur(j,4)-0.5d0*velsq/ur(j,1))
             if (press .le. 0.d0) inp = 0

             if(inp .eq. 0) then
                print *, "oops nonphysical in LLF"
                ur(j,:) = q(irow+1,j,:)
                qx(irow+1,j,:) = 0.d0
                qy(irow+1,j,:) = 0.d0
             endif

             inp = 1
             if (ul(j,1) .le. 0.d0) inp = 0
             velsq = ul(j,2)**2+ul(j,3)**2
             press = .4d0*(ul(j,4)-0.5d0*velsq/ul(j,1))
             if (press .le. 0.d0) inp = 0

             if(inp .eq. 0) then
                 print *, "oops nonphysical in LLF"
                 ul(j,:) = q(irow,j,:)
                 qx(irow,j,:) = 0.d0
                 qy(irow,j,:) = 0.d0
             endif
  611   continue

c
c store fluxes in ff
c
         call vrm_llf(ur,ul,ff,lwidth-2,mjtot-lwidth+3,xrp,msize)
c
         do 721  m = 1, meqn
         do 721 j = lwidth-2, mjtot-lwidth+3
            f(irow+1,j,m) = ff(j,m)
  721    continue
c
 900  continue



c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  multiply fluxes by mesh spacing.
c  this zeros out unused fluxes so solid vals dont get updated.
c
         do 580 m = 1, meqn
         do 580 i = 2, mitot
         do 580 j = 2, mjtot
            f(i,j,m) = f(i,j,m) * ffluxlen(i,j)
            g(i,j,m) = g(i,j,m) * gfluxlen(i,j)
 580     continue

c
c
c      # finite volume update
         c      = coeff(istage)
         ar(-1) = 1.d0   ! prevent zero divides for solid cells
         do 917 i = ix1-lwidth+istage, ixn+lwidth-istage
         do 917 j = iy1-lwidth+istage, iyn+lwidth-istage
         do 917 m = 1, meqn
            k = irr(i,j)
            dU(i,j,m) = - dtn/ar(k)*( f(i+1,j,m) -f(i,j,m)
     .      + g(i,j+1,m) -g(i,j,m)  - firreg(k,m))
 917     continue




      return
      end
