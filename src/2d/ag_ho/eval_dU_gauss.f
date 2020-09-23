      subroutine eval_dU_gauss(qin,dU,f_in,g_in,irr,mitot,mjtot,lwidth,
     &               dtn,dtnewn,lstgrd,dx,dy,flag,iorder,xlow,ylow,mptr,
     &                  vtime,steady,qx,qy,level,difmax,lastout,
     &                  meqn,time, ffluxlen, gfluxlen, istage)

      implicit double precision (a-h,o-z)
      parameter ( msize =5200 )
      dimension qin(mitot,mjtot,meqn)
      dimension q(mitot,mjtot,meqn)
      dimension f_in(mitot,mjtot,meqn),g_in(mitot,mjtot,meqn)
      dimension qx(mitot,mjtot,meqn),qy(mitot,mjtot,meqn)
      dimension qxx(mitot,mjtot,meqn),qyy(mitot,mjtot,meqn)
      dimension qxy(mitot,mjtot,meqn)
      dimension qxxx(mitot,mjtot,meqn),qxyy(mitot,mjtot,meqn)
      dimension qxxy(mitot,mjtot,meqn),qyyy(mitot,mjtot,meqn)
      dimension ur(msize,meqn),ul(msize,meqn)
      dimension res(mitot,mjtot,meqn)
      dimension ff(msize,meqn)
      dimension ffluxlen(msize,msize),gfluxlen(msize,msize)
      dimension xp(msize), yp(msize)
      dimension Uout(meqn)

      integer   irr(mitot,mjtot)
      logical IS_GHOST

      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)

      include "cirr.i"
      include "./quadrature.i"

      common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                  ismp,gradThreshold,ilts,ibc,imodal

      logical    flag, debug, vtime, steady, quad, nolimiter
      logical    inside, missing
      common /order2/ ssw, quad, nolimiter
      common/eqncount/ieqncount


      dimension  firreg(-1:irrsize,meqn)
      dimension  dU(mitot,mjtot,meqn),dU2(mitot,mjtot,meqn)
      dimension  dlimit(mitot,mjtot,meqn)
      dimension fakeState(4)

      integer    xrp, yrp
      data       debug/.false./
      data       xrp/1/, yrp/0/
      data       pi/3.14159265357989d0/


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

      ur(:,1) = fakestate(1)
      ur(:,2) = fakestate(2)
      ur(:,3) = fakestate(3)
      ur(:,4) = fakestate(4)

      ul(:,1) = fakestate(1)
      ul(:,2) = fakestate(2)
      ul(:,3) = fakestate(3)
      ul(:,4) = fakestate(4)



c
      ix1 = lwidth + 1
      ixn = mitot - lwidth
      iy1 = lwidth + 1
      iyn = mjtot - lwidth
      ar(-1) = 1.d0


c   :::::   rk with linear reconstruction follows ::::::
c
c  store primitive variables in f_in for now
c
      q = qin
      f_in = 0.d0
      g_in = 0.d0
      qx   = 0.d0
      qy   = 0.d0
      qxx  = 0.d0
      qxy  = 0.d0
      qyy  = 0.d0
      qxxx = 0.d0
      qxxy = 0.d0
      qxyy = 0.d0
      qyyy = 0.d0


      ! only in primitive variables for p = 1



c     ### call for exterior bcs at each stage so can use slopes
      xhigh= xlow + mitot*dx
      yhigh = ylow + mjtot*dy
      call pphysbdlin(xlow,xhigh,ylow,yhigh,level,mitot,mjtot,
     &                      meqn,q,time,dx,dy,qx,qy,irr,lstgrd)




      if(quad .and. (nolimiter .eqv. .true.) ) then ! yes in primitive variables, and without limiter (this is an accurate projection between the two
      call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlow,ylow,meqn)

      call vctoprm(mitot,mjtot,meqn,
     . dx, dy, lstgrd, xlow, ylow, irr,
     . q,q,qx,qy,qxx, qxy, qyy, qxxx, qxxy, qxyy, qyyy)

      call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlow,ylow,meqn)




      elseif(quad .and. (nolimiter .eqv. .false.) ) then ! yes in primitive variables, and with limiter
                                                      ! only works for p = 1

      call vctoprm(mitot,mjtot,meqn,
     . dx, dy, lstgrd, xlow, ylow, irr,
     . q,q,qx,qy,qxx, qxy, qyy, qxxx, qxxy, qxyy, qyyy)


      call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlow,ylow,meqn)

      ! limiter is only applied to the primitive variables.
!      if(imodal .eq. 0) then
            call limiter(q,  qx, qy,  mitot, mjtot, irr, meqn,dx, dy,
     .               lwidth, lstgrd, xlow, ylow)
!      else
!            call limiter_modal(q,  qx, qy,  qxx,qxy,qyy,
!     .             qxxx, qxxy, qxyy, qyyy,
!     .             mitot, mjtot, irr, meqn,dx, dy,
!     .             lwidth, lstgrd, xlow, ylow)
!      endif



      elseif (ssw .ne. 0.d0) then   ! conserved variables
         call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlow,ylow,meqn)
      endif







c
c
c  loop through rows of q calculating fluxes one row at time
c  vertical riemann problem first
c
      do 33 nn = 1,nlinequad
      do 800 jcol = lwidth-2, mjtot-lwidth+3
c
         do 511 i = lwidth-2, mitot-lwidth+3
            call getYface_gauss(i,jcol,xface,yface,irr,mitot,mjtot,
     .               xlow,ylow,dx,dy,lstgrd, nn, missing)
            if(missing) cycle

            call getCellCentroid(lstgrd,i,jcol+1,xcentp,ycentp,
     .                      xlow,ylow,dx,dy,irr(i,jcol+1))
            call getCellCentroid(lstgrd,i,jcol,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(i,jcol))

            kp = irr(i,jcol+1)
            k = irr(i,jcol)

            xp(i) = xface
            yp(i) = yface

            dxp = xface-xcentp
            dyp = yface-ycentp

            dxm = xface-xcent
            dym = yface-ycent

            do 512 m = 1, meqn
             if (gfluxlen(i,jcol+1) .ne. 0.d0) then  ! real face

       call evalU(Uout,dxp, dyp, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,jcol+1, kp, mitot, mjtot, meqn, m)
       ur(i,m) = Uout(m)

       if(ibc .eq. 1 .and. IS_GHOST(i,jcol+1)) then ! replace right state with exact
            call f(Uout, xface, yface, time, iprob)
            if(quad) then ! convert to primitive variables
            rho = Uout(1)
            u   = Uout(2)/Uout(1)
            v   = Uout(3)/Uout(1)
            p   = (gamma - 1.d0) * (Uout(4)
     .               -0.5d0*(Uout(2)**2 + Uout(3)**2)/ Uout(1) )
            Uout(1) = rho
            Uout(2) = u
            Uout(3) = v
            Uout(4) = p
            endif

            ur(i,m) = Uout(m)
       endif


       call evalU(Uout,dxm, dym, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,jcol, k, mitot, mjtot, meqn, m)
       ul(i,m) = Uout(m)

       if(ibc .eq. 1 .and. IS_GHOST(i,jcol)) then ! replace left state with exact
            call f(Uout, xface, yface, time, iprob)
            if(quad) then ! convert to primitive variables
            rho = Uout(1)
            u   = Uout(2)/Uout(1)
            v   = Uout(3)/Uout(1)
            p   = (gamma - 1.d0) * (Uout(4)
     .               -0.5d0*(Uout(2)**2 + Uout(3)**2)/ Uout(1) )
            Uout(1) = rho
            Uout(2) = u
            Uout(3) = v
            Uout(4) = p
            endif

            ul(i,m) = Uout(m)
!             ul(i,m) = ur(i,m) ! try this ...
       endif



             else
               ur(i,m) = fakeState(m)
               ul(i,m) = fakeState(m)
             endif
  512    continue

!             if( ul(i,4) < 0.d0 .or. ur(i,4) < 0.d0) then
!             print *,"here"
!             endif

  511    continue


c
c store fluxes in ff vector, copy into permanent flux array
c
         call vrm(ur,ul,ff,lwidth-2,mitot-lwidth+3,yrp,msize,xp,yp)
c
         do 720 i = lwidth-2, mitot-lwidth+3
         do 720 m = 1, meqn

            g_in(i,jcol+1,m) = g_in(i,jcol+1,m) + wline(nn)*ff(i,m)
  720    continue

c
c

 800  continue
  33  continue


c
c
c    Horizontal riemann problems next
c
      do 34 nn = 1,nlinequad
      do 900 irow = lwidth-2, mitot-lwidth+3
c
         xright = xlow + (dfloat(irow)-.5d0)*dx
         do 611 j = lwidth-2, mjtot-lwidth+3
            call getXface_gauss(irow,j,xface,yface,irr,mitot,mjtot,
     .               xlow,ylow,dx,dy,lstgrd,nn, missing)
            if(missing) cycle

            call getCellCentroid(lstgrd,irow+1,j,xcentp,ycentp,
     .                      xlow,ylow,dx,dy,irr(irow+1,j))
            call getCellCentroid(lstgrd,irow,j,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(irow,j))
            kp = irr(irow+1,j)
            k = irr(irow,j)

!            if( irow .eq. 5 .and. j .eq. 30) then
!            print *, here
!            endif


            xp(j) = xface
            yp(j) = yface

            dxp = xface-xcentp
            dyp = yface-ycentp

            dxm = xface-xcent
            dym = yface-ycent

            do 611 m = 1, ieqncount
             if (ffluxlen(irow+1,j) .ne. 0.) then ! real face

       call evalU(Uout,dxp, dyp, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            irow+1,j, kp, mitot, mjtot, meqn, m)
       ur(j,m) = Uout(m)

       if(ibc .eq. 1 .and. IS_GHOST(irow+1,j)) then ! replace right state with exact
            call f(Uout, xface, yface, time, iprob)

            if(quad) then ! convert to primitive variables
            rho = Uout(1)
            u   = Uout(2)/Uout(1)
            v   = Uout(3)/Uout(1)
            p   = (gamma - 1.d0) * (Uout(4)
     .               -0.5d0*(Uout(2)**2 + Uout(3)**2)/ Uout(1) )
            Uout(1) = rho
            Uout(2) = u
            Uout(3) = v
            Uout(4) = p
            endif

            ur(j,m) = Uout(m)
       endif

       call evalU(Uout,dxm, dym, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            irow,j, k, mitot, mjtot, meqn, m)
       ul(j,m) = Uout(m)

       if(ibc .eq. 1 .and. IS_GHOST(irow,j)) then ! replace left state with exact
            call f(Uout, xface, yface, time, iprob)
            if(quad) then ! convert to primitive variables
            rho = Uout(1)
            u   = Uout(2)/Uout(1)
            v   = Uout(3)/Uout(1)
            p   = (gamma - 1.d0) * (Uout(4)
     .               -0.5d0*(Uout(2)**2 + Uout(3)**2)/ Uout(1) )
            Uout(1) = rho
            Uout(2) = u
            Uout(3) = v
            Uout(4) = p
            endif

            ul(j,m) = Uout(m)
       endif




             else
               ur(j,m) = fakeState(m)
               ul(j,m) = fakeState(m)
             endif


  611   continue

c
c store fluxes in ff
c



         call vrm(ur,ul,ff,lwidth-2,mjtot-lwidth+3,xrp,msize, xp,yp)


         do 721  m = 1, ieqncount
         do 721 j = lwidth-2, mjtot-lwidth+3
            f_in(irow+1,j,m) = f_in(irow+1,j,m)+wline(nn)*ff(j,m)
  721    continue
c
c        if(irow+1 .eq. 54) print *, "f_in is ",ff(5,1), ffluxlen(54,5),
c     .  ffluxlen(54,6),ffluxlen(54,7)
 900  continue
 34   continue

c
c irregflux computes the cut cell bndry flux. since no flow
c  through bndry use eval pressure there.
          firreg = 0.d0
          call irregflux_gauss(q,firreg,irr,mitot,mjtot,dx,dy,lstgrd,
     .                   xlow,ylow,mptr,qx,qy,qxx,qxy,qyy,
     .                   qxxx, qxxy, qxyy, qyyy,
     .                   lwidth,msize,
     .                   time)
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  multiply fluxes by mesh spacing.
c  this zeros out unused fluxes so solid vals dont get updated.
c
         do 580 m = 1, ieqncount
         do 580 i = 2, mitot
         do 580 j = 2, mjtot
            f_in(i,j,m) = f_in(i,j,m) * ffluxlen(i,j)
            g_in(i,j,m) = g_in(i,j,m) * gfluxlen(i,j)
 580     continue

c
c
c      # finite volume update
         ar(-1) = 1.d0   ! prevent zero divides for solid cells
!         do 917 i = ix1-lwidth+istage, ixn+lwidth-istage
!         do 917 j = iy1-lwidth+istage, iyn+lwidth-istage
         do 917 i = ix1, ixn
         do 917 j = iy1, iyn
            k = irr(i,j)
            if(k .eq. -1) cycle

         do 917 m = 1, ieqncount
            if(ihob .eq. 0) then

            dU(i,j,m) = - dtn/ar(k)*( f_in(i+1,j,m) -f_in(i,j,m)
     .      + g_in(i,j+1,m) -g_in(i,j,m)  - firreg(k,m))

            if(ilts .eq. 1) then
            dU(i,j,m) = - dtn/(dx*dy)*( f_in(i+1,j,m) -f_in(i,j,m)
     .      + g_in(i,j+1,m) -g_in(i,j,m)  - firreg(k,m))
            endif

            else

            dU(i,j,m) = - dtn/ar_ho(k)*( f_in(i+1,j,m) -f_in(i,j,m)
     .      + g_in(i,j+1,m) -g_in(i,j,m)  - firreg(k,m))

            if(ilts .eq. 1) then
            dU(i,j,m) = - dtn/(dx*dy)*( f_in(i+1,j,m) -f_in(i,j,m)
     .      + g_in(i,j+1,m) -g_in(i,j,m)  - firreg(k,m))
            endif


            endif

 917     continue




      return
      end






      subroutine evalU(Uout,deltax, deltay, dx, dy,
     .                 q, qx, qy,
     .                 qxx, qxy, qyy,
     .                 qxxx, qxxy, qxyy, qyyy,
     .                 i, j, kirr, mitot, mjtot, nvar, m)
      implicit double precision (a-h,o-z)
      include "cirr.i"
      common /order2/ ssw, quad, nolimiter
      dimension Uout(nvar)
      dimension q(mitot,mjtot, nvar),qx(mitot,mjtot, nvar)
      dimension qy(mitot,mjtot, nvar)
      dimension qxx(mitot,mjtot, nvar),qxy(mitot,mjtot, nvar)
      dimension qyy(mitot,mjtot, nvar)
      dimension qxxx(mitot,mjtot, nvar),qxxy(mitot,mjtot, nvar)
      dimension qxyy(mitot,mjtot, nvar)
      dimension qyyy(mitot,mjtot, nvar)
      logical quad

      if(ssw .eq. -10) then
      Uout(m) = q(i,j,m) + deltax*qx(i,j,m)+ deltay*qy(i,j,m)
      return
      endif

      select case(ihob)
      case(0)
      Uout(m) = q(i,j,m) + deltax*qx(i,j,m)+ deltay*qy(i,j,m)
     .+ 0.5d0*qxx(i,j,m)*( deltax**2-(dx**2)*poly(8,1,kirr)    )
     .+       qxy(i,j,m)*( deltax*deltay-dx*dy*poly(9,1,kirr)  )
     .+ 0.5d0*qyy(i,j,m)*( deltay**2-(dy**2)*poly(10,1,kirr)   )
     .+(1.d0/6.d0)*qxxx(i,j,m)*( deltax**3-(dx**3)*dcubicshifts(1,kirr))
     .+(1.d0/2.d0)*qxxy(i,j,m)*( deltay*deltax**2
     .                                 -(dx**2)*dy*dcubicshifts(2,kirr))
     .+(1.d0/2.d0)*qxyy(i,j,m)*( deltax*deltay**2
     .                                 -(dy**2)*dx*dcubicshifts(3,kirr))
     .+(1.d0/6.d0)*qyyy(i,j,m)*( deltay**3-(dy**3)*dcubicshifts(4,kirr))

      case(1)
      Uout(m) = q(i,j,m) + deltax*qx(i,j,m)+ deltay*qy(i,j,m)
     .+ 0.5d0*qxx(i,j,m)*( deltax**2-(dx**2)*poly_ho(8,1,kirr)    )
     .+       qxy(i,j,m)*( deltax*deltay-dx*dy*poly_ho(9,1,kirr)  )
     .+ 0.5d0*qyy(i,j,m)*( deltay**2-(dy**2)*poly_ho(10,1,kirr)   )
     .+(1.d0/6.d0)*qxxx(i,j,m)*(deltax**3
     .-(dx**3)*dcubicshifts_ho(1,kirr))
     .+(1.d0/2.d0)*qxxy(i,j,m)*( deltay*deltax**2
     .                      -(dx**2)*dy*dcubicshifts_ho(2,kirr))
     .+(1.d0/2.d0)*qxyy(i,j,m)*( deltax*deltay**2
     .                       -(dy**2)*dx*dcubicshifts_ho(3,kirr))
     .+(1.d0/6.d0)*qyyy(i,j,m)*( deltay**3
     .                       -(dy**3)*dcubicshifts_ho(4,kirr))

      end select

!        if(i .eq. 5 .and. j .eq. 24) then
!            print *, "here"
!        endif

      return
      end subroutine
