c
c -----------------------------------------------------
c
      subroutine mymethod(q,qold,mitot,mjtot,lwidth,
     &                  dtn,dtnewn,dx,dy,nvar,xlow,ylow,mptr,maux,
     &                  aux,irr,lstgrd,ncount,numHoods,vtime,
     &                  istage,time,iir,jir)
c
c -----------------------------------------------------
c  compute one stage of an RK method. What comes in and what you do will vary
c  with stage. Set for 3 stage TVD RK now. 
c

      use amr_module
      implicit double precision (a-h,o-z)

      include "quadrature.i"
      dimension q(nvar,mitot,mjtot), qold(nvar,mitot,mjtot)
      dimension dU(nvar,mitot,mjtot)
      dimension f(nvar,mitot,mjtot),g(nvar,mitot,mjtot)
      dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
      dimension qxx(nvar,mitot,mjtot),qyy(nvar,mitot,mjtot)
      dimension qxy(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
      dimension ur(nvar,max(mitot,mjtot)),ul(nvar,max(mitot,mjtot))
      dimension ff(nvar,max(mitot,mjtot))
      dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)
      dimension Uout(nvar)
      dimension xp(max1d),yp(max1d)
      dimension iir(mitot,mjtot),jir(mitot,mjtot) 

      integer   irr(mitot,mjtot),ncount(mitot,mjtot)
      integer   numHoods(mitot,mjtot)

      common   /RKmethod/ coeff(5),mstage
      common   /order2/ ssw,quad,nolimiter
      include "cuserdt.i"

      logical   debug, vtime 
      logical   missing
      dimension firreg(nvar,-1:irrsize)
      dimension fakeStateCons(4), fakeStatePrim(4)
      data fakeStateCons/1.d0,0.d0,0.d0,2.5d0/
      data fakeStatePrim/1.d0,0.d0,0.d0,1.0d0/
      character ch

      integer    xrp, yrp
      data       debug/.false./
      data       xrp/1/, yrp/0/
      data       pi/3.14159265357989d0/

c
c xrp (= 1) to solve x riemann problem, yrp(= 0) for y riemann problem
c
c       Modified 2/3/89 to solve riemann problems one row at a time
c     # 
c
c cell (i,j) owns the fluxes to the left and bottom
c
      msize = max(mitot,mjtot)

       ! if this grid patch is all solid cells skip it
       ! this method only counts interior cells, not ghost cells.
       call countCellType(irr,mitot,mjtot,lwidth,numSolid,numCut,
     .                    numFull,lstgrd)
       nx = mitot-2*lwidth
       ny = mjtot-2*lwidth
       if (numSolid .eq. nx*ny) then
c        write(*,*) "grid ",mptr," all solid" 
         dtnewn = rinfinity     ! output var, need to return someting
         return
       endif
c      if (numFull .eq. nx*ny) write(*,*) "grid ",mptr," all flow" 
c      if (numCut+numSolid+numFull .ne. nx*ny) then
c         write(*,*)"count doesn't work for grid ",mptr
c      endif

c
      ar(-1) = 1.d0   ! fake value for solid area to avoid zero divides
c
c     # initializations
      firreg(:,-1)  = 0.d0
      f = 0.d0
      g = 0.d0
      qx = 0.d0
      qy = 0.d0
      qxx = 0.d0
      qxy = 0.d0
      qyy = 0.d0

      if (istage .eq. 3) then ! need to evaluate at qtemp, to q2 to make q3
         q = .75d0*qold + .25*q    ! q is q2  ! this is qtemp, called q to use below
      endif

c need routine to set face lengths and  midpoints
         call getirrlen(irr,mitot,mjtot,dtn,dx,dy,lstgrd,
     &               mptr,nvar,ffluxlen,gfluxlen)
c
c  turn on for debugging, or could move to regridding section of code
c        call cellsClose(ffluxlen,gfluxlen,mitot,mjtot,irr,lstgrd,
c    .                    lwidth)
 
c   ::::: 1 rk stage with linear reconstruction follows ::::::
c
c  store primitive variables in f for now
c
       !! q comes in as conserved variables. Computes slopes 
       call qslopes(q,qx,qy,qxx,qxy,qyy,
     &                mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &                 xlow,ylow,mptr,nvar,istage)

       !! now convert to pointwise primitive values
       !! for first test reconstruct in conserved vars`
       !call vctoprm(mitot,mjtot,nvar,
    !&               dx,dy,lstgrd,xlow,ylow,irr,
    !&               q,q,qx,qy,qxx,qxy,qyy)

       !! now compute slopes in primitive vars 
    !   call qslopes(q,qx,qy,qxx,qxy,qyy,
    ! &                mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
    ! &                 xlow,ylow,mptr,nvar,istage)


c
c  loop through rows of q calculating fluxes one row at time
c  vertical riemann problem first
c
      do 800 jcol = 2, mjtot-2 
c
         do 512  nn = 1, nlinequad
            do ii = 1, msize
              ur(:,ii) = fakeStatePrim(:) !initialize so ghost and missing vals still ok
              ul(:,ii) = fakeStatePrim(:) 
            end do
            do 511 i = 2, mitot-2
               call getYface_gauss(i,jcol,xface,yface,irr,mitot,mjtot,
     &                             xlow,ylow,dx,dy,lstgrd,nn,missing)
               if (missing) cycle

              call getCellCentroid(lstgrd,i,jcol+1,xcentp,ycentp,
     .                      xlow,ylow,dx,dy,irr(i,jcol+1))
              call getCellCentroid(lstgrd,i,jcol,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(i,jcol))
              kp = irr(i,jcol+1)
              k  = irr(i,jcol)
              xp(i) = xface
              yp(i) = yface

              dxp = xface - xcentp
              dyp = yface - ycentp

              dxm = xface - xcent
              dym = yface - ycent

              call evalU(Uout,dxp,dyp,dx,dy,q,qx,qy,qxx,qxy,qyy,
     &                  i,jcol+1,kp,mitot,mjtot,nvar)
              ur(:,i) = Uout(:)

              call evalU(Uout,dxm,dym,dx,dy,q,qx,qy,qxx,qxy,qyy,
     &                   i,jcol,k,mitot,mjtot,nvar)
              ul(:,i) = Uout(:)
c
c store fluxes in ff vector, copy into permanent flux array
c
  511      continue
           call vrmc(ur,ul,ff,2,mitot-2,yrp,msize)
c
           do 720 k = lwidth-2, mitot-lwidth+3
              g(:,k,jcol+1) = g(:,k,jcol+1) + wline(nn)*ff(:,k)
  720      continue
  512    continue ! end loop over # gauss pts
 800  continue
c
c
c    Horizontal riemann problems next
c
      do 900 irow = 2, mitot-2
c
        do 612  nn = 1, nlinequad
          do 611 j = 2, mjtot-2
            call getXface(irow,j,xface,yface,irr,mitot,mjtot,
     .               xlow,ylow,dx,dy,lstgrd,nn,missing)
            if (missing) then
               ur(:,j) = fakeState
               ul(:,j) = fakeState
               cycle
            endif

            call getCellCentroid(lstgrd,irow+1,j,xcentp,ycentp,
     .                      xlow,ylow,dx,dy,irr(irow+1,j))
            call getCellCentroid(lstgrd,irow,j,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(irow,j))

            kp = irr(irow+1,j)
            k = irr(irow,j)

            xp(j) = xface
            yp(j) = yface

            dxp = xface-xcentp
            dyp = yface-ycentp

            dxm = xface-xcent
            dym = yface-ycent

            call evalU(Uout,dxm,dym,dx,dy,q,qx,qy,qxx,qxy,qyy,irow,j,
     &                 k,mitot,mjtot,nvar)
            ur(:,j) = Uout(:) 

            call evalU(Uout,dx,dyp,dx,dy,q,qx,qy,qxx,qxy,qyy,irow+1,j,
     &                 kp,mitot,mjtot,nvar)
            ul(:,j) = Uout(:) 
  611   continue
  
c
c store fluxes in ff 
c
         call vrmc(ur,ul,ff,2,mjtot-2,xrp,msize)
c
         do 721 k = lwidth-2, mjtot-lwidth+3
            f(:,irow+1,k) = f(:,irow+1,k) + ff(:,k)
  721    continue
  612  continue  ! end loop over gauss pts
c
 900  continue
c
c irregflux computes the cut cell bndry flux. since no flow
c  through bndry use eval pressure there.
          call irregFlux_gauss(q,firreg,irr,mitot,mjtot,dx,dy,lstgrd,
     &                   xlow,ylow,mptr,qx,qy,qxx,qxy,qyy,
     &                   lwidth,nvar,time)
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  multiply fluxes by mesh spacing. 
c  this zeros out unused fluxes so solid vals dont get updated.
c
         do j = 1, mjtot
         do i = 1, mitot
            f(:,i,j) = f(:,i,j) * ffluxlen(i,j)
            g(:,i,j) = g(:,i,j) * gfluxlen(i,j)
         end do
         end do
c 
c      # finite volume update
c
         c      = coeff(istage)
         ar(-1) = 1.d0   ! prevent zero divides for solid cells
         do 918 j = iy1, iyn
         do 917 i = ix1, ixn

         k = irr(i,j)
           do m = 1, nvar
              dU(m,i,j) = -dtn/ar(k) * (f(m,i+1,j) - f(m,i,j) +
     &                                  g(m,i,j+1)-g(m,i,j)
     &                                - firreg(m,k))
           end do
 917     continue
 918     continue

         ! qold is always q(t_n), the old time
         if (istage .eq. 1) then
            q = q  + dU
         else if (istage .eq. 2) then  ! q1 came in, used in eval to make q2
            q = q  + dU
         else if (istage .eq. 3) then ! form and return q3
            ! done this way so SRD can be called on output q, then qfinal formed in advance
            q = q + dU    ! remember that q_temp was renamed to q
         endif
            

c       call checkPhys(q,irr,mitot,mjtot,mptr,istage,
c    .                  lstgrd,'from my_method')
c       call checkPhysInt(q,mitot,mjtot,mptr,istage,
c    .                    lwidth,'from my_method')
c       write(*,*)"done with physCheckInt for grid ",mptr

c     write(*,*)"calling symcheck after update stage ",istage,
c    .          " before srd"
c     call symcheck(q,irr,mitot,mjtot,nvar,lwidth)

c
c     # output irregular fluxes for debugging purposes:
      if (debug) then
         write(11,*)"grid ",mptr," flux:",k,i,j
         k = lstgrd
         do 811 ik = 1, irrsize
            k = iabs(nxtirr(k))
            if (k.eq.0) go to 822
               write(11,820) k,ixg(k),iyg(k),firreg(1,k)
               do 810 m = 2, nvar
                  write(11,821) firreg(m,k)
  810             continue
  811             continue
  822    continue
      endif
  820 format(3i4,2d22.10)
  821 format(12x,2d22.10)
c
c     # output fluxes for debugging purposes:
      if (debug) then
         write(11,*)"Inviscid: f             g          "
         do 834 i = lwidth+1, mitot-1
            do 833 j = lwidth+1, mjtot-1
               write(11,831) i,j,f(1,i,j),g(1,i,j)
               do 830 m = 2, nvar
                  write(11,832) f(m,i,j),g(m,i,j)
  830             continue
  833             continue
  834             continue
  831          format(2i4,4d20.10)
  832          format(8x, 4d20.10)
      endif

c
      return
      end
c
c --------------------------------------------------------------------
c
      subroutine symcheck(q,irr,mitot,mjtot,nvar,nghost)

      implicit real*8 (a-h,o-z)

      dimension q(nvar,mitot,mjtot)
      dimension irr(mitot,mjtot)

      jend = (mjtot-nghost)/2
      tol = 1.d-10

      do i = nghost+1,mitot-nghost
      do j = nghost+1,jend
         if (irr(i,j) .eq. -1) cycle
         j2 = mjtot - j + 1
         dif1 = abs(q(1,i,j) - q(1,i,j2))
         dif4 = abs(q(4,i,j) - q(4,i,j2))
         dif2 = abs(q(2,i,j) - q(2,i,j2))
         dif3 = abs(q(3,i,j) + q(3,i,j2))
         if (dif1 .gt. tol .or. dif2 .gt. tol .or.
     .       dif3 .gt. tol .or. dif4 .gt. tol) then
            write(*,100) i,j,dif1,dif2,dif3,dif4
 100        format(2i5," difs ",5e12.5)
         endif
      end do
      end do
      return
      end
