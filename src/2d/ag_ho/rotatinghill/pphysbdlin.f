


c
c ------------------------------------------------------------------
c
      subroutine pphysbdlin(xleft,xright,ybot,ytop,level,nrow,ncol,nvar,
     1                  val,time,hx,hy,qx,qy,irr,lstgrd)
 
c     This routine takes an (enlarged) grid (or grid patch)
c     with mesh widths hx,hy, and sets the values of any piece of
c     of the patch which extends outside the physical domain using the
c     values given by the boundary conditions. 
c
c 
c
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
      dimension val(nrow,ncol,nvar)
      dimension qx(nrow,ncol,nvar), qy(nrow,ncol,nvar)
      dimension irr(nrow,ncol)
      include "../cirr.i"


      hxmarg = hx*.01
      hymarg = hy*.01

c     left boundary (inflow)
 
      if (xleft .lt. -hxmarg) then
 
        nxl = (hxmarg-xleft)/hx
 
           do 400 i = 1,nxl
           do 400 j = 1,ncol
               call exact(xleft,xright,ybot,ytop,level,nrow,ncol,nvar,
     1                  val,time,hx,hy,qx,qy,irr,lstgrd, i, j)
400        continue
c
      endif
 
c     top boundary (extrap from interior)
 
      if (ytop .gt. yprob+hymarg) then
 
        nyt = (ytop - yprob + hymarg)/hy
        jbeg = max0(ncol-nyt+1, 1)
 
           do 100 j= jbeg,ncol
           do 100 i    = 1, nrow
             call exact(xleft,xright,ybot,ytop,level,nrow,ncol,nvar,
     1                  val,time,hx,hy,qx,qy,irr,lstgrd, i, j)

100        continue
c
      endif
 
c     right boundary. 
      if (xright .gt. xprob+hxmarg) then
 
        nxr = (xright - xprob + hxmarg)/hx
        nxr = nxr - 1
 
        ibeg = max0(nrow-nxr, 1)

c start extrap at bottom of grid, not including ghost cells
           if (ybot .lt. -hymarg) then
              jbeg = (hymarg-ybot)/hy + 1
           else
              jbeg = 1
           endif

           do 300 j = jbeg, ncol
           do 300 i = ibeg, nrow
            call exact(xleft,xright,ybot,ytop,level,nrow,ncol,nvar,
     1                  val,time,hx,hy,qx,qy,irr,lstgrd, i, j)
300        continue
c
      endif
 
 
c     bottom boundary. 
 
      if (ybot .lt. -hymarg) then
        nyb = (hymarg-ybot)/hy
           do 200 j = 1,nyb
           do 200 i = 1,nrow   ! the first lwidth set in left bc above
               call exact(xleft,xright,ybot,ytop,level,nrow,ncol,nvar,
     1                  val,time,hx,hy,qx,qy,irr,lstgrd, i, j)
200        continue
c
      endif
c
 99   return
      end

      subroutine exact(xleft,xright,ybot,ytop,level,nrow,ncol,
     1                      nvar,val,time,hx,hy,qx,qy,irr,lstgrd, i, j)

      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp

      dimension val(nrow,ncol,nvar)
      dimension qx(nrow,ncol,nvar), qy(nrow,ncol,nvar)
      dimension irr(nrow,ncol), outval(4)
      dimension state(nvar)
      include "../cirr.i"
      include "../quadrature.i"

      common  /order2/  ssw, quad, nolimiter
      logical quad
!       if(i .eq. 1 .and. j .eq. 189) then
!       print *, "here"
!       endif
      val(i,j,:) = 0.d0

          kirr = irr(i,j)
          if(kirr .eq. -1) return


          ! pointwise
          if(ssw .eq. -10 .or. ssw .eq. -1 .or. ssw .eq. 1) then
          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
             xcen = xcirr(kirr)
             ycen = ycirr(kirr)
          else if (kirr .eq. lstgrd) then
              xcen = xleft + (dfloat(i)-0.5d0)*hx
              ycen = ybot  + (dfloat(j)-0.5d0)*hy
          endif

          call f(state, xcen, ycen, 0.d0, iprob)
          val(i,j,:) = state(:)
          return
          endif



          if (kirr .ne. lstgrd) then
               ! cut cell projection



          ivert = 1
          do 20 while (poly(ivert+1,1,kirr) .ne. -11.)
            ivert = ivert + 1
  20     continue

          if(ihob .eq. 0) then
          itri = ivert - 3
          arr = ar(kirr)
          else
          itri = ivert - 4
          arr = ar_ho(kirr)
          endif

          idx1 = 1
          do 21 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,kirr)
            y1 = poly(idx1,2,kirr)

            x2 = poly(idx2,1,kirr)
            y2 = poly(idx2,2,kirr)

            x3 = poly(idx3,1,kirr)
            y3 = poly(idx3,2,kirr)

            artri = area(x1, x2, x3, y1, y2, y3)




            do 22 itq = 1,ntriquad

                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))


                  call f(state, xval, yval, time, iprob)

              val(i,j,1) = val(i,j,1) + (artri/arr)*wtri(itq) * state(1)
              val(i,j,2) = val(i,j,2) + (artri/arr)*wtri(itq) * state(2)
              val(i,j,3) = val(i,j,3) + (artri/arr)*wtri(itq) * state(3)
              val(i,j,4) = val(i,j,4) + (artri/arr)*wtri(itq) * state(4)
  22        continue ! for each quadrature point on each triangle
  21      continue ! for each triangle

      if(ssw .eq. 2 .and. ihob .eq. 1) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)
      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)
      x3 = bdry(3,1,kirr)
      y3 = bdry(3,2,kirr)
      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)

      do nq = 1,ntriquad_ho
      da = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)
      xc = rs2xy_2(x1,x2,x3,x4,nq)
      yc = rs2xy_2(y1,y2,y3,y4,nq)
      call f(state, xc, yc, time, iprob)
      val(i,j,:) = val(i,j,:) + (da/arr)*wtri_ho(nq) * state(:)
      enddo

      elseif(ssw .eq. 3  .and. ihob .eq. 1) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)
      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)
      x3 = bdry(4,1,kirr)
      y3 = bdry(4,2,kirr)
      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)
      x5 = bdry(3,1,kirr)
      y5 = bdry(3,2,kirr)

      do nq = 1,ntriquad_ho
      da = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)
      xc = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yc = rs2xy_3(y1,y2,y3,y4,y5,nq)
      call f(state, xc, yc, time, iprob)
      val(i,j,:) = val(i,j,:) + (da/arr)*wtri_ho(nq) * state(:)
      enddo
      endif




          else ! whole cell projection

          val(i,j,:) = 0.d0
          xcen1 = xleft + (dfloat(i)-1.d0)*hx
          ycen1 = ybot + (dfloat(j)-1.d0)*hy

          xcen3 = xleft + (dfloat(i)-0.d0)*hx
          ycen3 = ybot + (dfloat(j)-0.d0)*hy

        do 11 iq = 1, nquadquad
          xval = xcen1 * (1.d0 - rquad(iq))/2.
     .         + xcen3 * (1.d0 + rquad(iq))/2.
          yval = ycen1 * (1.d0 - squad(iq))/2.
     .         + ycen3 * (1.d0 + squad(iq))/2.

          call f(state, xval, yval, time, iprob)

          val(i,j,1) = val(i,j,1) + wquad(iq) * state(1)
          val(i,j,2) = val(i,j,2) + wquad(iq) * state(2)
          val(i,j,3) = val(i,j,3) + wquad(iq) * state(3)
          val(i,j,4) = val(i,j,4) + wquad(iq) * state(4)
 11   continue

          endif



      return
      end



      subroutine f(state, x, y, time, iprob)
        implicit double precision (a-h,o-z)
        dimension state(4), dir(2)
!        if(iprob .eq. 31) then
              pi = 3.14159265358979d0

            xp = (x-1.5d0) * cos(2*pi*time) + (y-1.5d0) *sin(2*pi*time)
            yp =-(x-1.5d0) * sin(2*pi*time) + (y-1.5d0) *cos(2*pi*time)




              theta = atan2(yp,xp) - pi/2.d0
              arg1 = (pi/6.d0 - theta)/dsqrt(4.d0/100.d0)
              arg2 = (pi/6.d0 + theta)/dsqrt(4.d0/100.d0)
              state(1) = 0.5d0*(erf(arg1) + erf(arg2))
              state(2) = 0.d0
              state(3) = 0.d0
              state(4) = 1.4d0

!              call getdir(dir, xp, yp)
!              state(1) = (x-dir(1)*time)**1
!               state(1) = (x-1.5d0)**2.d0 + (y-1.5d0)**2.d0

!          endif

      end


      function speed(q, i, j,mitot,mjtot,nvar, irr,lstgrd,
     .               xlow,ylow,dx,dy)
      implicit double precision (a-h,o-z)
        common /order2/ ssw, quad, nolimiter
        logical quad
        common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
        dimension q(mitot,mjtot,nvar)
        dimension dir(2), irr(mitot, mjtot)

        pi = 3.14159265358979d0
!        speed = 1.d0
!        speed = 2 * pi * 1.5d0

        kirr = irr(i,j)
        call getCellCentroid(lstgrd,i,j,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(i,j))

        call getdir(dir, xcent, ycent)
        speed = sqrt(dir(1)**2 + dir(2)**2)
        return
      end


      subroutine getbothspeeds(q, i, j, spx, spy,
     .               mitot,mjtot,meqn,irr,lstgrd,
     .               xlow,ylow,dx,dy)
      implicit double precision (a-h,o-z)

        common /order2/ ssw, quad, nolimiter
        logical quad
        common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
        dimension q(mitot,mjtot,nvar)
        dimension irr(mitot,mjtot)
        dimension dir(2)

        kirr = irr(i,j)
        call getCellCentroid(lstgrd,i,j,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(i,j))

        call getdir(dir, xcent, ycent)

        spx = dabs(dir(1))
        spy = dabs(dir(2))

      end


!      subroutine rarefpos(xleft,xright,ybot,ytop,level,nrow,ncol,
!     1                      nvar,val,time,hx,hy,qx,qy,irr,lstgrd,
!     1                      xcen,ycen)
!
!      implicit double precision (a-h,o-z)
!      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
!      dimension val(nvar)
!      dimension qx(nrow,ncol,nvar), qy(nrow,ncol,nvar)
!      dimension irr(nrow,ncol)
!      dimension dir(2)
!      include "../cirr.i"
!
!
!          pi = 3.14159265358979d0
!          dangle = pi/5.d0
!          dnx = cos(dangle)
!          dny = sin(dangle)
!
!          call getdir(dir)
!!          rhot = sin(dnx*(xcen-dir(1)*time) + dny*(ycen-dir(2)*time) )
!!          rhot = dnx*(xcen-dir(1)*time) + dny*(ycen-dir(2)*time)
!           rhot = xcen*xcen
!          pt = 1.d0
!
!
!          val(1) = rhot
!          val(2) = 0.d0
!          val(3) = 0.d0
!          val(4) = pt
!
!
!      return
!      end


