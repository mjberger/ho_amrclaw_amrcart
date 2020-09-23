


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
      include "../quadrature.i"
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
      dimension val(nrow,ncol,nvar)
      dimension qx(nrow,ncol,nvar), qy(nrow,ncol,nvar)
      dimension irr(nrow,ncol), outval(4)
      dimension state(nvar)
      include "../cirr.i"


      val(i,j,:) = 0.d0

          kirr = irr(i,j)
          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
               ! cut cell projection


!           if(kirr .eq. 7) then
!           print *, "here"
!           endif


          arr = ar(kirr)
          ivert = 1
          do 20 while (poly(ivert+1,1,kirr) .ne. -11.)
            ivert = ivert + 1
  20     continue

          itri = ivert - 3
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




            do 22 itq = 1,3

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


      ! always return primitive variables
      subroutine f(state, x, y, time, iprob)
        implicit double precision (a-h,o-z)
        dimension state(4), dir(2)
        common /order2/ ssw, quad, nolimiter
        logical quad
        gamma = 1.4d0
        pi = 3.14159265358979d0
!        if(iprob .eq. 21) then

        dspos = 1.d0/6.d0 + 10.d0 * time


        if(x < dspos) then
            ut =  8.25d0
            vt = 0.d0
            rhot = 8.d0
            pt = 116.5d0
        else
            ut = 0.d0
            vt = 0.d0
            rhot = gamma
            pt = 1.d0
        endif

         state(1) = rhot
         state(2) = rhot * ut
         state(3) = rhot * vt
         state(4) = pt/(gamma-1.d0) + 0.5d0*rhot*(ut**2 + vt**2)


!        endif

      end





      subroutine getbothspeeds(q, i, j, ux, uy, mitot,mjtot,nvar)
      implicit double precision (a-h,o-z)
        common /order2/ ssw, quad, nolimiter
        logical quad
        common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
        dimension q(mitot,mjtot,nvar)



        if(quad) then ! everything is in primitive variables
        c2 = gamma*q(i,j,4)/q(i,j,1)
        if (c2 .le. 0.d0) then
            print *, "non physical."
            speed = 0.d0
            return
        endif
        c = dsqrt(c2)
        u = q(i,j,2)
        v = q(i,j,3)
c               ::: dt using muscl cfl limit
        ux = max( abs(u+c), abs(u-c) )
        uy = max( abs(v+c), abs(v-c) )
!        speed = dsqrt(u*u+v*v) + c  !use max norm for muscl
        return

        else

        p = gamma1* (q(i,j,4)- .5d0* (q(i,j,2)*q(i,j,2)/
     &              q(i,j,1) + q(i,j,3)* q(i,j,3)/q(i,j,1)))
        c2 = gamma*p/q(i,j,1)
        if (c2 .le. 0.d0) then
            print *, "non physical."
            speed = 0.d0
            return
        endif
        c = dsqrt(c2)
        u = q(i,j,2)/q(i,j,1)
        v = q(i,j,3)/q(i,j,1)
c               ::: dt using muscl cfl limit
        ux = max( abs(u+c), abs(u-c) )
        uy = max( abs(v+c), abs(v-c) )

!        speed = dsqrt(u*u+v*v) + c  !use max norm for muscl
        return

        endif




      end






      function speed(q, i, j,mitot,mjtot,nvar)
      implicit double precision (a-h,o-z)
        common /order2/ ssw, quad, nolimiter
        logical quad
        common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
        dimension q(mitot,mjtot,nvar)


        if(quad) then ! everything is in primitive variables
        c2 = gamma*q(i,j,4)/q(i,j,1)
        if (c2 .le. 0.d0) then
            print *, "non physical."
            speed = 0.d0
            return
        endif
        c = dsqrt(c2)
        u = q(i,j,2)
        v = q(i,j,3)
c               ::: dt using muscl cfl limit
        speed = dsqrt(u*u+v*v) + c  !use max norm for muscl
        return

        else

        p = gamma1* (q(i,j,4)- .5d0* (q(i,j,2)*q(i,j,2)/
     &              q(i,j,1) + q(i,j,3)* q(i,j,3)/q(i,j,1)))
        c2 = gamma*p/q(i,j,1)
        if (c2 .le. 0.d0) then
            print *, "non physical."
            speed = 0.d0
            return
        endif
        c = dsqrt(c2)
        u = q(i,j,2)/q(i,j,1)
        v = q(i,j,3)/q(i,j,1)
c               ::: dt using muscl cfl limit
        speed = dsqrt(u*u+v*v) + c  !use max norm for muscl
        return

        endif
      end

!      subroutine raref(xleft,xright,ybot,ytop,level,nrow,ncol,
!     1                      nvar,val,time,hx,hy,qx,qy,irr,lstgrd, i, j)
!
!      implicit double precision (a-h,o-z)
!      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
!      dimension val(nrow,ncol,nvar)
!      dimension qx(nrow,ncol,nvar), qy(nrow,ncol,nvar)
!      dimension irr(nrow,ncol)
!      include "../cirr.i"
!
!
!          pi = 3.14159265358979d0
!          kirr = irr(i,j)
!          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
!             xcen = xcirr(kirr)
!             ycen = ycirr(kirr)
!          else
!             xcen = xleft + (dfloat(i)-0.5d0)*hx
!             ycen = ybot  + (dfloat(j)-0.5d0)*hy
!             kirr = lstgrd
!          endif
!
!c          dangle = pi/6.d0
!          dangle = pi/5.d0
!          dnx = cos(dangle)
!          dny = sin(dangle)
!          xi = (dnx*(xcen+0.25d0) + dny*(ycen+0.25d0))/(time+1.d0)
!
!
!          rhol = 8.d0
!          ul = 8.25d0
!          pl = 116.5d0
!          cl = dsqrt( gamma * pl / rhol )
!
!       solumag = ( (gamma - 1.d0)*ul + 2.d0*(cl + xi) ) / (gamma + 1.d0)
!          ut   = solumag * dnx
!          vt   = solumag * dny
!          rhot =( ( rhol**gamma * (solumag - xi)**2.d0 ) / (gamma*pl) )
!     .**( 1.d0/(gamma-1.d0) )
!          pt   = ( pl/(rhol**gamma) ) * rhot**gamma
!
!          val(i,j,1) = rhot
!          val(i,j,2) = ut
!          val(i,j,3) = vt
!          val(i,j,4) = pt
!
!
!
!      return
!      end
