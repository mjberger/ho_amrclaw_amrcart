c
c------------------------------------------------------------
c
       subroutine slopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       dimension dumax(nvar),dumin(nvar),phimin(nvar), graddot(nvar)
       dimension dalpha(nvar), recon(nvar), b(2)
       logical  regular, quad, nolimiter
       logical MC/.true./
       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
       common /order2/ ssw, quad, nolimiter
       logical IS_GHOST
       dimension a(2,2), rhs(2,nvar)
       include "cirr.i"
c      regular(i,j) = ((i.gt. lwidth).and.(i.le.mitot-lwidth).and.
c    &                 (j.gt. lwidth).and.(j.le.mjtot-lwidth))

c
c      # ssw = slope switch (1. for slopes, 0 for donor cell 0 slopes)
c      # now set in amrcart
c
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   compute slopes using muscl limiter
c   qp contains the primitive variables
c
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)
c
c  initialize all slopes to 0 first
c
      qx = 0.d0
      qy = 0.d0

      isloperecon = 3
      ! 0 is no slopes
      ! 1 is central
      ! 2 is MC
      ! 3 is LSQ + BJ

      if(isloperecon .eq. 0) then ! no slopes
       do 5 i = 1, mitot
       do 5 j = 1, mjtot
       do 5 m = 1, nvar 
          qx(i,j,m) = 0.d0
          qy(i,j,m) = 0.d0
 5     continue

      elseif(isloperecon .eq. 1) then ! central

       do 10 m = 1, nvar
       do 10 i = 2, mitot-1
       do 10 j = 1, mjtot
          if (irr(i,j) .ne. lstgrd) go to 10
          ducc = qp(i+1,j,m) - qp(i-1,j,m)
          qx(i,j,m) = .5d0*ducc/hx
 10    continue

       do 20 m = 1, nvar
       do 20 i = 1, mitot
       do 20 j = 2, mjtot-1
       if (irr(i,j) .ne. lstgrd) go to 20
          ducc = qp(i,j+1,m) - qp(i,j-1,m)
          qy(i,j,m) = .5d0*ducc/hy
 20    continue

      elseif(isloperecon .eq. 2) then ! MC

       do 30 m = 1, nvar
       do 30 i = 2, mitot-1
       do 30 j = 1, mjtot
          if (irr(i,j) .ne. lstgrd) go to 30
          ducc = qp(i+1,j,m) - qp(i-1,j,m)
          dupc = qp(i+1,j,m) - qp(i,j,m)
          dumc = qp(i,j,m)   - qp(i-1,j,m)
          du   = dmin1(dabs(dupc),dabs(dumc))
          du   = dmin1(2.d0*du, .5d0*dabs(ducc))
          fl = dmax1(0.d0,dsign(1.d0, dupc*dumc))*ssw
          qx(i,j,m) = du*dsign(1.d0,ducc)*fl/hx
 30    continue
       do 40 m = 1, nvar
       do 40 i = 1, mitot
       do 40 j = 2, mjtot-1
       if (irr(i,j) .ne. lstgrd) go to 40
          ducc = qp(i,j+1,m) - qp(i,j-1,m)
          dupc = qp(i,j+1,m) - qp(i,j,m)
          dumc = qp(i,j,m)   - qp(i,j-1,m)
          du   = dmin1(dabs(dupc),dabs(dumc))
          du   = dmin1(2.d0*du, .5d0*dabs(ducc))
          fl = dmax1(0.d0,dsign(1.d0, dupc*dumc))*ssw
          qy(i,j,m) = du*dsign(1.d0,ducc)*fl/hy
 40    continue

      elseif(isloperecon .eq. 3) then !LSQ + BJ

       do 34 m = 1, nvar
       do 34 i = 1, mitot
       do 34 j = 1, mjtot
            if (irr(i,j) .ne. lstgrd) go to 34
            if (IS_GHOST(i,j)) go to 34

            rhs = 0.d0 ! initialize for accumulation
            a = 0.d0
            do 22 joff = -1, 1
            do 22 ioff = -1, 1
                if (ioff .eq. 0 .and. joff .eq. 0) go to 22 ! no eqn to solve
                if (IS_GHOST(i+ioff,j+joff)) go to 22

                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 22

               deltax =  ioff*hx
               deltay =  joff*hy
               a(1,1) = a(1,1) + deltax*deltax
               a(1,2) = a(1,2) + deltax*deltay
               a(2,2) = a(2,2) + deltay*deltay
               rhs(1,:) = rhs(1,:) + deltax *
     .                    (qp(i+ioff,j+joff,:) - qp(i,j,:))
               rhs(2,:) = rhs(2,:) + deltay *
     .                    (qp(i+ioff,j+joff,:) - qp(i,j,:))
 22          continue

             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
             ! will have to add robustness checks
             c11 = sqrt(a(1,1))
             c12 = a(1,2) / c11
             c22 = sqrt(a(2,2) - c12**2)

             ! now back solve (C^t C = rhs of A^tdu ) to get x and y gradient for all variables
             do mm = 1, nvar
               b(1) = rhs(1,mm) / c11
               b(2) = (rhs(2,mm) - c12*b(1))/c22
               qy(i,j,mm) = b(2) / c22
               qx(i,j,mm) = (b(1) - c12*qy(i,j,mm))/c11
             end do
 34    continue

      endif


        do 50 j = lwidth+1, mjtot-lwidth
        do 50 i = lwidth+1, mitot-lwidth
            k = irr(i,j)
            if (k .eq. -1) go to 50 ! solid cells have no gradient

            ! find max and min needed for BJ limiting
            dumax = 0.d0
            dumin = 0.d0
            do 31 joff = -1, 1
            do 31 ioff = -1, 1
              if (ioff .eq. 0 .and. joff .eq. 0) go to 31
              if (IS_GHOST(i+ioff,j+joff)) go to 31
              koff = irr(i+ioff,j+joff)
              if (koff .eq. -1) go to 31
              dumax = max(dumax,qp(i+ioff,j+joff,:)-qp(i,j,:))
              dumin = min(dumin,qp(i+ioff,j+joff,:)-qp(i,j,:))
 31         continue

            phimin = 1.d0

            do 93 iface = 1,4
                if(iface .eq. 1) then
                    diffx = 0.5d0*hx
                    diffy = 0.d0
                else if(iface .eq. 2) then
                    diffx = 0.d0
                    diffy = 0.5d0*hy
                else if(iface .eq. 3) then
                    diffx = -0.5d0*hx
                    diffy = 0.d0
                else if(iface .eq. 4) then
                    diffx = 0.d0
                    diffy = -0.5d0*hy
                endif
                graddot  = qx(i,j,:)*diffx + qy(i,j,:)*diffy
                recon = qp(i,j,:) + graddot
                do m = 1,4
                   if (graddot(m) > 0.d0) then
                      dalpha(m) = min(1.d0, dumax(m)/graddot(m))
                   else if (graddot(m) < 0.d0) then
                      dalpha(m) = min(1.d0, dumin(m)/graddot(m))
                   else
                      dalpha(m) = 1.d0
                   endif
                end do
                ! one last check for positivity
                if (recon(1) .le. 0.d0) dalpha = 0.d0
                velsq = recon(2)**2+recon(3)**2
                press = .4d0*(recon(4)-0.5d0*velsq/recon(1))
                if (press .le. 0.d0) dalpha = 0.d0
                phimin = min(phimin, dalpha)
 93         continue

            qx(i,j,:) = qx(i,j,:)*phimin(:)
            qy(i,j,:) = qy(i,j,:)*phimin(:)
 50     continue








 99    continue
       return
       end
