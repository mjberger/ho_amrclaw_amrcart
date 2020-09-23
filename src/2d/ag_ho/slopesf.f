c
c------------------------------------------------------------
c
       subroutine slopesf(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar,dlimit)
       implicit double precision(a-h,o-z)
       dimension qp(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       logical  regular, quad, nolimiter
       logical MC/.true./
       dimension dlimitx(mitot,mjtot)
       dimension dlimity(mitot,mjtot)
       dimension dlimit(mitot,mjtot)
       dimension dumax(nvar),dumin(nvar),phimin(nvar), graddot(nvar)
       dimension dalpha(nvar), recon(nvar)

       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
       common /order2/ ssw, quad, nolimiter
       include "cirr.i"
       logical IS_GHOST
c      regular(i,j) = ((i.gt. lwidth).and.(i.le.mitot-lwidth).and.
c    &                 (j.gt. lwidth).and.(j.le.mjtot-lwidth))

      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)
c
c      # ssw = slope switch (1. for slopes, 0 for donor cell 0 slopes)
c      # now set in amrcart
c
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   compute slopes using muscl limiter
c   qp contains the primitive variables
c
      dlimitx = 1.d0
      dlimity = 1.d0
c
c  initialize all slopes to 0 first
c
       do 5 i = 1, mitot
       do 5 j = 1, mjtot
       do 5 m = 1, nvar 
          qx(i,j,m) = 0.d0
          qy(i,j,m) = 0.d0
 5     continue
       if (ssw .eq. 0.d0) go to 99
c
       do 10 m = 1, nvar
       do 10 i = 2, mitot-1
       do 10 j = 1, mjtot
          if (irr(i,j) .ne. lstgrd) go to 10
          ducc = qp(i+1,j,m) - qp(i-1,j,m)
          dupc = qp(i+1,j,m) - qp(i,j,m)
          dumc = qp(i,j,m)   - qp(i-1,j,m)

c         one sided at domain boundaries
c         x = xlow + (dfloat(i)-.5d0)*hx          ! cell center of this cell
c         if (x+hx .gt. xprob) ducc = 2.d0*dumc ! last pt has rt nbor outside domain
c         if (x-hx .lt. 0.d0)  ducc = 2.d0*dupc
c         turn off limiting if executing next 2 lines, with go to
          if (nolimiter) then
             qx(i,j,m) = .5d0*ducc/hx
c                # adjust for irregularity, if not done in qslopes
c$$$             if (irr(i+1,j) .ne. lstgrd) then
c$$$                 qx(i,j,m) = dumc/hx
c$$$             else if (irr(i-1,j) .ne. sltgrd) then
c$$$                 qx(i,j,m) = dupc/hx
c$$$             endif
             go to 10
          endif
         if(MC) then
         du   = dmin1(dabs(dupc),dabs(dumc))
         du   = dmin1(2.d0*du, .5d0*dabs(ducc))
         fl = dmax1(0.d0,dsign(1.d0, dupc*dumc))*ssw
         qx(i,j,m) = du*dsign(1.d0,ducc)*fl/hx
         else
         qx(i,j,m) = .5d0*ducc/hx
         endif
 10    continue
c
       do 20 m = 1, nvar
       do 20 i = 1, mitot
       do 20 j = 2, mjtot-1
       if (irr(i,j) .ne. lstgrd) go to 20
          ducc = qp(i,j+1,m) - qp(i,j-1,m)
          dupc = qp(i,j+1,m) - qp(i,j,m)
          dumc = qp(i,j,m)   - qp(i,j-1,m)
c         turn off limiting if execute next 2 lines, with go to
          if (nolimiter) then
             qy(i,j,m) = .5d0*ducc/hy
c                # adjust for irregularity, if not done in qslopes
c$$$             if (irr(i,j+1) .ne. lstgrd) then
c$$$                qy(i,j,m) = dumc/hy
c$$$             else if (irr(i,j-1) .ne. lstgrd) then
c$$$                qy(i,j,m) = dupc/hy
c$$$             endif
             go to 20
          endif

          if(MC) then
          du   = dmin1(dabs(dupc),dabs(dumc))
          du   = dmin1(2.d0*du, .5d0*dabs(ducc))
          fl = dmax1(0.d0,dsign(1.d0, dupc*dumc))*ssw
          qy(i,j,m) = du*dsign(1.d0,ducc)*fl/hy
          else
          qy(i,j,m) = .5d0*ducc/hy
          endif
 20    continue






      if(.not. MC) then
c        do 30 j = lwidth+1, mjtot-lwidth
c        do 30 i = lwidth+1, mitot-lwidth
         do 30 j = lwidth+2, mjtot-lwidth-1
        do 30 i = lwidth+2, mitot-lwidth-1
            k = irr(i,j)
            if (k .eq. -1) go to 30 ! solid cells have no gradient

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
 30     continue
      endif











       do 25 m = 1, nvar
       do 25 i = lwidth+2, mitot-lwidth-1
       do 25 j = lwidth+2, mjtot-lwidth-1
       if (irr(i,j) .ne. lstgrd) go to 25
         ducc = qp(i+1,j,m) - qp(i-1,j,m)
         if(  dsqrt(qx(i,j,1)*qx(i,j,1)
     .      + qy(i,j,1)*qy(i,j,1)) > 1e-10 .and. m .eq. 1) then
            dlimitx(i,j) = dabs(qx(i,j,m)/ (.5d0*ducc/hx))
         else if(dsqrt(qx(i,j,1)*qx(i,j,1)
     .      + qy(i,j,1)*qy(i,j,1)) .lt. 1e-10 .and. m .eq. 1) then
            dlimitx(i,j) = 1.d0
         endif

        ducc = qp(i,j+1,m) - qp(i,j-1,m)
         if(  dsqrt(qx(i,j,1)*qx(i,j,1)
     .      + qy(i,j,1)*qy(i,j,1)) > 1e-8  .and. m .eq. 1) then
            dlimity(i,j) = dabs(qy(i,j,m)/ (.5d0*ducc/hy))
         else if(dsqrt(qx(i,j,1)*qx(i,j,1)
     .      + qy(i,j,1)*qy(i,j,1)) .lt. 1e-8 .and. m .eq. 1) then
            dlimity(i,j) = 1.d0
         endif

 25    continue
        dlimitx = 1.d0 -dlimitx
        dlimity = 1.d0 -dlimity
        dlimit(:,:) = dmax1(dlimitx(:,:), dlimity(:,:))
        dlimit(:,:) = 0.d0
 99    continue
       return
       end
