      subroutine limiter_modal(q, qx, qy, qxx, qxy,
     .qyy, qxxx, qxxy, qxyy,
     .qyyy,
     .mitot, mjtot, irr, nvar, hx, hy,
     . lwidth, lstgrd, conserved)
       implicit double precision(a-h,o-z)
       include "./reconinfo.i"
       include "./quadrature.i"
       dimension q(mitot,mjtot,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
      dimension qxx(mitot,mjtot,4),qxy(mitot,mjtot,4)
      dimension qyy(mitot,mjtot,4)
      dimension qxxx(mitot,mjtot,4),qxxy(mitot,mjtot,4)
      dimension qxyy(mitot,mjtot,4),qyyy(mitot,mjtot,4)



      dimension iactive(mitot,mjtot,4, 3)


       common /order2/ ssw, quad, nolimiter
       dimension dumax(nvar),dumin(nvar),phimin(nvar), graddot(nvar)
       dimension dalpha(nvar), recon(nvar), reconlocal(nvar)
       logical quad, conserved

       if(quad .eqv. .false.) then
        print *, "PROBLEM: THE CURRENT LIMITER IMPLEMENTATION ONLY WORKS
     .            WHEN RECONSTRUCTING IN PRIMITIVE VARIABLES."
       endif


!      r3 = 2.d0
!      r2 = 1.5d0
!      r1 = 2.d0
      ! 1 means reconstruct to centroid
      ! 2 means reconstrut to edge

      r3 = 1.d0
      r2 = 1.d0
      r1 = 1.d0


      iactive = 0
      if(ssw .eq. 1) goto 111
      if(ssw .eq. 2) goto 222


      ! first do the quadratic coefficients
      do 40 i = 1, mitot
      do 40 j = 1, mjtot
      if (inuf(i,j) .eq. 1) then


        do nn = 1, nvar
        dup = (qxx(i+1,j,nn) - qxx(i,j,nn))/hx
        dum = (qxx(i,j,nn) - qxx(i-1,j,nn))/hx
        qxxx(i,j,nn) =  dminmod2( r3*dup, qxxx(i,j,nn) , r3*dum,
     .                        iactive(i,j,nn,3))

        dup = (qyy(i,j+1,nn) - qyy(i,j,nn))/hy
        dum = (qyy(i,j,nn) - qyy(i,j-1,nn))/hy
        qyyy(i,j,nn) =  dminmod2( r3*dup, qyyy(i,j,nn) , r3*dum,
     .                        iactive(i,j,nn,3))

!        if( iactive(i,j,nn,3) .eq. 0 ) cycle
!        iactive(i,j,nn,3) = 0 ! reset

        dup = (qxy(i+1,j,nn) - qxy(i,j,nn))/hx
        dum = (qxy(i,j,nn) - qxy(i-1,j,nn))/hx
        qxxy(i,j,nn) =  dminmod2( r3*dup, qxxy(i,j,nn) , r3*dum,
     .                        iactive(i,j,nn,3))

        dup = (qxx(i,j+1,nn) - qxx(i,j,nn))/hy
        dum = (qxx(i,j,nn) - qxx(i,j-1,nn))/hy
        qxxy(i,j,nn) =  dminmod2(r3*dup, qxxy(i,j,nn) , r3*dum,
     .                        iactive(i,j,nn,3))




        dup = (qyy(i+1,j,nn) - qyy(i,j,nn))/hx
        dum = (qyy(i,j,nn) - qyy(i-1,j,nn))/hx
        qxyy(i,j,nn) =  dminmod2(r3*dup, qxyy(i,j,nn) ,r3*dum,
     .                        iactive(i,j,nn,3))

        dup = (qxy(i,j+1,nn) - qxy(i,j,nn))/hy
        dum = (qxy(i,j,nn) - qxy(i,j-1,nn))/hy
        qxyy(i,j,nn) =  dminmod2(r3*dup, qxyy(i,j,nn) ,r3*dum,
     .                        iactive(i,j,nn,3))







        enddo
      else
          qxxx(i,j,:) = 0.d0
          qxxy(i,j,:) = 0.d0
          qxyy(i,j,:) = 0.d0
          qyyy(i,j,:) = 0.d0
          qxx(i,j,:) = 0.d0
          qxy(i,j,:) = 0.d0
          qyy(i,j,:) = 0.d0
          qx(i,j,:) = 0.d0
          qy(i,j,:) = 0.d0
      endif
 40     continue














 222  continue
      if(ssw .eq. 2) iactive(:,:,:,3) = 1
      ! first do the quadratic coefficients
      do 50 i = 1, mitot
      do 50 j = 1, mjtot



      if (inuf(i,j) .eq. 1) then
        do nn = 1, nvar
        if( iactive(i,j,nn,3) .eq. 1) then
        dup = (qx(i+1,j,nn) - qx(i,j,nn))/hx
        dum = (qx(i,j,nn) - qx(i-1,j,nn))/hx
        qxx(i,j,nn) =  dminmod2(r2*dup, qxx(i,j,nn) ,r2*dum,
     .                        iactive(i,j,nn,2))

        dup = (qy(i,j+1,nn) - qy(i,j,nn))/hy
        dum = (qy(i,j,nn) - qy(i,j-1,nn))/hy
        qyy(i,j,nn) =  dminmod2(r2*dup, qyy(i,j,nn) ,r2*dum,
     .                        iactive(i,j,nn,2))

!        if( iactive(i,j,nn,2) .eq. 0 ) cycle
!        iactive(i,j,nn,2) = 0 ! reset

        dup = (qy(i+1,j,nn) - qy(i,j,nn))/hx
        dum = (qy(i,j,nn) - qy(i-1,j,nn))/hx
        qxy(i,j,nn) =  dminmod2(r2*dup, qxy(i,j,nn) ,r2*dum,
     .                        iactive(i,j,nn,2))

        dup = (qx(i,j+1,nn) - qx(i,j,nn))/hy
        dum = (qx(i,j,nn) - qx(i,j-1,nn))/hy
        qxy(i,j,nn) =  dminmod2(r2*dup, qxy(i,j,nn) ,r2*dum,
     .                        iactive(i,j,nn,2))

        endif
        enddo
      else
          qxx(i,j,:) = 0.d0
          qxy(i,j,:) = 0.d0
          qyy(i,j,:) = 0.d0
          qx(i,j,:) = 0.d0
          qy(i,j,:) = 0.d0
      endif
 50     continue






 111  continue
      if(ssw .eq. 1) iactive(:,:,:,2) = 1

        ! next do the linear coefficients
        do 60 i = 1, mitot
        do 60 j = 1, mjtot
            if (inuf(i,j) .eq. 1) then
              do nn = 1, nvar
              if( iactive(i,j,nn,2) .eq. 1) then
              dup = (q(i+1,j,nn) - q(i,j,nn))/hx
              dum = (q(i,j,nn) - q(i-1,j,nn))/hx
              qx(i,j,nn) =  dminmod2(r1*dup, qx(i,j,nn) ,r1*dum,
     .                        iactive(i,j,nn,1))

              dup = (q(i,j+1,nn) - q(i,j,nn))/hy
              dum = (q(i,j,nn) - q(i,j-1,nn))/hy
              qy(i,j,nn) =  dminmod2(r1*dup, qy(i,j,nn) ,r1*dum,
     .                        iactive(i,j,nn,1))


              endif
              enddo
            else
                qx(i,j,:) = 0.d0
                qy(i,j,:) = 0.d0
            endif

 60     continue

!      call checkphysical(q, qx, qy, qxx, qxy, qyy, qxxx, qxxy, qxyy,
!     .                   qyyy, mitot, mjtot, irr, nvar, hx, hy,
!     .                   lwidth, lstgrd, conserved)






      return
      end


      function dminmod2(a, b, c, iactive)
        implicit double precision(a-h,o-z)

        dminmod2 = 0.d0
        if(a > 0.d0 .and. b > 0.d0 .and. c > 0.d0) then
            dminmod2 = min(min(a,b),c)
        elseif(a < 0.d0 .and. b < 0.d0 .and. c < 0.d0) then
            dminmod2 = max(max(a,b),c)
        endif

        if(abs(dminmod2 - b) > 1.d-9) iactive = 1


      return

      end function

!      subroutine checkphysical(q, qx, qy, qxx, qxy, qyy, qxxx, qxxy,
!     .qxyy,qyyy,mitot, mjtot, irr, nvar, hx, hy,lwidth, lstgrd,
!     . conserved)
!       implicit double precision(a-h,o-z)
!       include "./reconinfo.i"
!       include "./quadrature.i"
!       dimension q(mitot,mjtot,4),qx(mitot,mjtot,4),
!     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
!      dimension qxx(mitot,mjtot,4),qxy(mitot,mjtot,4)
!      dimension qyy(mitot,mjtot,4)
!      dimension qxxx(mitot,mjtot,4),qxxy(mitot,mjtot,4)
!      dimension qxyy(mitot,mjtot,4),qyyy(mitot,mjtot,4)
!       dimension qlocal(4),qxlocal(4),qylocal(4)
!      dimension qxxlocal(4),qxylocal(4),qyylocal(4)
!      dimension qxxxlocal(4),qxxylocal(4)
!      dimension qxyylocal(4),qyyylocal(4)
!      dimension recon(nvar), reconlocal(nvar)
!      logical quad, conserved
!      TOL = 1.d-10
!!      return
!
!        do 70 i = 1, mitot
!        do 70 j = 1, mjtot
!            k = irr(i,j)
!            if(k .eq. -1) cycle
!
!            do 934 iface = 1,4
!            do 934 nn = 1,nlinequad
!                if(iface .eq. 1) then
!                    diffx = 0.5d0*hx * 0.5d0*(1+rline(nn))
!     .                     -0.5d0*hx * 0.5d0*(1-rline(nn))
!                    diffy = -0.5d0*hy
!                else if(iface .eq. 2) then
!                    diffx =  0.5d0*hx
!                    diffy =  0.5d0*hy * 0.5d0*(1+rline(nn))
!     .                      -0.5d0*hy * 0.5d0*(1-rline(nn))
!                else if(iface .eq. 3) then
!                    diffx = 0.5d0*hx * 0.5d0*(1+rline(nn))
!     .                     -0.5d0*hx * 0.5d0*(1-rline(nn))
!                    diffy = 0.5d0*hy
!                else if(iface .eq. 4) then
!                    diffx = -0.5d0*hx
!                    diffy =  0.5d0*hy * 0.5d0*(1+rline(nn))
!     .                      -0.5d0*hy * 0.5d0*(1-rline(nn))
!                endif
!
!      do m = 1,4
!       call evalU(recon,diffx, diffy, hx, hy,
!     .            q, qx, qy,
!     .            qxx, qxy, qyy,
!     .            qxxx, qxxy, qxyy, qyyy,
!     .            i,j, k, mitot, mjtot, nvar, m)
!      enddo
!
!
!      ! limit density first
!      if (recon(1) .le. TOL) then
!          rho = recon(1)
!          avgrho = q(i,j,1)
!
!        thetarho = max(min( ( avgrho - TOL)/(avgrho-rho) , 1.d0),0.d0)
!
!          qx(i,j,  1) = qx(i,j,  1) * thetarho
!          qy(i,j,  1) = qy(i,j,  1) * thetarho
!          qxx(i,j, 1) = qxx(i,j, 1) * thetarho
!          qxy(i,j, 1) = qxy(i,j, 1) * thetarho
!          qyy(i,j, 1) = qyy(i,j, 1) * thetarho
!          qxxx(i,j,1) = qxxx(i,j,1) * thetarho
!          qxxy(i,j,1) = qxxy(i,j,1) * thetarho
!          qxyy(i,j,1) = qxyy(i,j,1) * thetarho
!          qyyy(i,j,1) = qyyy(i,j,1) * thetarho
!      endif
!
!      ! then limit pressure if we are in primitive variables
!      if((.not. conserved) .and. recon(4) .le. TOL) then
!            p = recon(4)
!            avgp = q(i,j,4)
!
!            thetap = max(min( ( avgp - TOL)/(avgp-p) , 1.d0),0.d0)
!
!
!            qx(i,j,  4) = qx(i,j,  4) * thetap
!            qy(i,j,  4) = qy(i,j,  4) * thetap
!            qxx(i,j, 4) = qxx(i,j, 4) * thetap
!            qxy(i,j, 4) = qxy(i,j, 4) * thetap
!            qyy(i,j, 4) = qyy(i,j, 4) * thetap
!            qxxx(i,j,4) = qxxx(i,j,4) * thetap
!            qxxy(i,j,4) = qxxy(i,j,4) * thetap
!            qxyy(i,j,4) = qxyy(i,j,4) * thetap
!            qyyy(i,j,4) = qyyy(i,j,4) * thetap
!      endif
! 934         continue
! 70     continue
!
!
!
!
!       ! now properly limit pressure
!       if(conserved) then
!        do 80 i = 1, mitot
!        do 80 j = 1, mjtot
!            k = irr(i,j)
!            if(k .eq. -1) cycle
!
!
!            tfinal = 1.d0
!
!            do 944 iface = 1,4
!            do 944 nn = 1,nlinequad
!                if(iface .eq. 1) then
!                    diffx = 0.5d0*hx * 0.5d0*(1+rline(nn))
!     .                     -0.5d0*hx * 0.5d0*(1-rline(nn))
!                    diffy = -0.5d0*hy
!                else if(iface .eq. 2) then
!                    diffx =  0.5d0*hx
!                    diffy =  0.5d0*hy * 0.5d0*(1+rline(nn))
!     .                      -0.5d0*hy * 0.5d0*(1-rline(nn))
!                else if(iface .eq. 3) then
!                    diffx = 0.5d0*hx * 0.5d0*(1+rline(nn))
!     .                     -0.5d0*hx * 0.5d0*(1-rline(nn))
!                    diffy = 0.5d0*hy
!                else if(iface .eq. 4) then
!                    diffx = -0.5d0*hx
!                    diffy =  0.5d0*hy * 0.5d0*(1+rline(nn))
!     .                      -0.5d0*hy * 0.5d0*(1-rline(nn))
!                endif
!
!      do m = 1,4
!       call evalU(recon,diffx, diffy, hx, hy,
!     .            q, qx, qy,
!     .            qxx, qxy, qyy,
!     .            qxxx, qxxy, qxyy, qyyy,
!     .            i,j, k, mitot, mjtot, nvar, m)
!      enddo
!
!
!      p1 = 0.4d0*(recon(4)-0.5d0*(recon(2)**2 + recon(3)**2) / recon(1))
!
!      if(p1 .gt. TOL) then
!        ttemp = 1.d0
!      else
!        ! solve the quadtratic
!        p0 = 0.4d0*(q(i,j,4)-0.5d0*(q(i,j,2)**2
!     .                      + q(i,j,3)**2) / q(i,j,1))
!        rhs0 = p0
!        rhs1 = p1
!        ttemp = p0/(p0-p1)
!
!!        ttemp = 0.d0
!      endif !not physical so had to compute quadratic
!      tfinal = min(tfinal, ttemp)
!
! 944         continue
!
!        qx(i,j,  :) = qx(i,j,  :) * tfinal
!        qy(i,j,  :) = qy(i,j,  :) * tfinal
!        qxx(i,j, :) = qxx(i,j, :) * tfinal
!        qxy(i,j, :) = qxy(i,j, :) * tfinal
!        qyy(i,j, :) = qyy(i,j, :) * tfinal
!        qxxx(i,j,:) = qxxx(i,j,:) * tfinal
!        qxxy(i,j,:) = qxxy(i,j,:) * tfinal
!        qxyy(i,j,:) = qxyy(i,j,:) * tfinal
!        qyyy(i,j,:) = qyyy(i,j,:) * tfinal
!
!
!      ! test the result
!        qlocal(:) =    q(i,j,:)
!        qxlocal(:) =   qx(i,j,  :)
!        qylocal(:) =   qy(i,j,  :)
!        qxxlocal(:) =  qxx(i,j, :)
!        qxylocal(:) =  qxy(i,j, :)
!        qyylocal(:) =  qyy(i,j, :)
!        qxxxlocal(:) = qxxx(i,j,:)
!        qxxylocal(:) = qxxy(i,j,:)
!        qxyylocal(:) = qxyy(i,j,:)
!        qyyylocal(:) = qyyy(i,j,:)
!
!            do 924 iface = 1,4
!            do 924 nn = 1,nlinequad
!                if(iface .eq. 1) then
!                    diffx = 0.5d0*hx * 0.5d0*(1+rline(nn))
!     .                     -0.5d0*hx * 0.5d0*(1-rline(nn))
!                    diffy = -0.5d0*hy
!                else if(iface .eq. 2) then
!                    diffx =  0.5d0*hx
!                    diffy =  0.5d0*hy * 0.5d0*(1+rline(nn))
!     .                      -0.5d0*hy * 0.5d0*(1-rline(nn))
!                else if(iface .eq. 3) then
!                    diffx = 0.5d0*hx * 0.5d0*(1+rline(nn))
!     .                     -0.5d0*hx * 0.5d0*(1-rline(nn))
!                    diffy = 0.5d0*hy
!                else if(iface .eq. 4) then
!                    diffx = -0.5d0*hx
!                    diffy =  0.5d0*hy * 0.5d0*(1+rline(nn))
!     .                      -0.5d0*hy * 0.5d0*(1-rline(nn))
!                endif
!
!       call evalUlocal(reconlocal,diffx, diffy, hx, hy,
!     .            qlocal, qxlocal, qylocal,
!     .            qxxlocal, qxylocal, qyylocal,
!     .            qxxxlocal, qxxylocal, qxyylocal, qyyylocal,nvar, k)
!
!      ptest = 0.4d0*(reconlocal(4)-0.5d0*(reconlocal(2)**2
!     . + reconlocal(3)**2) / reconlocal(1))
!      if(ptest < TOL .or. reconlocal(1) < TOL) then
!      print *,"here"
!      endif
! 924  continue






!
!
! 80     continue
!
!      endif
!
!
!      end subroutine




      !        qlocal(:) =            q(i,j,:)
!        qxlocal(:) =   0.5d0 * qx(i,j,  :)
!        qylocal(:) =   0.5d0 * qy(i,j,  :)
!        qxxlocal(:) =  0.5d0 * qxx(i,j, :)
!        qxylocal(:) =  0.5d0 * qxy(i,j, :)
!        qyylocal(:) =  0.5d0 * qyy(i,j, :)
!        qxxxlocal(:) = 0.5d0 * qxxx(i,j,:)
!        qxxylocal(:) = 0.5d0 * qxxy(i,j,:)
!        qxyylocal(:) = 0.5d0 * qxyy(i,j,:)
!        qyyylocal(:) = 0.5d0 * qyyy(i,j,:)
!
!       call evalUlocal(reconlocal,diffx, diffy, hx, hy,
!     .            qlocal, qxlocal, qylocal,
!     .            qxxlocal, qxylocal, qyylocal,
!     .            qxxxlocal, qxxylocal, qxyylocal, qyyylocal,nvar, k)
!      pmid = 0.4d0*(reconlocal(4)-0.5d0*(reconlocal(2)**2
!     . + reconlocal(3)**2) / reconlocal(1))
!      rhsmid = pmid - TOL
!
!
!      ! zeros of the quadratic are
!      zero1 = ((3 * rhs0) - (4 * rhsmid) + (rhs1) + sqrt((rhs0 ** 2
!     .- 2 * rhs0 * rhs1 - 8 * rhs0 * rhsmid + rhs1 ** 2
!     . - 8 * rhsmid * rhs1 + 16 * rhsmid ** 2))) / (rhs0
!     .- 2 * rhsmid + rhs1) / 0.4D1
!      zero2 = -(-(3 * rhs0) + (4 * rhsmid) - (rhs1) + sqrt((rhs0 ** 2
!     .- 2 * rhs0 * rhs1 - 8 * rhs0 * rhsmid + rhs1 ** 2
!     .- 8* rhsmid * rhs1 + 16 * rhsmid ** 2))) / (rhs0
!     .- 2 * rhsmid + rhs1) / 0.4D1
!
!      if(0 .lt. zero1 .and. zero1 .lt. 1.d0) then
!      ttemp = zero1
!      elseif(0 .lt. zero2 .and. zero2 .lt. 1.d0) then
!      ttemp = zero2
!      else
!      ttemp = -1.d0
!      print *, "NO S FOUND UH OH! ", ttemp, zero1, zero2
!      endif
!
!      ! test the result
!        qlocal(:) =    q(i,j,:)
!        qxlocal(:) =   qx(i,j,  :)
!        qylocal(:) =   qy(i,j,  :)
!        qxxlocal(:) =  qxx(i,j, :)
!        qxylocal(:) =  qxy(i,j, :)
!        qyylocal(:) =  qyy(i,j, :)
!        qxxxlocal(:) = qxxx(i,j,:)
!        qxxylocal(:) = qxxy(i,j,:)
!        qxyylocal(:) = qxyy(i,j,:)
!        qyyylocal(:) = qyyy(i,j,:)
!      call evalUlocal(reconlocal,diffx, diffy, hx, hy,
!     .            qlocal, qxlocal, qylocal,
!     .            qxxlocal, qxylocal, qyylocal,
!     .            qxxxlocal, qxxylocal, qxyylocal, qyyylocal,nvar, k)
!
!      ptest = 0.4d0*(reconlocal(4)-0.5d0*(reconlocal(2)**2
!     . + reconlocal(3)**2) / reconlocal(1))
!
!      ! test the result
!        qlocal(:) =    q(i,j,:)
!        qxlocal(:) =   ttemp * qx(i,j,  :)
!        qylocal(:) =   ttemp * qy(i,j,  :)
!        qxxlocal(:) =  ttemp * qxx(i,j, :)
!        qxylocal(:) =  ttemp * qxy(i,j, :)
!        qyylocal(:) =  ttemp * qyy(i,j, :)
!        qxxxlocal(:) = ttemp * qxxx(i,j,:)
!        qxxylocal(:) = ttemp * qxxy(i,j,:)
!        qxyylocal(:) = ttemp * qxyy(i,j,:)
!        qyyylocal(:) = ttemp * qyyy(i,j,:)
!      call evalUlocal(reconlocal,diffx, diffy, hx, hy,
!     .            qlocal, qxlocal, qylocal,
!     .            qxxlocal, qxylocal, qyylocal,
!     .            qxxxlocal, qxxylocal, qxyylocal, qyyylocal,nvar, k)
!
!      ptest1 = 0.4d0*(reconlocal(4)-0.5d0*(reconlocal(2)**2
!     . + reconlocal(3)**2) / reconlocal(1))
!
!      print *,"here"
!
!      ttemp2 = p0/(p0-p1)
!
!      ! test the result
!        qlocal(:) =    q(i,j,:)
!        qxlocal(:) =   ttemp2 * qx(i,j,  :)
!        qylocal(:) =   ttemp2 * qy(i,j,  :)
!        qxxlocal(:) =  ttemp2 * qxx(i,j, :)
!        qxylocal(:) =  ttemp2 * qxy(i,j, :)
!        qyylocal(:) =  ttemp2 * qyy(i,j, :)
!        qxxxlocal(:) = ttemp2 * qxxx(i,j,:)
!        qxxylocal(:) = ttemp2 * qxxy(i,j,:)
!        qxyylocal(:) = ttemp2 * qxyy(i,j,:)
!        qyyylocal(:) = ttemp2 * qyyy(i,j,:)
!      call evalUlocal(reconlocal,diffx, diffy, hx, hy,
!     .            qlocal, qxlocal, qylocal,
!     .            qxxlocal, qxylocal, qyylocal,
!     .            qxxxlocal, qxxylocal, qxyylocal, qyyylocal,nvar, k)
!
!      ptest2 = 0.4d0*(reconlocal(4)-0.5d0*(reconlocal(2)**2
!     . + reconlocal(3)**2) / reconlocal(1))
!
!                  print *,"here"

