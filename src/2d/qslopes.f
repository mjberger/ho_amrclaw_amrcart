c
c ---------------------------------------------------------------------
c
       subroutine qslopes(qp,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,lstgrd,
     &                    lwidth,hx,hy,xlow,ylow,mptr,nvar) 

      use amr_module

      implicit double precision(a-h,o-z)
      include "cuserdt.i"

      dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &          qy(nvar,mitot,mjtot),irr(mitot,mjtot)
      dimension qxx(nvar,mitot,mjtot), qyy(nvar,mitot,mjtot)
      dimension qxy(nvar,mitot,mjtot)

      ! driver routine to call the correct gradient routine

      ! first call slope routines for regular cells
      call reg_slopes(qp,qx,qy,qxx,qxy,qyy,mitot,mjtot,
     &                irr,lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar)


      ! (use for both cell gradients and merge nhood gradients)
      ! igradChoice =     
      !  0 = no gradient
      !  1 = first order accurate gradients  (2 terms in least squares)
      !  2 = pointwise quadratic gradients (5 terms in least squares)
      !  3 = cellwise average  quadratic gradients (5 terms in least squares)


      if (igradChoice .eq. 1 .or. igradChoice .eq. 2) then
         call qslopesPtQuad2(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                       hx,hy,xlow,ylow,mptr,nvar)
      else if (igradChoice .eq. 3) then
         call qslopesCellAvgQuad(qp,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                       lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar)
      else if (igradChoice .eq. 0) then
         write(*,*)"should not be here  should have set ssw = 0"
      else
         write(*,*)"unrecognized gradient choice",igradChoice
         write(*,*)"should test in setprob and stop there"
         stop
      endif

      if (limitTile .eq. 1) then ! BJ limiter
          call limitCellBJ(qp,qx,qy,mitot,mjtot,irr,nvar,hx,hy,
     &                     lwidth,lstgrd,xlow,ylow)
      endif

      return
      end
