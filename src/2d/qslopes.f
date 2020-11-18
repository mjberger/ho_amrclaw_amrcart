c
c ---------------------------------------------------------------------
c
       subroutine qslopes(qca,qpt,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                    lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar,
     &                    iir,jjr,istage) 

      use amr_module

      implicit double precision(a-h,o-z)
      include "cuserdt.i"

!    qca is for cell average
!    qpt is pointwise vals

      dimension qca(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &          qy(nvar,mitot,mjtot),irr(mitot,mjtot)
      dimension qxx(nvar,mitot,mjtot), qyy(nvar,mitot,mjtot)
      dimension qxy(nvar,mitot,mjtot)
      dimension qpt(nvar,mitot,mjtot)
      dimension iir(mitot,mjtot),jjr(mitot,mjtot)
      logical ag, orig

      ! driver routine to call the correct gradient routine

      ! first call slope routines for regular cells
      ag = .false.
      orig = .false.
      if (ag) then
         call reg_slopes2(qca,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar,istage)
      else if (orig) then
         call reg_slopes(qca,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar,istage)
      else ! higher order
         call reg_slopes3(qca,qpt,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar,istage)
      endif


      ! (use for both cell gradients and merge nhood gradients)
      ! igradChoice =     
      !  0 = no gradient
      !  1 = first order accurate gradients  (2 terms in least squares)
      !  2 = pointwise quadratic gradients (5 terms in least squares)
      !  3 = cellwise average  quadratic gradients (5 terms in least squares)
      !  4 = cellwise average  quads using Andrews iir approach             


      if (igradChoice .eq. 1 .or. igradChoice .eq. 2) then
         call qslopesPtQuad2(qca,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                       hx,hy,xlow,ylow,mptr,nvar)
      else if (igradChoice .eq. 3) then
         call qslopesCellAvgQuad(qca,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                       lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar,
     &                       istage)
      else if (igradChoice .eq. 4) then
         call qslopesQuadWithIIR(qca,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                       lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar,
     &                       iir,jjr,istage)
      else if (igradChoice .eq. 5) then
         call qslopes3(qca,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                       lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar,
     &                       iir,jjr,istage)
      endif

      if (limitTile .eq. 1) then ! BJ limiter
          call limitCellBJ(qca,qx,qy,mitot,mjtot,irr,nvar,hx,hy,
     &                     lwidth,lstgrd,xlow,ylow)
      endif

      return
      end
