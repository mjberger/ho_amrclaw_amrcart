c
c -----------------------------------------------------
c
      subroutine outgeom (lstgrd, mitot, mjtot, irr,xlow, ylow, dx, dy)

      implicit double precision (a-h,o-z)

      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      include "cirr.i"
      character*23  filename
      dimension irr(mitot,mjtot)
      dimension poly2(100,2,2)

      filename = 'outgeom.dat'
c      call makeNHood(volMerge,xcentMerge,ycentMerge,ncount,irr,numHoods,
c     .               mitot,mjtot,lwidth,lstgrd,xlow,ylow,dx,dy)
c
c      call getPerimeter(lstgrd,mitot,mjtot, lwidth, ncount, irr, dx, dy,
c     .                  xlow,ylow, perimeter)

      open(14,file=filename,status='unknown',form='formatted')
      write(*,*)" writing geometry file "


      xl   = rnode(1, mstart)
      yb   = rnode(2, mstart)
      xr   = rnode(7, mstart)
      yt   = rnode(4, mstart)

      write(14,*) mitot-2*lwidth, mjtot-2*lwidth
      write(14,*) xl, yb, xr, yt
      write(14,*) dx, dy

      do 15 i = lwidth+1, mitot-lwidth
           write(14,*)  (irr(i,j),j = lwidth+1, mjtot-lwidth)
  15  continue

      maxirr = -1
      do 16 j = lwidth+1, mjtot-lwidth
      do 16 i = lwidth+1, mitot-lwidth
        if(irr(i,j) > maxirr) maxirr = irr(i,j)


  16  continue


      do 17 k = 1,maxirr
!          write(14,*) "nco ", ncount_irr(k)
          kside = 1
          do 218 while (poly(kside+1,1,k) .ne. -11.)
            kside = kside + 1
  218       continue
          write(14,*) kside-1

          kside = 1
          do 28 while (poly(kside+1,1,k) .ne. -11.)
            write(14,*) poly(kside,1,k), poly(kside,2,k)
            kside = kside + 1
  28       continue
c           nco = ncount_irr(k)
c           if(i_irr(k) .eq. -1 .or. j_irr(k) .eq. -1) then
c               write(14,*) "nls ",0
c               goto 17
c           endif
c           call cellPerimeter(lstgrd, mitot, mjtot, lwidth, nco,
c     . irr,dx, dy, xlow, ylow, i_irr(k),j_irr(k), poly2)
c
c           ! count number of vertices in large polygon
c           kside = 1
c           do 219 while (poly2(kside,1,1) .ne. -11.)
c                kside = kside + 1
c  219       continue
c           write(14,*) "nls ", kside-1
c           kside = 1
c           do 29 while (poly2(kside,1,1) .ne. -11.)
c            write(14,*) poly2(kside,1,1), poly2(kside,2,1),
c     .                  poly2(kside,1,2), poly2(kside,2,2)
c            kside = kside + 1
c  29       continue


  17  continue


      close(14)

      return
      end





