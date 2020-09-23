c
c -------------------------------------------------------------------
c
      subroutine outval(val,maxip1,maxjp1,nvar,irr,mitot,mjtot,
     1                  mptr,outgrd,lstgrd)
c
      implicit double precision (a-h,o-z)
      dimension val(maxip1,maxjp1,nvar), irr(mitot,mjtot)
      dimension pout(20,2), valprim(4)
      logical    outgrd
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
c     common  /cirr2/  rlink(9,2,irrsize)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead

      include "cirr.i"

      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      data    spval/1.e10/
c
c if cell outside domain, dont print soln. value - meaningless.
c account for border of lwidth cells around grid in irr array.
c
      if (.not. outgrd) go to 99
c
      ar(-1) = 0.d0
      hx     = rnode(9,mptr)
      hy     = rnode(10,mptr)
      xlow   = rnode(1,mptr) - lwidth*hx
      ylow   = rnode(2,mptr) - lwidth*hy
      ar(lstgrd) = rnode(9,mptr)*rnode(10,mptr)
      do 20 i=1,maxip1-1
          write(6,*)
      do 20 j=1,maxjp1-1
         if (iprob .ne. 20 .and. iprob .ne. 25) then
          valprim(1) = val(i,j,1)
          valprim(2) = val(i,j,2)/val(i,j,1)
          valprim(3) = val(i,j,3)/val(i,j,1)
          valprim(4) = (val(i,j,4)-.5*valprim(1)*(valprim(2)*valprim(2)+
     .                 valprim(3)*valprim(3)))*gamma1
         else
          valprim(1) = val(i,j,1)
         endif
          x  = rnode(1,mptr) + hx*(dfloat(i)-1.5)
          y  = rnode(2,mptr) + hy*(dfloat(j)-1.5)
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          if (irr(iadj,jadj) .eq. -1) then
c           solid body cell
c           write(6,109) x,y,iadj,jadj,(val(i,j,ivar),ivar=1,nvar)
109         format(2h *,2hx=,f7.3,3h y=,f7.3,4h, i=,i3,4h, j=,i3,
     *             ' a = ',4(f9.4,1x))
           elseif (irr(iadj,jadj) .eq. lstgrd) then
c          regular cell
           write(6,107) x,y,iadj,jadj,
     .            (valprim(ivar),ivar=1,nvar)
c    .            (val(i,j,ivar),ivar=1,nvar)
c107        format(2x,2hx=,f7.3,3h y=,f7.3,4h, i=,i3,4h, j=,i3,' a = ',
107        format(2x,2hx=,f7.3,3h y=,f7.3,4h, i=,i6,4h, j=,i6,' a = ',
     *                     4(e14.6,1x))
c    *                     4(f9.4,1x)) . . . mjb debug
           else
c          irregular cell
           kirr = irr(iadj,jadj)
           x    = xcirr(kirr)
           y    = ycirr(kirr)
           write(6,106) x,y,iadj,jadj,(valprim(ivar),ivar=1,nvar)
106        format(2h +,2hx=,f7.3,3h y=,f7.3,4h, i=,i6,4h, j=,i6,' a = ',
     *                      4(e14.6,1x))
c    *                      4(f9.4,1x))
          endif
 20   continue
c
 99   return
      end
