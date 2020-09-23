c
c -------------------------------------------------------------------
c
      subroutine outbig(val,irreg,mitot,mjtot,
     1                  mptr,lstgrd)
c
      implicit double precision (a-h,o-z)
      dimension val(mitot,mjtot,4), irreg(mitot,mjtot)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
c
c if cell outside domain, dont print soln. value - meaningless.
c account for border of lwidth cells around grid in irr array.
c
      do 20 i=1,mitot
      do 20 j=1,mjtot
          x  = rnode(1,mptr) + rnode( 9,mptr)*(dfloat(i)-2.5)
          y  = rnode(2,mptr) + rnode(10,mptr)*(dfloat(j)-2.5)
          if (irreg(i,j) .eq. -1) then
            write(6,109) x,y,i,j,(val(i,j,ivar),ivar=1,4)
109         format(2h *,2hx=,f7.3,3h y=,f7.3,4h, i=,i3,4h, j=,i3,
     *             ' a = ',4(f9.4,1x))
	   else if (irreg(i,j) .eq. lstgrd) then
              write(6,107) x,y,i,j,(val(i,j,ivar),ivar=1,4)
107         format(2x,2hx=,f7.3,3h y=,f7.3,4h, i=,i3,4h, j=,i3,' a = ',
     *                     4(f9.4,1x))
	   else
           write(6,106) x,y,i,j,(val(i,j,ivar),ivar=1,4)
106        format(2h +,2hx=,f7.3,3h y=,f7.3,4h, i=,i3,4h, j=,i3,' a = ',
     *                      4(f9.4,1x))
          endif
 20   continue
c
 99   return
      end
