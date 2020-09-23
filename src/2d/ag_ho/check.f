c
c ---------------------------------------------------------
c
      subroutine check(nsteps,time,lfix)
c
      implicit double precision (a-h,o-z)
      logical            graf
      character  chkname*8
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "cirr.i"
      include "calloc.i"
      common   /space/  lfree(150,2),lenf,idimf
      common /userdt/cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .               ismp,gradThreshold
      common /stats/  evol, rvol, rvoll(maxlv), lentot, lenmax
      common  /cloops/  xloops(10),yloops(10),nloops
c
      if (.false.) then
      chkname = 'chkxxxxx'
      nstp = nsteps
      do 20 ipos = 8, 4, -1
         idigit = mod(nstp,10)
         chkname(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
 20   continue
      open(unit=10,file=chkname,form='unformatted')
c
      isize = allocsize
      write(10) lenmax, lentot, isize

c     output in chunks for solaris workaround
c     ist = 1
c     ichunk = 500000
c     numchunks = lenmax/ichunk + 1
c     write(*,*)" numchunks ",numchunks," lenmax ",lenmax
c     do  k = 1, numchunks
c        iend = min(ist+ichunk-1, lenmax)
c        write(*,*)" writing from",ist," to ",iend
c        write(10) (alloc(i),i=ist,iend)
c     ist = iend + 1
c     if (iend .ge. lenmax) go to 11
c     end do
c11   continue
c     write(6,*)" wrote checkpoint file ",chkname," in ",k,
c    .           " chunks of size ", ichunk
c     write(6,*) "numchunks = ", numchunks," lenmax ",lenmax
c     write(6,*) "free list length = ", lenf," of ",idimf

c     write(10)(alloc(i),i=1,lenmax)



      write(10) alloc
c
      write(10) hxposs,hyposs,possk,icheck,intcnt
      write(10) lfree,lenf,idimf
      write(10) rnode,node,lstart,newstl,listsp,tol,
     1  bzonex,bzoney,mstart,ndfree,lfine,iorder,mxnest,intrat,
     2        kcheck,maxcl, graf,lhead,nsteps,time
      write(10) cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp,lfix
      write(10) evol, rvol, rvoll
      write(10) poly,ar,xcirr,ycirr,points
      write(10) wt,ix,iy,nxtirr
      write(10) xloops,yloops,nloops
c
      close(10)
c     ierr = rename('chkfile',chkname)
c

        endif
      return
      end
