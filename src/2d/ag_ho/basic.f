c
c ------------------------------------------------------
c
      subroutine basic (time,lst,lend,iwrite)
c
      implicit double precision (a-h,o-z)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
c  basic = outputs basic information needed by every graphics
c          output routine (valout,drawrg). at the given time,
c          write the entire levellist, from level 1 to lfine,
c          and the tree structure from level lst to lend.
c
      write(iwrite,100) time
100   format(8h*TIME = ,f10.3)
      write(iwrite,101) lfine, (lstart(i),i=1,lfine), lwidth
101   format(10i6)
      write(iwrite,105) xprob,yprob,ismp
105   format(2e15.8,i5)
      write(iwrite,102) lst, lend
102   format(2i6)
c
      level = lst
 10   if (level .gt. lend) go to 99
          mptr = lstart(level)
 20       if (mptr .eq. 0) go to 30
              write(iwrite,103) mptr, (node(i,mptr),i=1,17)
              write(iwrite,104) (rnode(i,mptr),i=1,12)
103           format(10i7)
104           format(5e15.8)
              locirr = node(14,mptr)
              mitot  = node(5,mptr)-1+2*lwidth
              mjtot  = node(6,mptr)-1+2*lwidth
              call outirr(alloc(locirr),mitot,mjtot,node(17,mptr)) 
              mptr = node(10,mptr)
          go to 20
 30       level = level + 1
      go to 10
c
 99   return
      end
