c
c -----------------------------------------------------------
c
      subroutine flglvl(nvar,lcheck,nxypts,index,ignore,lbase,
     .                  steady)
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
      logical            fprint, steady
      data               fprint/.false./
c
c *********************************************************
c
c flglvl = controls the error estimation/flagging bad pts. for
c          an entire level of grids.  returns pointer into alloc
c          where the (x,y) coordinations of the flagged pts. are.
c input parameters:
c           lcheck = level to be flagged
c output parameters:
c           nxypts = no. of flagged pts. total
c           index  = starting index in alloc of the flagged pts.
c                    (which occupy 2*nxypts locations).
c
c ***************************************************************
c
c  reserve space (nbuff points) for overflowing buffer points. will later 
c  search for other grids to contain those points to flag.c
c
      nbuff = 50000
      lbuff = igetsp(3 * nbuff)
      nxypts = 0
      ignore = 0
      iflow  = 0
c
      numbad = 0
      call errest(nvar,numbad,lcheck,alloc(lbuff),nbuff,iflow,
     .            steady)
      nxypts = nxypts + numbad
c
c find neighboring grids (if any) to receive buffer points that should
c have been flagged that exceeded grid boundaries
c
      numbuf = 0
      if (iflow .ne. 0) call flgbuf(alloc(lbuff),iflow,numbuf,lcheck)
      call reclam(lbuff,3*nbuff)
      nxypts = nxypts + numbuf
c
c  project newly created fine grids on to some grid at this level
c  to insure level nesting.
c
      numpro  =  0
      if (lcheck+2 .le. maxlv) call projec(lcheck, numpro, nvar)
      nxypts = nxypts + numpro
c
c  put all flagged points on into one array of their (x,y) coords.
c
 30   if (fprint) write(6,100) nxypts, lcheck
 100  format(' ',i5,' points flagged for level ',i4)
      if (nxypts .le. 0) go to 99
          index = igetsp(2*nxypts)
          call colate(alloc(index),nxypts,lcheck,nvar,ignore,lbase)
c
 99   return
      end
