c
c ------------------------------------------------------------
c
      integer function nodget(dummy)
c
      implicit double precision (a-h,o-z)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
c
c nodget =  get first free node of the linked list kept in node
c            array. adjust pointers accordingly.
c
      if (ndfree .ne. 0) go to 10
          write(6,100) maxgr
100       format(' out of nodal space - allowed ',i5,' grids')
          stop
c
 10     nodget         = ndfree
        ndfree         = node(2,ndfree)
c
c  initialize nodal block
c
        do 20 i        = 1,17
        node(i,nodget) = 0
 20     continue
        do 30 i         = 1,12
        rnode(i,nodget) = 0.0
 30     continue
c
      return
      end
