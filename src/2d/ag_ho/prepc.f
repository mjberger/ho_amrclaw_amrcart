c
c ----------------------------------------------------------
c
      subroutine prepc(level,nvar)
c
      implicit double precision (a-h,o-z)
      logical            ovrlap
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
c
c this routine called because regridding just changed the fine grids.
c modify coarse grid boundary lists to store fluxes in appropriate
c fine grids.
c assume new fine grids have node(15) initialized to point to null
c
c  first compute max. possible number of list cells. allocate
c  initially so that one pass through is enough. at end, release
c  unused space
c
      maxsp  = 0
      mkid   = lstart(level+1)
 10   if (mkid .eq. 0) go to 20
         ikeep  = (node(5,mkid)-1)/intrat(level)
         jkeep  = (node(6,mkid)-1)/intrat(level)
         maxsp  = maxsp + 2*(ikeep+jkeep)
      mkid = node(10,mkid)
      go to 10
 20   listsp(level) = maxsp
      if (maxsp .eq. 0) go to 99
c
      mpar = lstart(level)
 30   if (mpar .eq. 0) go to 99
c
       ispot   = 0
       hxpar   = rnode(9, mpar)
       hypar   = rnode(10,mpar)
       mipar   = node(5,mpar)
       mjpar   = node(6,mpar)
       locbc   = igetsp(5*maxsp)
       node(15,mpar) = locbc
c
c  initialize list space to 0 (use 0 terminator to indicate end of bc
c  list processing).
c
       call ziplst(alloc(locbc),maxsp)
c
       mkid = lstart(level+1)
 40    if (mkid .eq. 0) go to 60
       if (.not. ovrlap(mpar,mkid)) go to 50
          mikid  = node(5,mkid)
          mjkid  = node(6,mkid)
          hxkid  = rnode(9,mkid)
          hykid  = rnode(10,mkid)
c
          ist  = (rnode(1,mkid)-rnode(1,mpar))/hxpar + 1.01
          ist  = max0(ist,1)
          xist = rnode(1,mpar)+dfloat(ist-1)*hxpar
          kst  = (xist - rnode(1,mkid))/hxkid + 1.01
          iend = (rnode(7,mkid)-rnode(1,mpar))/hxpar+1.01
          iend = min0(iend,mipar)
          jst  = (rnode(2,mkid)-rnode(2,mpar))/hypar + 1.01
          jst  = max0(jst,1)
          xjst = rnode(2,mpar) + dfloat(jst-1)*hypar
          lst  = (xjst - rnode(2,mkid))/hykid + 1.01
          jend = (rnode(4,mkid)-rnode(2,mpar))/hypar+1.01
          jend = min0(jend,mjpar)
c
          call setuse(alloc(locbc),maxsp,ist,iend,jst,jend,ispot,mkid,
     1          kst,lst,intrat(level),mikid,mjkid,mipar,mjpar)
c
 50     mkid = node(10,mkid)
        go to 40
c
c  done with subgrid cycle. if no cells would need fixing, all done
c  else cycle through again to set up list with info. for bc processing
c
 60     continue
c  for now, leave unused space allocated to the grid. alternative is
c  to return (maxsp-ispot) amount starting at loc. node(15,mpar)+ispot.
c
       mpar = node(10,mpar)
       go to 30
c
 99    return
       end
