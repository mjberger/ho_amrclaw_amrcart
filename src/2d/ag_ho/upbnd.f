c
c ------------------------------------------------------------
c
       subroutine upbnd(listbc,val,nvar,maxip1,maxjp1,
     1                  maxsp,iused,dtfine,irreg,mitot,mjtot)
 
      implicit double precision (a-h,o-z)
      parameter  (maxgr = 192, maxlv=12)
      logical    graf
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      include "cirr.i"
 
       dimension val(maxip1,maxjp1,nvar),listbc(5,maxsp),
     1           iused(maxip1,maxjp1),irreg(mitot,mjtot)
 
c We now correct the coarse grid with the flux differences stored
c with each of the fine grids. We use an array iused
c to store whether the flux has been updated or not for that zone.
c iused(i,j) = sum from (l=1,4) i(l)*2**(l-1), where i(l) = 1 if the
c flux for the  l-th side of the (i,j)-th cell has already been
c updated, and i(l) = 0 if not. low end bit is side 1
 
c  the fluxes are adjust by dtfine instead of dtcoarse because
c  when they were saved by fluxad/fluxsv into svdflx, a factor
c  of intrat (refrat)  was included in the coarse grid to avoid
c  dividing each time the fine grid fluxes..
c
 
      do 10 j=1,maxjp1
      do 10 i=1,maxip1
         iused(i,j) = 0.
 10   continue
 
      do 40 ispot = 1,maxsp
         icrse = listbc(1,ispot)
         if (icrse.eq.0) go to 99
 
         jcrse = listbc(2,ispot)
         iside = listbc(3,ispot)
         norm = 2**(iside-1)
         iflag =iused(icrse,jcrse)/norm
         if (mod(iflag,2).eq.1) go to 40
         mkid = listbc(4,ispot)
         lkid = listbc(5,ispot)
         sgnm = 1.
         if (mod(iside,4).gt.1) sgnm = -1.
         kidlst = node(16,mkid)
c        # only update coarse cell if in fact it's interior to domain
	 if (irreg(icrse+lwidth-1,jcrse+lwidth-1) .ne. -1) then
            dtnorm = dtfine/ar(irreg(icrse+lwidth-1,jcrse+lwidth-1))
c   no time adjustment for steady calcs with local time stepping
             dtnorm = 1.d0/ar(irreg(icrse+lwidth-1,jcrse+lwidth-1))
             do 20 ivar = 1,nvar
                val(icrse,jcrse,ivar) = val(icrse,jcrse,ivar) +
     1          sgnm*dtnorm* alloc(kidlst+nvar*(lkid-1)+ivar-1)
 20          continue
	 endif
         iused(icrse,jcrse) = iused(icrse,jcrse) + norm
 40   continue
c
 99   return
      end
