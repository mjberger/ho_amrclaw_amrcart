c
c -------------------------------------------------------------------
c
      subroutine outval(val,maxip1,maxjp1,nvar,irr,mitot,mjtot,
     1                  mptr,outgrd,lstgrd)
c
      implicit double precision (a-h,o-z)
      dimension val(maxip1,maxjp1,nvar), irr(mitot,mjtot)
      dimension pout(20,2)
      logical    outgrd
      logical            graf
      parameter  (maxgr = 192, maxlv = 6)
      common  /cirr2/  rlink(9,2,2933)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      common  /cirr/ poly(10,2,2933),hs(10,2,2933),ar(-1:2933),
     *               points(24,2,2933),wt(24,-1:2933),
     *               xcirr(2933),ycirr(2933),ixsave(2933),iysave(2933),
     *               nxtirr(2933)
      common   /userdt/  gamma,gamma1,xprob,yprob,iprob,ismp,cfl
      data    spval/1.e10/
c
c if cell outside domain, don't print soln. value - meaningless.
c account for border of lwidth cells around grid in irr array.
c
      if (.not. outgrd) go to 99
c
c     cyl quantites
      velu   = 0.d0
      velv   = 0.d0
      ubd    = 0.d0
      vbd    = 0.d0
      uflerr = 0.d0
      vflerr = 0.d0
      ufltot = 0.d0
      vfltot = 0.d0
c     boundary quantities
      wdenerr = 0.d0
      wden    = 0.d0
      pbd = 0.d0
      pr1 = 0.d0
      odenerr = 0.d0
      oden    = 0.d0
      opbd = 0.d0
      opr1 = 0.d0
c     flowfield quantities
      totden = 0.d0
      totpr  = 0.d0
      totderr = 0.d0
      totperr = 0.d0
      wrongden = 0.d0
      artot    = 0.d0

      ar(-1) = 0.d0
      hx     = rnode(9,mptr)
      hy     = rnode(10,mptr)
      xlow   = rnode(1,mptr) - lwidth*hx
      ylow   = rnode(2,mptr) - lwidth*hy
      ar(lstgrd) = rnode(9,mptr)*rnode(10,mptr)
      do 20 i=2,maxip1-1
      do 20 j=2,maxjp1-1
          x  = rnode(1,mptr) + hx*(dfloat(i)-1.5)
          y  = rnode(2,mptr) + hy*(dfloat(j)-1.5)
	  iadj = i + lwidth - 1
	  jadj = j + lwidth - 1
          if (irr(iadj,jadj) .eq. -1) then
c           solid body cell
	    err1 = spval
	    err2 = spval
	    err3 = spval
	    err4 = spval
            write(6,109) x,y,iadj,jadj,(val(i,j,ivar),ivar=1,nvar)
109         format(2h *,2hx=,f7.3,3h y=,f7.3,4h, i=,i3,4h, j=,i3,
     *             ' a = ',4(f9.4,1x))
	   else if (irr(iadj,jadj) .eq. lstgrd) then
c               regular cell
                write(6,107) x,y,iadj,jadj,(val(i,j,ivar),ivar=1,nvar)
	        if (iprob.ge.15 .and. iprob .le. 19) then
	           if (iprob .eq. 17) call p17tru(x,y,rho,u,v,p)
	           if (iprob .eq. 18) call p18tru(x,y,rho,u,v,p)
	           if (iprob .eq. 19) then
		      call makep(poly(1,1,lstgrd),iadj,jadj,xlow,ylow,
     *                           hx,hy)
                      rlink(2,1,lstgrd) = -11
c                     call preppoly(poly(1,1,lstgrd),rlink(1,1,lstgrd),
c    .                              pout)
c		      call p19tru(x,y,rho,u,v,p,pout,lstgrd)
 		      call p19tru(x,y,rho,u,v,p,poly(1,1,lstgrd),lstgrd)
		   endif
	           if (iprob .eq. 15) then
		      call makep(poly(1,1,lstgrd),iadj,jadj,xlow,ylow,
     *                           hx,hy)
 		      call p15tru(x,y,rho,u,v,p,poly(1,1,lstgrd))
		   endif
	           ru = rho*u
	           rv = rho*v
	           en = p/gamma1 + 0.5d0*rho*(u**2+v**2)
	           write(6,117) rho,ru,rv,en
	           err1 = val(i,j,1)-rho
	           err2 = val(i,j,2)-ru
	           err3 = val(i,j,3)-rv
	           err4 = val(i,j,4)-en
	           write(6,127) err1,err2,err3,err4
		   err2 = val(i,j,2)/val(i,j,1) - u
		   err3 = val(i,j,3)/val(i,j,1) - v
	         endif
107         format(2x,2hx=,f7.3,3h y=,f7.3,4h, i=,i3,4h, j=,i3,' a = ',
     *                     4(f9.4,1x))
117         format(34x,'true: ',4(f9.4,1x))
127         format(33x,'error: ',4(d9.3,1x))
122         format('.',2i3,2f6.3,4e12.4)
	      else
c                irregular cell
		 kirr = irr(iadj,jadj)
	         x    = xcirr(kirr)
	         y    = ycirr(kirr)
                 write(6,106) x,y,iadj,jadj,(val(i,j,ivar),ivar=1,nvar)
	         if (iprob.ge.15 .and. iprob .le. 19) then
	           if (iprob .eq. 17) call p17tru(x,y,rho,u,v,p)
	           if (iprob .eq. 18) call p18tru(x,y,rho,u,v,p)
	           if (iprob .eq. 15) then
 		      call p15tru(x,y,rho,u,v,p,poly(1,1,kirr))
		   endif
	           if (iprob .eq. 19) then
c                     call preppoly(poly(1,1,kirr),rlink(1,1,kirr),pout)
c		      call p19tru(x,y,rho,u,v,p,pout,kirr)
 		      call p19tru(x,y,rho,u,v,p,poly(1,1,kirr),kirr)
		   endif
	           ru = rho*u
	           rv = rho*v
	           en = p/gamma1 + 0.5d0*rho*(u**2+v**2)
	           write(6,117) rho,ru,rv,en
	           err1 = val(i,j,1)-rho
	           err2 = val(i,j,2)-ru
	           err3 = val(i,j,3)-rv
	           err4 = val(i,j,4)-en
	           write(6,127) err1,err2,err3,err4
		   err2 = val(i,j,2)/val(i,j,1) - u
		   err3 = val(i,j,3)/val(i,j,1) - v
		 endif
106        format(2h +,2hx=,f7.3,3h y=,f7.3,4h, i=,i3,4h, j=,i3,' a = ',
     *                      4(f9.4,1x))
123        format('+',2i3,2f6.3,4e12.4)
          endif
	  if (iprob .ge. 15 .and. iprob .le. 19) then
	     k = irr(iadj,jadj)
	     pcomp = gamma1*(val(i,j,4)-.5d0*(val(i,j,2)**2+
     .                      val(i,j,3)**2)/val(i,j,1))
	     perr = pcomp - p
	     if (k.ne.lstgrd .and. k.ne.-1) then
		 wdenerr = wdenerr + sqrt(ar(k))*dabs(err1)
		 wden    = wden +    sqrt(ar(k))*rho
		 pbd     = pbd +     sqrt(ar(k))*dabs(perr)
		 pr1     = pr1 +     sqrt(ar(k))*p
c                :: also use alternate weighting for boundary norm
		 do 21 kside = 1, 6
		    if (poly(kside+2,1,k).eq.-11) then
		       sidex2 = poly(kside,1,k)
		       sidey2 = poly(kside,2,k)
		       sidex1 = poly(kside+1,1,k)
		       sidey1 = poly(kside+1,2,k)
		       go to 22
		    endif
 21              continue
 22		 rlen = dsqrt((sidey1-sidey2)**2 + (sidex1-sidex2)**2)
		 odenerr = odenerr + rlen*dabs(err1)
		 oden    = oden +    rlen*rho
		 opbd     = opbd +   rlen*dabs(perr)
		 opr1     = opr1 +   rlen *p
c                cyl quantities
		 velu    = velu + sqrt(ar(k))*dabs(u)
		 velv    = velv + sqrt(ar(k))*dabs(v)
		 ubd     = ubd + sqrt(ar(k))*dabs(err2)
		 vbd     = vbd + sqrt(ar(k))*dabs(err3)
	     endif
	     totderr = totderr + ar(k)*dabs(err1)
	     totperr = totperr + ar(k)*dabs(perr)
	     totden = totden + ar(k)*rho
	     totpr  = totpr  + ar(k)*p
	     artot  = artot + ar(k)
	     if (k .ne. -1) then
	       wrongden = wrongden + ar(k)*dabs(err1)/rho
	     endif
c            cyl quantities
	     uflerr = uflerr + ar(k)*dabs(err2)
	     vflerr = vflerr + ar(k)*dabs(err3)
	     ufltot = ufltot + ar(k)*dabs(u)
	     vfltot = vfltot + ar(k)*dabs(v)
	 endif
 20   continue
c
 99   if (iprob .ge. 15 .and. iprob .le. 19)  then
	write(6,*)'  '
        write(6,*)' 1 norm of sboundary density error     = ', wdenerr
        write(6,*)' 1 norm of sboundary density           = ', wden
        write(6,*)' 1 norm of spressure error             = ', pbd
        write(6,*)' 1 norm of spressure                   = ', pr1
	write(6,*)'  '
        write(6,*)' 1 norm of rboundary density error     = ', odenerr
        write(6,*)' 1 norm of rboundary density           = ', oden
        write(6,*)' 1 norm of rpressure error             = ', opbd
        write(6,*)' 1 norm of rpressure                   = ', opr1
	write(6,*)'  '
	write(6,*)' 1 norm total ff density error        = ',totderr
	write(6,*)' reg 1 norm of ff density             = ',totden
	write(6,*)' 1 norm total ff pressure error       = ',totperr
	write(6,*)' reg 1 norm of ff pressure            = ',totpr 
	write(6,*)' wrong 1norm of ff density rel. err.  = ',wrongden
	write(6,*)' total area                           = ',artot

	if (iprob .eq.17) then
	  write(6,*)'  '
	  write(6,*)' 1 norm of u velocity bdry error      = ', ubd
	  write(6,*)' 1 norm of v velocity bdry error      = ', vbd
	  write(6,*)' 1 norm of u bdry velocity            = ', velu
	  write(6,*)' 1 norm of v bdry velocity            = ', velv

	  write(6,*)'  '
	  write(6,*)' 1 norm of u velocity ff error        = ', uflerr
	  write(6,*)' 1 norm of v velocity ff error        = ', vflerr
	  write(6,*)' 1 norm of u ff velocity              = ', ufltot
	  write(6,*)' 1 norm of v ff velocity              = ', vfltot
	endif
      endif

      return
      end
