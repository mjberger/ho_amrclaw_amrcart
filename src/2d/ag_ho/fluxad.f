c
c -------------------------------------------------------
c
      subroutine fluxad(xflux,yflux,svdflx,mptr,ndimx,ndimy,
     1                   maxi,maxj,nvar,lenbc,lratio)
c
      implicit double precision (a-h,o-z)
      dimension xflux(ndimx,ndimy,nvar), yflux(ndimx,ndimy,nvar)
      dimension svdflx(nvar,lenbc)
 
      maxjc = (maxj-1)/lratio
      maxic = (maxi-1)/lratio
 
c ::::: left side saved first
      lind = 0

      do 100 j=1,maxjc
         lind = lind + 1
         jfine = (j-1)*lratio
         do 110 ivar = 1, nvar
            do 120 l=1,lratio
               svdflx(ivar,lind) = svdflx(ivar,lind) +
     1				   xflux(1,jfine+l,ivar)
120         continue
110      continue
100   continue
 
c ::::: top side
      do 200 i=1,maxic
         lind = lind + 1
         ifine = (i-1)*lratio
         do 210 ivar = 1, nvar
            do 220 l=1,lratio
               svdflx(ivar,lind) = svdflx(ivar,lind) + 
     1				   yflux(ifine+l,maxj,ivar)
220         continue
210      continue
200   continue
 
c ::::: right side
      do 300 j=1,maxjc
         lind = lind + 1
         jfine = (j-1)*lratio
         do 310 ivar = 1, nvar
            do 320 l=1,lratio
               svdflx(ivar,lind) = svdflx(ivar,lind) + 
     1				   xflux(maxi,jfine+l,ivar)
320         continue
310      continue
300   continue
 
c ::::: bottom side
      do 400 i=1,maxic
         lind = lind + 1
         ifine = (i-1)*lratio
         do 410 ivar = 1, nvar
            do 420 l=1,lratio
               svdflx(ivar,lind) = svdflx(ivar,lind) +
     1				   yflux(ifine+l,1,ivar)
420         continue
410      continue
400   continue
 
      return
      end
