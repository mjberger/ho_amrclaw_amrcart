c
c ------------------------------------------------------------------
c
       subroutine setuse(listbc,maxsp,ist,iend,jst,jend,ispot,
     1                   mkid,kst,lst,lratio,mikid,mjkid,maxi,maxj)
c
c set up boundary list for coarse grid. loop around boundary of
c fine grids to do this.  each entry has
c     i, j, side #, fine grid #, loc in fine grid list for fluxes.
c  for example, side 1 of fine grid fixes side 3 of coarse grid,
c  so coarse grid list will store the # 3.
c
      implicit double precision (a-h,o-z)
       dimension listbc(5,maxsp)
c
c  left side
c
       if ((ist .eq. 1) .or. (jst .ge. jend)) go to 20
       lkid        = (lst-1)/lratio + 1
       do 10 j     = jst+1, jend
          ispot              = ispot + 1
          listbc(1,ispot)    = ist
          listbc(2,ispot)    = j
          listbc(3,ispot)    = 3
          listbc(4,ispot)    = mkid
          listbc(5,ispot)    = lkid
          lkid               = lkid + 1
 10    continue
c
c   top side
c
 20    if ((jend .ge. maxj) .or. (ist .ge. iend)) go to 40
       lkid       = (mjkid-1)/lratio + (kst-1)/lratio + 1
       do 30 i    = ist+1, iend
          ispot              = ispot + 1
          listbc(1,ispot)    = i
          listbc(2,ispot)    = jend + 1
          listbc(3,ispot)    = 4
          listbc(4,ispot)    = mkid
          listbc(5,ispot)    = lkid
          lkid               = lkid + 1
 30    continue
c
c  right side (numbered from bottom to top, so not continuous)
c
 40    if ((iend .ge. maxi) .or. (jst .ge. jend)) go to 60
       lkid       = (mikid+mjkid-2)/lratio + (lst-1)/lratio + 1
       do 50 j         = jst+1, jend
          ispot              = ispot + 1
          listbc(1,ispot)    = iend + 1
          listbc(2,ispot)    = j
          listbc(3,ispot)    = 1
          listbc(4,ispot)    = mkid
          listbc(5,ispot)    = lkid
          lkid   = lkid + 1
 50    continue
c
c  bottom side (numbered left to right, so not continuous)
c
 60    if ((jst .le. 1) .or. (ist .ge. iend)) go to 80
       lkid       = (2*mjkid+mikid-3)/lratio + (kst-1)/lratio + 1
       do 70 i          = ist+1, iend
          ispot              = ispot + 1
          listbc(1,ispot)    = i
          listbc(2,ispot)    = jst
          listbc(3,ispot)    = 2
          listbc(4,ispot)    = mkid
          listbc(5,ispot)    = lkid
          lkid   = lkid + 1
 70    continue
c
 80    continue
       return
       end
