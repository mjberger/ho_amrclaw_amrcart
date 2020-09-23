c
c  -----------------------------------------------------------
c
      subroutine reclam (index, nwords)
c
c  return of space. add to free list.
c  iplace points to next item on free list with larger index than
c  the item reclaiming, unless said item is greater then
c  everything on the list.
c
      implicit double precision (a-h,o-z)
      parameter  (maxgr = 192, maxlv=12)
      common   /space/ lfree(150,2),lenf,idimf
      common   /stats/ evol, rvol, rvoll(maxlv), lentot, lenmax
      logical          sprint
      data             sprint/.false./
c
      do 20 i = 1, lenf
      iplace  = i
      if (lfree(i,1) .gt. index) go to 30
 20   continue
         write(6,902)
 902     format(' no insertion pointer into freelist. error stop')
         stop
c
c  check previous segment for merging
c
 30      iprev = iplace - 1
         if (lfree(iprev,1)+lfree(iprev,2) .lt. index) go to 40
         lfree(iprev,2) = lfree(iprev,2) + nwords
         go to 50
c
c  check after segment - no previous merge case
c
 40   nexti = index + nwords
      if (lfree(iplace,1).ne. nexti) go to 70
         lfree(iplace,1) = index
         lfree(iplace,2) = lfree(iplace,2) + nwords
         go to 99
c
c  check following segment - yes previous merge case
c
 50   nexti = index + nwords
      if (lfree(iplace,1) .ne. nexti) go to 99
c
c forward merge as well, bump all down 1
c
      lfree(iprev,2) = lfree(iprev,2)+lfree(iplace,2)
      ipp1           = iplace + 1
         do 60 i = ipp1, lenf
         lfree(i-1,1) = lfree(i,1)
 60      lfree(i-1,2) = lfree(i,2)
         lenf = lenf - 1
         go to 99
c
c  no merges case - insert and bump future segments up to make room
c
 70   if (lenf .eq. idimf) go to 900
      do 80 ii = iplace, lenf
      i          = lenf + 1 - ii + iplace
      lfree(i,1) = lfree(i-1,1)
 80   lfree(i,2) = lfree(i-1,2)
      lenf            = lenf + 1
      lfree(iplace,1) = index
      lfree(iplace,2) = nwords
      go to 99
c
 900  write(6,901) idimf
 901  format('  free list full with ',i5,' items')
      stop
c
 99   lentot = lentot - nwords
      if (sprint) write(6,100) nwords, index, lentot
 100  format('     reclaiming ',i10,' words at loc. ',i10,
     .      ' lentot ',i10)
      return
      end
