c
c ----------------------------------------------------------
c
      function igetsp (nwords)
c
      implicit double precision (a-h,o-z)
      parameter  (maxgr = 192, maxlv=12)
      common   /space/  lfree(150,2),lenf,idimf
      common   /stats/  evol, rvol, rvoll(maxlv), lentot, lenmax
      logical        sprint
      data           sprint/.false./
c
c  allocate contiguous space of length nword in main storage array
c  alloc. user code (or pointer to the owner of this storage)
c  is  mptr.  lenf = current length of lfree list.
c
c
c  find first fit from free space list
c
      itake = 0
      do 20 i = 1, lenf
         if (lfree(i,2) .lt. nwords) go to 20
         itake = i
         go to 25
 20   continue
      go to 900
c
c  put remaining words back on list.
c
 25   left = lfree(itake,2) - nwords
      igetsp = lfree(itake,1)
c
c  the following code which is ignored for now adds the new
      if (left .le. 0) go to 30
      lfree(itake,2) = left
      lfree(itake,1) = lfree(itake,1) + nwords
      go to 99
c
c  item is totally removed.  move next items in list up one.
c
 30   lenf = lenf - 1
      do 40 i = itake, lenf
      lfree(i,1) = lfree(i+1,1)
 40   lfree(i,2) = lfree(i+1,2)
      go to 99
c
 900  write(6,901) nwords
 901  format('  require ',i10,' words - either none left or not big',
     1         '  enough space')
      write(6,902) ((lfree(i,j),j=1,2),i=1,lenf)
 902  format(' free list: ',//,2x,50(i10,4x,i10,/,2x))
      mstart = 1   ! dump entire grid tree (always starts at mptr 1)
      call outtre(mstart,.false.,nvar)
      stop
c
 99   lentot = lentot + nwords
      if (lenmax .lt. lentot) lenmax = lentot
      if (sprint) write(6,100) nwords, igetsp, lentot, lenmax
 100  format('  allocating ',i10,' words in location ',i10,
     1       ' lentot, lenmax ', 2i10)
      return
      end
