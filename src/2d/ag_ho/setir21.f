c
c ------------------------------------------------------
c
       subroutine setir21(irreg2,mi2tot,mj2tot,
     1                irreg,mitot,mjtot,lwidth,lstgrd)
c
       implicit double precision (a-h,o-z)
       dimension  irreg2(mi2tot,mj2tot),irreg(mitot,mjtot)
c
c make the irregular array for the (interior) coarsened grid using the 
c irregular array for the interior of the finer grid.
c irreg = 1  -> regular cell
c irreg = -1 -> outside the domain
c irreg = k  -> set index to "0".
c set the boundary cells to regular for now as well.
c
c$dir scalar
       do 10 i = lwidth+1, mi2tot-lwidth
       ifine = 2*(i-lwidth) + lwidth - 1
c$dir scalar
       do 10 j = lwidth+1, mj2tot-lwidth
       jfine = 2*(j-lwidth) + lwidth - 1
       if ((irreg(ifine,jfine) .eq. lstgrd) .and.
     1 	      (irreg(ifine+1,jfine) .eq. lstgrd) .and.
     2	      (irreg(ifine,jfine+1) .eq. lstgrd) .and.
     3	      (irreg(ifine+1,jfine+1) .eq. lstgrd)) then
	     irreg2(i,j) = lstgrd
	  elseif ((irreg(ifine,jfine) .eq. -1) .and.
     1 	      (irreg(ifine+1,jfine) .eq. -1) .and.
     2	      (irreg(ifine,jfine+1) .eq. -1) .and.
     3	      (irreg(ifine+1,jfine+1) .eq. -1)) then
	     irreg2(i,j) = -1
	  else
c	     irreg2(i,j) = 0
	     irreg2(i,j) = -1
	  endif
 10    continue
c
c set boundary cells to be the same, without coarsening
c
       do 20 i = lwidth+1, mi2tot-lwidth
        ifine = 2*(i-lwidth) + lwidth - 1

	do 21 j = 1, lwidth
	   if ((irreg(ifine,j) .eq. -1) 
     .               .or. (irreg(ifine+1,1).eq.-1)) then
	       irreg2(i,j) = -1
	   else
	       irreg2(i,j) = lstgrd
	   endif
	   if ((irreg(ifine,mjtot+1-j) .eq. -1) 
     .               .or. (irreg(ifine+1,mjtot+1-j).eq.-1)) then
	       irreg2(i,mj2tot+1-j) = -1
	   else
	       irreg2(i,mj2tot+1-j) = lstgrd
	   endif
 21      continue
 20    continue
c
       do 30 j = lwidth+1, mj2tot-lwidth
        jfine = 2*(j-lwidth) + lwidth - 1

	do 31 i = 1, lwidth
        if ((irreg(i,jfine).eq.-1) .or. (irreg(i,jfine+1).eq.-1)) then
	     irreg2(i,j) = -1
	  else 
	     irreg2(i,j) = lstgrd
	 endif
        if ((irreg(mitot+1-i,jfine).eq.-1) .or. 
     &                 (irreg(mitot+1-i,jfine+1).eq.-1)) then
	     irreg2(mi2tot+1-i,j) = -1
	  else 
	     irreg2(mi2tot+1-i,j) = lstgrd
	 endif
 31      continue
 30    continue
c
c now the corners - lwdith by lwidth squares
c
       do 41 i = 1, lwidth
       do 41 j = 1, lwidth
	 irreg2(i,j)          = irreg(i,j)
	 irreg2(i,mj2tot+1-j) = irreg(1,mjtot+1-j)
	 irreg2(mi2tot+1-i,j) = irreg(mitot+1-i,j)
	 irreg2(mi2tot+1-i,mj2tot+1-j) = irreg(mitot+1-i,mjtot+1-j)
 41    continue
c
       return
       end
