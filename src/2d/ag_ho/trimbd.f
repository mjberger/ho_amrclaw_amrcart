c
c ---------------------------------------------------------------
c
        subroutine trimbd(used,nx,ny,flag,il,ir,jb,jt)
c
c  if used array is completely set (=1.) then return flag=true, 
c  otherwise return false, and the dimensions of the smallest 
c  rectangles containing all unset points in il,ir,jb,jt.
c
        implicit double precision (a-h,o-z)
        dimension used(nx,ny)
        logical flag
 
        utot = 0.
        do 100 i = 1,nx
        do 100 j = 1,ny
100        utot = utot + used(i,j)
        if (utot .ge. dfloat(nx*ny)) then
                flag = .true.
                return
        endif
 
        flag = .false.
 
        uleft = 1.
        do 200 i = 1,nx
           do 220 j = 1,ny
              uleft = dmin1(uleft,used(i,j))
220        continue
           il = i
           if (uleft .eq. 0.) go to 230
200     continue

230     uright = 1.
        do 300 i = 1,nx
           do 320 j = 1,ny
              uright = dmin1(uright,used(nx - i + 1,j))
320        continue
           ir = nx - i + 1
           if (uright .eq. 0.) go to 330
300     continue

330     ubot = 1.
        do 400 j = 1,ny
           do 420 i = 1,nx
              ubot = dmin1(ubot,used(i,j))
420        continue
           jb = j
           if (ubot .eq. 0.) go to 430
400        continue
 
430     utop = 1.
        do 500 j = 1,ny
           do 520 i = 1,nx
              utop = dmin1(utop,used(i,ny - j + 1))
520        continue
           jt = ny - j + 1
           if (utop .eq. 0.) go to 530
500     continue
 
530     return
        end
