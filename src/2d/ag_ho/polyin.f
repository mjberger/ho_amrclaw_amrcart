c
c
c ---------------------------------------------------------------
c
      subroutine polyin(p1,p2,p3,area,xc,yc)
      implicit double precision (a-h,o-z)
      dimension p1(10,2),p2(10,2),p3(10,2)
c
c     # Find the intersection of polygons p1 and p2.  Usually one is a
c     # wave and the other a grid cell.
c     # The intersection polygon is stored in p3 (which may overwrite
c     # p1 or p2 if desired).  The area of the intersection is stored
c     # in area. 
c     #   The intersection is found by intersecting p2 with the half space
c     # determined by the first two vertices of p1.  This polygon is then
c     # intersected with the half space determined by the next two vertices
c     # of p1, and so on.
c     #   p2 should be a closed polygon, but p1 need not be.  In particular
c     # p1 may be a half space with only two vertices.  The routine is 
c     # used in this manner to find the intersection of a wave with the
c     # half space generated by a boundary segment, for example.
c
c     # modified to also return the centroid (xc,yc) of p3 ..  Jan. 12, 1991
c
      dimension x(20,3), y(20,3)
      logical in,inold
      common/smarea/ epsar
c
      if (p1(1,1).eq.-11 .or. p2(1,1).eq.-11) then
c        # empty polygon
	 p3(1,1) = -11
	 area = 0.d0
	 go to 80
	 endif
c
      do 5 i=1,10
	 x(i,1) = p1(i,1)
	 y(i,1) = p1(i,2)
	 x(i,2) = p2(i,1)
	 y(i,2) = p2(i,2)
    5    continue
      ip1 = 1
      ip2 = 2
      ip3 = 3
      do 30 i=1,19
         k = 0
         x(1,ip3) = -11
	 if (x(i+1,ip1) .eq. -11) go to 40
	 x1 = x(i,ip1)
	 y1 = y(i,ip1)
	 x2 = x(i+1,ip1)
	 y2 = y(i+1,ip1)
	 delx = x1 - x2
	 dely = y2 - y1
c	 write(6,*) 'ip1,ip2,ip3,x1,y1,x2,y2:',ip1,ip2,ip3,
c    &	       x1,y1,x2,y2
c
	 do 10 j=1,10
c	    write(6,*) 'j,x,y:',j,x(j,ip2),y(j,ip2)
	    if (x(j,ip2) .eq. -11) go to 20
	    s = dely*(x(j,ip2)-x1) + delx*(y(j,ip2)-y1)
	    in = (s.ge.-epsar)
	    if (j.eq.1 .and. in) go to 2
	    if (j.eq.1) go to 3
	    icase = 1
	    if (in .and. inold) icase = 2
	    if (.not. in .and. .not. inold) icase = 3
	    go to (1,2,3) icase
    1       continue
c           # calculate intersection of lines:
	    det = (x2-x1)*(y(j-1,ip2)-y(j,ip2)) - (x(j-1,ip2)-x(j,ip2))
     &            *(y2-y1)
	    if (det.eq.0.) go to 2
	    alf = ((x(j-1,ip2)-x1)*(y(j-1,ip2)-y(j,ip2)) - (x(j-1,ip2)-
     &            x(j,ip2))*(y(j-1,ip2)-y1)) / det
	    k = k+1
	    x(k,ip3) = x1 + alf*(x2-x1)
	    y(k,ip3) = y1 + alf*(y2-y1)
	    if (k.eq.1) go to 7
	    if (x(k,ip3).eq.x(k-1,ip3) .and. y(k,ip3).eq.y(k-1,ip3))
     &         k = k-1
    7       continue
	    if (inold) go to 3
c
    2       continue
	    k = k+1
	    x(k,ip3) = x(j,ip2)
	    y(k,ip3) = y(j,ip2)
c
    3       continue
	    inold = in
c	    write(6,*) 'j,in,icase,k,det:',j,in,icase,k,det
   10       continue
         write(6,*) '*** too many vertices'
   20    continue
         if (x(1,ip3) .ne. -11) then
	 if (x(k,ip3).eq.x(1,ip3) .and. y(k,ip3).eq.y(1,ip3)) go to 22
         endif
         k = k+1
         x(k,ip3) = x(1,ip3)
         y(k,ip3) = y(1,ip3)
   22    continue
         x(k+1,ip3) = -11
         ip2 = 5-ip2
         ip3 = 5-ip3
c	 write(6,*) 'after intersection with half-space:'
c	 do 25 j=1,k
c  25	    write(6,*) x(j,ip2),y(j,ip2)
   30    continue
   40 continue
c
c     # compute area:
c
      area = 0.d0
      y0   = y(1,ip2)
      do 50 i=1,19
	 if (x(i+1,ip2).eq.-11) go to 60
	 area = area + .5d0*((y(i,ip2)-y0)+(y(i+1,ip2)-y0))
     &          *(x(i+1,ip2)-x(i,ip2))
   50    continue
      write(6,*) '*** too many vertices'
   60 continue
      p3(1,1) = x(1,ip2)
      p3(1,2) = y(1,ip2)
      if (x(i,ip2).eq.-11) go to 80
      k = 1
      do 70 i=2,20
	 if (x(i,ip2).eq.p3(k,1) .and. y(i,ip2).eq.p3(k,2)) go to 70
	 k = k+1
         if (k.gt.10) then
            write(6,*) '**** too many vertices in p3'
            go to 75
            endif
	 p3(k,1) = x(i,ip2)
	 p3(k,2) = y(i,ip2)
	 if (x(i,ip2).eq.-11) go to 80
   70    continue
c
      write(6,*) '*** too many vertices in x'
   75 continue
         write(6,*) '   p1:'
         do 110 i=1,10
            write(6,610) p1(i,1),p1(i,2)
  610       format(e23.6,e16.6)
  110       continue
         write(6,*) '   p2:'
         do 120 i=1,10
            write(6,610) p2(i,1),p2(i,2)
  120       continue
         write(6,*) '   p3:'
         do 130 i=1,10
            write(6,610) p3(i,1),p3(i,2)
  130       continue
   80 continue
c
c     # compute centroid:
c     if (area.gt.0.d0) then
      if (area.gt.1.d-16) then
          call centrd(p3,xc,yc,ajunk)
	else
	  xc = p1(1,1)
	  yc = p1(1,2)
	endif
      return
      end
