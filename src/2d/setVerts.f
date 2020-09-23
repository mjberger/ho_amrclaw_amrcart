c
c -------------------------------------------------------------
c
       subroutine setVerts(poly,artri,x1,y1,x2,y2,x3,y3,ivert)

       implicit double precision (a-h,o-z)
       dimension poly(10,2)

       x1 = poly(1,1)
       y1 = poly(1,2)

       x2 = poly(ivert-1,1)
       y2 = poly(ivert-1,2)

       x3 = poly(ivert,1)
       y3 = poly(ivert,2)

       !artri = 0.5d0*((x2-x1) * (y3-y1) - (x3-x1) * (y2-y1) )
       artri = 0.5d0*((x2-x1)*(y1+y2)+(x3-x2)*(y2+y3)+(x1-x3)*(y3+y1))

       return
       end


