      function area(x1, x2, x3, y1, y2, y3)
        implicit double precision (a-h,o-z)

      area = dabs(0.5d0 * ( (x2 - x1) * (y3 - y1)
     .                    - (x3 - x1) * (y2 - y1) ))

      return
      end


