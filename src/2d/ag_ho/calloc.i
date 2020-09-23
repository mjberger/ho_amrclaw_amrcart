      integer allocSize
      parameter  (allocSize = 250000000)
      common   /calloc/  alloc(allocSize),hxposs(maxlv),hyposs(maxlv),
     *                   possk(maxlv),icheck(maxlv),intcnt(maxlv)

