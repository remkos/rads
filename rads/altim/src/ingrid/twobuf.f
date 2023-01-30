      subroutine twobuf(kb1,kb2,reuse,error,quit)
*
* Prompt for source and target buffer
* If reuse=.true., source and target may be the same
*
      include "ingrid.inc"
      logical error,quit,reuse
      integer kb1,kb2
*
      call oldbuf('Enter source buffer',kb1,reuse,error,quit)
      if (error .or. quit) return
      if (kb1.lt.0) then
         kb1=-kb1
         kb2=kb1
      else
         blkd(kb1)=(.not.reuse)
         call newbuf('Enter target buffer',
     .   kb2,nx(kb1),ny(kb1),error,quit)
      endif
      end
