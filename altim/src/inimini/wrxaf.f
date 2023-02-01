      subroutine wrxaf(iobs,ssh,ts)
      implicit none
      include "data.inc"
      include "buffer.inc"
      integer ibuff,ios,writef,seekf,iobs
      real    ssh,ts

      if (iobs.lt.minobs.or.iobs.gt.maxobs) 
     |call fin('inimini: error in WRXAF')

      ibuff=(iobs-minobs)*xaflen+1
      call emxaf(buffer(ibuff),buffer(ibuff),ssh,ts)

* Write old buffer

      if (iobs.eq.maxobs) then
         ios=seekf(fd,minobs*xaflen,0)
         ios=writef(fd,(maxobs-minobs+1)*xaflen,buffer)
      endif

      end

      subroutine emxaf(line2,line4,ssh,ts)
      integer*2 line2(*)
      integer*4 line4(*)
      include "init.inc"
      include "data.inc"
      real    ssh,ts

      if (short) then
         line4( 5)=nint(ssh*1e6)
         line2(13)=nint(ts*1e3)
      else
         line4( 6)=nint(ssh*1e6)
         line2(13)=nint(ts*1e3)
      endif

      end
