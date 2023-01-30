      subroutine wrxxf(iobs,ssh1,ssh2,ts1,ts2)
      implicit none
      include "data.inc"
      include "buffer.inc"
      integer ibuff,ios,writef,seekf,iobs
      real    ssh1,ssh2,ts1,ts2

      if (iobs.lt.minobs.or.iobs.gt.maxobs) 
     |call fin('inimini: error in WRXXF')

      ibuff=(iobs-minobs)*xxflen+1
      call emxxf(buffer(ibuff),buffer(ibuff),ssh1,ssh2,ts1,ts2)
      
*     Write old buffer

      if (iobs.eq.maxobs) then
	 ios=seekf(fd,minobs*xxflen,0)
         ios=writef(fd,(maxobs-minobs+1)*xxflen,buffer)
      endif

      end

      subroutine emxxf(line2,line4,ssh1,ssh2,ts1,ts2)

      integer*2 line2(*)
      integer*4 line4(*)
      include "data.inc"
      real    ssh1,ssh2,ts1,ts2

      if (short) then
         line4(8)=nint(ssh1*1e6)
         line4(9)=nint(ssh2*1e6)
         line2(23)=nint(ts1*1e3)
         line2(24)=nint(ts2*1e3)
      else
         line4(8)=nint(ssh1*1e6)
         line4(9)=nint(ssh2*1e6)
         line2(15+npar*4)=nint(ts1*1e3)   
         line2(16+npar*4)=nint(ts2*1e3)   
      endif

      end
