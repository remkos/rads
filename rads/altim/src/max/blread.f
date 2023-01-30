c
c     BLREAD is a subroutine that reads in the block memories
c     and converts the integers to double precision
c
      subroutine blread(ibuf,irec,rdata,adata)
      integer*4 ibuf,irec,i
      real*8    rdata(*)
      integer*2 adata(*)
      include "maxcom.inc"

      if (ladr.eq.24) then
         rdata(1)=bdata ( 1,irec,ibuf)+bdata(2,irec,ibuf)*1d-6
         rdata(2)=bdata ( 3,irec,ibuf)*murad
         rdata(3)=bdata ( 4,irec,ibuf)*murad
         rdata(4)=bdata ( 5,irec,ibuf)*1d-6
         rdata(5)=bdata2(11,irec,ibuf)*hfac
      else
         rdata(1)=bdatx ( 1,irec,ibuf)+bdatx(2,irec,ibuf)*1d-6
         rdata(2)=bdatx ( 3,irec,ibuf)*murad
         rdata(3)=bdatx ( 4,irec,ibuf)*murad
         rdata(4)=bdatx ( 5,irec,ibuf)*1d-6
         rdata(5)=bdatx2(11,irec,ibuf)*hfac
	 do i=1,naux
	    adata(i)=bdatx2(11+i,irec,ibuf)
	 enddo
      endif

      end
