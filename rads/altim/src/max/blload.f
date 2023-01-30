**BLLOAD -- Load a block of ADR records into memory
*-
      SUBROUTINE BLLOAD (IBUF, IBLK)
      implicit none
      integer*4 ibuf,iblk,irec0,irec,ifile,openf
      integer*4 fd,jblk(2)
      character*4 spec
      data jblk/2*0/
      save jblk
      include "maxcom.inc"

      if (iblk.eq.0) then
      else if (jblk(ibuf).eq.iblk) then
*	 write (*,*) 'blload: not loading',iblk
	 return
      else if (jblk(3-ibuf).eq.iblk) then
*	 write (*,*) 'blload: copy buffer ',3-ibuf,' to ',ibuf
	 if (ladr.eq.24) then
            call blcopy(ladr*blkobs(iblk),
     |		bdata(1,1,3-ibuf),bdata(1,1,ibuf))
	 else
            call blcopy(ladr*blkobs(iblk),
     |		bdatx(1,1,3-ibuf),bdatx(1,1,ibuf))
	 endif
	 jblk(ibuf)=iblk
	 return
      endif

      irec0=blkrec(iblk)
      irec =blkobs(iblk)
      ifile=blkfil(iblk)
*     write (*,*) 'blload:',ibuf,iblk,ifile,irec0,irec

      fd=openf(infile(ifile),'r')
      call readf(fd,4,spec)

      hfac=1d-3
      if (spec.eq.'@ADR') hfac=1d-2
      if (spec.eq.'xADR') then
	 ladr=22+3*2
	 call seekf(fd,irec0*ladr,0)
	 call readf(fd,irec *ladr,bdatx(1,1,ibuf))
      else
	 ladr=24
         call seekf(fd,irec0*ladr,0)
         call readf(fd,irec *ladr,bdata(1,1,ibuf))
      endif
      call closef(fd)
      jblk(ibuf)=iblk
      nrget=nrget+1
      end

      subroutine blcopy(n,b1,b2)
      implicit none
      integer*2 b1(*),b2(*)
      integer*4 n,i
      do i=1,n/2
         b2(i)=b1(i)
      enddo
      end
