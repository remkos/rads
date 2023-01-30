**WRIXAF -- Write one entire tracks to XAF file
*-
* 16-Oct-1994 - Old formats removed
*-----------------------------------------------------------------------
      subroutine wrixaf(iunit,isub)
      include "maxcom.inc"
      integer*4 i4data(7)
      integer*2 i2data(2)
      real*8    omega
      integer*4 i,k,irec,iunit,isub
      equivalence (i4data(7),i2data)
      save irec
      data irec/1/
*
      do k=1,isub
	 call grouplon(minlon,maxlon,track(3,k))
         i4data(1)=nint(track(1,k))
         i4data(2)=nint(track(2,k)/murad)
         i4data(3)=nint(track(3,k)/murad)
         omega=theta(itrknr)*(track(1,k)-time(3,itrknr))
         if (.not.asc(itrknr)) omega=omega+pi
         i4data(4)=nint(track(5,k)*1d6)
         i4data(5)=i4data(4)
         i4data(6)=nint(omega/murad)
	 i2data(1)=nint(altsig(satid(itrknr))*1d3)
         i2data(2)=itrknr
	 write (iunit,rec=irec+k) i4data
	 if (naux.gt.0) write (iunit+1,rec=irec+k+1)
     |		(atrack(i,k),i=1,naux)
      enddo
      irec=irec+isub
      end
