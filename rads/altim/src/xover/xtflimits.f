**XTFLIMITS - Determine approximate start and end record for XAF based on XTF

      SUBROUTINE XTFLIMITS (NAME, T0, T1, TRK0, TRK1, L0, L1)
      CHARACTER*(*) NAME
      INTEGER*4 T0, T1, TRK0, TRK1, L0, L1

      integer*4 i,l,nrec,npar,xtf4(8),lrec,lnblnk,
     |		ios,fd,openf,readf,seekf
      integer*2 xtf2(37)
      character spec*4,filenm*80
      equivalence (xtf2,xtf4)

      l=index(name,'.xaf')
      if (l.le.0) l=lnblnk(name)
      filenm=name
      filenm(l:l+3)='.xtf'
      fd=openf(filenm,'r')
      ios=readf(fd,4,spec)
      if (spec.ne.'@XTB' .or. ios.ne.4) goto 9999
      ios=readf(fd,4,nrec)
      ios=readf(fd,4,npar)
      if (ios.ne.4) goto 9999
      lrec=34+8*npar
      ios=seekf(fd,lrec,0)
      if (ios.ne.lrec) goto 9999

      l=1
      l0=2147483647
      l1=0

      do i=1,nrec
         ios=readf(fd,lrec,xtf4)
	 if (ios.ne.lrec) then
	    goto 9999
	 else if (xtf2(4).eq.0) then
	 else if (xtf2(1).lt.trk0 .or. xtf2(1).gt.trk1) then
	 else if (xtf4(7).lt.t0 .or. xtf4(7).gt.t1) then
	 else
	    l0=min(l0,l)
	    l1=max(l1,l+xtf2(4)-1)
	 endif
	 l=l+xtf2(4)
      enddo
      return

9999  call closef(fd)
      end
