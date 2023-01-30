**WRIXTF -- Write entire XTF file
*+
* WRIXTF is the subroutine that writes the necessary output to
* the track catalogue file. For this additional information is
* calculated, such as the Epoch of ascending node,etc.
*-
* 16-Oct-1994 - Old binary formats removed
*-----------------------------------------------------------------------
      subroutine wrixtf
      include "maxcom.inc"
      integer*4   i,k,minnrx,null/0/
      parameter   (minnrx=2)
      integer*4   xtelen,i4xtf(6),i4par(maxpar),i4sig(maxpar)
      integer*2   itn,satin,nxor,noalt,flags
      character*4 filspf
*
* Open Track Catalogue File: XTF
*
      xtelen=34+npar*8
      open (wrunit,file=tfile,status='new',form='unformatted',
     .   access='direct',recl=xtelen,err=1010)
      filspf='@XTB'
      write (wrunit,err=1020,rec=1) filspf,null
*
* Compute ascending node, time of node, and start argument of latitude
* for valid tracks only.
* Fill track file. For invalid tracks write only the significant data;
* set other fields to zero.
*
      do k=1,itrknr
*        if (ixpt(k).lt.minnrx) good(k)=.false.
         if (good(k)) then
            inc(k)=abs(inc(k))/rad
            if (.not.asc(k)) then
               u0(k)=pi+u0(k)
               node(k)=node(k)-pi+pi/theta(k)*rotate
               time(3,k)=time(3,k)-pi/theta(k)
            endif
            node(k)=dmod(node(k)+6*pi,2*pi)
            i4xtf(1)=nint(inc(k)*1d6)
            i4xtf(2)=nint(u0(k)/murad)
            i4xtf(3)=nint(time(3,k))
            i4xtf(4)=nint(node(k)/murad)
         else
            i4xtf(1)=0
            i4xtf(2)=0
            i4xtf(3)=0
            i4xtf(4)=0
         endif
         i4xtf(5)=nint(time(1,k))
         i4xtf(6)=nint(time(2,k))
         flags=0
         if (asc(k))   flags=flags+1
         if (short(k)) flags=flags+2
         if (good(k))  flags=flags+256
         itn=k
         satin=satid(k)
         nxor=ixpt(k)
         noalt=iobs(k)
         do i=1,npar
            if (short(k) .and. .not.use(i)) then
               i4sig(i)=-nint(orberr(satin)*1d6)
            else
               i4sig(i)=nint(orberr(satin)*1d6)
            endif
	    i4par(i)=0
         enddo
	 i4par(1)=nint(-abias(satid(k))*1d6)
         write (wrunit,err=1020,rec=k+1) itn,satin,nxor,noalt,
     &         i4xtf,(i4par(i),i=1,npar),(i4sig(i),i=1,npar),flags
      enddo
*
* Write XTF header
*
      write (wrunit,err=1020,rec=1) filspf,itrknr,npar,bound,0
      close (wrunit)
      return
 1010 write (6,*) 'ERROR in wrixtf : cannot open XTF or XTE file'
      stop
 1020 write (6,*) 'ERROR in wrixtf : cannot write in XTF or XTE file'
      stop
      end
