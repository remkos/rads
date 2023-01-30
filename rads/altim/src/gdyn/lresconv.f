      PROGRAM lresconv

* Convert Geodyn residual file to (small) integer binary format
*-
* 18-Mar-2001 - Avoid overflow residuals
* 18-Nov-2002 - Ensure Big Endian output on Little Endian machines
*-----------------------------------------------------------------------

      implicit none
      integer maxobs,maxpas
      parameter(maxobs=200,maxpas=2500)
      real*8 t1985
      parameter(t1985=(46066-30000)*86400d0)
      real*8 times(maxobs,maxpas),resid(maxobs,maxpas),
     |sigma(maxobs,maxpas),deriv(maxobs,maxpas),elev(maxobs,maxpas)
      character*8 cpass(maxpas)
      integer*4 ipass(4,maxpas),npass,ip,i,j,irec,iarg,iargc,kpass,nobs
      integer*4 irec4(2),zero/0/
      integer*2 irec2(4)
      character*80 filenm1,filenm2,arg
      logical altim/.true./,ltlend,swap
      data filenm1/'fort.19'/,filenm2/' '/

* Read command line arguments

      do iarg=1,iargc()
	 call getarg(iarg,arg)
	 if (arg(1:2).eq.'-a') then
	    altim=.false.
	 else if (filenm2.eq.' ') then
	    filenm2=arg
	 else
	    filenm1=filenm2
	    filenm2=arg
	 endif
      enddo
      swap=ltlend()

* Load residuals from Geodyn-native file

      call resread(filenm1,maxobs,maxpas,ipass,cpass,
     |	times,resid,sigma,deriv,elev)

* Count number of passes

      npass=0
      do ip=1,maxpas
	 if (ipass(3,ip).gt.0) npass=npass+1
      enddo

* Open output file

      open (20,file=filenm2,form='unformatted',status='new',
     |access='direct',recl=16)
      irec=1
      kpass=0
      do ip=1,npass
	 if (altim .or. ipass(4,ip).lt.99) then
	 kpass=kpass+1
         irec=irec+1
	 nobs=ipass(3,ip)
	 if (swap) call i4swap(4,ipass(1,ip))
         write (20,rec=irec) (ipass(i,ip),i=1,4)
         irec=irec+1
         write (20,rec=irec) cpass(ip),cpass(ip)
         do i=1,nobs
	    j=nint(resid(i,ip)*1d3)
	    j=max(-32768,min(j,32767))
            irec4(1)=int(times(i,ip)-t1985)
            irec4(2)=nint((times(i,ip)-t1985-irec4(1))*1d6)
            irec2(1)=j
            irec2(2)=nint(sigma(i,ip)*1d3)
            irec2(3)=nint(deriv(i,ip)    )
            irec2(4)=nint(elev (i,ip)*1d2)
            irec=irec+1
	    if (swap) then
	       call i4swap(2,irec4)
	       call i2swap(4,irec2)
	    endif
            write (20,rec=irec) irec4,irec2
         enddo
	 endif
      enddo
      if (swap) call i4swap(1,kpass)
      write (20,rec=1) '@RES',kpass,zero,zero
      close (20)
      end
