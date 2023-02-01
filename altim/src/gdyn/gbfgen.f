**GBFGEN -- Generate a GBF file
*
* GBFGEN is a program to generate a GBF file (Geodyn Binary
* Formatted file). It uses station dependent sigmas from a sigma
* file.
*
* Syntax:
* gbfgen [options] outfile
* where
* outfile   : output GBF file
* and [options] are:
* t=t0,t1,dt: select period t0,t1 (sec85,mjd,yymmdd,yyddd)
*             and interval dt (seconds)
* sig=file  : read SIGMA cards from "file"
* range=rng : use fixed range (m)
*-
* 20-Jun-1997 -- Created by Remko Scharroo from gbfedit
* 12-Apr-2000 -- Added DATEARG. Removed unneeded variables
* 16-May-2000 -- Adjusted manual
*-----------------------------------------------------------------------
      program gbgen
      implicit none

* Arguments

      character*80 arg,sigfile/' '/,filenm/' '/
      integer iargc
      real*8 t0/-1d30/,t1/1d30/,dt/15d0/
      logical datearg

* Sigma cards

      integer nsta
      real*4  sigma(1000)
      integer ista(1000)

* GBF variables

      real*4 trop,cion,cmass
      integer*2 ityp,indtim
      integer*4 isat,ipre,mjd,idum1,idum2,imet,idum3,t
      real*8 fract,range

* General

      integer i

* Initialise

       trop=0
       cion=0
       cmass=0
       range=0
       idum1=0
       idum2=0
       idum3=0
       ipre=371982336
       imet=303219672
       isat=9105001
       ityp=20

* Scan arguments

      do i=1,iargc()
         call getarg(i,arg)
         if (arg(:4).eq.'sig=') then
            sigfile=arg(5:)
	 else if (datearg(arg,t0,t1,dt)) then
	 else if (arg(:6).eq.'satid=') then
	    read (arg(7:),*) isat
	 else if (arg(:6).eq.'range=') then
	    read (arg(7:),*) range
         else
            filenm=arg
         endif
      enddo
      if (filenm.eq.' ') then
         write (0,10)
         goto 9999
      endif
10    format (
     |'GBFGEN -- Generate GBF file for all stations on SIGMA cards'//
     |'syntax: gbfgen [options] outfile'//
     |'where:'/
     |'outfile   : output GBF file'//
     |'and [options] are:'/
     |'t=t0,t1,dt: select period t0,t1 (sec85,mjd,yymmdd,yyddd)'/
     |'            and stepsize dt (seconds, def:15)'/
     |'          ... or use mjd=, doy=, ymd=, sec='/
     |'sig=file  : read SIGMA cards from "file" (required)'/
     |'satid=id  : use satellite ID "id" (def: 9105001)'/
     |'range=rng : set ranges to "rng" (m, m/s) (default: 0)')

* Read sigma cards

      nsta=0
      if (sigfile.eq.' ') then
         call fin("sig=file required")
      else
         open (10,file=sigfile,status='old')
 110     format (a)
 200     read (10,110,end=299) arg
         nsta=nsta+1
         read (arg,220) ista(nsta),sigma(nsta)
 220     format (10x,i4,16x,f20.7)
         goto 200
 299     continue
         close (10)
      endif

* Open output GBF

      open (20,file=filenm,form='unformatted')

* Create GBF record for each station mentioned in SIGMA file

      do t=nint(t0),nint(t1),nint(dt)
	 fract=t/86400d0
	 mjd=int(fract)
	 fract=fract-mjd
	 mjd=mjd+46066
         do i=1,nsta

* Write GBF line to output file

*	 write (*,*) ista(i),t
         if (ista(i)/1000.eq.4) then
	    ityp=92		! DORIS
	 else
	    ityp=20		! SLR
	 endif
         write(20) isat,ityp,indtim,ista(i),ipre,mjd,fract,range,
     |      idum1,idum2,sigma(i),trop,imet,cion,idum3,cmass

*	 call chrdat(nint((t-46066)*86400),arg)
*	 write (*,480) ista,sig,arg,t
*480	 format (i4,f8.3,2x,a15,f20.12)

	 enddo
      enddo

      close (20)

9999  end
