      program ladr
*
* This program makes a simple list of the contents of Abridged data
* records.
*-
*  7-Feb-1994 -- Created by Remko Scharroo
*  3-Nov-1999 -- Added DATEARG function
*  6-Mar-2002 -- Ported to Linux
*----------------------------------------------------------------------
      character filenm*80/' '/,arg*80,date*15,spec*4,satel*8
      integer*2 lon0,lon1,lat0,lat1,ssh,sigssh,swh,baro,orbdot,
     |adr2_2(2),adr2_4(4)
      integer*4 isec,msec,lat,lon,orb,adr4(5)
      integer*4 nrec,irec,iargc,iarg,l0,l1,ios
      logical*4 datearg	!,ltlend
      real*8 hfac,t0/-1d40/,t1/1d40/,dum
      data l0/0/,l1/0/
      equivalence (adr2_2(1),ssh),(adr2_2(2),sigssh)
      equivalence (adr2_4(1),ssh),(adr2_4(2),swh)
      equivalence (adr2_4(3),baro),(adr2_4(4),orbdot)
      equivalence (adr4(1),isec),(adr4(2),msec),(adr4(3),lat)
      equivalence (adr4(4),lon),(adr4(5),orb)

* Scan arguments

      do iarg=1,iargc()
         call getarg(iarg,arg)
         if (arg(1:2).eq.'l=') then
	    read(arg(3:),*,iostat=ios) l0,l1
	 else if (datearg(arg,t0,t1,dum)) then
	 else if (arg(1:1).ne.'-') then
	    filenm=arg
	 endif
      enddo

      if (filenm.eq.' ') goto 1300

* Open ADR file

      open (10,file=filenm,status='old',form='unformatted',
     |	access='direct',recl=24,err=1310)
      read (10,rec=1) spec,satel,lon0,lon1,lat0,lat1,nrec
*     if (ltlend()) then
*        call i2swap(1,lon0)
*        call i2swap(1,lon1)
*        call i2swap(1,lat0)
*        call i2swap(1,lat1)
*        call i4swap(1,nrec)
*     endif
      write (*,600) spec,satel,nrec,lon0,lon1,lat0,lat1

* Print header

      if (spec.eq.'@ADR') then
         hfac=1d2
	 write (*,550) ' SigSSH'
      else if (spec.eq.'aADR') then
         hfac=1d3
	 write (*,550) ' SigSSH'
      else if (spec.eq.'xADR') then
	 hfac=1d3
	 close (10)
         open (10,file=filenm,status='old',form='unformatted',
     |	access='direct',recl=28,err=1310)
	 write (*,550) '    SWH Invbar  Orbdot'
      else
         goto 1320
      endif
550   format (a)
600   format('Filetype : ',a4,t30,'Satellite: ',a8,t60,'Records:',i10
     |/      'Longitude: ',2i5,t30,'Latitude : ',2i5/
     |'YYMMDD HH:MM:SS.SSSSSS    Latitude    Longitude     ',
     |'Altitude RelSSH',$)

* Determine line interval
* - If l0 or l1 is negative, count it from the end
* - If l0 and l1 are zero, print all
* - If l1 is zero, print till the end

      if (l0.lt.0) l0=nrec+l0+1
      if (l1.lt.0) l1=nrec+l1+1
      if (l0.eq.0) then
	 l0=1
	 l1=nrec
      else if (l1.lt.l0) then
	 l1=nrec
      endif
      l0=max(1,l0)
      l1=min(nrec,l1)

* Print out records

      if (spec.eq.'xADR') then
      do irec=l0+1,l1+1
         read (10,rec=irec) isec,msec,lat,lon,orb,ssh,
     |		swh,baro,orbdot
*	 if (ltlend()) then
*	    call i4swap(5,adr4)
*	    call i2swap(4,adr2_4)
*	 endif
	 if (isec.ge.t0 .and. isec.le.t1) then
	    call strf1985(date,'%y%m%d %H:%M:%S',isec)
	    write (*,610) date,msec,lat/1d6,lon/1d6,orb/1d3,ssh/hfac,
     |		swh/1d3,baro/1d3,orbdot/1d3
         endif
      enddo
      else
      do irec=l0+1,l1+1
         read (10,rec=irec) isec,msec,lat,lon,orb,ssh,sigssh
*	 if (ltlend()) then
*	    call i4swap(5,adr4)
*	    call i2swap(2,adr2_2)
*	 endif
	 call strf1985(date,'%y%m%d %H:%M:%S',isec)
	 write (*,610) date,msec,lat/1d6,lon/1d6,orb/1d3,ssh/hfac,
     |		sigssh/1d3
610	 format (a15,'.',i6.6,f12.6,f13.6,f13.3,3f7.3,f8.3)
      enddo
      endif
      goto 9999

* Some error messages

1300  write (0,1301)
1301  format ('ladr: list contents of ADR file'//
     |'usage: ladr [ options ] filenm'//
     |'required:'/
     |'  filenm  : ADR file name (aADR and xADR headers allowed)'/
     |'options:'/
     |'  t=t0,t1 : Specify time interval ([yy]yymmdd[hhmmss],mjd)'/
     |'        ... or use mjd=, doy=, ymd=, sec='/
     |'  l=n     : Print records starting with record n'/
     |'  l=n,m   : Print records n through m'/
     |'          (use negative n or m to start counting from the end)')
      goto 9999

*1310  call perror('ladr')
1310  goto 9999

1320  write (0,1321) spec
1321  format ('ladr: unknown file type: ',a4)
9999  end
