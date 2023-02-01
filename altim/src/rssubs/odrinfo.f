**ODRINFO -- Open ODR file and return some information
*+
      FUNCTION ODRINFO (UNIT, FILENM, SATEL, REP, ARC, VER, NREC,
     |			  TIME0, TIME1, TSTEP, BEGIN, END, REV)
      INTEGER*4 ODRINFO, UNIT, REP, ARC, VER, NREC,
     |		TIME0, TIME1, TSTEP, BEGIN, END
      REAL*8    REV
      CHARACTER FILENM*(*), SATEL*8
*
* This routine opens an ODR (Orbital Data Record) with filename FILENM,
* on unit UNIT and returns some information on the contents of the file
* Upon return the ODR file will NOT be closed, unless an error occurs.
*
* This routine can also be called as a function, in that case, ODRINFO
* will hold the value 0 when the old ODR format (specifier: @ODR) is
* read, and 1 when the new ODR format is encountered (specifier: xODR).
*
* Arguments:
*   UNIT   (output): File unit number.
*   FILENM  (input): Name of ODR file.
*   SATEL  (output): Satellite name.
*   REP    (output): Length of the repeat cycle in 10^-3 days.
*   ARC    (output): Arc number.
*   VER    (output): Orbit version number.
*   NREC   (output): Number of records (0 if file not found).
*   TIME0  (output): Time in UTC seconds of first record.
*   TIME1  (output): Time in UTC seconds of last record.
*   TSTEP  (output): Time-step in seconds.
*   BEGIN  (output): Begin of precise part of arc (in UTC seconds).
*   END    (output): End of precise part of arc (in UTC seconds).
*   REV    (output): Length of a revolution (in seconds).
*   ODRINFO(output): 0 = old format (@ODR), 1 = new format (xODR),
*                    -1 = can not find file, -2 = format error,
*                    -3 = unknown satellite
*-
* $Log: odrinfo.f,v $
* Revision 1.6  2006/01/10 14:19:15  rads
* - Removed FASTIO; reverted to FORTRAN I/O
*
* 13-Aug-2004 - Added longer default length of arc
*  5-Aug-2004 - Removed "unknown satellite"
*  6-Mar-2003 - Added CHAMP to the known satellites
* 20-Dec-2001 - Added extra argument to IOCONST because of new FASTIO
*  7-Apr-2000 - Prepared for CRYOSAT
* 22-Jul-1999 - Non-fatal exits.
* 20-Dec-1996 - Change in manual only.
* 12-Nov-1996 - Adjusted for new ODR version. Prepared for ENVISAT and JASON
*  6-Jul-1995 - ERS-2 168-day implemented. Handling of non-extended orbits included.
* 11-Mar-1993 - Created.
*-----------------------------------------------------------------------
      real*8	day,arclength,extension,endmax
      integer*4	nrev,l,lnblnk,freeunit,ios
      character*4 spec
      parameter (day=86400d0)
      logical*4	swap,ltlend
      save

* Open ODR and read headers

      nrec=0
      l=lnblnk(filenm)
      unit=freeunit()
      open (unit,file=filenm,status='old',access='direct',
     |		recl=16,iostat=ios)
      if (ios.ne.0) then
         write (0,550) 'odrinfo: can not open file ',filenm(:l)
	 odrinfo=-1
         return
      endif
      swap=ltlend()
      read (unit,rec=1) spec,satel,begin
      read (unit,rec=2) rep,arc,nrec,ver
      read (unit,rec=3) time0
      if (spec.eq.'@ODR') then
         odrinfo=0
      else if (spec.eq.'xODR') then
         odrinfo=1
      else
         write (0,550) 'odrinfo: wrong file type'
	 odrinfo=-2
	 close (unit)
	 return
      endif
      if (swap) then
         call i4swap(1,begin)
         call i4swap(1,rep)
         call i4swap(1,arc)
         call i4swap(1,nrec)
         call i4swap(1,ver)
         call i4swap(1,time0)
      endif
      read (unit,rec=nrec+2) time1
      if (swap) call i4swap(1,time1)
      tstep=nint(dble(time1-time0)/dble(nrec-1))

* Determine end of precise part of ODR

      if (begin.eq.0) then
	 begin=time0
	 end=time1
	 return
      else if (satel(1:4).eq.'ERS-' .or. satel(1:7).eq.'ENVISAT') then
         arclength=5.5d0
         extension=1d0
	 if (rep/1000.eq.168) then
	    rev=168*day/2357
         else if (rep/1000.eq.3) then
            rev= 3*day/43
         else
            rev=35*day/501
         endif
      else if (satel(1:3).eq.'T/P' .or. satel(1:5).eq.'JASON'
     |		.or. satel.eq.'TOPEX') then
         arclength=10.915636574d0
         extension=0.5d0
         rev=6745.759d0
      else if (satel.eq.'GEOSAT' .or. satel(1:3).eq.'GFO'
     |		.or. satel.eq.'SEASAT') then
         arclength=17.05
         extension=0d0
         rev=17.05*day/244
      else if (satel.eq.'CRYOSAT') then
         arclength=5.0d0
         extension=0.5d0
         rev=3*day/43
      else if (satel.eq.'CHAMP') then
         arclength=1.0d0
         extension=0.0d0
         rev=3*day/43
      else
	 rev=6035d0
	 arclength=10d0
	 extension=0d0
      endif

* ENDMAX is the latest moment the precise part stops. END is the actual end.

      endmax=time0+(arclength-extension)*day+rev/2
      nrev=int((endmax-begin)/rev)
      end=nint(begin+nrev*rev)

550   format (a,a)
      end
