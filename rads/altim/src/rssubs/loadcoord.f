**LOADCOORD -- Read station coordinates from file
*+
      SUBROUTINE LOADCOORD (FILENM, MJDREF, N, STA, POS, SIG)
      CHARACTER*(*) FILENM
      INTEGER MJDREF, N, STA(*)
      REAL*8  POS(6,*),SIG(6,*)

* This routine will load station coordinates from a variety of file
* formats, such as
* - XYZ format used by DUT/DEOS
* - LST format used by UT/CSR
* - SSC format used by ITRF
* - SINEX format
*
* Apart from the coordinates, the routine also attempts to read the
* velocities and the standard deviations of the station positions and the
* velocties.
* The X, Y, Z coordinates are stored in POS(1..3,I) and the
* X, Y, Z velocties are stored in POS(4..6,I), while STA(I) is the
* station number. The respective standard deviations (when available)
* are stored in SIG(1..6,ISTA).
*
* The XYZ format may contain the line 'NUVEL'. In this case the NUVEL1A
* velocities are calculated.
*
* Arguments:
*  FILENM  (input) : Name of the file in which coordinates are stored
*                    Use "-" for standard input.
*  MJDREF (output) : Reference epoch of the coordinates in MJD
*  N      (output) : Number of stations in the file
*  STA    (output) : Array containing the station numbers
*  POS    (putput) : Array containing the XYZ coordinates (in meters)
*                    and the XYZ velocities (in meters/year)
*  SIG    (output) : Standard deviations of the XYZ coordinates
*                    (in meters) and of the XYZ velocities (in
*                    meters/year), when available (otherwise 0)
*-
*  5-Apr-2001 - Created by Remko Scharroo
* 13-Apr-2001 - Allow reading from standard input
*-----------------------------------------------------------------------
      character	line*256,mmm*3,dum*32
      integer*4	i,yy,mm,dd,yymmdd,loadco1,format,unit,freeunit,mdate,
     |		plate
      logical	loadco2
      real*8	year,t
      parameter (year=86400*365.25d0)
      character months*36/'JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'/
      character string*6
      common	/cloadcoord/ unit

* Initialise

      n=0
      mjdref=0
      if (filenm.eq.'-') then
         unit=5
      else
         unit=freeunit()
         open (unit,file=filenm,status='old',err=900)
      endif

* Read first line and determine format

      if (loadco2(line)) goto 900
      if (line(:5).eq.'%=SNX') then
         format=1	! SINEX format
      else if (line(:13).eq.'  SLR STATION') then
         format=2	! LST format
      else if (index(line,'ITRF').ne.0) then
         format=3	! SSC format
      else if (line(:5).eq.'EPOCH') then
         format=4	! XYZ format
      else
         write (*,550) 'LOADCOORD: Sorry, unknown format'
	 return
      endif

* Read the station coordinates and the reference epoch from the
* different file formats

* SINEX format

      if (format.eq.1) then
	 i=loadco1(line,'+SOLUTION/ESTIMATE')
10       if (loadco2(line)) goto 900
         if (line.eq.'-SOLUTION/ESTIMATE') goto 900
	 i=index(line,':')-15
	 if (line(8:11).eq.'STAX') then
	    n=n+1
	    read (line(i:),510) sta(n),yy,dd,pos(1,n),sig(1,n)
	 else if (line(8:11).eq.'STAY') then
	    read (line(i:),510) sta(n),yy,dd,pos(2,n),sig(2,n)
	 else if (line(8:11).eq.'STAZ') then
	    read (line(i:),510) sta(n),yy,dd,pos(3,n),sig(3,n)
	 else if (line(8:11).eq.'VELX') then
	    read (line(i:),510) sta(n),yy,dd,pos(4,n),sig(4,n)
	 else if (line(8:11).eq.'VELY') then
	    read (line(i:),510) sta(n),yy,dd,pos(5,n),sig(5,n)
	 else if (line(8:11).eq.'VELZ') then
	    read (line(i:),510) sta(n),yy,dd,pos(6,n),sig(6,n)
	 endif
	 mjdref=mdate(2,yy*10000+0101)+(dd-1)	! Store reference epoch
         goto 10
510   format (i4,9x,i2,1x,i3,13x,d22.14,d12.5)
550   format (a)

* LST format

      else if (format.eq.2) then
	 i=loadco1(line,'EPOCH')		! Search for reference epoch
	 if (i.le.0) goto 900
	 read (line(i+6:),520) dd,mmm,yy
	 mm=index(months,mmm)+2/3
	 mjdref=mdate(2,yy*10000+mm*100+dd)

         if (loadco1(line,'X-DOT   Y-DOT   Z-DOT').le.0) goto 900
20       if (loadco2(line)) goto 900
	 if (line.eq.' ') goto 900
	 n=n+1
	 read (line,521) sta(n),(pos(i,n),i=1,6)
	 do i=4,6
	    pos(i,n)=pos(i,n)/1d3
	 enddo
	 do i=1,6
	    sig(i,n)=0d0
	 enddo
         goto 20
520   format (i2,1x,a3,1x,i4)
521   format (t19,i4,3f15.4,t81,3f8.2)

* SSC format

      else if (format.eq.3) then
         i=index(line,'EPOCH')		! Search for reference epoch
	 read (line(i+6:),'(i4)') yy
	 mjdref=mdate(2,yy*10000+0101)

	 if (loadco1(line,'--------').le.0) goto 900
	 if (loadco1(line,'--------').le.0) goto 900
30       if (loadco2(line)) goto 900
	 n=n+1
	 read (line,530) sta(n),(pos(i,n),i=1,3),(sig(i,n),i=1,3)
	 if (loadco2(line)) goto 900
	 read (line,531)        (pos(i,n),i=4,6),(sig(i,n),i=4,6)
	 goto 30
530   format (t34,i4,3f13.3,3f6.3)
531   format (t38,   3f13.4,3f6.4)

* XYZ format

      else if (format.eq.4) then
	 read (line,540) yymmdd
	 mjdref=mdate(2,yymmdd)
40       if (loadco2(line)) goto 900
	 n=n+1
	 read (line,541) string,sta(n),(pos(i,n),i=1,3)
	 if (loadco2(line)) goto 900
	 read (line,541) string,yy,(pos(i,n),i=4,6)
	 if (string.eq.'NUVEL ') then
	    call statinfo(1d30,sta(n),dum,dum,dum,dum,dum,plate,dum)
	    call nuvel1a(plate,pos(1,n),pos(4,n))
	    do i=4,6
	       pos(i,n)=pos(i,n)*year
	    enddo
	 endif
	 if (yy.ne.0) then
	    t=(mdate(2,yy*10000+0101)-mjdref)/365.25d0
	    do i=1,3
	       pos(i,n)=pos(i,n)-t*pos(i+3,n)
	    enddo
	 endif
	 goto 40
540   format (t15,i6)
541   format (a6,i4,3f15.4)

      endif

* End up here when all coordinates are read.

900   continue

* Close file

      if (unit.ne.5) close (unit)

* When EPOCH not found, use default value

      if (mjdref.eq.0) then
         write (*,550)
     |'LOADCOORD: Reference epoch not found, using 930101'
	 mjdref=mdate(2,19930101)
      endif

* Warning when no coordinates are found

      if (n.eq.0) write (*,550) 'LOADCOORD: No coordinates found'
      end

*--------- Supporting subroutines -----------------
* Search for the first line that contains string
      function loadco1(line,string)
      character line*(*),string*(*)
      logical loadco2
      integer loadco1
      loadco1=-1
10    if (loadco2(line)) return
      loadco1=index(line,string)
      if (loadco1.gt.0) return
      goto 10
      end
* Read line, but skip comments
      function loadco2(line)
      logical loadco2
      character*(*) line
      integer unit
      common /cloadcoord/ unit
      loadco2=.true.
10    read (unit,'(a)',end=9999) line
      if (line(:1).eq.'*' .or. line(:1).eq.'#') goto 10
      loadco2=.false.
9999  end
