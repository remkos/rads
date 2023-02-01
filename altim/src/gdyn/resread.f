**RESREAD -- Subroutine to read Geodyn II Residual File
*+
      SUBROUTINE RESREAD (FILENM, MAXOBS, MAXPAS,
     .   IPASS, CPASS, TIMES, RESID, SIGMA, DERIV, ELEV)
      CHARACTER*80 FILENM
      INTEGER*4   MAXOBS, MAXPAS, IPASS(4,MAXPAS)
      CHARACTER*8 CPASS(MAXPAS)
      REAL*8      TIMES(MAXOBS,MAXPAS), RESID(MAXOBS,MAXPAS),
     .   SIGMA(MAXOBS,MAXPAS), DERIV(MAXOBS,MAXPAS),
     .   ELEV(MAXOBS,MAXPAS)

* This routine loads the contents of a Geodyn II Residual File
* into a number of buffers.
*
* Arguments:
* FILENM (input) : Name of residual file
* MAXOBS (input) : Maximum number of observations per pass
* MAXPAS (input) : Maximum number of passes
* IPASS (output) : Pass information:
*                  IPASS(1,I) Station ID for pass I
*                  IPASS(2,I) Satellite ID for pass I
*                  IPASS(3,I) Number of observations for pass I
*                  IPASS(4,I) Measurement type for pass I
* CPASS (output) : CPASS(I)   Stations name for pass I
* TIMES (output) : TIMES(J,I) Time of observation J for pass I
* RESID (output) : RESID(J,I) Observation residual J for pass I
* SIGMA (output) : SIGMA(J,I) Observation input sigma
* DERIV (output) : DERIV(J,I) Observation time derivative J for pass I
* ELEV  (output) : ELEV(J,I)  Elevation of observation J for pass I
*-
*  7-Jul-1995 -- New manual
*  8-Feb-1996 -- Time derivative included
*  3-Mar-1999 -- Adjusted to work with GIIE 9606 too.
*  8-May-2002 -- Eelco: corrected giie.lt.9600 bug for 0104 compatibility
* 17-Jun-2003 -- Removed sysdep.h
*-----------------------------------------------------------------------
      integer idsta(100)
      real*8  rec(20),resrec(6000)
      character chrsta(100)*8,stanm*8
      integer ip,nc,nsta,i,ifloc,nloc,nm,mtype,ista,isat,nobs,narc
      integer k1,k2,k,unit,freeunit
      real*8  t0,x,giie
      logical swap,ltlend

* For integer file format

      character*4 spec
      integer*4 irec4(4),npass,irec
      integer*2 irec2(8)
      equivalence (irec4,irec2)
      equivalence (irec4(3),stanm)
      real*8 t1985
      parameter (t1985=(46066-30000)*86400d0)

* Clear pass information

      do ip=1,maxpas
	 ipass(1,ip)=0
	 ipass(2,ip)=0
	 ipass(3,ip)=0
	 ipass(4,ip)=0
      enddo
      swap=ltlend()

* Open residual file and check filetype

      unit=freeunit()
      open (unit,file=filenm,form='unformatted',status='old',err=1301,
     |access='direct',recl=16)
      read (unit,rec=1) spec
      if (spec.eq.'@RES') goto 1000
      close (unit)
      open (unit,file=filenm,form='unformatted',status='old',err=1301)

* Read global header

      read (unit,end=1300) rec
      nc=nint(rec(1))
      narc=nint(rec(2))
      nsta=nint(rec(10))
      giie=rec(20)
      if (nsta.gt.100) call fin('resread: too many stations')

* Skip input card deck

      do i=1,nc
         read (unit,end=1300)
      enddo

* Store station abbrevations, add ALTIM and XOVER

      do i=1,nsta
         read (unit,end=1300) chrsta(i),x
         idsta(i)=nint(x)
      enddo
      idsta (nsta+1)=99
      chrsta(nsta+1)='ALTIM'
      idsta (nsta+2)=100
      chrsta(nsta+2)='XOVER'

* Arc Header

   90 read (unit,end=1300) rec
      if (rec(1).ne.-1d12) goto 1310
      ifloc=nint(rec(7))
  100 if (nint(rec(8)).ne.0) goto 1320

* Lengths record

  120 read (unit,end=1300) rec
      t0=rec(1)
      if (t0.eq.-1d12) goto 100
      if (t0.eq.+1d12) then
	 narc=narc-1
	 if (narc.le.0) return
	 goto 90
      endif
      if (nint(rec(5)).ne.1) goto 1320
      if (nint(rec(9)+rec(10)).gt.2) goto 1320
      nloc=nint(rec(6))
      nm=nint(rec(8))
      mtype=nint(rec(4))
      ista=nint(rec(11))
      isat=nint(rec(12))

      if (mtype.eq.99 .or. mtype.eq.100) ista=mtype

* Look for identical pass (not for altimeter data)

      do ip=1,maxpas
	 nobs=ipass(3,ip)
	 if (nobs.eq.0) then
	    goto 140
         else if (ipass(1,ip).eq.ista .and. ipass(2,ip).eq.isat
     |		.and. ipass(4,ip).eq.mtype
     |       .and. abs(t0+rec(2)-times(1,ip)).lt.3600) then
	    goto 200
         endif
      enddo
      goto 1330

* Store the information for the (new) pass

  140 continue
      do i=1,nsta+2
       	 if (idsta(i).eq.ista) cpass(ip)=chrsta(i)
      enddo
      ipass(1,ip)=ista
      ipass(2,ip)=isat
      ipass(4,ip)=mtype

* Location Data Records

  200 continue
*      write (*,*) 'ista,mtype,nm,nloc:',ista,mtype,nm,nloc

* Location data records are given when:
* 1)  IFLOC > 0   AND
* 2a) GIIE version > 9600  OR
* 2b) GIIE version < 9600 and measurement is not xover/altim
*
* The location data records are multiple,
* depending on the entries on the location data header record.
* E.g: for SLR the location data header record says: 11.0 12.0 11.0 0.0 0.0 ...
* then there are 3 location data records following:
* - 1) Station location at observation start: time in elapsed seconds from
*      pass starting MJDS; cosine and sine of right ascension of Greenwich;
*      inertial X, Y, Z coordinate; inertial X, Y, Z velocity
* - 2) Satellite location at observation reflection: time in elapsed seconds
*      from pass starting MJDS; inertial X, Y, Z coordinate; inertial
*      X, Y, Z velocity
* - 3) Station location at observation arrival: time in elapsed seconds from
*      pass starting MJDS; cosine and sine of right ascension of Greenwich;
*      inertial X, Y, Z coordinate; inertial X, Y, Z velocity
* For altimetry (GIIE version > 9600 only) the header record usually says: 2.0
* 0.0 ...., giving satellite location information only.

      if (ifloc.le.0) then
      else if ((mtype.eq.99 .or. mtype.eq.100) .and. giie.lt.9600
     |         .and. giie.gt.0900) then
      else
         read (unit,end=1300) rec
*	 write (*,*) 'lengths:',rec
	 do i=1,20
	    k1=nint(rec(i))
	    do k2=1,nloc
	       if (k1.ne.0) read (unit,end=1300)
*	       if (k1.eq.0) then
*	       else if (k1.eq.11) then
*	          read (unit,end=1300) (resrec(k),k=1,9*nm)
**	          write (*,*)'locdat:',(resrec(k),k=1,9*nm)
*	       else if (k1.eq.12 .or. k1.eq.2) then
*	          read (unit,end=1300) (resrec(k),k=1,7*nm)
**	          write (*,*)'locdat:',(resrec(k),k=1,7*nm)
*	       else
*	          read (unit,end=1300)
*	       endif
	    enddo
         enddo
      endif

* For crossovers we have to read another record.
* If this record contains a small number, it is the number of xovers for
* which we have residuals coming. If it is a large number we have another
* lengths record, so we go back to line 120.

      if (mtype.eq.100) then
         read(unit) rec(1)
*	 write (*,*) 'nm(xover):',rec(1)
         nm=nint(rec(1))
         if (nm.gt.1000) then
            backspace(unit)
            goto 120
         endif
      endif

* Residual data records

      read (unit,end=1300) (resrec(k),k=1,6*nm)
*      write (*,*) 'residuals:',(resrec(k),k=1,6*nm)
      nobs=ipass(3,ip)
      if (nobs.gt.maxobs) then
         write (0,*) 'resread: Too many observations per pass'
	 write (0,*) '         ',nobs,' >',maxobs
         return
      endif
      do i=1,nm
	 times(nobs+i,ip)=t0+resrec(i)
	 resid(nobs+i,ip)=resrec(nm+i)
	 sigma(nobs+i,ip)=resrec(2*nm+i)
	 deriv(nobs+i,ip)=resrec(3*nm+i)
	 elev(nobs+i,ip)=resrec(5*nm+i)
      enddo
      ipass(3,ip)=nobs+nm
      goto 120

* END reading GEODYN II residual file

1000  continue

* Process short integer binary residual file
* Read file header

      read (unit,rec=1) irec4
      if (swap) call i4swap(1,irec4(2))
      npass=irec4(2)
      if (npass.gt.maxpas) goto 1330

      irec=1
      do ip=1,npass

* Read pass information

         irec=irec+1
         read (unit,rec=irec,err=1300) irec4
	 if (swap) call i4swap(4,irec4)
	 ipass(1,ip)=irec4(1)
	 ipass(2,ip)=irec4(2)
	 ipass(3,ip)=irec4(3)
	 ipass(4,ip)=irec4(4)
	 irec=irec+1
	 read (unit,rec=irec,err=1300) irec4
	 cpass(ip)=stanm

* Read observations

         do i=1,ipass(3,ip)
            irec=irec+1
            read (unit,rec=irec,err=1300) irec4
	    if (swap) then
               call i4swap(2,irec4(1))
               call i2swap(4,irec2(5))
	    endif
            times(i,ip)=irec4(1)+irec4(2)/1d6+t1985
            resid(i,ip)=irec2(5)/1d3
            sigma(i,ip)=irec2(6)/1d3
            deriv(i,ip)=irec2(7)
            elev (i,ip)=irec2(8)/1d2
         enddo
      enddo
      goto 9999

 1300 write (0,*) 'resread: Premature EOF'
      stop

 1301 write (0,*) 'resread: Unable to open file'
      stop

 1310 write (0,*) 'resread: Unexpected record'
      stop

 1320 write (0,*) 'resread: Residual file too complicated'
      stop

 1330 write (0,*) 'resread: Too many passes'
      stop
 9999 close (unit)
      end
