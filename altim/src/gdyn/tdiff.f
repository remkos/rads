      FUNCTION TDIFF (TIME, DELTAT1, DELTAT2)
      REAL*8   TIME, DELTAT1, DELTAT2
      INTEGER  TDIFF

* This routine returns for a given TIME the correcponding time offsets between
* UTC(BIPM), UTC(USNO), and UTC(GPS). The routine initialises upon the first
* call. It need access to file times.data
*
* Arguments:
*  TIME     (input) : Time in MJD
*  DELTAT1 (output) : UTC(BIPM)-UTC(USNO) in seconds
*  DELTAT2 (output) : UTC(BIPM)-UTC(GPS) = -(leap sec.) + Co in seconds
*
* Exit code:
*  TDIFF   (output) : 0=no error, 1=before start BIPM values,
*                     2=before start GPS values, 3=after end BIPM/GPS values
*-
*  4-Aug-1998 - Remko Scharroo
*-----------------------------------------------------------------------
      integer*4 mutc,nutc/0/,ngps/0/,i,iold,j(2),freeunit
      parameter (mutc=1000)
      integer*4 tutc(mutc)
      real*8	xutc(mutc),xgps(mutc),x(3),frac
      character*80 name
      
      save xutc,xgps,tutc,iold,nutc,ngps
      
      if (nutc.eq.0) then
         iold=freeunit()
         call checkenv('ALTIM',name,i)
         name(i+1:)='/data/tables/times.data'
         open (iold,file=name,status='old')
         rewind (iold)
	 do i=1,3
	    read (iold,*,end=999)
	 enddo
10       read (iold,*,end=100) j,x
	 if (abs(x(1)).lt.1d-10.and.abs(x(2)).lt.1d-10) goto 10
	 nutc=nutc+1 
	 if (nutc.gt.mutc)
     |		call fin('tdiff: too many entries in times.data')
	 tutc(nutc)=dble(j(2))
	 xutc(nutc)=x(1)*1d-6
	 xgps(nutc)=x(3)*1d-6
	 if (abs(x(3)).gt.1d-10 .and. ngps.eq.0) ngps=nutc
         goto 10
100      close (iold)
         iold=1
      endif

* When before first or after last use first or last values.
* This is better than the original, because that used lineair extrapolation.

      if (time.lt.tutc(1)) then
	 deltat1=xutc(1)
	 deltat2=xgps(1)
	 iold=1
	 tdiff=1
	 return
      else if (time.gt.tutc(nutc)) then
	 deltat1=xutc(nutc)
	 deltat2=xgps(nutc)
	 iold=1
	 tdiff=3
	 return
      endif
	 
      do i=iold,nutc-1
         if (time.ge.tutc(i) .and. time.le.tutc(i+1)) goto 900
      enddo
      do i=1,iold-1
         if (time.ge.tutc(i) .and. time.le.tutc(i+1)) goto 900
      enddo

900   iold=i
      frac=(time-tutc(i))/(tutc(i+1)-tutc(i))
         deltat1=xutc(i)+frac*(xutc(i+1)-xutc(i))
      if (i.ge.ngps) then
         deltat2=xgps(i)+frac*(xgps(i+1)-xgps(i))
         tdiff=0	 
      else
	 deltat2=0d0
         tdiff=2
      endif
      return

999   stop 'tdiff: incorrect header in times.data'
      end
