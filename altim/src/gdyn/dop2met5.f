c
c     conversion of prare data to geodyn II binary metric format
c
c     Created by Pieter Visser
c     020307 - Adapted for altim by Eelco Doornbos
c
      real*8 secday,freqtr0,dfrac,rmul,timsec,vlight,difft,alpha,
     |		azim,cexcal,cion2,cmasc2,cmcgs2,cphsz2,offsat,
     |		press,rdop,dopptr,drrdn,elev,wett,ctro2,temp,
     |		sec,tinteg,zqmes,sdmes,freqtr
      integer*4 i,icur,ipre,nmulti,nobs,nonpra,nutc,mdate,itag,
     |		iono,jahr,inpt,ndop,nsatel,numsta,icamp,mjd
      parameter(secday=86400d0)
      parameter (freqtr0=8489d6,alpha=749d0/880)
      integer*4 if(9),it(2),ipreh(31),imout(12),mjd0(0:99)
      character reckey*8,if10*1,output*60
      real*8 xutc(2000),tutc(2000),rmout(19)
      real*8 rmntro,rmnion,rmncom,rmnphs,rmncal,rmngst
      real*8 rmstro,rmsion,rmscom,rmsphs,rmscal,rmsgst

      common/utc/nutc,icur,xutc,tutc

      call getarg(1,output)

      open(unit=20,status='unknown',form='unformatted',file=output)
 
c     initialize statistics of corrections  
      nobs=0
      rmntro=0d0
      rmnion=0d0
      rmncom=0d0
      rmnphs=0d0
      rmngst=0d0
      rmncal=0d0
      rmstro=0d0
      rmsion=0d0
      rmscom=0d0
      rmsphs=0d0
      rmsgst=0d0
      rmscal=0d0

c     initialize
      do 10 i=1,12
      imout(i)=0
   10 continue
      do 20 i=1,19
      rmout(i)=0d0
   20 continue

      icur=1
      nonpra=0

c initialize
c - utc(bih)-utc(usno) table difft
c - table with mjd for jan 0.0 for 1960 to 2059 (mjd0)
      call initdt
      do 30 i=0,99
      mjd0(i)=mdate(2,10000*i+0101)-1
   30 continue

c PREPRO word
      ipre=0
      do 33 i=1,31
      ipreh(i)=0
   33 continue
c Flags relating to c.g. corrections
      ipreh(7)=0
      ipreh(8)=0
      ipreh(9)=0
c Flags relating to antenna axis displacement corrections
      ipreh(10)=0
      ipreh(11)=0
      ipreh(12)=0
c Flags relating to troposheric refraction corrections
      ipreh(13)=1
      ipreh(14)=0
      ipreh(15)=0
c Flags relating to ionospheric refraction corrections
      ipreh(16)=1
      ipreh(17)=0
      ipreh(18)=0
c Flags relating to transponder delay corrections
      ipreh(19)=0
      ipreh(20)=0
      ipreh(21)=0
c Flags relating to relativistic corrections
      ipreh(22)=0
      ipreh(23)=0
c Flages relating to Doppler
      ipreh(30)=0
      ipreh(31)=0
      do 36 i=1,31
      if (ipreh(i).ne.0) then
      ipre=ipre+2**(i-1)
      endif
   36 continue
      write(6,*) 'PREPRO word: ',ipre
      imout(11)=ipre

c     measurement type/time system indicator
      imout(7)=46000300

c speed of light
      vlight=2.99792458d8
c     rmout(17)=vlight-299792458d0
c fictitious wavelength (microns)
c      rmout(18)=8d9
      rmout(18)=vlight/freqtr0*1d6

c
c Initialize multipath counter
      nmulti=0

c
c read praredop record
   60 read (5,100,end=999) reckey,nsatel,jahr,itag,sec,it,
     .numsta,icamp,azim,elev,tinteg,zqmes,offsat,inpt,ndop,
     .press,temp,wett,ctro2,iono,cion2,cmasc2,cphsz2,cmcgs2,
     .cexcal,if,if10,sdmes
  100 format (a8,i7,i2,i3,f12.7,2i1,i5,i4,2f5.2,f6.3,f12.3,f9.3,
     .i1,i3,f5.1,2f4.1,f8.3,i7,f7.3,f6.3,3f5.3,9i1,a1,f7.3)

c     satellite id
      imout(1)=nsatel

c Check for multipath
      if (if(8).eq.1) then
      nmulti=nmulti+1
      goto 60
      endif

c Check whether PRARE doppler measurement
      if (reckey.ne.'PRAREDOP') goto 60

c Data edit flag
      if (if(9).eq.1) goto 60

c Station id.
      imout(2)=numsta

c correct time to end of integration interval
      if (it(1).eq.1) sec=sec+tinteg
c Observation in integral seconds after MJD 39855
      mjd=mjd0(jahr)+itag
      dfrac=sec/secday
      if (it(2).eq.5) then
      dfrac=dfrac-difft(mjd+dfrac)/1d6/secday
      endif
      timsec=(mjd+dfrac-39855)*secday
      imout(12)=int(timsec)
      rmout(1)=timsec-imout(12)

c Doppler counting interval
      rmout(4)=tinteg

c     apply tropospheric and ionospheric correction
      if (mod(if(1),2).eq.1) zqmes=zqmes+ctro2
      if (mod(if(2),2).eq.1) zqmes=zqmes+cion2
c     apply centre of mass correction
      if (mod(if(3),2).eq.1) zqmes=zqmes+cmasc2
c     apply onboard prare antenna phase centre correction
      if (mod(if(4),2).eq.1) zqmes=zqmes+cphsz2
c     apply ground station range rate correction
      if (mod(if(5),2).eq.1) zqmes=zqmes+cmcgs2
c     apply external calibration correction
      if (mod(if(6),2).eq.1) zqmes=zqmes+cexcal
c except ionospheric/tropospheric corrections
c      if (mod(if(1),2).ne.1) zqmes=zqmes+ctro2
c      if (mod(if(2),2).ne.1) zqmes=zqmes+cion2

c Correct freqtr0
      freqtr=freqtr0+offsat
      dopptr=freqtr*tinteg
      rdop=vlight*zqmes/(alpha*freqtr)/tinteg/2d0
      rmout(3)=rdop

c write corrections (in mm/s)
      rmul=vlight/(alpha*freqtr)/tinteg/2d0*1d3
      nobs=nobs+1
      rmntro=rmntro+ctro2*rmul
      rmnion=rmnion+cion2*rmul
      rmncom=rmncom+cmasc2*rmul
      rmnphs=rmnphs+cphsz2*rmul
      rmngst=rmngst+cmcgs2*rmul
      rmncal=rmncal+cexcal*rmul
      rmstro=rmstro+(ctro2*rmul)**2
      rmsion=rmsion+(cion2*rmul)**2
      rmscom=rmscom+(cmasc2*rmul)**2
      rmsphs=rmsphs+(cphsz2*rmul)**2
      rmsgst=rmsgst+(cmcgs2*rmul)**2
      rmscal=rmscal+(cexcal*rmul)**2
c      write(99,199) ctro2*rmul,cion2*rmul,cmasc2*rmul,cphsz2*rmul,
c     c  cmcgs2*rmul,cexcal*rmul,press,temp
c  199 format(8f10.4)

c Sum of dry and wet tropospheric refraction correction provided 
c SIGN OF TROPO AND IONO CORRECTION REVERSED SINCE 970807
c      rmout(12)=ctro2*vlight/(alpha*freqtr)/tinteg/2d0
      rmout(12)=-ctro2*vlight/(alpha*freqtr)/tinteg/2d0
c Ionopheric refraction correction provided
c      rmout(13)=cion2*vlight/(alpha*freqtr)/tinteg/2d0
      rmout(13)=-cion2*vlight/(alpha*freqtr)/tinteg/2d0
c Standard deviation
      drrdn =vlight/2/alpha/dopptr
      rmout(19)=abs(drrdn*sdmes)/2d0
 
      write(20) imout,rmout

      goto 60

  999 continue

c     write mean of corrections
      write(6,*) 'mean of tropo. corr. (mm/s): ',rmntro/nobs
      write(6,*) 'mean of ionos. corr. (mm/s): ',rmnion/nobs
      write(6,*) 'mean of c.o.m. corr. (mm/s): ',rmncom/nobs
      write(6,*) 'mean of phase  corr. (mm/s): ',rmnphs/nobs
      write(6,*) 'mean of grnd.  corr. (mm/s): ',rmngst/nobs
      write(6,*) 'mean of calib. corr. (mm/s): ',rmncal/nobs
      write(6,*) 'rms of tropo. corr. (mm/s): ',dsqrt(rmstro/nobs)
      write(6,*) 'rms of ionos. corr. (mm/s): ',dsqrt(rmsion/nobs)
      write(6,*) 'rms of c.o.m. corr. (mm/s): ',dsqrt(rmscom/nobs)
      write(6,*) 'rms of phase  corr. (mm/s): ',dsqrt(rmsphs/nobs)
      write(6,*) 'rms of grnd.  corr. (mm/s): ',dsqrt(rmsgst/nobs)
      write(6,*) 'rms of calib. corr. (mm/s): ',dsqrt(rmscal/nobs)
      write(6,*)
      write(6,*) 'Nr. of obs. edited due to multipath: ',nmulti

      end

      real*8 function difft(date)
c
c interpolate utc(bih)-utc(usno) from table
c
      integer*4 nutc,icur,i
      real*8 xutc(2000),tutc(2000),date
      common/utc/nutc,icur,xutc,tutc

      if (date.gt.tutc(icur+1)) then
         do 10 i=icur+1,nutc-1
   10       if (date.le.tutc(i+1)) goto 90
         difft=xutc(nutc)
         icur=nutc-1
         return
      else if (date.lt.tutc(icur)) then
         do 20 i=icur-1,1,-1
   20       if (date.ge.tutc(i)) goto 90
         difft=xutc(1)
         icur=1
         return
      else
         goto 100
      endif

  90  icur=i
 100  difft=xutc(icur)+(xutc(icur+1)-xutc(icur))*(date-tutc(icur))/
     .  (tutc(icur+1)-tutc(icur))
      if (icur+1.gt.nutc) icur=icur-1
      end


      subroutine initdt
c
c read times.data file and store utc(bih)-utc(usno) values in a table
c
      integer*4 i(2),nutc,icur,j
      real*8 x(2),xutc(2000),tutc(2000)
      common/utc/nutc,icur,xutc,tutc

      open (80,file='/u2a/geodyn/bin/Source/UTCTAB',
     c  status='old',form='formatted')
      rewind (80)
      nutc=0
      do 20 j=1,3
         read (80,10,end=1300)
  10     format (a80)
  20  continue

  30  read (80,*,end=100) i,x
      if (dabs(x(1)).lt.1.0d-10.and.dabs(x(2)).lt.1.0d-10) goto 30
      nutc=nutc+1
      if (nutc.gt.1000) then
         write (6,*) 'the number of utc(bih)-utc(usno) values exceeds',
     .    ' 2000.'
         stop 99
      end if
      tutc(nutc)=i(2)
      xutc(nutc)=x(1)
      goto 30
  100 close (80)
      return

 1300 write (6,*) 'improper utc(bih)-utc(usno) file.'
      stop 99
      close (80)
      end

**MDATE -- Convert MJD to YYMMDD or v.v.
*+
      FUNCTION MDATE (I, J)
      INTEGER  MDATE, I, J
*
* This function converts Modified Julian Dates (MJD) to Year-Month-Day
* in the format YYMMDD, or vice versa.
* Dates in th 21st century are given in 6-digits, similar to dates in
* the 20th century. So 13 October 2001 is written as 011013.
* The algorithm is good for the years 1960 till 2059 only.
*
* Arguments:
*  I      (input): I=1, convert MJD to YYMMDD
*                  I=2, convert YYMMDD to MJD
*  J      (input): Input value, either MJD or YYMMDD (depending on I).
*  MDATE (output): Output value, either YYMMDD or MJD (depending on I).
*-
* 15-Sep-1988 - Created: Remko Scharroo, DUT/SOM (c), DUT/KOSG
* 14-Jan-1992 - Revised version
*-----------------------------------------------------------------------
      integer   t1901,cal(13,0:1),t,y,m,d,leap
      parameter (t1901=15384)

* T1901 is the MJD of January 0.0, 1901

      save cal
      data cal /0,31,59,90,120,151,181,212,243,273,304,334,365,
     |          0,31,60,91,121,152,182,213,244,274,305,335,366/
      leap=0
      if (i.eq.1.and.j.lt.73459.and.j.gt.36934) then
          t=j-t1901-1
          y=idint((t+8d-1)/365.25d0)
          t=t-dint(y*365.25d0)
          y=y+1
          if (y.ge.100) y=y-100
          if (mod(y,4).eq.0) leap=1

* Estimate and correct month-number, T=0 for january 1

          m=t/31+1
          if (t.ge.cal(m+1,leap)) m=m+1
          d=t-cal(m,leap)+1
          mdate=y*10000+m*100+d
      else if (i.eq.2) then
          d=mod(j,100)
          t=(j-d)/100
          if (j.lt.600000) t=t+10000
          m=mod(t,100)
          y=(t-m)/100
          if (mod(y,4).eq.0) leap=1
          mdate=t1901+dint((y-1)*365.25d0)+cal(m,leap)+d
      else
          mdate=j
      endif

      return
      end
