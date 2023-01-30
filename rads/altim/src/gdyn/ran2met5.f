c
c     conversion of prare data to geodyn II binary metric format
c
c Created by Pieter Visser
c Adapted for altim by Eelco Doornbos
c
      implicit double precision (a-h,o-z)
      parameter(secday=86400d0)
      parameter (freqtr0=8489d6)
      integer phasz2,grouz2,gropz2,if(10),it(2),ipreh(31)
      character reckey*8,if11*1,input*60,output*60
      dimension xutc(2000),tutc(2000)
      dimension mjd0(0:99)
      dimension rmout(19),imout(12)

      common/utc/nutc,icur,xutc,tutc

      call getarg(1,output)

      open(unit=20,status='unknown',form='unformatted',file=output)
c      open(unit=30,status='unknown',form='formatted',file='RANGE.ASC')
 
c     initialize statistics of corrections  
      nobs=0
      rmntro=0d0
      rmnion=0d0
      rmncom=0d0
      rmnphs=0d0
      rmngro=0d0
      rmngrp=0d0
      rmstro=0d0
      rmsion=0d0
      rmscom=0d0
      rmsphs=0d0
      rmsgro=0d0
      rmsgrp=0d0

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
      imout(7)=45000300

c speed of light
      vlight=2.99792458d8
c     rmout(17)=vlight-299792458d0
c fictitious wavelength (microns)
c      rmout(18)=8d9
      rmout(18)=vlight/freqtr0*1d6

c
c initialize multipath counter
      nmulti=0

c
c read prareran record
   60 read (5,100,end=999) reckey,nsatel,jahr,itag,sec,it,
     .numsta,icamp,inpt,nrang,azim,elev,range2,press,temp,wett,
     .ictro2,iono,icion2,massz2,phasz2,grouz2,gropz2,i91val,
     .if,if11,isdmes
  100 format (a8,i7,i2,i3,f12.7,2i1,i5,i4,i1,i3,2f5.2,f12.0,f5.1,2f4.1,
     .3i7,i6,2i4,i6,i5,10i1,a1,i7)

c     satellite id
      imout(1)=nsatel

c Check for Multipath and eliminate measurements if true
      if (if(8).eq.1) then
      nmulti=nmulti+1
      goto 60
      endif

c Check whether PRARE range measurement
      if (reckey.ne.'PRARERAN') goto 60

c Data edit flag
      if (if(9).eq.1) goto 60

c Station id.
c      imout(2)=numsta+2000
      imout(2)=numsta

c Correct time to satellite receiving time
      if (it(1).eq.1) sec=sec+range2*1d-12

c Observation in integral seconds after MJD 39855
      mjd=mjd0(jahr)+itag
      dfrac=sec/secday
      if (it(2).eq.5) then
      dfrac=dfrac-difft(mjd+dfrac)/1d6/secday
      endif
      timsec=(mjd+dfrac-39855)*secday
      imout(12)=int(timsec)
      rmout(1)=timsec-imout(12)

c Apply all possible corrections to PRARE range measurement
      if (mod(if(1),2).eq.1) range2=range2-ictro2  
      if (mod(if(2),2).eq.1) range2=range2-icion2
      if (mod(if(3),2).eq.1) range2=range2+massz2
      if (mod(if(4),2).eq.1) range2=range2+phasz2
      if (mod(if(5),2).eq.1) range2=range2+grouz2
      if (mod(if(6),2).eq.1) range2=range2+gropz2
c except ionospheric/tropospheric corrections
c      if (mod(if(1),2).ne.1) range2=range2+ictro2  
c      if (mod(if(2),2).ne.1) range2=range2+icion2
      range2=range2*1d-12*vlight/2d0
c Observation value
      rmout(3)=range2

c Update correction statistics
      rmul=1d-12*vlight/2d0*1d2
      nobs=nobs+1
      rmntro=rmntro+ictro2*rmul
      rmnion=rmnion+icion2*rmul
      rmncom=rmncom+massz2*rmul
      rmnphs=rmnphs+phasz2*rmul
      rmngro=rmngro+grouz2*rmul
      rmngrp=rmngrp+gropz2*rmul
      rmstro=rmstro+(ictro2*rmul)**2
      rmsion=rmsion+(icion2*rmul)**2
      rmscom=rmscom+(massz2*rmul)**2
      rmsphs=rmsphs+(phasz2*rmul)**2
      rmsgro=rmsgro+(grouz2*rmul)**2
      rmsgrp=rmsgrp+(gropz2*rmul)**2
c      write(99,199) ictro2*rmul,icion2*rmul,massz2*rmul,
c     c  phasz2*rmul,grouz2*rmul,gropz2*rmul
c  199 format(6f10.4)

c Sum of the troposheric refraction corrections
      rmout(12)=ictro2*1d-12*vlight/2d0
c Sum of the ionospheric refraction corrections
      rmout(13)=icion2*1d-12*vlight/2d0
c Standard deviation
      rmout(19)=isdmes*1d-12*vlight/2d0
 
      write(20) imout,rmout
c      write(30,'(2i5,2d15.8)') jahr,itag,sec,range2

      goto 60

  999 continue

c     write mean of correction 
      write(6,*) 'mean of tropo. corr. (cm)  : ',rmntro/nobs
      write(6,*) 'mean of ionos. corr. (cm)  : ',rmnion/nobs
      write(6,*) 'mean of c.o.m. corr. (cm)  : ',rmncom/nobs
      write(6,*) 'mean of phase  corr. (cm)  : ',rmnphs/nobs
      write(6,*) 'mean of grnd.  corr. (cm)  : ',rmngro/nobs
      write(6,*) 'mean of calib. corr. (cm)  : ',rmngrp/nobs
      write(6,*) 'rms of tropo. corr. (cm)  : ',dsqrt(rmstro/nobs)
      write(6,*) 'rms of ionos. corr. (cm)  : ',dsqrt(rmsion/nobs)
      write(6,*) 'rms of c.o.m. corr. (cm)  : ',dsqrt(rmscom/nobs)
      write(6,*) 'rms of phase  corr. (cm)  : ',dsqrt(rmsphs/nobs)
      write(6,*) 'rms of grnd.  corr. (cm)  : ',dsqrt(rmsgro/nobs)
      write(6,*) 'rms of calib. corr. (cm)  : ',dsqrt(rmsgrp/nobs)
      write(6,*)
      write(6,*) 'Nr. of obs. edited due to multipath: ',nmulti

      end

      double precision function difft(date)
c
c interpolate utc(bih)-utc(usno) from table
c
      implicit double precision (a-h,o-z)
      common/utc/nutc,icur,xutc(2000),tutc(2000)

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
      implicit double precision (a-h,o-z)
      dimension i(2),x(2)
      common/utc/nutc,icur,xutc(2000),tutc(2000)

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
