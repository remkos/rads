      program xxo2alt
c
c rewrite XXO file to ASCII file, time sorted
c time in MJD, simulate altimeter measurement
c (opposed to Xovers in 'xxo2asc')
c
c usage: xxo2alt XXO-file ASCII-file Begin-(MJD) End-(MJD)
c
      implicit real*8(a-h,o-z)

      parameter(mjd85=46066,day=86400d0)

      integer*4 lat,lon,it1,ifract1,it2,ifract2
      integer*4 ssh1,ssh2,iorb1,iorb2,nrxo
      integer*2 trnr1,trnr2

      character*4 spec
      character*60 input,output,chartim1,chartim2
      
      call getarg(1,input)
      call getarg(2,output)
      call getarg(3,chartim1)
      call getarg(4,chartim2)
      
      if (chartim2.eq.' ') call fin
     |(' usage: xxo2alt XXO-file ASCII-file Begin-(MJD) End-(MJD)')

      iout=0
      if (output.eq.'-') iout=1

      read(chartim1,*) rmjd1
      read(chartim2,*) rmjd2

      open(unit=10,status='old',form='unformatted',
     c  access='direct',recl=44,file=input)
      if (iout.eq.0) then
      open(unit=20,status='unknown',form='formatted',file=output)
      endif
      
      rms=0d0
      nobs=0

      read(10,rec=1) spec,nrxo
      
      do 10 ixo=1,nrxo
      read(10,rec=ixo+1) lat,lon,it1,ifract1,it2,ifract2,
     c trnr1,trnr2,ssh1,ssh2,iorb1,iorb2
      dh=(ssh1-ssh2)*1d-4
      
c     compute orbital height (m)
      r1=iorb1/1d3-ssh1/1d6
      r2=iorb2/1d3-ssh2/1d6

c     compute longitude and latitude (degrees)
      dlon=lon*1d-6
      dlat=lat*1d-6

c     specify dummy geoid value (m)
      geoid=0d0

c     compute time in MJD and fraction of day in seconds
      t1=it1+ifract1*1d-6
      mjd1=t1/day
      sec1=t1-mjd1*day
      t2=it2+ifract2*1d-6
      mjd2=t2/day
      sec2=t2-mjd2*day

      day1=mjd85+mjd1+sec1/day
      day2=mjd85+mjd2+sec2/day

      nobs=nobs+1
      rms=rms+dh**2
      
c     write altimeter measurements 
      idummy=0
      if (day1.gt.rmjd1. and. day1.lt.rmjd2) then
      if (iout.eq.0) then
      write(20,100) mjd1+mjd85,sec1,dlat,dlon,r1,geoid,dh,idummy
      else
      write(6,100) mjd1+mjd85,sec1,dlat,dlon,r1,geoid,dh,idummy
      endif
      endif
      if (day2.gt.rmjd1. and. day2.lt.rmjd2) then
      if (iout.eq.0) then
      write(20,100) mjd2+mjd85,sec2,dlat,dlon,r2,geoid,dh,idummy
      else
      write(6,100) mjd2+mjd85,sec2,dlat,dlon,r2,geoid,dh,idummy
      endif
      endif
      
  100 format(i5,f12.5,2f12.6,f12.3,2f9.3,i7)

   10 continue

      if (iout.eq.0) then
      write(6,*) 'number of crossovers: ',nobs
      write(6,*) 'rms (cm): ',dsqrt(rms/nobs)
      endif
      
      end
