      program earthfit

* This program reads XGF or ADR altimeter data products that should be
* corrected for tides and mean sea surface topography. It then computes
* four parameters that constitute an altimeter bias and coordinate
* shift. These prameters are:
* - Delta R (radius correction), positive when sea surface is too high
* - Delta X, Delta Y, Delta Z (geocenter correction) moves geocenter
*   towards the sea surface height data
*-
* 28-Jul-2000 - Created by Remko Scharroo
* 29-Jul-2000 - SGI Complib routines replaced by LAPACK routines
* 22-Nov-2004 - Removed all ocean tide routines
*-----------------------------------------------------------------------
      integer*4 npar
      parameter (npar=4)
      real*8 h0/-1d30/,h1/1d30/,t0/-1d30/,t1/1d30/,s0/-1d30/,s1/1d30/
      real*8 lat0/-1d30/,lat1/1d30/,lon0/-1d30/,lon1/1d30/
      real*8 t,h,s,lat,lon,dum,w,pi,rad
      real*8 atwa(npar*(npar+1)/2),p(npar),atwb(npar)
      integer*4 hi(npar),ih,it1,it2,orb,iargc,narg,i,j,k,mode,iarg,
     |	ilat,ilon,info,irec,nrec,mrec/0/
      integer*2 js,jh,bounds(4)
      character spec*4,arg*80,satel*8,par(npar)*2/'dR','dX','dY','dZ'/
      logical datearg
      
* Initialise

      call matsy1(npar,hi)
      pi=4*atan(1d0)
      rad=pi/180
      do i=1,npar*(npar+1)/2
         atwa(i)=0
      enddo
      do i=1,npar
         atwb(i)=0
      enddo

* Check number of arguments. Should be at least one, otherwise error
* message is printed.

      narg=iargc()
      if (narg.lt.1) then
         write (*,600)
	 goto 9999
      endif
600   format ('earthfit -- Fit coordinate shift to XGF or ADR data'//
     |'syntax: earthfit [options] <file(s)>'//
     |'where'/
     |' file(s)      : name(s) of the ADR and/or XGF files'//
     |'and options are:'/
     |' lon=lon0,lon1: specify longitude boundaries'/
     |' lat=lat0,lat1: specify latitude  boundaries'/
     |' t=t0,t1      : specify time interval',
     |' (real:[yy]yymmdd[hhmmss],mjd)'/
     |'            ... or use mjd=, doy=, ymd=, sec='/
     |' h=h0,h1      : specify height interval (cm)'/
     |' s=s0,s1      : specify sigma interval (cm)')

* Walk through all arguments

      do iarg=1,narg
        call getarg(iarg,arg)
        if (arg(:2).eq.'h=') then
           read (arg(3:),*) h0,h1
	   h0=h0/100
	   h1=h1/100
        else if (arg(:2).eq.'s=') then
           read (arg(3:),*) s0,s1
	   s0=s0/100
	   s1=s1/100
        else if (arg(:4).eq.'lon=') then
           read (arg(5:),*) lon0,lon1
        else if (arg(:4).eq.'lat=') then
           read (arg(5:),*) lat0,lat1
        else if (datearg(arg,t0,t1,dum)) then
	else

* Open XGF and ADR input files one by one
* Verify file specification to decide wether input is XGF (mode=1) or
* ADR (mode=2)

	  open (10,file=arg,form='unformatted',status='old',
     |		recl=18,access='direct')
          read (10,rec=1) spec,nrec
	  if (spec.eq.'@XGF') then
	     mode=1
	  else if (spec.eq.'aADR') then
	     mode=2
	     close(10)
	     open (10,file=arg,form='unformatted',status='old',
     |		recl=24,access='direct')
             read (10,rec=1) spec,satel,bounds,nrec
	  else if (spec.eq.'xADR') then
	     mode=2
	     close(10)
	     open (10,file=arg,form='unformatted',status='old',
     |		recl=28,access='direct')
             read (10,rec=1) spec,satel,bounds,nrec
	  else
	     write (*,550) 'earthfit: unknown file type in file '//arg
	     goto 120
	  endif

* Start reading data

	  do irec=1,nrec
	     if (mode.eq.1) then
	        read (10,rec=irec+1) it1,ilat,ilon,ih,js
	        h=ih/1d6
		s=js/1d3
	     else
	        read (10,rec=irec+1) it1,it2,ilat,ilon,orb,jh,js
	        h=jh/1d3
		s=js/1d3
	     endif
	     s=1
	     w=1/(s*s)
	     t=it1+it2/1d6
	     lat=ilat/1d6
	     lon=ilon/1d6

* Skip data beyond selected interval of time, weight, height,
* latitude, or longitute

	     if (h.lt.h0 .or. h.gt.h1) goto 100
	     if (s.lt.s0 .or. s.gt.s1) goto 100
	     if (t.lt.t0 .or. t.gt.t1) goto 100
	     if (lon.lt.lon0 .or. lon.gt.lon1) goto 100
	     if (lat.lt.lat0 .or. lat.gt.lat1) goto 100

* Compute partials for this data point:
* p(1)=dh/dr, p(2)=dh/dx, p(3)=dh/dy, p(4)=dh/dz

	     p(1)=1
	     p(2)=cos(lat*rad)*cos(lon*rad)
	     p(3)=cos(lat*rad)*sin(lon*rad)
	     p(4)=sin(lat*rad)

* Add particals to normal matrix

	     do i=1,npar
	        do j=1,i
		   k=hi(i)+j
		   atwa(k)=atwa(k)+p(i)*w*p(j)
		enddo
		atwb(i)=atwb(i)+p(i)*w*h
	     enddo
	     mrec=mrec+1
100          continue
   	  enddo
	endif

* All measurements for this file are read, continue with next file

120     continue
      enddo

* AtWA matrix and right-hand side vector are filled.
* Invert the matrix and solve the parameters.

      call dpptrf('U',npar,atwa,info)			! factorisation
      if (info.ne.0) write (*,550) 'earthfit: factorisation impossible'
      call dpptri('U',npar,atwa,info)			! inverse
      if (info.ne.0) write (*,550) 'earthfit: inversion impossible'
      call dspmv('U',npar,1d0,atwa,atwb,1,0d0,p,1)	! solve

* Compute sigma and correlations

      do i=1,npar
	 k=hi(i)+i
         atwa(k)=sqrt(atwa(k))
      enddo
      do i=1,npar
         do j=1,i-1
	    k=hi(i)+j
	    atwa(k)=atwa(k)/atwa(hi(i)+i)/atwa(hi(j)+j)*100
	 enddo
      enddo

* Write results

      write (*,610) mrec,(par(i),p(i),atwa(hi(i)+i),i=1,npar),
     |(par(i),i=2,npar),(par(i),(atwa(hi(i)+j),j=i+1,npar),i=1,npar-1)
610   format (
     |'Number of measurements selected:',i10//
     |'Parameter    Value (m)   Sigma (m)'/
     |4(a9,2f12.6/)/
     |'Correl(%)',3a9/a9,3f9.3/a9,9x,2f9.3/a9,18x,f9.3)

550   format (a)
9999  end
