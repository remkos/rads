      program xxo2asc
*
* Convert XXO data set to ASCII format, which can be used as input
* to asc2gbx
*
* Remko Scharroo - 30 November 1993
* 03-Nov-1996 - Also (pseudo-) dual xovers
* 03-Nov-1999 - Added DATEARG function
* 10-May-2002 - Added satellite IDs for tp, jason1 and envisat (Eelco)
* 10-May-2002 - Reimplemented dual-satellite xovers (Remko)
*  5-Aug-2002 - More initialisations
*
      implicit none
      character*4 spec
      integer*4 itime(4),ssh(2),orb(2),nrec,lat,lon,naux
      integer*2 itrk(2)
      real*8    h,hmin/-1d40/,hmax/1d40/,dt,dtmin/-1d40/,dtmax/1d40/,
     |		tmin/-1d40/,tmax/1d40/,tlo/0d0/,thi/0d0/,sigma,sigma0/1d-1/
      integer*4 iarg,iargc,nargs,isat(2)/2*0/
      logical ocean/.false./,sland,t1in,t2in,t1ou,t2ou,datearg

      real*8    day,sec1,sec2,r1,r2,s1,s2,sref,t1,t2,dlat,dlon
      integer*4 mjd1,mjd2,mjd85,i,j,irec,lnblnk
      character*80 filenm
      parameter (mjd85=46066,day=86400d0)

      integer*4	nx,ny,maxgrd,grdidx
      parameter (maxgrd=360*180)
      real*4	var,ssvar(maxgrd),x0,x1,y0,y1,z0,z1,fx/1e0/,fy/1e0/,
     |		ssscale/1e0/
      logical	addvar/.false./

      integer*4	maxsat
      parameter (maxsat=40)
      character*16 satnm(maxsat)
      integer*4 satid(maxsat)
      real*8	rmean,rrms,rmin,rmax,dum
      integer*4 n
      common /cstat/ rmean,rrms,rmin,rmax,n
      namelist /xxo2asc_nml/ satid,satnm

* Initialise statistics

      n=0
      rrms=0
      rmean=0
      rmin=1d50
      rmax=-1d50
      do i=1,maxsat
	 satid(i)=0
      enddo

* Load namelist

      filenm="/user/altim"
      call checkenv('ALTIM',filenm,i)
      filenm(i+1:)='/nml/xxo2asc.nml'
      open (10,file=filenm,status='old')
      read (10,nml=xxo2asc_nml)
      close (10)

* Scan command line arguments

      nargs=iargc()
      if (nargs.lt.1) goto 1300
      do iarg=1,nargs
	 call getarg(iarg,filenm)
         if (filenm(1:3).eq.'dt=') then
            read (filenm(4:),*) dtmin,dtmax
         else if (filenm(1:2).eq.'h=') then
            read (filenm(3:),*) hmin,hmax
	    hmin=hmin/1d2
	    hmax=hmax/1d2
	 else if (filenm(1:1).ge.'A' .and. filenm(:1).le.'Z'
     |	    .and. datearg(filenm,tlo,thi,dum)) then
	 else if (datearg(filenm,tmin,tmax,dum)) then
	 else if (filenm(1:2).eq.'-o') then
	    ocean=.true.
	 else if (filenm(1:7).eq.'sigma0=') then
	    read (filenm(8:),*) sigma0
	    sigma0=sigma0/100
	 else if (filenm(1:6).eq.'ssvar=') then
	    i=index(filenm,'*')
	    if (i.gt.0) then
	       read (filenm(7:i-1),*) ssscale
	    else
	       i=6
	    endif
	    nx=0
	    ny=maxgrd
	    call gridrd4(filenm(i+1:),nx,ny,ssvar,x0,x1,y0,y1,z0,z1)
	    fx=(x1-x0)/(nx-1)
	    fy=(y1-y0)/(ny-1)
	    addvar=.true.
         else

* Open XXO file

      open (10,file=filenm,status='old',form='unformatted',
     |access='direct',recl=44)
      read (10,rec=1) spec,nrec,naux
      if (spec.eq.'@XXO') then
      else if (spec.eq.'xXXO') then
         close (10)
         open (10,file=filenm,recl=44+naux*4,access='direct',
     |          status='old',form='unformatted')
      else
         call fin('xxo2asc: input file not XXO')
      endif

* Determine satellite IDs.
* First satellite is the one with which the filename starts.
* Second satellite (in case of duals) follows after _ or -.
* Matches in other positions of the filename are ignored.

      isat(1)=0
      isat(2)=0
      do i=1,maxsat
	 if (satid(i).ne.0) then
	    j=index(filenm,satnm(i)(:lnblnk(satnm(i))))
	    if (j.le.0) then
	    else if (j.eq.1) then
	       isat(1)=satid(i)
	    else if (filenm(j-1:j-1).eq.'/') then
	       isat(1)=satid(i)
	    else if (filenm(j-1:j-1).eq.'_') then
	       isat(2)=satid(i)
	    else if (filenm(j-1:j-1).eq.'-') then
	       isat(2)=satid(i)
	    endif
         endif
      enddo
      if (isat(1).eq.0) then
         write (*,*) 'No satellite names recognised. File skipped'
	 goto 101
      endif
      if (isat(2).eq.0) isat(2)=isat(1)

* Process records

      do irec=1,nrec
         read (10,rec=irec+1) lat,lon,itime,itrk,ssh,orb
         dt=abs(itime(1)-itime(3))/day
         if (dt.lt.dtmin.or.dt.gt.dtmax) goto 100
	 h=(ssh(1)-ssh(2))/1d6
         if (h.lt.hmin.or.h.gt.hmax) goto 100
	 dlat=lat/1d6
	 dlon=lon/1d6
	 if (ocean.and.sland(sngl(dlat),sngl(dlon))) goto 100

	 t1=itime(1)+itime(2)/1d6
         mjd1=itime(1)/86400
	 sec1=t1-mjd1*day
         r1=orb(1)/1d3
	 s1=ssh(1)/1d6

	 t2=itime(3)+itime(4)/1d6
         mjd2=itime(3)/86400
	 sec2=t2-mjd2*day
         r2=orb(2)/1d3
	 s2=ssh(2)/1d6

	 t1in=(t1.ge.tmin .and. t1.le.tmax)
	 t2in=(t2.ge.tmin .and. t2.le.tmax)
	 if (tlo.gt.0) then

* For pseudo-duals:
* - Exclude data that are actually dual xovers
* - Take only data with one component in tmin-tmax and
*   one in tlo-thi

	    if (t1in.and.t2in) goto 100
	    t1ou=(t1.ge.tlo .and. t1.le.thi)
	    t2ou=(t2.ge.tlo .and. t2.le.thi)
	    if (.not.((t1in.and.t2ou).or.(t2in.and.t1ou))) goto 100
	 else

* For singles and duals
* - Exclude data outside tmin-tmax

	    if (.not.(t1in.and.t2in)) goto 100
	 endif

* Use emperical sigma function: sigma = 10 - 2 * 1.2^(-dt)
* (sigma in cm, dt in days)
* Or sigma = sqrt ( sigma0**2 + var**2 )

	 if (addvar) then
	    grdidx=nint((dlat-y0)/fy)*nx+nint((dlon-x0)/fx)+1
	    var=ssvar(grdidx)*ssscale
	    if (var.ge.1e20) goto 100
	    sigma=sqrt(sigma0**2+var**2)
	    if (ssscale.lt.0 .or. sigma0.lt.0) sigma=-sigma
	 else
	    sigma= sigma0 * (1d0 - 0.2d0 * (1.2d0)**(-dt))
	 endif

* Give the 2 tracks of dual satellite xovers the proper satid.
* j is the pointer to the j-th satellite ID for the first pass.

	 if (itrk(1).lt.itrk(2)) then
	    j=1
	 else
	    j=2
	 endif

* Write out the ASCII output

	 if (t1in.and.t2in) then
         write (*,1601)
     |		mjd1+mjd85,sec1,dlat,dlon,r1-s1,2*n+1,sigma,isat(j)
         write (*,1601)
     |		mjd2+mjd85,sec2,dlat,dlon,r2-s2,2*n+2,sigma,isat(3-j)
	 call addstat(s1-s2)
	 else if (t1in) then
	 sref=s2
         write (*,1601)
     |		mjd1+mjd85,sec1,dlat,dlon,r1-s1+sref,0,sigma,isat(j)
	 call addstat(s1-sref)
	 else
	 sref=s1
         write (*,1601)
     |		mjd2+mjd85,sec2,dlat,dlon,r2-s2+sref,0,sigma,isat(3-j)
	 call addstat(s2-sref)
	 endif
 1601 format (i5,f13.6,2f12.6,f12.3,i9,f6.3,i8)

  100 continue
      enddo
  101 close (10)

* End processing XXO file

	 endif
      enddo
      
      goto 900
 1300 write (0,1301)
 1301 format ('usage: xxo2asc [options] xxo-file(s)'//
     |'where [options] are:'/
     |'-o             : select deep ocean xovers only'/
     |'dt=dtmin,dtmax : select xovers with dtmin < dt < dtmax'/
     |'t=tmin,tmax    : select xovers for period (tmin,tmax)'/
     |'             ... or use mjd=, doy=, ymd=, sec='/
     |'T=tmin,tmax    : period for pseudo-dual xovers (tmin,tmax)'/
     |'             ... or use MJD=, DOY=, YMD=, SEC='/
     |'h=hmin,hmax    : specify range xover difference (cm)'/
     |'sigma0=sigma0  : base sigma for xovers (cm)'/
     |'ssvar=gridname : specify sea surface variability grid')
      goto 999

  900 rmean=rmean/n
      rrms=sqrt(rrms/n-rmean**2)
      write (0,910) n,rmin,rmax,rmean,rrms
  910 format (
     |i6,' records converted xxo -> asc'/
     |7x,'Min =',f8.3,'  Max =',f8.3,'  Mean =',f8.3,'  Var =',f8.3)
  999 end

      subroutine addstat(r)
      real*8 r
      real*8 rmean,rrms,rmin,rmax
      integer n
      common /cstat/ rmean,rrms,rmin,rmax,n
      rmin=min(r,rmin)
      rmax=max(r,rmax)
      rmean=rmean+r
      rrms=rrms+r*r
      n=n+1
      end

