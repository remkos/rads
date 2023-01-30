      program odrdif3

* Three dimensional differences between individual ODRs or sets of
* ODRs.
*-
* 20-Dec-1996 - Created by Remko Scharroo
* 16-Jan-1997 - Add more checkings.
*  9-May-2000 - Introduced DATEARG
*  7-Feb-2000 - Allow dt smaller than step size in orbit file
*  5-Sep-2002 - Optional plot of ascii orbit difference output file
*----------------------------------------------------------------------

      implicit none
      character*80 file(2)/2*' '/,arg,dev/' '/,title/'-'/,xgf/' '/,
     |		table/' '/,asc/' '/
      character	satel*8,date0*15,date1*15
      integer*4	narg,iarg,iargc,fd,ios,i,j,n,jmax,nout,pgbeg,type/4/
      integer*4 rep,arc,ver,nrec,t10,t11,t12,t13,dt1,nstep/0/
      real*8	t,t0,t1,dt/0d0/,dtmax/0d0/,rev,ymdhms
      integer*4 itime, idate
      real*8    ymdtime
      integer*4 mrec,lnblnk
      integer*2 sig/10/
      parameter (mrec=800000)
      real*8    lon(0:mrec),lat(0:mrec),orb(0:mrec),lon2,lat2,orb2,
     |		xyz(0:mrec,8),acr(0:mrec,4),ti(0:mrec),
     |		wk1(mrec),wk2(mrec),x(0:mrec),
     |		axyz(3),cxyz(3),dxyz(3),rxyz(3),r2,rad,scaprd,
     |		hmean,hrms,hsigma,ofac,hifac,hifreq,prob,amp,dmax/1d40/,
     |		rmax/1d40/
      real*4	h(mrec),f(mrec),hmax/0d0/,amax/-1./,step/0./
      logical	verbose/.false./,exist,datearg

      equivalence (x(0),xyz(0,1)),(wk1,lat(0)),(wk2,lon(0))
      equivalence (h,xyz(0,2)),(f,xyz(mrec/2,2))
      equivalence (xyz(0,4),orb(0)),(xyz(0,5),acr(0,1))

      integer*4     minint4,maxint4
      parameter     (maxint4=2147483647,minint4=-maxint4-1)

* Initialize

      narg=iargc()
      rad=atan(1d0)/45
      ofac=4
      hifreq=3.0
      rev=35d0*86400d0/501
      t0=minint4
      t1=maxint4

      do iarg=1,narg
	 call getarg(iarg,arg)
	 if (arg(:4).eq.'dev=') then
	    dev=arg(5:)
	 else if (datearg(arg,t0,t1,dt)) then
	 else if (arg(:2).eq.'-v') then
	    verbose=.true.
	 else if (arg(:5).eq.'ofac=') then
	    read (arg(6:),*) ofac
	 else if (arg(:7).eq.'hifreq=') then
	    read (arg(8:),*) hifreq
	 else if (arg(:6).eq.'type=z') then
	    type=3
	 else if (arg(:6).eq.'type=h') then
	    type=4
	 else if (arg(:6).eq.'type=a') then
	    type=5
	 else if (arg(:6).eq.'type=c') then
	    type=6
	 else if (arg(:6).eq.'type=r') then
	    type=7
	 else if (arg(:6).eq.'type=t') then
	    type=8
	 else if (arg(:6).eq.'title=') then
	    title=arg(7:)
	 else if (arg(:5).eq.'amax=') then
	    read (arg(6:),*,iostat=ios) amax,step,nstep
	 else if (arg(:5).eq.'dmax=') then
	    read (arg(6:),*) dmax
	    dmax=dmax/100
	 else if (arg(:5).eq.'rmax=') then
	    read (arg(6:),*) rmax
	    rmax=rmax/100
	 else if (arg(:4).eq.'xgf=') then
	    xgf=arg(5:)
	 else if (arg(:4).eq.'asc=') then
	    asc=arg(5:)
	 else if (arg(:4).eq.'tab=') then
	    table=arg(5:)
	 else if (file(1).eq.' ') then
	    file(1)=arg
	 else
	    file(2)=arg
	 endif
      enddo

* Print usage when not sufficient arguments

      if (file(2).eq.' ') goto 1300

* Check if given names might be files

      do i=1,2
	 arg=file(i)
	 if (arg(:1).ne.'+') then
            j=lnblnk(arg)
	    arg(j+1:)='/arclist'
	    inquire(exist=exist,file=arg)
	    if (.not.exist) then
	       arg(j+1:)='/ODR.001'
	       inquire(exist=exist,file=arg)
	    endif
	    if (.not.exist .and. index(file(i),'ODR.').gt.0) then
	       arg='+'
	       arg(2:)=file(i)
	       file(i)=arg
	    endif
	 endif
      enddo

* Limit time range when ODR files are used

      do i=1,2
      if (file(i)(:1).eq.'+') then
	 call odrinfo(fd,file(i)(2:),satel,rep,arc,ver,nrec,t10,t11,
     |		dt1,t12,t13,rev)
	 dtmax=max(dtmax,dble(dt1))
	 t0=max(t0,dble(t10))
	 t1=min(t1,dble(t11))
      endif
      enddo

* If dt not specified, use largest step size of orbit files

      if (dtmax.eq.0) dtmax=60
      if (dt.eq.0) dt=dtmax

* When time range not specified:

      if (t1-t0.gt.1d40) then
	 write (0,"(
     |'odrdif3: Use t= to indicate time range.',
     |' Or use + to indicate separate ODR files.')")
	 goto 9999
      endif

* Check number of steps

      nrec=nint((t1-t0)/dt)
      dt=(t1-t0)/nrec
      if (nrec.gt.mrec) then
	 write (0,"(
     |'odrdif3: Too many time steps. Increase step size.')")
	 goto 9999
      endif

* Check memory size for periodogram

      hifac=hifreq*(2*dt/rev)
      if ((table.ne.' ' .or. dev.ne.' ')
     |		.and. mrec.lt.nint(ofac*hifac*nrec*4)) then
	 write (0,"(
     |'odrdif3: Memory too low. Shorten time interval or',
     |' ofac or hifreq.')")
	 goto 9999
      endif

* Read first ODR or directory of ODRs

      do i=0,nrec
         t=t0+dt*i
	 call getorb(t,lat(i),lon(i),orb(i),file(1),verbose)
	 call geoxyz(lat(i)*rad,lon(i)*rad,orb(i),rxyz,acr(i,3))
*	 call polcar(lat(i)*rad,lon(i)*rad,orb(i),rxyz)
	 do j=1,3
	    xyz(i,j)=rxyz(j)
	 enddo
      enddo

* Open XGF or ASC file for output

      if (xgf.ne.' ') then
	 open (10,file=xgf,status='new',form='unformatted',
     |		access='direct',recl=18)
      endif
      if (asc.ne.' ') then
	 open (11,file=asc,form='formatted')
      endif

* Read second ODR or directory of ODRs and compute differences

      n=0
      do i=0,nrec
         t=t0+dt*i
	 call getorb(t,lat2,lon2,orb2,file(2),verbose)
	 call geoxyz(lat2*rad,lon2*rad,orb2,rxyz,r2)
*	 call polcar(lat2*rad,lon2*rad,orb2,rxyz)

* Determine crosstrack direction

	 if (i.lt.nrec) then
	    do j=1,3
	       dxyz(j)=xyz(i+1,j)-xyz(i,j)
	    enddo
	 endif
	 call vecprd(rxyz,dxyz,cxyz)
	 call vecnrm(cxyz,cxyz)

* Determine alongtrack direction

	 call vecprd(cxyz,rxyz,axyz)
	 call vecnrm(axyz,axyz)

	 do j=1,3
	    rxyz(j)=xyz(i,j)-rxyz(j)
	    xyz(n,j)=rxyz(j)
	 enddo
	 ti(n)=dt*i
	 lat(n)=lat(i)-lat2
	 lon(n)=lon(i)-lon2
	 orb(n)=orb(i)-orb2
	 acr(n,1)=scaprd(rxyz,axyz)
	 acr(n,2)=scaprd(rxyz,cxyz)
	 acr(n,3)=acr(i,3)-r2
	 acr(n,4)=sqrt(xyz(n,1)**2+xyz(n,2)**2+xyz(n,3)**2)
	 if (acr(n,4).le.dmax .and. abs(orb(n)).le.rmax) then
	    if (xgf.ne.' ') write (10,rec=n+2)
     |		nint(t),nint(lat2*1d6),nint(lon2*1d6),nint(orb(n)*1d6),sig
            ymdtime = ymdhms(t)
            idate = int(ymdtime/1d6)
            itime = nint(ymdtime-1d6*idate)
            if (asc.ne.' ') write (11,'(i6.6,i6.6,1x,f7.2,1x,
     |                                f7.2,1x,3f10.5)') 
     |		idate,itime,lat2,lon2,acr(n,1),acr(n,2),acr(n,3)
	    n=n+1
	 else
*           call strf1985(date0,'%y%m%d %H:%M:%S',nint(t))
*	    write
*    |		(*,'(a,5f12.6)') date0,lat(n),lon(n),orb(n),lat2,lon2
	 endif
      enddo

      if (xgf.ne.' ') then
	 write (10,rec=1) '@XGF',n
         close (10)
      endif
      if (asc.ne.' ') then
         close (11)
      endif

      call strf1985(date0,'%y%m%d %H:%M:%S',nint(t0))
      call strf1985(date1,'%y%m%d %H:%M:%S',nint(t1))
      write (*,606)
      write (*,605) nrec,n,date0,date1
605   format ('Period         ',i12,i12,3x,a15,' - ',a15)
606   format (75('-'))
      write (*,606)
      call prstat(n,lat(0)  ,'Latitude  (deg)')
      call prstat(n,lon(0)  ,'Longitude (deg)')
      call prstat(n,orb(0)  ,'Altitude    (m)')
      write (*,606)
      call prstat(n,acr(0,1),'Along-track (m)')
      call prstat(n,acr(0,2),'Cross-track (m)')
      call prstat(n,acr(0,3),'Radial      (m)')
      call prstat(n,acr(0,4),'Total       (m)')
      write (*,606)
      call prstat(n,xyz(0,1),'X (ECF)     (m)')
      call prstat(n,xyz(0,2),'Y (ECF)     (m)')
      call prstat(n,xyz(0,3),'Z (ECF)     (m)')
      write (*,606)

      if (dev.eq.' ' .and. table.eq.' ') goto 9999

* If plot or table is requested, make periodogram

      call statis(n,xyz(0,type),hmean,hrms,hsigma)
      call spfper(ti(0),xyz(0,type),n,hmean,hsigma,
     |ofac,hifac,wk1,wk2,mrec,nout,jmax,prob)

      do i=1,nout
         f(i)=wk1(i)*rev
         amp=2*hsigma*sqrt(wk2(i)/n)
         h(i)=amp*100
         hmax=max(hmax,h(i))
      enddo

      if (table.ne.' ') then
         open(20,file=table)
	 write (20,610)
	 do i=1,nout
	    write (20,'(2f7.4)') f(i),h(i)
	 enddo
	 close (20)
      endif

      if (dev.eq.' ') goto 9999

* If plot is requested, run output through pgplot

      if (pgbeg(0,dev,1,1).ne.1)
     |stop "odrdiff: error opening plot device"
      if (title.eq.'-') then
	 i=index(file(1),' ')-1
	 title=file(1)(:i)//' - '//file(2)
      endif
      if (amax.lt.0) amax=hmax*1.05
      call pgswin(0.,sngl(hifreq),0.0,amax)
      call pgbox('ABCNSTI',0.0,0,'ABCNSTI',step,nstep)
      call pglab('frequency (cycl/rev)','amplitude (cm)',
     |title)
      do i=1,nout
	 if (h(i).gt.amax/500) then
	    call pgmove(f(i),0.0)
	    call pgdraw(f(i),h(i))
	 endif
      enddo
      call pgend

      goto 9999

1300  write (0,600)
600   format('odrdif3 - ODR difference in three directions'//
     |'syntax: odrdif3 [ options ] orbit1 orbit2'//
     |'where [ options ] are:'/
     |'t=t0,t1,dt  : Specify time interval (MJD,YYMMDD,SEC85)',
     |' and step size'/
     |'              (sec, def: 60)'/
     |'          ... or use mjd=, doy=, ymd=, sec='/
     |'-v          : verbose'/
     |'xgf=filename: Specify name for optional XGF output file'/
     |'asc=filename: Specify name for optional ASC output file'/
     |'tab=filename: Specify name for optional periodogram table'/
     |'dev=dev     : Specify plot device'/
     |'title=title : Specify plot title'/
     |'hifreq=f    : Highest frequency in plot (def: 3.0)'/
     |'ofac=ofac   : Oversampling factor (def: 4.0)'/
     |'amax=amax   : Highest amplitude in plot (cm)'/
     |'dmax=dmax   : Maximum total orbit difference (cm)'/
     |'rmax=rmax   : Maximum radial orbit difference (cm)'/
     |'type=x      : Specify ploting of component x, where x is'/
     |'              z (z)',
     |' a (along), c (cross), r (radial), t (total), h (alt, def)')
610   format('# Orbit spectrum'/
     |'# Column 1: frequency (cycles/rev)'/
     |'# Column 2: amplitude (cm)')
9999  end

      subroutine prstat(n,x,a)
      integer*4 n,i
      real*8 x(n),xmin,xmax,xmean,xrms,xsigma
      character*(*) a

      xmin=1d40
      xmax=-1d40
      xmean=0
      xrms=0
      xsigma=0

      do i=1,n
	 xmin=min(xmin,x(i))
	 xmax=max(xmax,x(i))
	 xmean=xmean+x(i)
	 xrms=xrms+x(i)**2
      enddo
      xmean=xmean/n
      xrms=sqrt(xrms/n)
      xsigma=sqrt(xrms**2-xmean**2)
      write (*,610) a,xmin,xmax,xmean,xrms,xsigma
610   format (a,5f12.6)
      end
