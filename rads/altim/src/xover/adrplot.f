      program adrplot

* This is a simple plotting program for ADR files, particularly of daily
* ADR files. Within one plot there are to graphs:
* 1) On top a map with the location of the measurements
* 2) On the bottom a map with the sea level anomalies
* Each file is assigned a different colour. Maximum eight files can be
* plotted.
*-
* 13-Jul-2000 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      implicit none
      integer*4 nfile/0/,mfile,i,j,narg,iargc,mrec,l,lnblnk,mdate
      integer*4 isec,msec,lat,lon,orb
      integer*2 ssh,bound(4)
      parameter (mfile=6,mrec=86400)
      integer*4 nrec(mfile),ymd(mfile)
      character*80 dev/'/xs'/,file(mfile),arg,text
      character spec*4,satel(mfile)*8
      real*4 t(mrec,mfile),x(mrec),y(mrec),z(mrec,mfile)
      logical exist

* Read arguments first

      narg=iargc()

      do i=1,narg
	 call getarg (i,arg)
	 if (arg(:4).eq.'dev=') then
	    dev=arg(5:)
	 else
	    nfile=nfile+1
	    if (nfile.gt.mfile) call fin('Too many files')
	    file(nfile)=arg
	 endif
      enddo

* If arguments are missing, write syntax

      if (nfile.eq.0) then
	 write (0,1300)
	 goto 9999
      endif
1300  format ('adrplot: plot sea level anomalies from (daily)',
     |' ADR file'//
     |'syntax: adrplot [ options ] ADRfile ...'//
     |'where:'/
     |'  ADRfile: name(s) of ADR files to be plotted'/
     |'and [options] are:'/
     |'  dev=dev: specify PGPLOT device as "dev"')

* Open plot device and initialise colours

      call pgbeg(0,dev,1,2)
      call pgsvp(0.08,0.98,0.09,0.98)
      call pgsch(1.5)
      call pgscr(0,1.,1.,1.)
      call pgscr(1,0.,0.,0.)
      call pgscr(2,1.,0.,0.)
      call pgscr(3,0.,1.,0.)
      call pgscr(4,0.,0.,1.)
      call pgscr(5,0.,1.,1.)
      call pgscr(6,1.,0.,1.)
      call pgscr(7,1.,1.,0.)

* Plot map

      call pgpage
      call pmdef(1,0.,0.,0.)
      call pmswin(-180.,180.,-82.,82.)
      call pmbox('BCNST',0.0,0,'BCNST',0.0,0)
      call pmwdb('0.lnd',1,0)

* Open plot files, ignore non-existing ones, and store
* time, location and sea level anomaly. Plot locations.

      do i=1,nfile
	 nrec(i)=0
	 inquire(file=file(i),exist=exist)
	 if (.not.exist) goto 100
	 open (10,file=file(i),status='old',form='unformatted',
     |		recl=24,access='direct')
	 read (10,rec=1) spec,satel(i),bound,nrec(i)
	 if (spec.eq.'aADR') then
	 else if (spec.eq.'xADR') then
	    close (10)
	    open (10,file=file(i),status='old',form='unformatted',
     |		recl=28,access='direct')
	 else
	    call fin('Unrecognised format of input file')
	 endif
	 if (nrec(i).gt.mrec) call fin('Too many records in file')
	 do j=1,nrec(i)
	    read (10,rec=j+1) isec,msec,lat,lon,orb,ssh
	    x(j)=lon/1e6
	    y(j)=lat/1e6
	    z(j,i)=ssh/1e3
	    t(j,i)=(mod(isec,86400)+msec/1e6)/36e2
	 enddo
	 ymd(i)=mdate(1,isec/86400+46066)
	 call pgsci(i+1)
	 call pmconv(nrec(i),x,y)
	 call pgpt(nrec(i),x,y,-1)
100      continue
      enddo

* Plot lower graph with sea level anomalies

      call pgsci(1)
      call pgpage
      call pgswin(0.,24.,-0.4,1.2)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
      call pglab('time (hours)','sea level anomaly (m)',
     |'Sea level anomalies with respect to OSU MSS95')
      do i=1,nfile
	 if (nrec(nfile).gt.0) then
	    call pgsci(i+1)
	    l=lnblnk(satel(i))
	    write (text,'(a,": ",i6.6)') satel(i)(:l),ymd(i)
	    call pgmtxt('T',-1.-i,0.95,1.0,text)
	    call pgpt(nrec(i),t(1,i),z(1,i),-1)
	 endif
      enddo
      call pgend

* End

9999  end
