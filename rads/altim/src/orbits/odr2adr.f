      program odr2adr

      integer mx,my,nx,ny
      parameter (mx=2161,my=1081)
      integer getorb,t
      real*8 lat,lon,orb,h/0d0/
      real*8 xmin,xmax,ymin,ymax,zmin,zmax,x,y,dx/1d0/,dy/1d0/,
     |		grid(mx,my),gr3int
      real*8 t0,t1,dt/1d0/,lon0/-180d0/,lon1/180d0/,lat0/-90d0/,
     |		lat1/90d0/
      integer irec,narg,i,iargc
      character*80 geoname/' '/,adrname/' '/,arg,
     |		orbname/'/user/altim/data/ODR.ERS-2/dgm-e04'/
      integer*4 adr4(5)
      integer*2 adr2(2),bound(4)/-180,180,-90,90/
      logical geo,datearg

      narg=iargc()

      do i=1,narg
	 call getarg(i,arg)
	 if (arg(:4).eq.'orb=') then
	    orbname=arg(5:)
	 else if (arg(:4).eq.'geo=') then
	    geoname=arg(5:)
	 else if (arg(:4).eq.'lon=') then
	    read (arg(5:),*) lon0,lon1
	 else if (arg(:4).eq.'lat=') then
	    read (arg(5:),*) lat0,lat1
	 else if (datearg(arg,t0,t1,dt)) then
	 else
	    adrname=arg
	 endif
      enddo

      geo=geoname.ne.' '
      if (geo) then
	 nx=mx
	 ny=my
	 call gridrd8(geoname,nx,ny,grid,xmin,xmax,ymin,ymax,zmin,zmax)
	 dx=(xmax-xmin)/(nx-1)
	 dy=(ymax-ymin)/(ny-1)
      endif

      open (10,file=adrname,status='new',form='unformatted',
     |recl=24,access='direct')
      irec=1

      do t=nint(t0),nint(t1),nint(dt)
	 i=getorb(dble(t),lat,lon,orb,orbname,.true.)
	 if (i.gt.0) stop 'error in getorb'
	 if (lat.lt.lat0 .or. lat.gt.lat1) cycle
	 if (lon.lt.lon0) lon=lon+360d0
	 if (lon.gt.lon1) lon=lon-360d0
	 if (lon.lt.lon0) cycle
	 if (geo) then
	    x=(lon-xmin)/dx+1
	    y=(lat-ymin)/dy+1
	    h=gr3int(grid,x,y,nx,ny,mx)
	    if (abs(h).gt.1d20) h=0
	 endif
	 irec=irec+1
	 adr4(1)=t
	 adr4(2)=0
	 adr4(3)=nint(lat*1d6)
	 adr4(4)=nint(lon*1d6)
	 adr4(5)=nint(orb*1d3)
	 adr2(1)=nint(h*1d3)
	 write (10,rec=irec) adr4,adr2
      enddo

      write (10,rec=1) 'aADRERS-2   ',bound,irec-1
      end
