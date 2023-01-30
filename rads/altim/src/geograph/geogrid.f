      program geogrid
      implicit none

* Program to create grid of gravity induced radial orbit errors

      include "geograph.inc"
      integer i,k,l,m,ip
      integer mx,kx,ky,nx,ny,narg,iargc
      parameter (mx=361*181)
      real*8 mean(mx),var(mx),a(npar),b(npar),d,h0,lat,lon,
     |	x0,x1,y0,y1,dx,dy,x,y,cycles,period,zmin,zmax,zmean,zrms
      character*80 filenm1/' '/,filenm2/' '/,gridnm/' '/,arg
      logical debug/.false./
      namelist /geograph_nml/ incl,period,cycles,deadband,lmax,test,
     |  x0,x1,y0,y1,dx,dy

* Initialise

      lmax=ndeg
      k=0

* Read namelists

      arg='/user/altim'
      call checkenv('ALTIM',arg,l)
      arg(l+1:)='/nml/geograph.nml'
      open (7,file=arg,status='old')
      read (7,geograph_nml)
      close (7)
      open (7,file="geograph.nml",status='old',err=5)
      read (7,geograph_nml)
      close (7)
5     continue

* Scan argument list

      narg=iargc()
      i=0
10    i=i+1
      call getarg(i,arg)
      if (arg(1:4).eq.'nml=') then
	 open (7,file=arg(5:))
	 read (7,geograph_nml)
	 close (7)
      else if (arg.eq.'-') then
	 i=i+1
	 call getarg(i,filenm2)
      else if (arg(1:4).eq.'clm=') then
	 read (arg(5:),*) l,m
	 k=1
	 lmax=l
      else if (arg(1:4).eq.'slm=') then
	 read (arg(5:),*) l,m
	 k=2
	 lmax=l
      else if (arg.eq.'-v' .or. arg.eq.'-debug') then
	 debug=.true.
      else if (arg.eq.'-test') then
	 test=.true.
      else if (arg(1:5).eq.'lmax=') then
	 read (arg(6:),*) lmax
      else if (filenm1.eq.' ') then
	 call getarg(i,filenm1)
      else
         call getarg(i,filenm2)
      endif
      if (i.lt.narg-1) goto 10
      call getarg(narg,gridnm)

* Read gravity parameters

      if (gridnm.eq.' '.or.filenm1.eq.' ') then
	 write (*,1300)
	 goto 9999
      endif
      call gravrd(1.0,filenm1)
      if (filenm2.ne.' ') call gravrd(-1.0,filenm2)
      if (debug) write (0,*) 'reading gravity completed'

* Do some checks

      nx=nint((x1-x0)/dx)+1
      ny=nint((y1-y0)/dy)+1
      if (nx*ny.gt.mx) stop 'too many gridpoints'

* Convert satellite orbit information

      incl=incl*rad
      n0=(2*pi)/(period*86400d0/cycles)
      wmdot=n0
      ogdot=-nint(period)*n0/cycles
      a0=(gm/n0**2)**(1d0/3d0)
      h0=a0-ae  
      
      if (debug) write (0,*) 'lmax =',lmax
      if (debug) write (0,*) 'Initialize Flmp and Dlmp'
      call d_lmp

      if (k.ne.0) then
         do ip=1,ipmax
	    if (ideg(ip).eq.l .and. iord(ip).eq.m .and. ics(ip).eq.k)
     |		goto 12
         enddo
	 stop 'coefficient not found.'
12       continue
         if (debug) write (0,*) 'ip,l,m,ics,cs = ',
     |		ip,ideg(ip),iord(ip),ics(ip),cs(ip)
      else
         ip=0
      endif
      k=0
      do ky=1,ny
	 lat=y0+(ky-1)*dy
	 call geocen(lat*rad,h0,y,d)		! convert to geocentric
         do kx=1,nx
	    k=k+1
	    lon=x0+(kx-1)*dx
	    x=lon*rad
	    if (ip.eq.0) then
               call geograph(y,x,mean(k),var(k))
	    else
               call geogrcmp(y,x,a,b)
	       mean(k)=a(ip)
	       var (k)=b(ip)
	    endif
	    if (debug) write (*,610) lon,lat,mean(k),var(k)
	 enddo
      enddo
      l=index(gridnm,' ')-1
      call gridwr8(gridnm(1:l)//'-c.grd',nx,ny,mean,nx,x0,x1,y0,y1)
      call gridrms(nx,ny,nx,mean,m,zmin,zmax,zmean,zrms)
      write (*,600) nx,ny,gridnm(1:l)//'-c.grd',zmin,zmax,zmean,zrms
      call gridwr8(gridnm(1:l)//'-s.grd',nx,ny,var ,nx,x0,x1,y0,y1)
      call gridrms(nx,ny,nx,var,m,zmin,zmax,zmean,zrms)
      write (*,600) nx,ny,gridnm(1:l)//'-s.grd',zmin,zmax,zmean,zrms
  600 format ('#',2i4,a15,4f12.4)
  610 format (2f8.3,2f12.4)
 1300 format(
     |'geogrid -- make grid of geographically correlated radial',
     |' orbit error'//
     |'usage: geogrid [options] <model1> [- <model2>] <prefix>'//
     |'where [ options ] are:'/
     |'lmax=lmax  maximum degree/order'/
     |'clm=l,m    only error due to Clm'/
     |'slm=l,m    only error due to Slm'/
     |'nml=name   Use namelist (on top of geograph.nml)'//
     |'and:'//
     |'model1     gravity model1 (total or covariances)'/
     |'model2     gravity model2 (use difference with model1)'/
     |'prefix     prefix for gridnames (extended with -c/s.grd')

 9999 end
