      program makessh

* Program to generate an instantaneous T/P derived mean sea surface from
* the mean and annual cycles. When -ers option is provided, use
* ERS-2 derived model. Grids are stored in $ALTIM/data/world.
*
* Syntax: makessh date gridname
*
* 28-Oct-1996 - Created by Remko Scharroo
*  3-Nov-1999 - Added DATEARG function
* 17-Nov-2000 - Changed pathnames of grid files. Using checkenv.
* 21-Dec-2000 - Introduced ERS-2 model (default)
************************************************************************
      implicit none
      character*80 arg,nmean,nanc,nans,nsmc,nsms,gridnm,altimdir
      real*8 year,pi,rad,t0/-1d0/,t1/-1d0/,dum
      integer kx,ky,nx,ny,mx,my,iarg,iargc,l
      parameter (year=365.25d0*86400d0,mx=360,my=180)
      real*4 gmean(mx,my)
      real*4 ganc(mx,my),gans(mx,my),gsmc(mx,my),gsms(mx,my)
      real*4 grid(-2:mx+3,my)
      real*4 xmin,xmax,ymin,ymax,zmin,zmax
      real*4 sin1t0,cos1t0,sin2t0,cos2t0
      real*4 sin1t1,cos1t1,sin2t1,cos2t1,dt
      logical nomean/.false./,datearg

* Assign grid names

      altimdir='/user/altim'
      call checkenv('ALTIM',altimdir,l)
      nmean=altimdir(:l)//'/data/world/ers2/mean_mss95.grd'
      nanc =altimdir(:l)//'/data/world/ers2/annual_c.grd'
      nans =altimdir(:l)//'/data/world/ers2/annual_s.grd'
      nsmc =altimdir(:l)//'/data/world/ers2/semian_c.grd'
      nsms =altimdir(:l)//'/data/world/ers2/semian_s.grd'
      pi=4*atan(1d0)
      rad=pi/180

* Scan arguments

      do iarg=1,iargc()
	 call getarg(iarg,arg)
	 if (datearg(arg,t0,t1,dum)) then
	 else if (arg(:4).eq.'-nom') then
	    nomean=.true.
	 else if (arg(:3).eq.'-tp') then
      nmean=altimdir(:l)//'/data/world/topex/mean_egm96.grd'
      nanc =altimdir(:l)//'/data/world/topex/annual_c.grd'
      nans =altimdir(:l)//'/data/world/topex/annual_s.grd'
      nsmc =altimdir(:l)//'/data/world/topex/semian_c.grd'
      nsms =altimdir(:l)//'/data/world/topex/semian_s.grd'
	 else if (arg(:5).eq.'-ers1') then
      nmean=altimdir(:l)//'/data/world/ers1/mean_mss95.grd'
      nanc =altimdir(:l)//'/data/world/ers1/annual_c.grd'
      nans =altimdir(:l)//'/data/world/ers1/annual_s.grd'
      nsmc =altimdir(:l)//'/data/world/ers1/semian_c.grd'
      nsms =altimdir(:l)//'/data/world/ers1/semian_s.grd'
	 else if (arg(:4).eq.'-ers') then
      nmean=altimdir(:l)//'/data/world/ers2/mean_mss95.grd'
      nanc =altimdir(:l)//'/data/world/ers2/annual_c.grd'
      nans =altimdir(:l)//'/data/world/ers2/annual_s.grd'
      nsmc =altimdir(:l)//'/data/world/ers2/semian_c.grd'
      nsms =altimdir(:l)//'/data/world/ers2/semian_s.grd'
	 else
	    gridnm=arg
	 endif
      enddo

      if (gridnm.eq.' ' .or. t0.lt.0) then
	 write (0,600)
	 stop
      endif
600   format(
     |'makessh: make SSH grid for a given epoch, based on'/
     |'         annual and semi-annual cycle observed from altimetry'//
     |'Syntax: makessh [ options ] t=t0[,t1] gridname'//
     |'where:'/
     |'gridname : name of the output grid of SSH'/
     |'t=t0     : compute SSH at date t0'/
     |'t=t0,t1  : compute integral of SSH over period t0,t1'/
     |'       ... or use mjd=, doy=, ymd=, sec='//
     |'options:'/
     |'-tp      : use model computed from T/P, reference is EGM96+SST'/
     |'-ers1    : use model computed from ERS-1, reference is OSU MSS95'/
     |'-ers2    : use model computed from ERS-2, reference is OSU MSS95',
     |' (def)'/
     |'-nomean  : do not include the mean')

* Read grids

      nx=mx
      ny=my
      call gridrd4(nmean,nx,ny,gmean,xmin,xmax,ymin,ymax,zmin,zmax)
      call gridrd4(nanc ,nx,ny,ganc ,xmin,xmax,ymin,ymax,zmin,zmax)
      call gridrd4(nans ,nx,ny,gans ,xmin,xmax,ymin,ymax,zmin,zmax)
      call gridrd4(nsmc ,nx,ny,gsmc ,xmin,xmax,ymin,ymax,zmin,zmax)
      call gridrd4(nsms ,nx,ny,gsms ,xmin,xmax,ymin,ymax,zmin,zmax)

* Blank out mean if requested

      if (nomean) then
      do ky=1,ny
	 do kx=1,nx
	    if (gmean(kx,ky).lt.1d20) gmean(kx,ky)=0
	 enddo
      enddo
      endif

* Transform t0 to SEC85 and then to PHASE in the year

      t0=mod(t0/year,1d0)*2*pi
      cos1t0=cos(  t0)
      sin1t0=sin(  t0)
      cos2t0=cos(2*t0)
      sin2t0=sin(2*t0)

* If T1 is given, transform it to SEC85 and PHASE in the year

      if (t1.gt.0) then
         t1=mod(t1/year,1d0)*2*pi
         cos1t1=cos(  t1)
         sin1t1=sin(  t1)
         cos2t1=cos(2*t1)
         sin2t1=sin(2*t1)
	 dt=t1-t0
	 
* Compute integral

      do ky=1,ny
	 do kx=1,nx
	    grid(kx,ky)=gmean(kx,ky)
     |		+ganc(kx,ky)*(sin1t1-sin1t0)/dt
     |		-gans(kx,ky)*(cos1t1-cos1t0)/dt
     |		+gsmc(kx,ky)*(sin2t1-sin2t0)/dt/2
     |		-gsms(kx,ky)*(cos2t1-cos2t0)/dt/2
	 enddo
      enddo

* If T1 is not given, compute SSH at given moment T0

      else

      do ky=1,ny
	 do kx=1,nx
	    grid(kx,ky)=gmean(kx,ky)
     |		+ganc(kx,ky)*cos1t0+gans(kx,ky)*sin1t0
     |		+gsmc(kx,ky)*cos2t0+gsms(kx,ky)*sin2t0
	 enddo
      enddo

      endif

* Fill edges of grid

      do ky=1,ny
	 do kx=1,3
	    grid(kx+nx,ky)=grid(kx,ky)
	 enddo
	 do kx=nx-2,nx
	    grid(kx-nx,ky)=grid(kx,ky)
	 enddo
      enddo

* Write grid

      call gridwr4(gridnm,nx+6,ny,grid,mx+6,xmin-3.0,xmax+3.0,ymin,ymax)
      end
