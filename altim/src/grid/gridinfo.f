**GRIDINFO -- Compute and print statistics of DEOS grids
*+
      program gridinfo
*
* This is a small program to get some information on the statistics
* per grid. This includes minimum, maximum, mean and rms.
*
*-
*        1995 - Created
*        1999 - Several upgrades
* 19-Apr-2000 - Improved manual
*  5-Jan-2004 - Bug in first COG output fixed
*-----------------------------------------------------------------------
      implicit none
      integer*4 mpnt,nx,ny,l,nv,narg,iarg,iargc,mx,mode
      integer*4	kx,kx0,kx1,ky,ky0,ky1,k
      parameter (mpnt=2881*1441,mx=50000)
      real*4    grid(mpnt),x0,x1,dx,y0,y1,dy,z0,z1,xa,xb,ya,yb
      real*4	lat0/90./,lat1/-90./,lon0/720./,lon1/-720./
      real*4    zmean,zrms,zsigma,xcog(mx+1),ycog(mx+1)
      character filenm*80

* Initialisation and checks

      mode=0
      narg=iargc()
      if (narg.eq.0) then
	 write (0,699)
	 goto 9999
      endif

699   format ('gridinfo - compute statistics of gridfile(s)'//
     |'usage: gridinfo [options] gridfile(s)'//
     |'with [options]:'/
     |' lon=x0,x1: specify longitude range (deg)'/
     |' lat=y0,y1: specify latitude range (deg)'/
     |'  -s[e][w]: single line statistics'/
     |'  -s      : no editing, no weighting'/
     |'  -se     : with editing, no weighting'/
     |'  -sw     : no editing, with weighting'/
     |'  -sew    : with editing, with weighting')

* Scan arguments, process grids directly

      do iarg=1,narg
      call getarg(iarg,filenm)
      if (filenm.eq.'-s') then
         mode=100
      else if (filenm.eq.'-se') then
         mode=110
      else if (filenm.eq.'-sw') then
         mode=101
      else if (filenm.eq.'-swe' .or. filenm.eq.'-sew') then
         mode=111
      else if (filenm(:4).eq.'lat=') then
         read (filenm(5:),*) lat0,lat1
      else if (filenm(:4).eq.'lon=') then
         read (filenm(5:),*) lon0,lon1
      else
      nx=0
      ny=mpnt
      call gridrd4(filenm,nx,ny,grid,x0,x1,y0,y1,z0,z1)
      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)
      if (lon0.gt.lon1) then
         kx0=1
	 kx1=nx
      else
         kx0=max(1,nint((lon0-x0)/dx+1))
         kx1=min(nx,nint((lon1-x0)/dx+1))
      endif
      if (lat0.gt.lat1) then
         ky0=1
	 ky1=ny
      else
         ky0=max(1,nint((lat0-y0)/dy+1))
         ky1=min(ny,nint((lat1-y0)/dy+1))
      endif
      kx=kx1-kx0+1
      ky=ky1-ky0+1
      k=kx0+(ky0-1)*nx
      xa=x0+(kx0-1)*dx
      xb=x0+(kx1-1)*dx
      ya=y0+(ky0-1)*dy
      yb=y0+(ky1-1)*dy

      l=index(filenm,' ')-1
      if (mode.eq.0) write (*,600)
     |	filenm(:l),kx,xa,xb,dx,ky,ya,yb,dy

* Compute statistics without editing, first without weighting

      if (mode.eq.0 .or. mode.eq.100) then
         call grstat4
     |	    (grid(k),kx,ky,nx,0.,0.,-1.,nv,z0,z1,zmean,zrms,zsigma)
         if (mode.eq.0) write (*,602) kx*ky,z0,z1
         if (mode.eq.0) write (*,605) nv,zmean,zrms,zsigma
         call gridcog4(grid(k),kx,ky,nx,-1e30,1e30,xcog,ycog)
         if (mode.eq.0) write (*,620)
     |	    xa+(xcog(ky+1)-1)*dx,ya+(ycog(kx+1)-1)*dy
         if (mode.eq.100) write (*,630)
     |      kx*ky,nv,z0,z1,zmean,zrms,zsigma,filenm(:l)
      endif

* Compute statistics with cos(latitude) weighting

      if (mode.eq.0 .or. mode.eq.101) then
         call grstat4
     |	    (grid(k),kx,ky,nx,ya,yb,-1.,nv,z0,z1,zmean,zrms,zsigma)
         if (mode.eq.0) write (*,610) zmean,zrms,zsigma
         if (mode.eq.101) write (*,630)
     |      kx*ky,nv,z0,z1,zmean,zrms,zsigma,filenm(:l)
      endif

* Compute statistics with editing, first without weighting

      if (mode.eq.0 .or. mode.eq.110) then
         call grstat4
     |	    (grid(k),kx,ky,nx,0.,0.,3.5,nv,z0,z1,zmean,zrms,zsigma)
         if (mode.eq.0) write (*,612) 3.5,z0,z1
         if (mode.eq.0) write (*,605) nv,zmean,zrms,zsigma
         call gridcog4(grid(k),kx,ky,nx,
     |	    zmean-3.5*zsigma,zmean+3.5*zsigma,xcog,ycog)
         if (mode.eq.0) write (*,620)
     |	    xa+(xcog(ky+1)-1)*dx,ya+(ycog(kx+1)-1)*dy
         if (mode.eq.110) write (*,630)
     |	    kx*ky,nv,z0,z1,zmean,zrms,zsigma,filenm(:l)
      endif

* Compute statistics with cos(latitude) weighting

      if (mode.eq.0 .or. mode.eq.111) then
         call grstat4
     |	    (grid(k),kx,ky,nx,ya,yb,3.5,nv,z0,z1,zmean,zrms,zsigma)
         if (mode.eq.0) write (*,610) zmean,zrms,zsigma
         if (mode.eq.111) write (*,630)
     |	    kx*ky,nv,z0,z1,zmean,zrms,zsigma,filenm(:l)
      endif

      endif
      enddo

600   format (/32('*'),' GRIDINFO ',33('*')/'    Grid : ',a/
     |'X-points : ',i10,t35,' X-range : ',3f10.3/
     |'Y-points : ',i10,t35,' Y-range : ',3f10.3)
602   format (/'   Total : ',i10,t35,' Min/Max : ',2f10.3)
605   format ('   Valid : ',i10,t29,'Mean/RMS/Sigma : ',3f10.3)
610   format (t20,'Weighted Mean/RMS/Sigma : ',3f10.3)
612   format (/'Edit-sgm : ',f10.3,t35,' Min/Max : ',2f10.3)
620   format (t35,' COG X/Y : ',2f10.3)
630   format (2i10,5f10.3,1x,a)

9999  end
