      program hovmuller

* This program reads several grids and makes a longitudinal or latitudinal
* slice of the data, creating either a latitude-time or longitude-time
* grid.
*-
*  3-Nov-1999 - Added DATEARG function
*----------------------------------------------------------------------

      integer i,j,iargc,narg,nx,ny,mx,my,ngrid,mgrid,n,sel
      parameter (mgrid=1000000)
      real*4 lon0,lon1,lat0,lat1,z0,z1,lon,lat
      real*4 slice(mgrid),grid(mgrid)
      real*8 t0,t1,dum
      logical datearg
      character*80 arg

* Initialize

      narg=iargc()
      t0=1
      t1=0
      lat=0
      lon=1e30
      n=0
      sel=0
      ngrid=0

      if (narg.lt.2) then
         write (0,600)
         goto 9999
      endif

600   format ('HOVMULLER',
     |' -- Create longitude-time or latitude-time grids'//
     |'Syntax: hovmuller [ options ] inputgrids outputgrid'//
     |'Required fields:'/
     |' inputgrids : names of input longitude-latitude grids',
     |' (must have equal size)'/
     |' outputgrid : name of output grid'//
     |'Optional fields:'/
     |' lon=lon    : specify longitude along which to slice'/
     |' lat=lat    : specify latitude  along which to slice'/
     |' t=t0,t1    : specify time dimensions in output grid'/
     |'          ... or use mjd=, doy=, ymd=, sec='//
     |'Defaults:'/
     |' lon=0 t=0,NrOfGrids-1')

* Process command line arguments

      do i=1,narg-1
         call getarg(i,arg)
	 if (arg(:4).eq.'lon=') then
	    read (arg(5:),*) lon
	 else if (arg(:4).eq.'lat=') then
	    read (arg(5:),*) lat
	 else if (datearg(arg,t0,t1,dum)) then
	 else
	    ngrid=ngrid+1
	    nx=0
	    ny=mgrid
	    call gridrd4(arg,nx,ny,grid,lon0,lon1,lat0,lat1,z0,z1)
               if (nx*ny.eq.0)
     |		call fin("Error reading input grid")
	    if (ngrid.eq.1) then
	       mx=nx
	       my=ny
	       if (lon.lt.1e20) then
	          sel=nint((lon-lon0)/(lon1-lon0)*(nx-1))
		  sel=sel-nx+1
	       else
	          sel=nint((lat-lat0)/(lat1-lat0)*(ny-1))
		  sel=sel*nx
	       endif
	    else if (nx.ne.mx .or. ny.ne.my) then
	       write (*,*) arg
               write (*,*) "mx,my,nx,ny:",mx,my,nx,ny
      	       call fin("Input grids not of equal size")
            endif

* Copy slice to output grid

	    if (lon.lt.1e20) then
	       do j=1,ny
	          slice(n+j)=grid(sel+j*nx)
	       enddo
	       n=n+ny
	       if (n+ny.gt.mgrid)
     |		stop "hovmuller: too much input data"
	    else
	       do j=1,nx
	          slice(n+j)=grid(sel+j)
	       enddo
	       n=n+nx
	       if (n+nx.gt.mgrid)
     |		stop "hovmuller: too much input data"
            endif
	 endif
      enddo

* Dump output grid

      call getarg(narg,arg)
      if (t0.gt.t1) then
         z0=0
	 z1=ngrid-1
      else
	 z0=t0/86400d0/365.25d0+85
	 z1=t1/86400d0/365.25d0+85
      endif
      if (lon.lt.1e20) then
         call gridwr4(arg,my,ngrid,slice,my,lat0,lat1,z0,z1)
      else
         call gridwr4(arg,mx,ngrid,slice,mx,lon0,lon1,z0,z1)
      endif

9999  end
