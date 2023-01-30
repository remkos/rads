**GRID2GRID -- Program to transform grids to new dimensions
*+
      program grid2grid

* This program reads a DEOS grid file and reinterpolates it at
* a user defined raster. It handles wrap-around margins too.
*
*-
* 10-Feb-2003 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      implicit none
      character*80 gridin/' '/,gridout/' '/,arg
      integer*4 nx,ny,maxgrd,i,iargc,ios,k,kx,ky,grid1,gridbuff,inter/2/
      parameter (maxgrd=2000000)
      real*8	grid2(maxgrd)
      real*8    x0/-180d0/,x1/180d0/,y0/-90d0/,y1/90d0/,
     |		dx/1d0/,dy/1d0/,x,y,gridbinq,gridbint,gridbspl,lon0,lon1
      logical	wrap/.true./

* Read arguments

      do i=1,iargc()
	 call getarg(i,arg)
	 if (arg(:4).eq.'lon=') then
	    read (arg(5:),*,iostat=ios) x0,x1,dx
	 else if (arg(:4).eq.'lat=') then
	    read (arg(5:),*,iostat=ios) y0,y1,dy
	 else if (arg(:2).eq.'-q') then
	    inter=0
	 else if (arg(:2).eq.'-l') then
	    inter=1
	 else if (arg(:2).eq.'-s') then
	    inter=2
	 else if (arg(:2).eq.'-w') then
	    wrap=.true.
	 else if (arg(:2).eq.'-n') then
	    wrap=.false.
	 else if (gridin.eq.' ') then
	    gridin=arg
	 else
	    gridout=arg
	 endif
      enddo

* If output grid is not specified, print usage message

      if (gridout.eq.' ') then
         write (*,1300)
	 goto 9999
      endif

* Load input grid

      if (gridbuff(gridin,grid1).ne.0)
     |		call fin ('Error reading input grid')

* Check gridsize

      call gridbinf(grid1,nx,ny,lon0,lon1,x,x,x,x)
      nx=nint((x1-x0)/dx+1)
      ny=nint((y1-y0)/dy+1)
      if (nx*ny.gt.maxgrd) call fin ('Output grid too large')
      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)

* Cycle through grid

      k=0
      do ky=1,ny
         y=y0+(ky-1)*dy
	 do kx=1,nx
	    x=x0+(kx-1)*dx
	    k=k+1
	    if (wrap) then
	       if (x.gt.lon1) x=x-360d0
	       if (x.lt.lon0) x=x+360d0
	    endif
	    if (inter.eq.2) then
	       grid2(k)=gridbspl(grid1,x,y)
	    else if (inter.eq.1) then
	       grid2(k)=gridbint(grid1,x,y)
	    else
	       grid2(k)=gridbinq(grid1,x,y)
	    endif
	 enddo
      enddo

* Write output grid

      call gridwr8(gridout,nx,ny,grid2,nx,x0,x1,y0,y1)
      call dallocf(grid1)

* Formats

1300  format('grid2grid - Grid manipulation program'//
     |'syntax: grid2grid [options] inputgrid outputgrid'//
     |'where:'//
     |' inputgrid    : name of the source grid'/
     |' outputgrid   : name of the target grid'//
     |'and [options] are:'//
     |' lon=x0,x1,dx : select longitude range and step size',
     |' (def: -180,180,1)'/
     |' lat=y0,y1,dy : select latitude range and step size',
     |' (def: -90,90,1)'/
     |' -w           : allow wrapping of the longitude coordinate',
     |' (default)'/
     |' -n           : do not wrap the longitude coordinate'/
     |' -s           : use spline interpolation (default)'/
     |' -l           : use linear interpolation'/
     |' -q           : use no interpolation (closest point)')
9999  end
