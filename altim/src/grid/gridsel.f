      program gridsel

* Program to select a subsection of a grid (without changing the scale).
*-
*  5-Feb-2002 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      character*80 line,filein/' '/,fileout/' '/
      integer*4 fdin,fdout,ios,openf,readf,writef,seekf,
     |		nx,mx,ix0,ix1,ny,my,iy0,iy1,i,maxbuf,iargc,iwrap/0/
      parameter (maxbuf=10000)
      integer*2 buffer(maxbuf)
      logical	wrap/.false./
      real*8	zmid,fact,xmin,xmax,ymin,ymax,dx,dy,
     |		lon0/1d0/,lon1/-1d0/,lat0/1d0/,lat1/-1d0/

* Read arguments

      do i=1,iargc()
         call getarg(i,line)
	 if (line(:4).eq.'lon=') then
	    read (line(5:),*) lon0,lon1
	 else if (line(:4).eq.'lat=') then
	    read (line(5:),*) lat0,lat1
	 else if (line(:2).eq.'-w') then
	    wrap=.true.
	 else if (filein.eq.' ') then
	    filein=line
	 else
	    fileout=line
	 endif
      enddo

* If argument list is incomplete, write usage information

      if (fileout.eq.' ') then
         write (*,1300)
	 goto 9999
      endif
1300  format ('gridsel -- Select portion of grid'//
     |'usage: gridsel lon=lon0,lon1 lat=lat0,lat1 gridin gridout')
	 
* Open grid files. Check if file exists.

      fdin=openf(filein,'r')
      if (fdin.lt.0) goto 1310
      fdout=openf(fileout,'w')
      if (fdout.lt.0) goto 1311

* Read header

      ios=readf(fdin,80,line)
      if (ios.ne.80) goto 1320

* Interpret grid header

      if (line(1:5).ne.'@GRID') goto 1320
      read (line,600,err=1320) zmid,fact,nx,ny,xmin,xmax,ymin,ymax
600   format (5x,2d14.7,2i7,1x,4f8.3)
601   format ('@GRID',2d14.7,2i7,1x,4f8.3)
550   format (a)

* Use default boundaries when they are not specified

      if (lon0.gt.lon1) then
         lon0=xmin
	 lon1=xmax
      endif
      if (lat0.gt.lat1) then
         lat0=ymin
	 lat1=ymax
      endif

* Determine latitude index boundaries

      dy=(ymax-ymin)/(ny-1)
      iy0=nint((lat0-ymin)/dy+1)
      iy1=nint((lat1-ymin)/dy+1)
      if (iy0.lt.1 .or. iy1.gt.ny) goto 1330
      my=iy1-iy0+1

* Determine longitude index boundaries

      dx=(xmax-xmin)/(nx-1)
      ix0=nint((lon0-xmin)/dx+1)
      ix1=nint((lon1-xmin)/dx+1)
      if (ix0.lt.1 .and. wrap) then
	 iwrap=nint(360d0/dx)
      else if (ix1.gt.nx .and. wrap) then
         iwrap=-nint(360d0/dx)
      endif
      if (ix0+iwrap.lt.1 .or. ix1+iwrap.gt.nx) goto 1331

      mx=ix1-ix0+1
      if (mx.gt.maxbuf) goto 1340

* Round of boundaries

      lon0=xmin+(ix0-1)*dx
      lon1=xmin+(ix1-1)*dx
      lat0=ymin+(iy0-1)*dy
      lat1=ymin+(iy1-1)*dy

* Write new grid header

      write (line,601) zmid,fact,mx,my,lon0,lon1,lat0,lat1
      ios=writef(fdout,80,line)
      if (ios.lt.80) goto 1321

* Jump to first grid line

      ios=seekf(fdin,((iy0-1)*nx+(ix0+iwrap-1))*2,1)
      if (ios.lt.0) goto 1320

* Read all required grid lines and write out selected part

      do i=iy0,iy1
         ios=readf(fdin,mx*2,buffer)
	 if (ios.ne.mx*2) goto 1320
	 ios=writef(fdout,mx*2,buffer)
	 if (ios.ne.mx*2) goto 1320
	 if (i.lt.iy1 .and. mx.lt.nx) then
	    ios=seekf(fdin,(nx-mx)*2,1)
	    if (ios.lt.0) goto 1320
         endif
      enddo

* Close files

      call closef(fdin)
      call closef(fdout)
      goto 9999
1310  write (*,550) 'Error opening input file'
      goto 9999
1311  write (*,550) 'Error opening output file'
      goto 9999
1320  write (*,550) 'Error reading input grid'
      goto 9999
1321  write (*,550) 'Error writing to output grid file'
      goto 9999
1330  write (*,550) 'Latitude boundaries beyond grid boundaries'
      goto 9999
1331  write (*,550) 'Longitude boundaries beyond grid boundaries'
      goto 9999
1340  write (*,550) 'Grid line too long'
9999  end
