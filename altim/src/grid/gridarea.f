      program gridinfo
*
* This is a small program to get some information on the contents of a
* grid.
*
      integer*4 nx,ny,i,iargc,nargs
      real*8    x0,x1,y0,y1,zmid,fact
      real*8    lon0,lon1,lat0,lat1
      character filenm*80,line*80,fspec*5

      nargs=iargc()
      if (nargs.eq.0) then
	 write (6,1300)
1300     format ('gridarea: change grid boundaries'//
     .'usage: gridarea lat=lat0,lat1 lon=lon0,lon1 gridnames')
	 goto 9999
      endif

      do i=1,nargs
         call getarg(i,filenm)
	 if (filenm(1:4).eq.'lat=') then
	    read (filenm(5:),*) lat0,lat1
	 else if (filenm(1:4).eq.'lon=') then
	    read (filenm(5:),*) lon0,lon1
	 else
         open (10,file=filenm,status='old',form='unformatted',
     .		recl=80,access='direct')
         read (10,rec=1) line
         if (line(1:5).ne.'@GRID') goto 10
	 read (line,600) fspec,zmid,fact,nx,ny,x0,x1,y0,y1
600      format (a5,2e14.7,2i7,1x,4f8.3)
	 write (6,610) filenm(1:40),x0,x1,y0,y1
610      format (a40,' : ',4f8.3)
         if (lon0.eq.lon1 .or. lat0.eq.lat1) then
	    write (line,600) fspec,zmid,fact,nx,ny,lon0,lon1,lat0,lat1
         else if
     &		(x0.ne.lon0.or.x1.ne.lon1.or.lat0.ne.y0.or.lat1.ne.y1) then
	    write (line,600) fspec,zmid,fact,nx,ny,lon0,lon1,lat0,lat1
	    write (6,610) 'Changed to',lon0,lon1,lat0,lat1
            write (10,rec=1) line
	 endif
         close (10)
	 endif
10	 continue
      enddo
9999  end
