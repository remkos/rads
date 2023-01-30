      real*8 x(3),lat,lon,height,rad,r
      integer iarg,iargc
      character*80 arg
      logical lonlat

      data lat/1d30/,lon/1d30/,lonlat/.false./

      rad=atan(1d0)/45

      do iarg=1,iargc()
         call getarg(iarg,arg)
	 if (arg.eq.'-latlon') then
	    lonlat=.false.
	 else if (arg.eq.'-lonlat') then
	    lonlat=.true.
	 else if (lat.ge.1d20) then
	    read (arg,*) lat
	 else if (lon.ge.1d20) then
	    read (arg,*) lon
	 else
	    read (arg,*) height
	 endif
      enddo

      if (lon.lt.1d20) then
         if (lonlat) then
	    call geoxyz(lon*rad,lat*rad,height,x,r)
	 else
	    call geoxyz(lat*rad,lon*rad,height,x,r)
	 endif
	 write (*,100) x
      else
10	 read (5,*,end=200) lat,lon,height
         if (lonlat) then
	    call geoxyz(lon*rad,lat*rad,height,x,r)
	 else
	    call geoxyz(lat*rad,lon*rad,height,x,r)
	 endif
	 write (*,100) x
	 goto 10
      endif

100   format (3f15.6)
200   end
