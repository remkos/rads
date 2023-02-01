      real*8 x(3),lat,lon,height,rad,r
      integer iarg,iargc
      logical lonlat/.false./
      character*80 arg

      data x/1d30,1d30,1d30/

      rad=atan(1d0)/45

      do iarg=1,iargc()
         call getarg(iarg,arg)
	 if (arg.eq.'-latlon') then
	    lonlat=.false.
	 else if (arg.eq.'-lonlat') then
	    lonlat=.true.
	 else if (x(1).ge.1d20) then
	    read (arg,*) x(1)
	 else if (x(2).ge.1d20) then
	    read (arg,*) x(2)
	 else
	    read (arg,*) x(3)
	 endif
      enddo

      if (x(3).lt.1d20) then
	 call xyzgeo(x,r,lat,lon,height)
         if (lonlat) then
	    write (*,100) lon/rad,lat/rad,height
	 else
	    write (*,100) lat/rad,lon/rad,height
	 endif
      else
10       read (5,*,end=200) x
	 call xyzgeo(x,r,lat,lon,height)
         if (lonlat) then
	    write (*,100) lon/rad,lat/rad,height
	 else
	    write (*,100) lat/rad,lon/rad,height
	 endif
	 goto 10
      endif

100   format (2f15.9,f15.3)
200   end
