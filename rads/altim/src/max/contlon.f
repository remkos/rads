      function contlon(prog,lon0,lon1)
      logical contlon,prog
      real*8 lon0,lon1
      REAL*8 PI,TWOPI
      PARAMETER (PI=3.14159265358979D0,TWOPI=2*PI)

      if (prog) then
         if (lon1.lt.lon0) then
            lon1=lon1+twopi
            contlon=(lon1.lt.lon0)
	 else
	    contlon=.false.
         endif
      else
         if (lon1.gt.lon0) then
            lon1=lon1-twopi
            contlon=(lon1.gt.lon0)
	 else
	    contlon=.false.
         endif
      endif
      end
