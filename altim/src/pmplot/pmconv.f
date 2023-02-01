**PMCONV -- convert X and Y coordinates of several points.
*+
      FUNCTION PMCONV (N, X, Y)
      INTEGER PMCONV, N
      REAL    X(N),Y(N)
*
* Convert real-world coordinates (longitude,latitude) to map coordinates (x,y).
* The real-world coordinates are the position coordinates in degrees,
* the map coordinates are mapped linearly in the viewport and
* depend on the projection type used. For equi-rectangular projection
* real-world and map coordinates are identical, and makes this conversion
* optional. For non-linear projections real-world coordinates must be
* changed into map coordinates by means of PMCONV before using PGLINE, PGPTXT,
* etc.
* Upon return PMCONV contains the number of points that do not have
* invalid values, occurring when a point is not convertable.
*
* Arguments:
*  N    (input) : Number of points of which the coordinates must be converted.
*  X,Y  (input) : Longitude and latitude of the points (degrees).
*      (output) : Map coordinates.
*  PMCONV (out) : Number of valid points after conversion.
*--
*  9-Jan-1991 - created [Remko Scharroo]
* 14-Jan-1991 - conic projections included
* 16-Jan-1991 - azimuthal projections included
* 30-Oct-1991 - more projections included.
* 13-Jan-1992 - Standardize PMPLOT.
*  1-Apr-1993 - Include tilted rectangular projection.
* 16-Nov-1993 - Take care of high latitudes in Molweide.
* 30-Jun-1994 - Include polar projections.
* 14-Aug-1997 - Add ERS projection
*-----------------------------------------------------------------------
      include 'pmplot.inc'
      real a,b,f,er,rc,xy,z,zmin,cosa,cosb,sinb
      integer i,j
*
      PMCONV=n
      if (ptype.eq.2) then
	 do i=1,n
            y(i)=factk*sin(y(i)*rad)
         enddo
      else if (ptype.eq.3 .or. ptype.eq.4) then
	 do i=1,n
            y(i)=factk*log(tan(qpi+y(i)*factl))
         enddo
      else if (ptype.eq.5) then
	 do i=1,n
            y(i)=factk*tan(y(i)*rad2)
         enddo
      else if (ptype.eq.6) then
         do i=1,n
	    if (abs(y(i)).gt.82) then
	       y(i)=1e35
	    else
*	       y(i)=asin(tan(y(i)*rad)/factl)*factk
	       y(i)=sin(tan(y(i)*rad)/factl)*factk
	    endif
	 enddo
      else if (ptype.eq.11) then
	 do i=1,n
	    a=(x(i)-loncen)*rad
	    b=y(i)*rad
	    cosa=cos(a)
	    sinb=sin(b)
	    cosb=cos(b)
	    z=sinb0*sinb+cosb0*cosb*cosa
	    x(i)=cosb*sin(a)
	    y(i)=sinb*cosb0-cosb*sinb0*cosa
	    if (z.lt.0) then
	       call pgnorm(1e20,x(i),y(i))
	       pmconv=pmconv-1
	    endif
         enddo
      else if (ptype.eq.12) then
	 rc=rview/rmean
	 zmin=1/rc
	 do i=1,n
	    a=(x(i)-loncen)*rad
	    b=y(i)*rad
	    cosa=cos(a)
	    sinb=sin(b)
	    cosb=cos(b)
	    z=sinb0*sinb+cosb0*cosb*cosa
	    xy=(rc-1)/(rc-z)
	    x(i)=xy*cosb*sin(a)
	    y(i)=xy*(sinb*cosb0-cosb*sinb0*cosa)
	    if (z.lt.zmin) then
	       call pgnorm(1e20,x(i),y(i))
	       pmconv=pmconv-1
	    endif
         enddo
      else if (ptype.eq.13) then
	 do i=1,n
	    a=(x(i)-loncen)*rad
	    b=y(i)*rad
	    cosa=cos(a)
	    sinb=sin(b)
	    cosb=cos(b)
	    z=sinb0*sinb+cosb0*cosb*cosa
	    xy=sqrt(2/(1+z))
	    x(i)=xy*cosb*sin(a)
	    y(i)=xy*(sinb*cosb0-cosb*sinb0*cosa)
	    if (z.lt.0) then
	       call pgnorm(1e20,x(i),y(i))
	       pmconv=pmconv-1
	    endif
         enddo
      else if (ptype.eq.14) then
	 zmin=cos(rview/rmean)
	 do i=1,n
	    a=(x(i)-loncen)*rad
	    b=y(i)*rad
	    cosa=cos(a)
	    sinb=sin(b)
	    cosb=cos(b)
	    z=sinb0*sinb+cosb0*cosb*cosa
	    er=acos(z)
	    xy=er/sin(er)
	    x(i)=xy*cosb*sin(a)
	    y(i)=xy*(sinb*cosb0-cosb*sinb0*cosa)
	    if (z.lt.zmin) then
	       call pgnorm(1e20,x(i),y(i))
	       pmconv=pmconv-1
	    endif
         enddo
      else if (ptype.eq.21 .or. ptype.eq.22) then
	 do i=1,n
	    a=factk-y(i)
	    f=factl*(x(i)-loncen)
	    x(i)=a*sin(f)
            y(i)=-a*cos(f)
         enddo
      else if (ptype.eq.31) then
	 do i=1,n
            x(i)=(x(i)-loncen)*cos(y(i)*rad)
         enddo
      else if (ptype.eq.32) then
	 do i=1,n
	    if (abs(y(i)).gt.89.5) then
	       x(i)=0.
	    else
	       a=y(i)*rad
	       f=pi*sin(a)
	       do j=1,10
                  a=f-sin(a)
               enddo
	       y(i)=sin(a/2)*90.
               x(i)=(x(i)-loncen)*sqrt(1-(y(i)/90.)**2)
	    endif
	 enddo
      else if (ptype.eq.33) then
	 do i=1,n
	    x(i)=(x(i)-loncen)*cos(factl*y(i))
            y(i)=factk*tan(y(i)*rad2)
	 enddo
      else if (ptype.eq.34) then
	 do i=1,n
	    y(i)=y(i)*factk
	    x(i)=x(i)+y(i)*factl
	 enddo
      else if (ptype.eq.41) then
         do i=1,n
            f=(90.-y(i))*factk
            a=(x(i)-loncen)*rad
            x(i)= f*sin(a)
            y(i)=-f*cos(a)
         enddo
      else if (ptype.eq.42) then
         do i=1,n
            f=(y(i)+90.)*factk
            a=(x(i)-loncen)*rad
            x(i)=f*sin(a)
            y(i)=f*cos(a)
         enddo
      endif
      end
