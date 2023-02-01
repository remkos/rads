**PMCINV -- convert X and Y coordinates of several points (inverse).
*+
      SUBROUTINE PMCINV (N, X, Y)
      INTEGER N
      REAL X(N), Y(N)
*
* Convert map coordinates (x,y) to real-world coordinates (longitude,latitude).
* The real-world coordinates are the position coordinates in degrees,
* the map coordinates are mapped linearly in the viewport and
* depend on the projection type used. For equi-rectangular projection
* real-world and map coordinates are identical, and makes this conversion
* optional.
*
* Arguments:
*  N    (input) : Number of points of which the coordinates must be converted.
*  X,Y  (input) : Map coordinates.
*      (output) : Longitude and latitude of the points (degrees).
*--
* 14-Jan-1991 - created [Remko Scharroo]
* 16-Jan-1991 - support azimuthal projections.
* 30-Oct-1991 - support more projections.
* 19-Mar-1992 - Standardize PMPLOT.
*  2-Apr-1993 - Include tilted rectangular projection.
* 16-Nov-1993 - Take care of high latitudes in Molweide.
* 30-Jun-1994 - Include polar projections (41 & 42)
* 28-Jun-1995 - Include orthographic projection (11)
* 14-Aug-1997 - Add ERS projection
*-----------------------------------------------------------------------
      include 'pmplot.inc'
      real a,f
      integer i
*
      if (ptype.eq.2) then
	 do i=1,n
	    f=min(1.,max(-1.,y(i)/factk))
            y(i)=asin(f)/rad
	 enddo
      else if (ptype.eq.3 .or. ptype.eq.4) then
	 do i=1,n
            y(i)=(atan(exp(y(i)/factk))-qpi)/factl
	 enddo
      else if (ptype.eq.5) then
	 do i=1,n
            y(i)=atan(y(i)/factk)/rad2
	 enddo
      else if (ptype.eq.6) then
         do i=1,n
*	    y(i)=atan(sin(y(i)/factk)*factl)/rad
	    y(i)=atan(asin(y(i)/factk)*factl)/rad
	 enddo
      else if (ptype.eq.11) then
	 do i=1,n
	    f=1-x(i)**2-y(i)**2
	    if (f.gt.0.) then
*	    write (6,'(2f8.4,$)') x(i),y(i)
	       f=sqrt(f)
	       a=-y(i)*sinb0+f*cosb0
	       f= y(i)*cosb0+f*sinb0
	       x(i)=loncen+atan2(x(i),a)/rad
	       y(i)=asin(f)/rad
*	    write (6,'(6f10.4)') loncen,latcen,a,f,x(i),y(i)
	    else
	       x(i)=1e30
	       y(i)=1e30
	    endif
	 enddo
      else if (azmtal) then
	 do i=1,n
	    x(i)=1e30
            y(i)=1e30
	 enddo
      else if (ptype.eq.21 .or. ptype.eq.22) then
	 do i=1,n
	    a=sign(sqrt(x(i)**2+y(i)**2),factk)
	    f=atan2(x(i)/a,-y(i)/a)
	    x(i)=f/factl+loncen
	    y(i)=factk-a
	 enddo
      else if (ptype.eq.31) then
	 do i=1,n
	    f=cos(y(i)*rad)
	    if (f.eq.0) then
	       x(i)=loncen
	    else
	       x(i)=x(i)/f+loncen
	    endif
	 enddo
      else if (ptype.eq.32) then
	 do i=1,n
	    if (abs(y(i)).gt.89.5) then
	       x(i)=1e30
	       y(i)=1e30
	    else
	       x(i)=x(i)/sqrt(1-(y(i)/90.)**2)+loncen
	       a=2*asin(y(i)/90.)
               y(i)=asin((a+sin(a))/pi)/rad
	    endif
	 enddo
      else if (ptype.eq.33) then
	 do i=1,n
	    y(i)=atan(y(i)/factk)/rad2
            x(i)=x(i)/cos(factl*y(i))+loncen
	 enddo
      else if (ptype.eq.34) then
	 do i=1,n
	    x(i)=x(i)-y(i)*factl
	    y(i)=y(i)/factk
	 enddo
      else if (ptype.eq.41) then
	 do i=1,n
	    f=sqrt(x(i)**2+y(i)**2)
	    a=atan2(x(i),-y(i))
	    x(i)=loncen+a/rad
	    y(i)=90.-f/factk
	 enddo
      else if (ptype.eq.42) then
	 do i=1,n
	    f=sqrt(x(i)**2+y(i)**2)
	    a=atan2(x(i),y(i))
	    x(i)=loncen+a/rad
	    y(i)=f/factk-90.
	 enddo
      endif
      end
