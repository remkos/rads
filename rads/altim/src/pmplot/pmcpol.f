**PMCPOL -- convert X and Y coordinates of a polygon.
*+
      FUNCTION PMCPOL (N, X, Y)
      INTEGER PMCPOL, N
      REAL X(N),Y(N)
*
* Convert real-world coordinates (longitude,latitude) to map coordinates (x,y).
* The real-world coordinates are the position coordinates in degrees,
* the map coordinates are mapped linearly in the viewport and
* depend on the projection type used.
* For map projections other than azimuthal, PMCPOL is equivalent to
* PMCONV. When an azimuthal projection is chosen, PMCPOL maps
* out-of-area coordinates to the circular edge.
*
* Arguments:
*  N       (input) : Number of points in the polygon.
*  X,Y     (input) : Longitude and latitude of the points (degrees).
*         (output) : Map coordinates.
*  PMCPOL (output) : Number of points inside window.
*--
* 18-Mar-1992 - Created [Remko Scharroo].
*  7-Jul-1994 - Define all variables
* 17-Jul-1996 - Complete revision. Renamed from PMCPOLY.
*-----------------------------------------------------------------------
      INCLUDE 'pgplot.inc'
      INTEGER I, PMCONV
*
      PMCPOL=PMCONV(N,X,Y)
      IF (PMCPOL.NE.N) THEN
      DO I=1,N
	 IF (X(I)**2+Y(I)**2.GT.1E20) THEN
	    X(I)=X(I)/1E20
	    Y(I)=Y(I)/1E20
	 ENDIF
      ENDDO
      ENDIF
      END
