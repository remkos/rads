**PMCPOLY -- convert X and Y coordinates of a polygon.
*+
      FUNCTION PMCPOLY (N, X, Y)
      INTEGER PMCPOLY, N
      REAL X(N),Y(N)
*
* This function is now called PMCPOL. See there for a full description.
*-
* 19-Jul-1996 - Replaced by PMCPOL
*-----------------------------------------------------------------------
      INTEGER PMCPOL
      PMCPOLY=PMCPOL(N,X,Y)
      END
