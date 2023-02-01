**POLCAR -- Convert Polar to Cartesian
*+
      SUBROUTINE POLCAR (LAT, LON, R, VECTOR)
      REAL*8 LAT, LON, R, VECTOR(3)
*
* Convert polar coordinates (LAT,LON,R) to Cartesian VECTOR() = (X,Y,Z)
*
* Arguments:
* LAT     (input) : Latitude in radians
*                   = angle from XY-plane towards +Z-axis
* LON     (input) : Longitude in radians
*                   = angle in XY-plane measured from +X-axis towards +Y-axis
* R       (input) : Radius
* VECTOR (output) : Vector of X, Y, and Z coordinates.
*-
* 28-Apr-1991 - Created: Remko Scharroo, DUT/SOM (c)
*-----------------------------------------------------------------------
      VECTOR(1)=R*COS(LAT)*COS(LON)
      VECTOR(2)=R*COS(LAT)*SIN(LON)
      VECTOR(3)=R*SIN(LAT)
      END
