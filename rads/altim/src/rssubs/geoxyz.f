**GEOXYZ -- Convert geodetic latitude, longitude and height to ECF coordinates.
*+
      SUBROUTINE GEOXYZ (LAT, LON, HEIGHT, XYZ, R)
      REAL*8 LAT, LON, HEIGHT, XYZ(3), R
*
* This subroutine converts geodetic latitude, longitude and height above the
* GRS80 reference ellipsoid to Earth Centered Fixed coordinates X, Y and Z.
*
* Arguments:
* LAT     (input) : Geodetic latitude (rad).
* LON     (input) : Geodetic longitude (rad).
* HEIGHT  (input) : Height above the reference ellipsoid (m).
* XYZ(3) (output) : Earth Centered Fixed coordinates X, Y, and Z (m).
* R      (output) : Distance to geocenter (m).
*-
* 19-Nov-1991: Created - Remko Scharroo
*-----------------------------------------------------------------------
      REAL*8 LATC
      CALL GEOCEN (LAT, HEIGHT, LATC, R)
      CALL POLCAR (LATC, LON, R, XYZ)
      END
