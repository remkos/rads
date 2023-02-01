**GEOCEN -- Convert geodetic coordinates to geocentric coordinates
*+
      SUBROUTINE GEOCEN (LAT, HEIGHT, LATC, R)
      REAL*8 LAT, HEIGHT, LATC, R
*
* This subroutine converts the geodetic coordinates (latitude and height)
* on the GRS80 reference ellipsoid to geocentric coordinates (geocentric
* latitude and geocentric distance).
*
* Arguments:
* LAT    (input) : Geodetic latitude (rad).
* HEIGHT (input) : Height above the reference ellipsoid (m).
* LATC  (output) : Geocentric latitude (rad).
* R     (output) : Distance to geocenter (m).
*-
*  5-Mar-1991: Created, Remko Scharroo
* 18-Nov-1991: Revised
* 16-Jan-2002: Version with adjustable ellipsoid coefficients
*-----------------------------------------------------------------------
      REAL*8 RS, LATS, FLAT, FFACT, AE, FINV

      CALL GETEARTH(AE,FINV)
      FLAT=1D0/FINV
      FFACT=FLAT*(FLAT-2D0)
      LATS=DATAN((FFACT+1D0)*DTAN(LAT))
      RS=AE*(1D0-FLAT)/DSQRT(1D0+FFACT*DCOS(LATS)**2)
      R=DSQRT(HEIGHT**2+RS**2+2D0*HEIGHT*RS*DCOS(LAT-LATS))
      LATC=LATS+DASIN(HEIGHT*DSIN(LAT-LATS)/R)
      END
