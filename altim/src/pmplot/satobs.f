**SATOBS -- compute sub-satellite point from range, elevation and azimuth.
*+
      SUBROUTINE SATOBS (OBSLON, OBSLAT, RANGE, ELEV, AZIM, LON, LAT)
      REAL OBSLON, OBSLAT, RANGE, ELEV, AZIM, LON, LAT
*
* Compute the coordinates of a the sub-satellite point from range,
* elevation and azimuth observation of a satellite from a given station.
*
* Arguments:
*  OBSLON (input) : Longitude of the observer (degrees).
*  OBSLAT (input) : Latitude of the observer (degrees).
*  RANGE  (input) : Distance between observer and satellite (kilometers).
*  ELEV   (input) : Elevation of the range observation (degrees).
*  AZIM   (input) : Azimuth of the range observation (degrees).
*  LON   (output) : Longitude of the sub-satellite point (degrees).
*  LAT   (output) : Latitude of the sub-satellite point (degrees).
*--
* 13-Mar-1991 - created [Remko Scharroo]
*  7-Jul-1994 - All variables defined
*-----------------------------------------------------------------------
      REAL RAD, EPS, SINE, CECF, CESF, SECF, SESF

      RAD=ATAN(1.)/45
      EPS=ATAN2(COS(ELEV*RAD),SIN(ELEV*RAD)+6371./RANGE)
      SINE=SIN(EPS)
      CECF=COS(EPS)*COS(OBSLAT*RAD)
      CESF=COS(EPS)*SIN(OBSLAT*RAD)
      SECF=SINE*COS(OBSLAT*RAD)
      SESF=SINE*SIN(OBSLAT*RAD)
      LON=ATAN2(SINE*SIN(AZIM*RAD),CECF-SESF*COS(AZIM*RAD))/RAD+OBSLON
      LAT=ASIN(CESF+SECF*COS(AZIM*RAD))/RAD
      END
