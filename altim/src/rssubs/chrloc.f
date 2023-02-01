**CHRLOC -- Write latitude and longitude in character format
*+
      SUBROUTINE CHRLOC (LAT, LON, POSITION)
      INTEGER*4     LAT, LON
      CHARACTER*(*) POSITION
*
* Arguments:
* LAT       (input) : Latitude in microdegrees
* LON       (input) : Longitude in microdegrees
* POSITION (output) : Position in format '72.0N 179.9E'
*-
* 13-May-1991: created - Remko Scharroo
*-----------------------------------------------------------------------
      CHARACTER A1*1,A2*1
      REAL*8 XLAT,XLON
*
      XLAT=LAT/1D6
      XLON=LON/1D6
      IF (XLAT.LT.0) THEN
         XLAT=-XLAT
         A1='S'
      ELSE
         A1='N'
      ENDIF
      IF (XLON.GT.180) THEN
         XLON=360-XLON
         A2='W'
      ELSE
         A2='E'
      ENDIF
      WRITE (POSITION,10) XLAT,A1,XLON,A2
*
   10 FORMAT (F4.1,A1,F6.1,A1)
      END
