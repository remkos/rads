**GEODET -- Convert geocentric coordinates to geodetic coordinates
*+
      SUBROUTINE GEODET (LATC, HEIGHT, LAT, R)
      REAL*8 LATC, HEIGHT, LAT, R
*
* This subroutine converts the geocentric latitude (LATC) and
* height (HEIGHT) above the reference ellipsoid GRS80 to the
* geodetic latitude (LAT) and distance to the geocenter (R).
* The method is iterative and stops after 100 iterations or
* when LAT is converged to 1.D-11 radian and R up to 1.D-4 meter.
*
* Arguments:
* LATC   (input) : Geocentric latitude (rad).
* HEIGHT (input) : Height above the reference ellipsoid (m).
* LAT   (output) : Geodetic latitude (rad).
* R     (output) : Distance to geocenter (m).
*-
*  5-Mar-1991: Created, Remko Scharroo
* 18-Nov-1991: Revised
* 19-May-1999: Implemented fixed number of iterations.
* 16-Jan-2002: Version with adjustable ellipsoid coefficients
*-----------------------------------------------------------------------
      REAL*8 LAT0,LATS,R0,RS,FLAT,FFACT,AE,FINV
      INTEGER ITER

      CALL GETEARTH(AE,FINV)
      FLAT=1D0/FINV
      FFACT=FLAT*(FLAT-2D0)
      LATS=LATC
      LAT0=LATC
      R0=AE+HEIGHT
      DO ITER=1,100
         LAT=DATAN(DTAN(LATS)/(FFACT+1D0))
         RS=AE*(1D0-FLAT)/DSQRT(1D0+FFACT*DCOS(LATS)**2)
         R=DSQRT(HEIGHT**2+RS**2+2D0*HEIGHT*RS*DCOS(LAT-LATS))
         IF (DABS(LAT-LAT0).LT.1D-11 .AND. DABS(R-R0).LT.1D-4) RETURN
         LAT0=LAT
         R0=R
         LATS=LATC-DASIN(HEIGHT*DSIN(LAT-LATS)/R)
      ENDDO
      write (0,'(a)') 'geodet: more than 100 iterations'
      END
