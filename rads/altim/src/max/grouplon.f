**GROUPLON -- 'Groups' longitude within boundaries by adding +/- 360 deg
*+
      SUBROUTINE GROUPLON (LON0, LON1, LON)
      REAL*8 LON0, LON1, LON
*
* This routine tries to convert LON to be within boundaries LON0 and
* LON1 by adding or subtracting 2pi.
* If LON is already within LON0 and LON1, LON remains unchanged.
*-
* 13-Oct-1994 - Created. Usage suspect.
*-----------------------------------------------------------------------
      REAL*8 PI,TWOPI
      PARAMETER (PI=3.14159265358979D0,TWOPI=2*PI)

      if (lon.lt.lon0) then
	 lon=lon+twopi
      else if (lon.gt.lon1) then
	 lon=lon-twopi
      endif
      end
