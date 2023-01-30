**MCW -- Convert SIGMA0 to wind speed using Modified Chelton-Wentz model
*+
      FUNCTION MCW (SIGMA0, U10)
      REAL*8    SIGMA0, U10
      INTEGER*4 MCW
*
* This function computes wind speed (U10, referenced to 10 m above
* sea level) from altimeter backscatter coefficient (SIGMA0) based
* on the Modified Chelton-Wentz algorithm [Witter and Chelton, 1991].
*
* The MCW model function provides wind speed estimates for SIGMA0 values
* ranging from 19.6 dB to 7.0 dB at intervals of 0.2 dB, corresponding
* to 10 m wind speeds between 0 and 20.2 m/s. This subroutine
* interpolates in that table, found in Witter and Chelton [1991].
*
* In case SIGMA0 exceeds the upper limit of the table (19.6 dB), the
* wind speed is assumed to be zero. At the same time MCW is set to +1.
* Wind speeds for SIGMA0 less than the lower limit of the table (7.0 dB)
* are determined by linear extrapolation of the first two table values.
* In that case MCW is set to -1.
* When SIGMA0 is VOID (> 1d20) the value returned for U10 is VOID (1d30)
* and the function value MCW is set to 2.
* When SIGMA0 is NaN, U10 will be NaN and MCW is set to 2.
*
* Arguments:
*  SIGMA0 (input): Backscatter coefficient in dB
*  U10   (output): Wind speed (in m/s) referenced to 10 m above sea level
*  MCW      (out): Function value indicating success of table
*                  interpolation:
*                    0 = success, SIGMA0 value within bounds
*                   -1 = SIGMA0 value below minimum; linear extrapolation
*                   +1 = SIGMA0 value above maximum; wind speed set to 0
*                    2 = SIGMA0 > 1d20
*
* Reference:
*
* Witter D. L. and D. B. Chelton, A Geosat altimeter wind speed
* algorithm and a method for wind speed algorithm development,
* J. Geophys. Res., 96(C5), pp 8853-8860, 1991.
*-
* $Log: mcw.f,v $
* Revision 1.5  2004/07/18 20:10:14  remko
* Added checking for NaN
*
*
*  3-May-2003 - Improved precision in table by using integers
* 21-Feb-2001 - Coded by Remko Scharroo
*-----------------------------------------------------------------------
      real*8	sigmin,sigmax,dsig,x
      parameter (sigmin=7.0d0,sigmax=19.6d0,dsig=0.2d0)
      integer*4 i,tab(64)
      logical	isnan
      data tab /20154,19597,19038,18463,17877,17277,16655,16011,
     |		15348,14669,13976,13273,12557,11830,11092,10345,
     |		 9590, 8827, 8059, 7298, 6577, 5921, 5321, 4763,
     |		 4252, 3792, 3378, 3014, 2708, 2447, 2208, 1992,
     |		 1817, 1676, 1547, 1419, 1292, 1167, 1056,  972,
     |		  915,  873,  833,  794,  755,  716,  677,  637,
     |		  599,  559,  520,  481,  442,  403,  363,  324,
     |		  285,  246,  207,  167,  128,   89,   50,   11/
      save tab

      if (isnan(sigma0)) then
         ! If SIGMA0 is NaN, U10=NaN and MCW=2
	 u10=sigma0
	 mcw=2
	 return
      else if (sigma0.gt.1d20) then
	 ! If SIGMA0 is VOID, U10=VOID and MCW=2
	 u10=1d30
	 mcw=2
	 return
      else if (sigma0.gt.sigmax) then
	 ! If SIGMA0 is larger than maximum table entry, U10=0 and MCW=1
         u10=0d0
	 mcw=1
	 return
      else if (sigma0.lt.sigmin) then
	 ! MCW is set to -1 when SIGMA0 < sigmin
         mcw=-1
	 ! Extrapolate using first two table entries
	 x=(sigma0-sigmin)/dsig
	 i=1
      else
	 ! Proper interpolation within table: MCW=0
         mcw=0
	 ! Interpolate using table entries i and i+1
	 x=(sigma0-sigmin)/dsig+1d0
	 i=int(x)
	 x=x-i
      endif
      ! Linear inter/extrapolation of table, also when SIGMA0 < sigmin
      u10=(tab(i)+(tab(i+1)-tab(i))*x)*1d-3
      return
      end
