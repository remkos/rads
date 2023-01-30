**SPRAND -- Generate random number
*+
      REAL*8 FUNCTION SPRAND()
*
* This routine generates a random number in the range [0,1]. The expected
* distribution is uniform.
*
* Arguments:
*   SPRAND (output): The random number in the range [0,1]
*-
*  2-Jul-1993: Created by Remko Scharroo.
*-----------------------------------------------------------------------
      INTEGER ISEED
      SAVE ISEED

      ISEED=2045*ISEED+1
      ISEED=ISEED-(ISEED/1048576)*1048576
      SPRAND=(ISEED+1)/1048577D0
      END
