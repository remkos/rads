**NOISE -- Create Gaussian noise
*+
        FUNCTION NOISE ()
        REAL*8 NOISE
*
* This function returns a value from a Gaussian distribution with
* expectation = 0 and variance = 1
* To give the "seed value" to the randomizer, use the standard Fortran
* routine SRAND
*
* Example:
*       REAL*8 R, NOISE
*       INTEGER SEED
*
*       CALL SRAND(SEED)
*       R = NOISE()
*
* Argument:
* NOISE (output): Random value from Gaussian distribution
*-
*  1-May-1995 - Created (Remko Scharroo and Pieter Visser)
* 10-May-2000 - Updated for real*8 version
*-----------------------------------------------------------------------
        integer*4 k
	real*8 rand

        noise=-6d0
        do k=1,12
           noise=noise+rand()
        enddo
        end
