**PMRND -- find the smallest "round" number greater than X (in degrees)
*+
      REAL FUNCTION PMRND (X, NSUB)
      REAL X
      INTEGER NSUB
*
* Routine to find the smallest "round" number larger than X to be used in
* maps and the appropriate sub intervals. Round numbers are:
* 180 (4), 90 (3), 60 (4), 30 (3), 10 (5), 5 (5), 2 (4), 1 (4),
* 30' (3), 10' (5), 5' (5), 1' (4).
* This routine is used by PMBOX for choosing  tick intervals.
*
* Arguments:
*
*  PMRND (output): the "round" number.
*  X    (input)  : the number to be rounded.
*  NSUB (output) : a suitable number of subdivisions
*--
*  19-Mar-1992 - Created [Remko Scharroo]
*   2-Apr-1993 - 60 degrees included
*   7-Jul-1994 - Define all variables
*-----------------------------------------------------------------------
      INTEGER N, I
      REAL PGRND
      PARAMETER (N=12)
      REAL RND(N)
      INTEGER SUB(N)
      SAVE RND, SUB
      DATA RND
     ./180.,90.,60.,30.,10.,5.,2.,1.,.5,.166667,.083333,.016667/
      DATA SUB
     ./  4 , 3 , 4 , 3 , 5 ,5 ,4 ,4 , 3,     5 ,     5 ,     4 /

      DO I=N,1,-1
         IF (X.LT.RND(I)) EXIT
      ENDDO
      IF (I.EQ.0) THEN
	       PMRND=PGRND (X, NSUB)
      ELSE
	       I=MIN(N,I+1)
	       NSUB=SUB(I)
	       PMRND=RND(I)
      ENDIF
      END
