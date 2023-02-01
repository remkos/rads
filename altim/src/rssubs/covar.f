**COVAR -- Determine undulation covariance of two points.
*+
      FUNCTION COVAR (PSI)
      REAL*8 PSI, COVAR
*
* Compute the covariance between the geoid height values in two points
* as a function of their spherical distance in spherical approximation.
* WARNING: this function can only be used after
*          `CALL COVINI (FILENM, ITYPE, MINDEG, MAXDEG)'
*
* Arguments:
*  PSI    (input): Spherical distance between the two points (rad).
*  COVAR (output): Geoid undulation covariance of the two points (m**2).
*-
*  1-Feb-1991: Created. Pieter Visser and Remko Scharroo.
*-----------------------------------------------------------------------
      INCLUDE 'covar.inc'
      INTEGER*4 I
      REAL*8 S(0:MX),T

      S(NMAX+2)=0
      S(NMAX+1)=0
      T=DCOS(PSI)

* Recurrent relation for the series:
*     nmax
* N =  SUM  l  * P  (t)
*      n=0   n    n

      DO I=NMAX,0,-1
         S(I)=LN(I)-AN(I+1)*S(I+1)*T-BN(I+2)*S(I+2)
      ENDDO

      COVAR=S(0)
      END


**COVAR2 -- Determine gravity covariance of two points.
*+
      FUNCTION COVAR2 (Q1, Q2, PSI)
      INTEGER*4 Q1, Q2
      REAL*8 PSI, COVAR2
*
* Compute the covariance between the gravity values in two points,
* or some of its derived quantities, as a function of the spherical
* distance between the points in spherical approximation.
* WARNING: this function can only be used after
*          `CALL COVINI (FILENM, ITYPE, MINDEG, MAXDEG)'
*
* Arguments:
*  Q1, Q2  (input): Integers denoting the quantities that are considered:
*                   1 = Gravity potential (units: m**2/s**2)
*                   2 = Gravity anomaly   (units: mgal)
*                   3 = Geoid undulation  (units: m)
*  PSI     (input): Spherical distance between the two points (rad).
*  COVAR2 (output): Covariance of the two points.
*-
*  1-Feb-1991: Created. Remko Scharroo.
*-----------------------------------------------------------------------
      INCLUDE 'covar.inc'
      INTEGER*4 I
      REAL*8 S(0:MX),T,FACOVAR
*
      S(NMAX+2)=0
      S(NMAX+1)=0
      T=DCOS(PSI)

* Recurrent relation for the series:
*     nmax
* N =  SUM  l  * fac * P  (t)
*      n=0   n          n

      DO I=NMAX,0,-1
         S(I)=LN(I)*FACOVAR(Q1,I)*FACOVAR(Q2,I)
     |            -AN(I+1)*S(I+1)*T-BN(I+2)*S(I+2)
      ENDDO
      COVAR2=S(0)
      END

      function facovar(q,n)
      real*8 facovar
      integer*4 q,n
      include 'covar.inc'
      real*8 gamma
      parameter (gamma=gm/ae**2)
      if (q.eq.1) then
	 facovar=gamma
      else if (q.eq.2) then
	 facovar=gamma/ae*(n-1)*1d5
      else
	 facovar=1
      endif
      end


**COVINI -- Initiate covariance function.
*+
      SUBROUTINE COVINI (FILENM, ITYPE, MINDEG, MAXDEG)
      CHARACTER*(*) FILENM
      INTEGER*4 ITYPE, MINDEG, MAXDEG
*
* This subroutine initializes the geoid undulation degree variances
* (ln) taking the differences between two geoid models complete to
* degree `ndeg' and a model for higher degrees.
* This subroutine must be executed before using the function `covar'.
*
* Arguments:
*  FILENM (input): Name of the file containing the degree variances for low
*                  degrees.
*  ITYPE  (input): Type of model used:
*                  0 = all higher orders set to zero.
*                  1 = Rapp (1979).
*                  2 = Kaula.
*                  3 = Kaula (to fit OSU89b order 36-360).
*                  4 = Kaula (to fit OSU89b-OSU86f order 36-360).
*  MINDEG (input): Minimum degree to be included in covariance function.
*  MAXDEG (input): Maximum degree to be included in covariance function.
*-
*  5-Feb-1991. Pieter Visser and Remko Scharroo
* 18-Sep-1991. Message added if file not found.
*-----------------------------------------------------------------------
      INTEGER*4 B1,B2,NDEG,I
      REAL*8 DEGVAR,S1,S2,A,A1,A2,B,FACTOR
      INCLUDE 'covar.inc'

* Check if array dimensions are not exceeded

      IF (MAXDEG.GT.MX-2) THEN
         WRITE (0,"('COVINI: model truncated at degree ',I6)") MX-2
         MAXDEG=MX-2
      ENDIF
      NMAX=MAXDEG

* Load degree variances of difference between two models and convert
* them to undulation degree variances

      NDEG=0
      OPEN (10,STATUS='OLD',FILE=FILENM,ERR=15)
      REWIND(10)
   10 READ (10,*,END=20) NDEG,DEGVAR
      LN(NDEG)=DEGVAR*AE**2
      GOTO 10

* Add model for higher orders

   15 WRITE (0,"('COVINI: model unknown; only extension used')")
   20 IF (ITYPE.EQ.0) THEN
         DO I=NDEG+1,MAXDEG
            LN(I)=0
         ENDDO
      ELSE IF (ITYPE.EQ.1) THEN
         A1=3.4050D-10
         A2=140.03D-10
         S1=.998006
         S2=.914232
         B1=1
         B2=2
         FACTOR=AE**6/GM**2
         DO I=NDEG+1,MAXDEG
            LN(I)=FACTOR/(I-1)
     |            *(A1*S1**(I+2)/(I+B1)+A2*S2**(I+2)/(I-2)/(I+B2))
         ENDDO
      ELSE IF (ITYPE.EQ.2) THEN
         A=10
         B=2
         A=(AE*A*1D-6)**2
         B=2*B
         DO I=NDEG+1,MAXDEG
            LN(I)=(2*I+1)*A/I**B
         ENDDO
      ELSE IF (ITYPE.EQ.3) THEN
         A=23.53D0
         B=2.187D0
         A=(AE*A*1D-6)**2
         B=2*B
         DO I=NDEG+1,MAXDEG
            LN(I)=(2*I+1)*A/I**B
         ENDDO
      ELSE IF (ITYPE.EQ.4) THEN
         A=.8616D0
         B=1.656D0
         A=(AE*A*1D-6)**2
         B=2*B
         DO I=NDEG+1,MAXDEG
            LN(I)=(2*I+1)*A/I**B
         ENDDO
      ENDIF

* Clear orders that are not to be included

      DO I=0,MINDEG
         LN(I)=0
      ENDDO

      DO I=1,MAXDEG+2
         AN(I)=DFLOAT(1-2*I)/I
         BN(I)=DFLOAT(I-1)/I
      ENDDO
      END
