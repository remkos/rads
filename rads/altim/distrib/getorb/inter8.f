**INTER8 -- Interpolate vector using 8th order Legendre interpolation
*+
      SUBROUTINE INTER8 (N, T, V, VT)
      INTEGER*4 N
      REAL*8  T, V(N,9), VT(N)
*
* This subroutine interpolates a N dimensional vector V at time T given
* 9 equi-distant N dimensional vectors V. These vectors refer to:
*  T = 0  ->  V(1,1) ... V(N,1)
*  T = 1  ->  V(1,2) ... V(N,2)
*   ...
*  T = 8  ->  V(1,9) ... V(N,9)
*
* Arguments:
* N   (input) : Number of dimensions.
* T   (input) : Time at which VT must be interpolated. E.g. T=3.5 refers
*               to interpolation between V(.,4) and V(.,5).
* V   (input) : 9 consecutive vectors of dimension N.
* VT (output) : Interpolated vector at time T.
*-
* 28-Apr-1991 - Created: Remko Scharroo, DUT/SOM (c)
* 22-Feb-1996 - SAVE before DATA for HP compatibility
*-----------------------------------------------------------------------
* For computational reasons T is changed to X relative to one step
* prior to the interpolation interval
*
      REAL*8  X,TELLER,COEFF
      INTEGER*4 NOEMER(9),KX,I
      SAVE NOEMER
      DATA NOEMER /40320,-5040,1440,-720,576,-720,1440,-5040,40320/

      X=T+1
      TELLER=(X-1)*(X-2)*(X-3)*(X-4)*(X-5)*(X-6)*(X-7)*(X-8)*(X-9)

      IF (TELLER.EQ.0) THEN
         KX=NINT(X)
         DO I=1,N
            VT(I)=V(I,KX)
         ENDDO
         RETURN
      ENDIF

      DO I=1,N
         VT(I)=0
      ENDDO
      DO KX=1,9
         COEFF=TELLER/NOEMER(KX)/(X-KX)
         DO I=1,N
            VT(I)=VT(I)+COEFF*V(I,KX)
         ENDDO
      ENDDO
      END
