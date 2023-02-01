C    ****************************************************************** 
C                                                      55. NULTWE
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE ROOTS OF A QUADRATIC EQUATION
C                                                                       
C    USAGE                                                              
C    CALL NULTWE(B,C,E)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    B,C  I COEFFICIENTS OF THE QUADRATIC EQUATION
C            X*X-2*B*X+C = 0
C         O ROOTS, SEE IE
C    IE   O EXITCODE
C           IE=0: REAL ROOTS B AND C, WHERE B>C
C           IE=1: COMPLEX ROOTS B+IC AND B-IC
C
C    REFERENCES
C    RC-TWA-81002 ROOTS OF QUADRATIC AND CUBIC EQUATIONS
C
C    *******************************************************************
      SUBROUTINE NULTWE(B,C,IE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA Z,U/0.D0,1.D0/
      IE=0
      IF (C.EQ.Z) THEN
      T=2.D0*B
      IF (B.GT.Z) THEN
      B=T
      ELSE
      C=T
      B=Z
      ENDIF
      ELSE
      R=B*B-C
      T=DSQRT(DABS(R))
      IF (R.LT.Z) THEN
      IE=1
      C=T
      ELSE
      IF (B.GT.Z) THEN
      B=B+T
      C=C/B
      ELSE
      R=C
      C=B-T
      B=R/C
      ENDIF
      ENDIF
      ENDIF
      RETURN
      END
