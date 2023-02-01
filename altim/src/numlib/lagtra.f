C    ****************************************************************** 
C                                                      31. LAGTRA
C                                                                       
C    PURPOSE                                                            
C    TO TRANSFORM COEFFICIENTS CI OF THE LAGUERRE SERIES
C    C0*L0(X)+.......+CN*LN(X) INTO COEFFICIENTS CFI OF
C    THE SERIES CF0*U**0+.........+CFN*U**N, WHERE X=A*U+B
C                                                                       
C    USAGE                                                              
C    CALL LAGTRA(N1,A,B,C,CF,P)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N1   I MAXIMUM INDEX +1 IN THE SERIES
C    A    I COEFFICIENT IN X=A*U+B
C    B    I COEFFICIENT IN X=A*U+B
C    C    I ARRAY C(N1), WHERE C(I+1)=CI, WHERE N1=N+1
C    CF   O ARRAY CF(N1), WHERE CF(I+1)=CFI
C           CONTAINS THE RESULTANT COEFFICIENTS
C    P      ARRAY P(M), WHERE M IS AT LEAST 2*N+3, WORKING SPACE
C
C    REMARKS
C    1. FOR N<0, THE SUBROUTINE HAS NO EFFECT
C
C    *******************************************************************
      SUBROUTINE LAGTRA(N1,A,B,C,CF,P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(N1),CF(N1),P(*)
      DATA Z,U,T/0.D0,1.D0,2.D0/
      IF (N1.LT.1) GOTO 40
      N2=N1+1
      CF(1)=C(1)
      IF (N1.EQ.1) GOTO 40
      P(N2)=Z
      P(N2+1)=U-B
      P(N2+2)=-A
      P(1)=U
      P(2)=Z
      CF(1)=C(1)+C(2)*P(N2+1)
      CF(2)=C(2)*P(N2+2)
      IF (N1.GE.3) THEN
      DO 20 I=3,N1
      CF(I)=Z
      P(I+N2)=Z
      P(I)=Z
      AA=-A/(I-1)
      BB=T-(B+U)/(I-1)
      CC=U-U/(I-1)
      DO 10 J=I,1,-1
      PN=AA*P(J+N1)+BB*P(J+N2)-CC*P(J)
      CF(J)=CF(J)+C(I)*PN
      P(J)=P(J+N2)
  10  P(J+N2)=PN
  20  CONTINUE
      ENDIF
  40  CONTINUE
      RETURN
      END
