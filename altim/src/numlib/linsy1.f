C    ****************************************************************** 
C                                                                       
C    PURPOSE                                                            
C    AUXILLIARY SUBROUTINE TO SUBROUTINE LINSYD
C                                                                       
C    *******************************************************************
      SUBROUTINE LINSY1(S,L,J,N,A,H)
      DOUBLE PRECISION A,T
      INTEGER S,H
      DIMENSION H(N),A(*)
      IF (S.LT.J) THEN
      M=S+H(J)
      ELSE
      M=J+H(S)
      ENDIF
      IF (L.LT.J) THEN
      I=L+H(J)
      ELSE
      I=J+H(L)
      ENDIF
      T=A(M)
      A(M)=A(I)
      A(I)=T
      RETURN
      END
