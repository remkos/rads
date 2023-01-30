C    ****************************************************************** 
C                                                      50. MATSY1
C                                                                       
C    PURPOSE                                                            
C    TO FILL AN INTEGER VECTOR H, USED IN SUBROUTINES
C    CONCERNING SYMMETRIC MATRICES
C                                                                       
C    USAGE                                                              
C    CALL MATSY1(N,H)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I NUMBER OF ELEMENTS IN H
C    H    O RESULTANT INTEGER VECTOR
C                                                                       
C    REMARKS
C    H(I)=I*(I-1)/2
c
C    *******************************************************************
      SUBROUTINE MATSY1(N,H)
      INTEGER H(N)
      H(1)=0
      DO 10 I=2,N
   10    H(I)=H(I-1)+I-1
      END
