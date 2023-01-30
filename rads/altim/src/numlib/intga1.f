C    ****************************************************************** 
C                                                      20. INTGA1
C                                                                       
C    PURPOSE                                                            
C    TO EVALUATE A 5-GAUSS-LEGENDRE QUADRATURE FORMULA
C                                                                       
C    USAGE                                                              
C    VAR=INTGA1(X,H,F)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    X    I LOWER BOUND OF THE INTERVAL
C    H    I LENGTH OF THE INTERVAL
C    F    I NAME OF DOUBLE PRECISION FUNCTION TO BE INTEGRATED
C                                                                       
C    REMARKS
C    FUNCTION INTGA1 IS USED IN SUBROUTINE INTGAU
C                                                                       
C    *******************************************************************
      DOUBLE PRECISION FUNCTION INTGA1(X,H,F)
      DOUBLE PRECISION X,Y,H,R,T,F
      EXTERNAL F
      Y=X+H/2.D0
      R=0.453089922969332D0*H
      T=0.1184634425280945D0*(F(Y+R)+F(Y-R))
      R=0.2692346550528415D0*H
      T=T+0.239314335249683D0*(F(Y+R)+F(Y-R))
      INTGA1=H*(T+0.2844444444444444D0*F(Y))
      RETURN
      END
