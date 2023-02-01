C    ******************************************************************
C                                                      14. GAMMAF
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE GAMMA FUNCTION GAMMA(X)
C                                                                       
C    USAGE                                                              
C    CALL GAMMAF(X,IC,G)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    X    I ARGUMENT VALUE
C    IC   O RETURN CODE:
C           X < -50            IC=-2
C           -50 <= X < 1       IC=-1
C           1 <= X <= 50       IC=1
C           50 < X             IC=2
C           FOR X < 1 AND ABS(X-I) <= 1.D-63 THEN IC=0,
C           WHERE I IS AN INTEGER NUMBER
C           FOR RC=-2          G= 1.D-63
C           FOR RC= 0          G= 1.D63
C           FOR RC= 2          G= 1.D63
C    G    O THE VALUE OF GAMMA(X)
C
C    REMARKS
C    FOR X > 50 GAMMAF(X) WOULD BE LARGER THAN THE MAXIMUM NUMBER.
C    IN SUCH CASE LOG GAMMA(X) SHOULD BE COMPUTED.
C
C    METHOD                                                             
C    REFER TO : C.W. CLENSHAW, MATHEMATICAL TABLES, VOL 5 TABLE 8
C
C    *******************************************************************
      SUBROUTINE GAMMAF(X,IC,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(22)
      DATA A/0.9417855977954947D0,0.44153813248410D-2,0.568504368159936D
     1-1,-0.42198353964186D-2,0.13268081812125D-2,-0.1893024529799D-3,
     20.360692532744D-4,-0.60567619045D-5,0.10558295463D-5,
     3-0.1811967366D-6,0.311772496D-7,-0.53542196D-8,0.9193276D-9,
     4-0.1577941D-9,0.270798D-10,-0.46468D-11,0.7973D-12,-0.1386D-12,
     50.253D-13,-0.40D-14,0.7D-15,-0.1D-15/
      DATA Z,U,GM/0.0D0,1.0D0,1.D-63/
      IF (DABS(X).LE.50.D0) THEN
      IF (X.GE.U) THEN
      IC=1
      ELSE
      IC=-1
      ENDIF
      ELSE
      IF (X.GT.Z) THEN
      IC=2
      ELSE
      IC=-2
      ENDIF
      ENDIF
      GOTO (10,20,30,40,50),IC+3
  10  GA=GM
      GOTO 200
  20  N=X
      IF (X.LT.Z) N=N-1
      T=X-N
      IF ((T.LT.GM).OR.((U-T).LE.GM)) THEN
      GOTO 30
      ELSE
      S=T+T-U
      CALL ERROR1(22,S,A,GA)
      DO 24 I=0,N,-1
  24  GA=GA/(T+I)
      ENDIF
      GOTO 200
  40  N=X
      T=X-N
      S=T+T-U
      CALL ERROR1(22,S,A,GA)
      DO 44 I=N-1,1,-1
  44  GA=GA*(T+I)
      GOTO 200
  30  IC=0
  50  GA=1.D63
 200  G=GA
      RETURN
      END
