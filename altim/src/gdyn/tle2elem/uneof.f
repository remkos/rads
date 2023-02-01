C================================================================               
      SUBROUTINE UNEOF ( IFILE , CSUB , N )                                     
C                                                                               
C----------------------------------------------------------------               
C  UNEOF generates error message for unexpected end-of-file                     
C----------------------------------------------------------------               
C                                                                               
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )                               
      CHARACTER * 6 CSUB                                                        
C                                                                               
      WRITE  ( IFILE , 20 ) CSUB , N                                            
  20  FORMAT ( 1X , 1A6 , ': unexpected end-of-file on unit' , 1I3 ,            
     .         '. stop.' )                                                      
      STOP 99                                                                   
C                                                                               
      END                                                                       
