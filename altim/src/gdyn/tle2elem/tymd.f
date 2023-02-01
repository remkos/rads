C***********************************************************************        
      SUBROUTINE TYMD(T,IYMD,DHMS)                                              
C                                                                               
C  subroutine TYMD converts T (in mjd.dd) to IYMD (in yymmdd)                   
C  and DHMS (in hhmmss.ss)                                                      
C                                                                               
      IMPLICIT REAL * 8 (A-H,O-Z)                                               
C                                                                               
      MJD=T                                                                     
      CALL MJDATE(1,MJD,IYMD,IFLUT1,IFLUT2,IFLUT3)                              
      TFRAC=T-DFLOAT(MJD)                                                       
      IH=TFRAC*24.0D0                                                           
      IM=TFRAC*1440.0D0-DFLOAT(IH)*60.0D0                                       
      DSEC=TFRAC*86400.0D0-DFLOAT(IM)*60.0D0-DFLOAT(IH)*3600.0D0                
      DHMS=DFLOAT(IH)*10000.0D0+DFLOAT(IM)*100.0D0+DSEC                         
C                                                                               
      RETURN                                                                    
      END                                                                       
