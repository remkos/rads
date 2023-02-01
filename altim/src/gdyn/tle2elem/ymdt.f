C***********************************************************************        
      SUBROUTINE YMDT(IYMD,DHMS,T)                                              
C                                                                               
C  subroutine YMDT converts IYMD (in yymmdd) and DHMS (in hhmmss.ss)            
C  to T (in mjd.dd).                                                            
C                                                                               
      IMPLICIT REAL * 8 (A-H,O-Z)                                               
C                                                                               
      CALL MJDATE(2,MJD,IYMD,IFLUT1,IFLUT2,IFLUT3)                              
      IH=DHMS/1.0D4                                                             
      IM=DHMS/1.0D2                                                             
      IM=IM-100*IH                                                              
      DSEC=DHMS-DFLOAT(IH)*1.0D4-DFLOAT(IM)*1.0D2                               
      T=DFLOAT(MJD)+DFLOAT(IH)/24.0D0+DFLOAT(IM)/1440.0D0                       
     .  +DSEC/86400.0D0                                                         
C                                                                               
      RETURN                                                                    
      END                                                                       
