      program df2gbf

**DF2GBF -- Conversion of DORIS doppler data into GBF

* February 1990  : created - Dirk Kuijper 
* March 1995     : stored - Pieter Visser
* September 2001 : adapted - Eelco Doornbos

C PURPOSE          : This program converts the DORIS data format to the         
C                    GEODYN binary format.                                      
C                                                                               
C SR CALLED        :                                                            
C                                                                               
C REMARKS          :                                                            
C                                                                               
C AUTHORS          : D.C.KUIJPER                                                
C                                                                               
C CREATED          : 07.02.1990                                                 
C                                                                               
C ADJUSTED         : 13.10.1992 BY D.C. KUIJPER FOR THE CONVEX                  C                                                                               
C COPYRIGHT        : DELFT UNIVERSITY OF TECHNOLOGY                             
C                    FACULTY OF AEROSPACE ENGINEERING                           
C                    SECTION ORBITAL MECHANICS                                  
C                    KLUYVERWEG 1                                               
C                    2629 HS  DELFT                                             
C                    THE NETHERLANDS                                            
C                                                                               
C PARAMETER SPECIFICATION:                                                      
C                                                                               
C MTYP              : Measurement Type                         (I2)             
C TIMIND            : Time System Indicator                    (I2)             
C NTSI              : N-value of Time System Indicator         (I4)             
C MTSI              : M-value of Time System Indicator         (I4)             
C                                                                               
C GEODYN BINARY FORMAT                                                          
C --------------------                                                          
C ISAT              : Satellite ID                             (I4)             
C ISTAT             : Station Number                           (I4)             
C ICOUNT            : Counting Interval (10E-8 sec) (d-type D) (I4)             
C IPRE              : Preprocessing Indicator/Report           (I4)             
C MJD               : Modified Julian Date of observation      (I4)             
C MJCENT            : Modified Julian Century                  (I4)             
C OSD               : Standard Deviation of Observations (m/s) (R4)             
C TROPRC            : Tropospheric Refraction Correction (m/s) (R4)             
C IONRC             : Ionospheric Refraction Correction  (m/s) (R4)             
C AOC               : Antenna Offset Correction            (m) (R4)             
C BOC               : Beacon Offset Correction             (m) (R4)             
C FRACT             : Fraction of day past midnight (GMT)      (R8)             
C DRRATE            : Observation values (meters/sec)          (R8)             
C IDT               : Doppler type (destruct type -> IDT = 1)  (I4)             
C IC                : Speed of Light (2.99792458*10**8, IC = 3)(I4)             
C                                                                               
C DORIS FORMAT                                                                  
C ------------                                                                  
C SATID             : Satellite ID                             (A6)             
C STATID            : Station Number                           (A6)             
C IYEAR             : Year of observation                      (I4)             
C IDAYS             : Day of year                              (I4)             
C ISEC              : Seconds after midnight                   (I4)             
C IFRACT            : Fraction of day past midnight (GMT)      (I4)             
C PREIND            : Preprocessing Indicator/Report           (A3)             
C DCOUNT            : Counting Interval (0.1 microsec)         (R8)             
C IRRATE            : Observation values (micrometers/sec)     (I4)             
C IPRES             : Surface Pressure (millibars)             (I4)             
C ITEMP             : Surface Temperature (degrees kelvin)     (I4)             
C IRHUM             : Relative Humidity (%)                    (I4)             
C IOSD              : Observation Standard Deviation (m/s)     (I4)             
C IRC               : Ionospheric Refraction Correction  (m/s) (I4)             
C ITRC              : Tropospheric Refraction Correction (m/s) (I4)             
C IMETDS            : Meteorological Data Source               (I4)             
C IPCE              : Phase Center Effect (micrometers/sec)    (I4)             
C ICA               : Ionosphere "correction Applied           (I4)             
C ITCA              : Troposphere Correction Applied           (I4)             
C IPS               : Point status                             (I4)             
C ICHANNEL          : Channel (for dual-channel receivers)     (I4)
C                                                                               
C  Input:                                                                       
C  - FT10 : the DORIS formatted input file                                      
C  - FT20 : the station coordinates input file                                  
C                                                                               
C  Output:                                                                      
C  - FT15 : the resulting GEODYN BINARY formatted file.                         
C  - FT25 : the skip output file containing all the records skipped             
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                               
C  The output (GBF data) is written to: 'gbffile'
C  Skipped records are written to: 'skipfile'                                   
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      PARAMETER (M8 = 2**8)                                                     
      PARAMETER (M16 = 2**16)                                                   
      PARAMETER (ISTMAX = 1000)                                                
      INTEGER*2 MTYP, TIMIND                                                    
      INTEGER IUNIT, OUNIT, DCOUNT, STATI
      INTEGER BDATE, EDATE, DBEG, DEND 
      REAL OSD, TROPRC, IONRC, AOC, BOC
      CHARACTER*3 PREIND                                                        
      CHARACTER*5 STATID, STATOLD, STNAME(ISTMAX)
      CHARACTER*6 SATID                                                         
      CHARACTER*12 STATFILE
      CHARACTER*96 LINE                                                         
      CHARACTER*96 systemdata
      INTEGER*4 len
      REAL*8 t0, t1, t2, tsec, sec85
      LOGICAL FIRST, STATUS, ADJUST, DATEARG, CHECKDATES
      DIMENSION IDST(ISTMAX)                                                    
      DATA IT,IH,ID,ITT,IDT,IC,MJCENT / 10,100,1000,5,1,3,15019 /               
      DATA ITOI / 1 /
      DATA DMICRO, SPD, DT / 1.D-06, 86400.D0, 10.D0 /                          
C                                                                               
C Initialization                                                                
C                                                                               
      checkdates = .false.
      NREC = 0                                                                  
      NMES = 0                                                                  
      NSKIP = 0                                                                 
      NSKTM = 0
      NSKED = 0
      NSKST = 0
      I = 1                                                                     
      IUNIT = 5                                                                
      OUNIT = 15                                                                
      STATI = 20                                                                
      MUNIT = 6                                                                
      FIRST = .TRUE.
      STATUS = .TRUE.
      ISUMI = 0
      ISUMT = 0
      ISUMP = 0
      AOC = 0.D0
      BOC = 0.D0

* Check arguments

      do i=1,iargc()
        call getarg(i,line)
        if (datearg(line,t0,t1,t2)) then
          checkdates = .true.          
        end if
        if (line(:2).eq.'-h') then
          write(0,600)
600       format (
     |'df2gbf - Convert DORIS data to GBF'/
     |'syntax: df2gbf t=t0,t1',
     |'[outputfile] < inputfile'//
     |'where: '/
     |'  t0,t1      : time boundaries'/
     |'  inputfile  : input file in DORIS data format'/
     |'  outputfile : output file in GBF format'//
     |'Further print output comes to standard output')
          goto 9999
        else
          OPEN(OUNIT,FILE=line,FORM='unformatted')
        endif
      enddo
C
      DBEG=0
      DEND=366
C
C Open input and output files
C
      CALL GETENV("ALTIM",systemdata)
      lngth = LNBLNK(systemdata)
      IF (systemdata .eq. "") THEN
        systemdata = '/uwa/altim/data/tables/system.data'
      ELSE 
        systemdata = systemdata(1:lngth) // '/data/tables/system.data'
      END IF
      OPEN (STATI,FILE=systemdata,FORM='formatted')
C                                                                               
C Read station info file                                                        
C                                                                               
 10   READ(STATI,'(A96)',END=20)LINE
      IF(LINE(1:1).eq.'4') THEN
        READ(LINE,15) IDST(I),STNAME(I) 
 15     FORMAT(I4,6X,A5)                                                       
        I = I + 1                                                               
      END IF
C
      GOTO 10                                                                   
 20   CONTINUE                                                                  
      IMAX = I
C                                                                               
C Write header data status file                                                 
C                                                                               
      IF (STATUS) WRITE(MUNIT,25)
 25   FORMAT ('SATELLITE',' STATION',' YEAR',' MM/DD','  START TIME',
     . ' END TIME','  NUM OBS',' IRMS','    TRMS','    PRMS',/)
C                                                                               
C Read a DORIS FORMATTED file                                                   
C                                                                               
 40   READ (IUNIT,'(A96)',END=100)LINE                                          
      READ (LINE,50)SATID,MTYP,TIMIND,STATID,IYEAR,IDAYS,ISEC,                  
     .              IFRACT,PREIND,DCOUNT,RRATE,IPRES,ITEMP,IRHUM,               
     .              IOSD,IRC,ITRC,IBEACON,IMETDS,ICHANNEL,IPCE  
 50   FORMAT(1X,A6,2I2,A5,I2,I3,I5,I6,A3,I10,                                  
     .          F11.0,I4,2I3,I6,I8,I7,I1,I1,I1,I6)                             
      IF (IDAYS.LT.DBEG) GOTO 40
      IF (IDAYS.GT.DEND) GOTO 40
      NREC = NREC + 1                                                           
      IF (MTYP.NE.39) GOTO 40                                                   
      READ(PREIND,'(3I1)')ICA,ITCA,IPS                                          
      IF (IPS .GT. 0) THEN                                                      
         NSKIP = NSKIP + 1                                                      
         NSKED = NSKED + 1
         GOTO 40                                                                
      END IF                                                                    
C                                                                               
      IF (FIRST) THEN
	 STATOLD = STATID
	 ISTOLD = ISTAT
	 NP=0
      END IF
C                                                                               
C Preprocessing the data                                                        
C                                                                               
      IF (SATID.EQ.'SPOT2 ') ISAT = 9000501
      IF (SATID.EQ.'TOPEX1') ISAT = 9205201
      IF (SATID.EQ.'JASON1') ISAT = 0105501
      IF (SATID.EQ.'ENVIS1') ISAT = 0200901
      MTYP = 92                                                                 
      NTSI = TIMIND/IT                                                          
      MTSI = TIMIND - NTSI*IT                                                   
      IF (NTSI .EQ. 3) THEN                                                     
         NTSI = 1                                                               
      END IF                                                                    
      TIMIND = IT*NTSI + MTSI                                                   
      DO 60 N = 1, IMAX
         IF (STATID .EQ. STNAME(N)) THEN                                        
            ISTAT = IDST(N)                                                     
            GOTO 70                                                             
         END IF                                                                 
 60   CONTINUE                                                                  
      NSKIP = NSKIP + 1                                                         
      NSKST = NSKST + 1
      WRITE(MUNIT,61) STATID
 61   FORMAT('Station ',A5,' unknown. Please check system.data.')
      GOTO 40                                                                   
 70   CONTINUE                                                                  
C                                                                               
      NMES = NMES + 1                                                           
C  Ionospheric correction applied
      IF (ICA .EQ. 1) THEN                                                      
         RRATE = RRATE + DBLE(IRC)
         ICA = 0                                                                
      END IF                                                                    
      IONRC = - DBLE(IRC) * DMICRO
C  Tropospheric correction applied
      IF (ITCA.EQ.1) THEN                                                     
         IF (ITRC.NE.0.AND..NOT.ADJUST) THEN
            RRATE = RRATE + DBLE(ITRC)
            ITCA = 0                                                         
         ELSE                                                                   
            IF (ITEMP .NE. 0) THEN                                              
               ITCA = 5                                                         
            END IF                                                              
         END IF                                                                 
      END IF                                                                    
C  Center of mass correction applied
      IF (.NOT.ADJUST) THEN
         IF (IPCE.NE.0) THEN
            RRATE = RRATE + DBLE(IPCE)
            ITOI = 0
         END IF
      END IF
c     AOC = - DBLE(IPCE) * DMICRO
      IPRE = ((16*ITOI+8*IDT+2*IC)*M8+(8*ITCA+4*ICA))*M16
C  Fix the millennium bug
      IF(IYEAR.LT.85) IYEAR = IYEAR + 100
C
      MJD = MJCENT + IDAYS + (1461*IYEAR - 1)/4                                 
C  Time tag of observation is at end of count interval                          
      FRACT=((DBLE(IFRACT)+DBLE(DCOUNT)/DT)*DMICRO+DBLE(ISEC))/SPD
      IF (FRACT.GE.1.D0) THEN                                                   
         MJD=MJD+1                                                              
         FRACT=FRACT-1.D0                                                       
      END IF                                                                    
      DRRATE = RRATE * DMICRO                                                   
      OSD = DBLE(IOSD) * DMICRO
      ICOUNT = NINT(DBLE(DCOUNT)*DT)  ! TYPE 42
c     ICOUNT = INT(DBLE(DCOUNT)*DT)
      if (checkdates) then
        tsec = sec85(1,dble(mjd))
        if (tsec.lt.t0 .or. tsec.gt.t1) then
         NSKIP = NSKIP + 1                                                      
         NSKTM = NSKTM + 1
         GOTO 40                                                                
        end if
      end if
C                                                                               
C Write the data to a GEODYN BINARY FORMATED file                               
C                                                                               
      IF (ITCA.EQ.5) THEN                                    
         ITROPC = IPRES + (16*ITEMP + M16*IRHUM)*M8                             
         WRITE(OUNIT)ISAT,MTYP,TIMIND,ISTAT,IPRE,MJD,FRACT,DRRATE,              
     .               IDUM,IDUM,OSD,ICOUNT,ITROPC,IONRC,AOC,BOC                
      ELSE                                                                      
         TROPRC = - DBLE(ITRC) * DMICRO
         WRITE (OUNIT)ISAT,MTYP,TIMIND,ISTAT,IPRE,MJD,FRACT,DRRATE,
     .                IDUM,IDUM,OSD,ICOUNT,TROPRC,IONRC,AOC,BOC
      END IF                                                                    
C                                                                               
C Calculate the statistics for each pass
C                                                                               
      IF (.NOT.STATUS) GOTO 40
      IF (STATID.EQ.STATOLD) THEN
	 IF (FIRST) THEN
	   CALL MJDATE(1,MJD,IDUM,IYY,IMM,IDD)
	   DSEC=DBLE(ISEC)
	   IH1=INT(DSEC/3600.D0)
	   IM1=INT(DSEC/60.D0-IH1*60)
  	   IS1=ISEC-(IH1*60+IM1)*60  
	   FIRST=.FALSE.
	 END IF
	 RSUMI=RSUMI+(IRC*IRC)*1.D-06
	 RSUMT=RSUMT+(ITRC*ITRC)*1.D-06
	 RSUMP=RSUMP+(IPCE*IPCE)*1.D-06
	 DSOLD=DBLE(ISEC)
	 NP=NP+1
      ELSE
	 RSUMI=SQRT(RSUMI/NP)
	 RSUMT=SQRT(RSUMT/NP)
	 RSUMP=SQRT(RSUMP/NP)
	 IH2=INT(DSOLD/3600.D0)
	 IM2=INT(DSOLD/60.D0-IH2*60)
  	 IS2=INT(DSOLD)-(IH2*60+IM2)*60  
  	 WRITE(MUNIT,90)SATID,STATOLD,IYY,IMM,IDD,IH1,IM1,IS1,IH2,IM2,IS2,   
     .			NP,RSUMI,RSUMT,RSUMP
 90	 FORMAT(1X,A7,A8,I5,I4,'/',I2,1X,3I3,2X,3I3,I7,3F8.3)
	 IF (NP.LE.9) THEN
 91	    FORMAT(1X,A7,A8,I5,I4,'/',I2,1X,3I3,2X,3I3,I7)
	    IM1 = IM1 - 2
	    IF (IM1.LT.0) THEN
	       IM1=IM1+60
	       IH1=IH1-1
	       IF (IH1.LT.0) THEN
		  IH2=IH1+24
		  IDD=IDD-1
	       END IF
	    END IF
	    IM2 = IM2 + 2
	    IF (IM2.GE.60) THEN
	       IM2=IM2-60
	       IH2=IH2+1
	       IF (IH2.GE.24) THEN
		  IH2=IH2-24
		  IDD=IDD+1
	       END IF
	    END IF
	 END IF
	 STATOLD = STATID
	 ISTOLD =ISTAT
	 CALL MJDATE(1,MJD,IDUM,IYY,IMM,IDD)
	 DSEC=DBLE(ISEC)
	 IH1=INT(DSEC/3600.D0)
	 IM1=INT(DSEC/60.D0-IH1*60)
  	 IS1=ISEC-(IH1*60+IM1)*60  
	 RSUMI=(IRC*IRC)*1.D-12
	 RSUMT=(ITRC*ITRC)*1.D-12
	 RSUMP=(IPCE*IPCE)*1.D-12
	 NP=1
      END IF 
C
      GOTO 40                                                                   
C                                                                               
 100  CONTINUE                                                                  
C
C  Write statistics of last pass 
C
      IF (STATUS) THEN
	 RSUMI=SQRT(RSUMI/NP)
	 RSUMT=SQRT(RSUMT/NP)
	 RSUMP=SQRT(RSUMP/NP)
	 IH2=INT(DSOLD/3600.D0)
	 IM2=INT(DSOLD/60.D0-IH2*60)
  	 IS2=INT(DSOLD)-(IH2*60+IM2)*60  
  	 WRITE(MUNIT,90)SATID,STATOLD,IYY,IMM,IDD,IH1,IM1,IS1,IH2,IM2,IS2,   
     .   		NP,RSUMI,RSUMT,RSUMP
      END IF
C                                                                               
      WRITE (6,110) IMAX, NREC, NSKTM, NSKED, NSKST, NSKIP, NMES
 110  FORMAT(' ','Total number of stations     : ',I6/,
     .       ' ','Total number of records      : ',I6/,                       
     .       ' ',' - not in time interval      : ',I6/,
     .       ' ',' - edited by CNES            : ',I6/,                         
     .       ' ',' - station unknown           : ',I6/,                         
     .       ' ','Total number of skiped lines : ',I6/,
     .       ' ','Total number of measurements : ',I6)                        
C                                                                               
9999  END                                                                       
