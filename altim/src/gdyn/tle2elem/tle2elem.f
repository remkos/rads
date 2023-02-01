C      DRIVER                                                 3 NOV 80
C
C   WGS-72 PHYSICAL AND GEOPOTENTIAL CONSTANTS
C          CK2= .5*J2*AE**2     CK4=-.375*J4*AE**4
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  This program generates a state-vector for Geodyn II, based on NORAD
C  two line elements.
C
C  This is modified version of the original Fortran IV code of SPACETRACK.
C  Modifications have been made in terms of Fortran coding only;
C  - elimination of entries
C  - addition of common blocks
C  - addition of block data
C  - turning off of checks for initialization
C  - all parameters 64 bits
C  No substantial alterations were made.
C  After these modifications, the test results as mentioned in the document
C  "Spacetrack Report No. 3" (state vector after integration over 24 minutes)
C  could be reproduced at the level of 1-7 m (equivalent to 1 part in 7)
C  for each of the 5 integrators.
C
C  input:
C
C  output:
C
C  R. Noomen. April 9, 1999
C  E. Doornbos. May 27, 2002 : Input using command-line parameters. 
C                              Additional satids for altimetry.
C                              Output to stdout, errors to stderr. Other suppressed.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
      CHARACTER * 4 CTYPE
C
      PARAMETER ( MAXSAT = 12 )
C
      CHARACTER*8 satnam(MAXSAT)
      CHARACTER*255 filename, altim
      integer*4 length
C
      DIMENSION ISATNO(MAXSAT) , ISATCO(MAXSAT) ,
     .          AREA(MAXSAT) , DMASS(MAXSAT)
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
      DATA QO,SO,XJ2,XJ4/
     2            120.0,78.0,1.082616E-3,-1.65597E-6/
C
C  data for 
C  Starlette, Lageos-1, EGP, Lageos-2,
C  Stella, GFZ-1, Westpac,
C  ERS-1, ERS-2, Envisat,
C  TOPEX/POSEIDON, Jason-1
      DATA satnam / "starlette", "lageos-1" , "egp" ,    "lageos-2" ,
     .              "stella",    "gfz-1",     "westpac",
     .              "ers-1",     "ers-2",     "envisat",
     .              "topex",     "jason-1" /
C??? NORAD catalog numbers
      DATA ISATNO / 7646 , 8820 , 16908 , 22195 ,
     .              22824 , 23558 , 25398 ,
     .              21574 , 23560 , 27386 , 
     .              22076 , 26997 /
C??? cospar id for GFZ-1 ????
      DATA ISATCO / 7501001 , 7603901 , 8606101 , 9207002 ,
     .              9306102 , 8601799 , 9804301 ,
     .              9105001 , 9502101 , 0200901 ,
     .              9205201 , 0105501 /
C??? area for satellites ....???
      DATA AREA   / 0.0D0 , 0.283D0 , 0.0D0 , 0.283D0 ,
     .              0.0D0 , 0.0D0   , 0.04714D0 ,
     .              21.4D0, 21.4D0  , 40.0D0 ,
     .              0.0D0 , 0.0D0 /
C??? mass for satellites ....???
      DATA DMASS  / 0.0D0    , 407.0D0 , 0.0D0 , 407.0D0 ,
     .              0.0D0    , 0.0D0 , 23.5D0 ,
     .              2312.0D0 , 2496.0D0, 8078.070D0,
     .              0.0D0    , 481.0D0 /
C
      DATA DKM2M / 1.0D3 /
C
C     SELECT EPHEMERIS TYPE AND OUTPUT TIMES
C
      CK2=.5*XJ2*AE**2
      CK4=-.375*XJ4*AE**4
      QOMS2T=((QO-SO)*AE/XKMPER)**4
      S=AE*(1.+SO/XKMPER)
C
C  read input parameters
C      CALL READ5 ( ISATREQ , IYMD , IHMS , CTYPE )
      CALL READARG ( ISATREQ , IYMD , IHMS, CTYPE , IYMD2, IHMS2 )
      do i=1,maxsat
        if (isatco(i).eq.isatreq) then 
          isatreq=isatno(i)
          j=lnblnk(satnam(i))
          altim = "/uwa/altim"
          call checkenv("ALTIM",altim,length)
          filename = altim(1:length) //
     |              "/data/tle/"//satnam(i)(:j)//".txt"
        end if
      enddo
C
C  read state-vector which matches required epoch best
      open (unit=10,file=filename)
      CALL READ10 ( IYMD , IHMS , ISATREQ , EPOCHREQ , EPOCHTLE )
      close(10)
C
      IDEEP=0
C
      IF(XNO.LE.0.) STOP
      XNDD6O=XNDD6O*(10.**IEXP)
      XNODEO=XNODEO*DE2RA
      OMEGAO=OMEGAO*DE2RA
      XMO=XMO*DE2RA
      XINCL=XINCL*DE2RA
      TEMP=TWOPI/XMNPDA/XMNPDA
      XNO=XNO*TEMP*XMNPDA
      XNDT2O=XNDT2O*TEMP
      XNDD6O=XNDD6O*TEMP/XMNPDA
C
C     INPUT CHECK FOR PERIOD VS EPHEMERIS SELECTED
C     PERIOD GE 225 MINUTES  IS DEEP SPACE
C
      A1=(XKE/XNO)**TOTHRD
      TEMP=1.5*CK2*(3.*COS(XINCL)**2-1.)/(1.-EO*EO)**1.5
      DEL1=TEMP/(A1*A1)
      AO=A1*(1.-DEL1*(.5*TOTHRD+DEL1*(1.+134./81.*DEL1)))
      DELO=TEMP/(AO*AO)
      XNODP=XNO/(1.+DELO)
      IF((TWOPI/XNODP/XMNPDA) .GE. .15625) IDEEP=1
C
      BSTAR=BSTAR*(10.**IBEXP)/AE
C  TSINCE is difference w.r.t. epoch of TLE [min]
      TSINCE = EPOCHREQ - EPOCHTLE
      TSINCE = TSINCE * 1440.0D0
      IFLAG=1
      IF(IDEEP .EQ. 1 .AND. (CTYPE .EQ. 'SGP ' .OR. CTYPE .EQ. 'SGP4'
     1           .OR. CTYPE .EQ. 'SGP8')) THEN
        WRITE(0,930)
  930   FORMAT("SHOULD USE DEEP SPACE EPHEMERIS")
        STOP 99
      END IF
      IF (IDEEP .EQ. 0 .AND.
     .   (CTYPE .EQ. 'SDP4' .OR. CTYPE .EQ. 'SDP8')) THEN
        WRITE(0,940)
  940   FORMAT("SHOULD USE NEAR EARTH EPHEMERIS")
        STOP 99
      END IF
      IF ( CTYPE .EQ. 'SGP ' ) CALL SGP(IFLAG,TSINCE)
      IF ( CTYPE .EQ. 'SGP4' ) CALL SGP4(IFLAG,TSINCE)
      IF ( CTYPE .EQ. 'SGP8' ) CALL SGP8(IFLAG,TSINCE)
      IF ( CTYPE .EQ. 'SDP4' ) CALL SDP4(IFLAG,TSINCE)
      IF ( CTYPE .EQ. 'SDP8' ) CALL SDP8(IFLAG,TSINCE)
      X=X*XKMPER/AE
      Y=Y*XKMPER/AE
      Z=Z*XKMPER/AE
      XDOT=XDOT*XKMPER/AE*XMNPDA/86400.
      YDOT=YDOT*XKMPER/AE*XMNPDA/86400.
      ZDOT=ZDOT*XKMPER/AE*XMNPDA/86400.
C  convert to [m] and [m/s]
      X    = X    * DKM2M
      Y    = Y    * DKM2M
      Z    = Z    * DKM2M
      XDOT = XDOT * DKM2M
      YDOT = YDOT * DKM2M
      ZDOT = ZDOT * DKM2M
C  print resulting state-vector in Geodyn II format
      ISATP  = 9999999
      AREAP  = 9999.0D0
      DMASSP = 9999.0D0
      DO I = 1 , MAXSAT
        IF ( ISATREQ .EQ. ISATNO(I) ) THEN
          ISATP  = ISATCO(I)
          AREAP  = AREA(I)
          DMASSP = DMASS(I)
        END IF
      END DO
      WRITE  ( 6 , 801 ) ISATP , AREAP , DMASSP
 801  FORMAT ( 'SATPAR0' , 10X , 1I7.7 , 1F20.5 , 1F15.5 )
      DHMS = DBLE ( IHMS )
      DHMS2 = DBLE ( IHMS2 )
      WRITE  ( 6 , 802 ) IYMD , DHMS , IYMD , DHMS , IYMD2 , DHMS2
 802  FORMAT ( 'EPOCH' , 15X , 3 ( 1I6.6 , 1F14.7 ) )
      WRITE  ( 6 , 803 ) X , Y , Z
 803  FORMAT ( 'ELEMS110   300' , 6X , 3D20.14 )
      WRITE  ( 6 , 804 ) XDOT , YDOT , ZDOT
 804  FORMAT ( 'ELEMS2        ' , 6X , 3D20.14 )
C
      STOP
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      BLOCK DATA
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
      DATA DE2RA,E6A,PI,PIO2,TOTHRD,TWOPI,X3PIO2,XJ3,
     1            XKE,XKMPER,XMNPDA,AE/.174532925E-1,1.E-6,
     2            3.14159265,1.57079633,.66666667,
     4            6.2831853,4.71238898,-.253881E-5,
     5            .743669161E-1,6378.135,1440.,1./
C
      DATA EPS / 1.0D-10 /
C
      DATA              ZNS,           C1SS,          ZES/
     A                  1.19459E-5,    2.9864797E-6, .01675/
      DATA              ZNL,           C1L,           ZEL/
     A                  1.5835218E-4,  4.7968065E-7,  .05490/
      DATA              ZCOSIS,        ZSINIS,        ZSINGS/
     A                  .91744867,     .39785416,     -.98088458/
      DATA              ZCOSGS,        ZCOSHS,        ZSINHS/
     A                  .1945905,      1.0,           0.0/
C
      DATA Q22,Q31,Q33/1.7891679E-6,2.1460748E-6,2.2123015E-7/
      DATA G22,G32/5.7686396,0.95240898/
      DATA G44,G52/1.8014998,1.0508330/
      DATA G54/4.4108898/
      DATA ROOT22,ROOT32/1.7891679E-6,3.7393792E-7/
      DATA ROOT44,ROOT52/7.3636953E-9,1.1428639E-7/
      DATA ROOT54/2.1765803E-9/
      DATA THDT/4.3752691E-3/
      DATA RHO/.15696615/
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE READARG (ISATREQ,IYMD,IHMS,CTYPE,IYMD2,IHMS2)
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
      INTEGER*4 iargc, iarg
      REAL*8 IYMDHMS, IYMDHMS2
      CHARACTER * 4 CTYPE
      CHARACTER*255 arg, program

      call getarg(0,program)
      l=index(program,' ')-1
      if (iargc().lt.1) then
         write (*,600) program(:l)
         stop
      endif

600   format (/
     |'syntax: ',a,' SGPn satid=xxxxxx ymd=yymmddhhmmss ')

      do iarg=1,iargc()
      call getarg(iarg,arg)
      if (arg(:3).eq.'SGP') then
         CTYPE = arg(:4)
      else if (arg(:6).eq.'satid=') then
         read (arg(7:),*) ISATREQ
      else if (arg(:4).eq.'ymd=') then
         read (arg(5:),*,iostat=ios) IYMDHMS,IYMDHMS2
         IYMD = INT(IYMDHMS/1000000.0D0)
         IHMS = INT(IYMDHMS - 1000000.0D0 * IYMD)
         IYMD2 = INT(IYMDHMS2/1000000.0D0)
         IHMS2 = INT(IYMDHMS2 - 1000000.0D0 * IYMD2)
      endif
      enddo 

      IF ( CTYPE .NE. 'SGP ' .AND.
     .     CTYPE .NE. 'SGP4' .AND.
     .     CTYPE .NE. 'SGP8' .AND.
     .     CTYPE .NE. 'SDP4' .AND.
     .     CTYPE .NE. 'SDP8' ) THEN
        WRITE ( 0 , * ) 'READ5: incorrect value for CTYPE: ' , CTYPE
        STOP 99
      END IF

      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE READ5 ( ISATREQ , IYMD , IHMS , CTYPE )
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
      CHARACTER * 4 CTYPE
C
      OPEN ( 5 , FILE = 'fort.5' , BLANK = 'ZERO' )
C
C  read orbit integrator
      READ ( 5 , * , END = 999 ) CTYPE
      IF ( CTYPE .NE. 'SGP ' .AND.
     .     CTYPE .NE. 'SGP4' .AND.
     .     CTYPE .NE. 'SGP8' .AND.
     .     CTYPE .NE. 'SDP4' .AND.
     .     CTYPE .NE. 'SDP8' ) THEN
        WRITE ( 0 , * ) 'READ5: incorrect value for CTYPE: ' , CTYPE
        STOP 99
      END IF
C  read satellite id
      READ ( 5 , * , END = 999 ) ISATREQ
C  read epoch of required state-vector
      READ ( 5 , * , END = 999 ) IYMD , IHMS
C
      RETURN
C
 999  CALL UNEOF ( 6 , 'READ5 ' , 5 )
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE READ10 ( IYMD , IHMS , ISATREQ , EPOCHREQ , EPOCHTLE )
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
      CHARACTER * 80 C80
      LOGICAL L1ST
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
C
*      WRITE ( 0 , * ) 'READ10: requested epoch = ' , IYMD , IHMS
      DHMS = DBLE ( IHMS )
      CALL YMDT ( IYMD , DHMS , EPOCHREQ )
C
      L1ST = .TRUE.
C  read state-vectors
  10  READ   ( 10 , 20 , END = 100 ) C80
  20  FORMAT ( 1A80 )
      IF ( C80(1:2) .EQ. '1 ' ) THEN
        READ   ( C80 , 30 ) ISATTLE
  30    FORMAT ( 2X , 1I5 )
        IF ( ISATTLE .EQ. ISATREQ ) THEN
          READ   ( C80 , 40 ) IY , IDAY , TFRAC
  40      FORMAT ( 18X , 1I2 , 1I3 , 1F9.8 )
          IYMDTEST = IY * 10000 + 101
          CALL MJDATE ( 2 , MJDTEST , IYMDTEST , IYY , IMM , IDD )
          EPOCHTEST = DBLE ( MJDTEST + IDAY - 1 ) + TFRAC
          TDIFF = ABS ( EPOCHTEST - EPOCHREQ )
          IF ( L1ST ) THEN
            L1ST = .FALSE.
            TDIFFMIN = TDIFF
            EPOCHTLE = EPOCHTEST
C     READ IN MEAN ELEMENTS FROM 2 CARD T(TRANS) OR G(INTERN) FORMAT
            IF ( C80(80:80) .EQ. 'G' ) THEN
              READ  ( C80,701) XMO,XNODEO,OMEGAO
  701         FORMAT(43X,1X,3F8.4 )
              READ   ( 10 , 20 , END = 999 ) C80
              READ   ( C80 , 30 ) ISATTLE
              IF ( C80(1:2) .NE. '2 ' .OR. ISATTLE .NE. ISATREQ ) THEN
               WRITE ( 0 , * ) 'READ10: error in TLEs:'
               WRITE ( 0 , * ) C80
                STOP 99
	      END IF
              READ  ( C80,711) EO,XINCL,XNO,XNDT2O,
     1                XNDD6O,IEXP,BSTAR,IBEXP
  711         FORMAT(6X,F8.7,F8.4,1X,2F11.9,1X,F6.5,I2,
     1             4X,F8.7,I2)
            ELSE
              READ   ( C80,702) XNDT2O,XNDD6O,IEXP,BSTAR,IBEXP
  702         FORMAT(33X,F10.8,2(1X,F6.5,I2) )
              READ   ( 10 , 20 , END = 999 ) C80
              READ   ( C80 , 30 ) ISATTLE
              IF ( C80(1:2) .NE. '2 ' .OR. ISATTLE .NE. ISATREQ ) THEN
               WRITE ( 0 , * ) 'READ10: error in TLEs:'
               WRITE ( 0 , * ) C80
                STOP 99
	      END IF
              READ   ( C80,712) XINCL, XNODEO,EO,OMEGAO,XMO,XNO
  712         FORMAT(7X,2(1X,F8.4),1X,
     1             F7.7,2(1X,F8.4),1X,F11.8)
            END IF
            GOTO 10
          ELSE
            IF ( TDIFF .LT. TDIFFMIN ) THEN
              TDIFFMIN = TDIFF
              EPOCHTLE = EPOCHTEST
C     READ IN MEAN ELEMENTS FROM 2 CARD T(TRANS) OR G(INTERN) FORMAT
              IF ( C80(80:80) .EQ. 'G' ) THEN
                READ  ( C80,701) XMO,XNODEO,OMEGAO
                READ   ( 10 , 20 , END = 999 ) C80
                READ   ( C80 , 30 ) ISATTLE
                IF ( C80(1:2) .NE. '2 ' .OR. ISATTLE .NE. ISATREQ ) THEN
                 WRITE ( 0 , * ) 'READ10: error in TLEs:'
                 WRITE ( 0 , * ) C80
                  STOP 99
	        END IF
                READ  ( C80,711) EO,XINCL,XNO,XNDT2O,
     1                XNDD6O,IEXP,BSTAR,IBEXP
              ELSE
                READ   ( C80,702) XNDT2O,XNDD6O,IEXP,BSTAR,IBEXP
                READ   ( 10 , 20 , END = 999 ) C80
                READ   ( C80 , 30 ) ISATTLE
                IF ( C80(1:2) .NE. '2 ' .OR. ISATTLE .NE. ISATREQ ) THEN
                 WRITE ( 0 , * ) 'READ10: error in TLEs:'
                 WRITE ( 0 , * ) C80
                  STOP 99
  	        END IF
                READ   ( C80,712) XINCL, XNODEO,EO,OMEGAO,XMO,XNO
              END IF
            ELSE
C  skip next record
              READ ( 10 , 20 , END = 999 )
            END IF
            GOTO 10
          END IF
        ELSE
          READ ( 10 , 20 ) C80
        END IF
      END IF
      GOTO 10
C
 100  CONTINUE
C
      CALL TYMD ( EPOCHTLE , IYMDTLE , DHMSTLE )
*      WRITE ( 0 , * ) 'READ10: closest epoch   = ' , IYMDTLE , DHMSTLE
*      WRITE ( 0 , * ) '(propagated over difference)'
      IF ( L1ST ) THEN
       WRITE ( 0 , * ) 'READ10: no matching TLE could be found ' ,
     .                  'for satellite ' , ISATREQ
        STOP 99
      END IF
C
      RETURN
C
 999  CALL UNEOF ( 6 , 'READ10' , 10 )
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SDP4(IFLAG,TSINCE)
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
C  turned off. R. Noomen, March 5, 1999
C     IF  (IFLAG .EQ. 0) GO TO 100
C
C      RECOVER ORIGINAL MEAN MOTION (XNODP) AND SEMIMAJOR AXIS (AODP)
C      FROM INPUT ELEMENTS
C
      A1=(XKE/XNO)**TOTHRD
      COSIO=COS(XINCL)
      THETA2=COSIO*COSIO
      X3THM1=3.*THETA2-1.
      EOSQ=EO*EO
      BETAO2=1.-EOSQ
      BETAO=SQRT(BETAO2)
      DEL1=1.5*CK2*X3THM1/(A1*A1*BETAO*BETAO2)
      AO=A1*(1.-DEL1*(.5*TOTHRD+DEL1*(1.+134./81.*DEL1)))
      DELO=1.5*CK2*X3THM1/(AO*AO*BETAO*BETAO2)
      XNODP=XNO/(1.+DELO)
      AODP=AO/(1.-DELO)
C
C      INITIALIZATION
C
C      FOR PERIGEE BELOW 156 KM, THE VALUES OF
C      S AND QOMS2T ARE ALTERED
C
      S4=S
      QOMS24=QOMS2T
      PERIGE=(AODP*(1.-EO)-AE)*XKMPER
      IF(PERIGE .GE. 156.) GO TO 10
      S4=PERIGE-78.
      IF(PERIGE .GT. 98.) GO TO 9
      S4=20.
    9 QOMS24=((120.-S4)*AE/XKMPER)**4
      S4=S4/XKMPER+AE
   10 PINVSQ=1./(AODP*AODP*BETAO2*BETAO2)
      SING=SIN(OMEGAO)
      COSG=COS(OMEGAO)
      TSI=1./(AODP-S4)
      ETA=AODP*EO*TSI
      ETASQ=ETA*ETA
      EETA=EO*ETA
      PSISQ=ABS(1.-ETASQ)
      COEF=QOMS24*TSI**4
      COEF1=COEF/PSISQ**3.5
      C2=COEF1*XNODP*(AODP*(1.+1.5*ETASQ+EETA*(4.+ETASQ))+.75*
     1         CK2*TSI/PSISQ*X3THM1*(8.+3.*ETASQ*(8.+ETASQ)))
      C1=BSTAR*C2
      SINIO=SIN(XINCL)
      A3OVK2=-XJ3/CK2*AE**3
      X1MTH2=1.-THETA2
      C4=2.*XNODP*COEF1*AODP*BETAO2*(ETA*
     1         (2.+.5*ETASQ)+EO*(.5+2.*ETASQ)-2.*CK2*TSI/
     2         (AODP*PSISQ)*(-3.*X3THM1*(1.-2.*EETA+ETASQ*
     3         (1.5-.5*EETA))+.75*X1MTH2*(2.*ETASQ-EETA*
     4         (1.+ETASQ))*COS(2.*OMEGAO)))
      THETA4=THETA2*THETA2
      TEMP1=3.*CK2*PINVSQ*XNODP
      TEMP2=TEMP1*CK2*PINVSQ
      TEMP3=1.25*CK4*PINVSQ*PINVSQ*XNODP
      XMDOT=XNODP+.5*TEMP1*BETAO*X3THM1+.0625*TEMP2*BETAO*
     1         (13.-78.*THETA2+137.*THETA4)
      X1M5TH=1.-5.*THETA2
      OMGDOT=-.5*TEMP1*X1M5TH+.0625*TEMP2*(7.-114.*THETA2+
     1         395.*THETA4)+TEMP3*(3.-36.*THETA2+49.*THETA4)
      XHDOT1=-TEMP1*COSIO
      XNODOT=XHDOT1+(.5*TEMP2*(4.-19.*THETA2)+2.*TEMP3*(3.-
     1         7.*THETA2))*COSIO
      XNODCF=3.5*BETAO2*XHDOT1*C1
      T2COF=1.5*C1
      XLCOF=.125*A3OVK2*SINIO*(3.+5.*COSIO)/(1.+COSIO)
      AYCOF=.25*A3OVK2*SINIO
      X7THM1=7.*THETA2-1.
   90 IFLAG=0
      CALL DPINIT(EOSQ,SINIO,COSIO,BETAO,AODP,THETA2,
     1         SING,COSG,BETAO2,XMDOT,OMGDOT,XNODOT,XNODP)
C
C      UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
C
  100 XMDF=XMO+XMDOT*TSINCE
      OMGADF=OMEGAO+OMGDOT*TSINCE
      XNODDF=XNODEO+XNODOT*TSINCE
      TSQ=TSINCE*TSINCE
      XNODE=XNODDF+XNODCF*TSQ
      TEMPA=1.-C1*TSINCE
      TEMPE=BSTAR*C4*TSINCE
      TEMPL=T2COF*TSQ
      XN=XNODP
C     CALL DPSEC(XMDF,OMGADF,XNODE,EM,XINC,XN,TSINCE)
      CALL DPSEC(XMDF,OMGADF,XNODE,EM,XINC,XN,TSINCE,OMGDOT)
      A=(XKE/XN)**TOTHRD*TEMPA**2
      E=EM-TEMPE
      XMAM=XMDF+XNODP*TEMPL
C     CALL DPPER(E,XINC,OMGADF,XNODE,XMAM)
      CALL DPPER(E,XINC,OMGADF,XNODE,XMAM,TSINCE,COSIO,SINIO)
      XL=XMAM+OMGADF+XNODE
      BETA=SQRT(1.-E*E)
      XN=XKE/A**1.5
C
C      LONG PERIOD PERIODICS
C
      AXN=E*COS(OMGADF)
      TEMP=1./(A*BETA*BETA)
      XLL=TEMP*XLCOF*AXN
      AYNL=TEMP*AYCOF
      XLT=XL+XLL
      AYN=E*SIN(OMGADF)+AYNL
C
C      SOLVE KEPLERS EQUATION
C
      CAPU=MOD(XLT-XNODE , TWOPI )
      TEMP2=CAPU
      DO 130 I=1,10
      SINEPW=SIN(TEMP2)
      COSEPW=COS(TEMP2)
      TEMP3=AXN*SINEPW
      TEMP4=AYN*COSEPW
      TEMP5=AXN*COSEPW
      TEMP6=AYN*SINEPW
      EPW=(CAPU-TEMP4+TEMP3-TEMP2)/(1.-TEMP5-TEMP6)+TEMP2
      IF(ABS(EPW-TEMP2) .LE. E6A) GO TO 140
  130 TEMP2=EPW
C
C      SHORT PERIOD PRELIMINARY QUANTITIES
C
  140 ECOSE=TEMP5+TEMP6
      ESINE=TEMP3-TEMP4
      ECCANOMRON=ATAN2(ESINE,ECOSE)
      SMALLOMEG=EPW-ECCANOMRON
      ELSQ=AXN*AXN+AYN*AYN
      TEMP=1.-ELSQ
      PL=A*TEMP
      R=A*(1.-ECOSE)
      TEMP1=1./R
      RDOT=XKE*SQRT(A)*ESINE*TEMP1
      RFDOT=XKE*SQRT(PL)*TEMP1
      TEMP2=A*TEMP1
      BETAL=SQRT(TEMP)
      TEMP3=1./(1.+BETAL)
      COSU=TEMP2*(COSEPW-AXN+AYN*ESINE*TEMP3)
      SINU=TEMP2*(SINEPW-AYN-AXN*ESINE*TEMP3)
      U=ATAN2(SINU,COSU)
      SIN2U=2.*SINU*COSU
      COS2U=2.*COSU*COSU-1.
      TEMP=1./PL
      TEMP1=CK2*TEMP
      TEMP2=TEMP1*TEMP
C
C      UPDATE FOR SHORT PERIODICS
C
      RK=R*(1.-1.5*TEMP2*BETAL*X3THM1)+.5*TEMP1*X1MTH2*COS2U
      UK=U-.25*TEMP2*X7THM1*SIN2U
      XNODEK=XNODE+1.5*TEMP2*COSIO*SIN2U
      XINCK=XINC+1.5*TEMP2*COSIO*SINIO*COS2U
      RDOTK=RDOT-XN*TEMP1*X1MTH2*SIN2U
      RFDOTK=RFDOT+XN*TEMP1*(X1MTH2*COS2U+1.5*X3THM1)
C
C      ORIENTATION VECTORS
C
      SINUK=SIN(UK)
      COSUK=COS(UK)
      SINIK=SIN(XINCK)
      COSIK=COS(XINCK)
      SINNOK=SIN(XNODEK)
      COSNOK=COS(XNODEK)
      XMX=-SINNOK*COSIK
      XMY=COSNOK*COSIK
      UX=XMX*SINUK+COSNOK*COSUK
      UY=XMY*SINUK+SINNOK*COSUK
      UZ=SINIK*SINUK
      VX=XMX*COSUK-COSNOK*SINUK
      VY=XMY*COSUK-SINNOK*SINUK
      VZ=SINIK*COSUK
C
C      POSITION AND VELOCITY
C
      X=RK*UX
      Y=RK*UY
      Z=RK*UZ
      XDOT=RDOTK*UX+RFDOTK*VX
      YDOT=RDOTK*UY+RFDOTK*VY
      ZDOT=RDOTK*UZ+RFDOTK*VZ
C
      RETURN
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SDP8(IFLAG,TSINCE)
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
C  turned off. R. Noomen, March 5, 1999
C     IF  (IFLAG .EQ. 0) GO TO 100
C
C      RECOVER ORIGINAL MEAN MOTION (XNODP) AND SEMIMAJOR AXIS (AODP)
C      FROM INPUT ELEMENTS --------- CALCULATE BALLISTIC COEFFICIENT
C      (B TERM) FROM INPUT B* DRAG TERM
C
      A1=(XKE/XNO)**TOTHRD
      COSI=COS(XINCL)
      THETA2=COSI*COSI
      TTHMUN=3.*THETA2-1.
      EOSQ=EO*EO
      BETAO2=1.-EOSQ
      BETAO=SQRT(BETAO2)
      DEL1=1.5*CK2*TTHMUN/(A1*A1*BETAO*BETAO2)
      AO=A1*(1.-DEL1*(.5*TOTHRD+DEL1*(1.+134./81.*DEL1)))
      DELO=1.5*CK2*TTHMUN/(AO*AO*BETAO*BETAO2)
      AODP=AO/(1.-DELO)
      XNODP=XNO/(1.+DELO)
      B=2.*BSTAR/RHO
C
C      INITIALIZATION
C
      PO=AODP*BETAO2
      POM2=1./(PO*PO)
      SINI=SIN(XINCL)
      SING=SIN(OMEGAO)
      COSG=COS(OMEGAO)
      TEMP=.5*XINCL
      SINIO2=SIN(TEMP)
      COSIO2=COS(TEMP)
      THETA4=THETA2**2
      UNM5TH=1.-5.*THETA2
      UNMTH2=1.-THETA2
      A3COF=-XJ3/CK2*AE**3
      PARDT1=3.*CK2*POM2*XNODP
      PARDT2=PARDT1*CK2*POM2
      PARDT4=1.25*CK4*POM2*POM2*XNODP
      XMDT1=.5*PARDT1*BETAO*TTHMUN
      XGDT1=-.5*PARDT1*UNM5TH
      XHDT1=-PARDT1*COSI
      XLLDOT=XNODP+XMDT1+
     2           .0625*PARDT2*BETAO*(13.-78.*THETA2+137.*THETA4)
      OMGDT=XGDT1+
     1      .0625*PARDT2*(7.-114.*THETA2+395.*THETA4)+PARDT4*(3.-36.*
     2         THETA2+49.*THETA4)
      XNODOT=XHDT1+
     1       (.5*PARDT2*(4.-19.*THETA2)+2.*PARDT4*(3.-7.*THETA2))*COSI
      TSI=1./(PO-S)
      ETA=EO*S*TSI
      ETA2=ETA**2
      PSIM2=ABS(1./(1.-ETA2))
      ALPHA2=1.+EOSQ
      EETA=EO*ETA
      COS2G=2.*COSG**2-1.
      D5=TSI*PSIM2
      D1=D5/PO
      D2=12.+ETA2*(36.+4.5*ETA2)
      D3=ETA2*(15.+2.5*ETA2)
      D4=ETA*(5.+3.75*ETA2)
      B1=CK2*TTHMUN
      B2=-CK2*UNMTH2
      B3=A3COF*SINI
      C0=.5*B*RHO*QOMS2T*XNODP*AODP*TSI**4*PSIM2**3.5/SQRT(ALPHA2)
      C1=1.5*XNODP*ALPHA2**2*C0
      C4=D1*D3*B2
      C5=D5*D4*B3
      XNDT=C1*(
     1  (2.+ETA2*(3.+34.*EOSQ)+5.*EETA*(4.+ETA2)+8.5*EOSQ)+
     1   D1*D2*B1+   C4*COS2G+C5*SING)
      XNDTN=XNDT/XNODP
      EDOT=-TOTHRD*XNDTN*(1.-EO)
      IFLAG=0
      CALL DPINIT(EOSQ,SINI,COSI,BETAO,AODP,THETA2,SING,COSG,
     1          BETAO2,XLLDOT,OMGDT,XNODOT,XNODP)
C
C      UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
C
  100 Z1=.5*XNDT*TSINCE*TSINCE
      Z7=3.5*TOTHRD*Z1/XNODP
      XMAMDF=XMO+XLLDOT*TSINCE
      OMGASM=OMEGAO+OMGDT*TSINCE+Z7*XGDT1
      XNODES=XNODEO+XNODOT*TSINCE+Z7*XHDT1
      XN=XNODP
C     CALL DPSEC(XMAMDF,OMGASM,XNODES,EM,XINC,XN,TSINCE)
      CALL DPSEC(XMAMDF,OMGASM,XNODES,EM,XINC,XN,TSINCE,OMGDT)
      XN=XN+XNDT*TSINCE
      EM=EM+EDOT*TSINCE
      XMAM=XMAMDF+Z1+Z7*XMDT1
C     CALL DPPER(EM,XINC,OMGASM,XNODES,XMAM)
      CALL DPPER(EM,XINC,OMGASM,XNODES,XMAM,TSINCE,COSI,SINI)
      XMAM=MOD(XMAM , TWOPI )
C
C      SOLVE KEPLERS EQUATION
C
      ZC2=XMAM+EM*SIN(XMAM)*(1.+EM*COS(XMAM))
      DO 130 I=1,10
      SINE=SIN(ZC2)
      COSE=COS(ZC2)
      ZC5=1./(1.-EM*COSE)
      CAPE=(XMAM+EM*SINE-ZC2)*
     1   ZC5+ZC2
      IF(ABS(CAPE-ZC2) .LE. E6A) GO TO 140
  130 ZC2=CAPE
C
C      SHORT PERIOD PRELIMINARY QUANTITIES
C
  140 AM=(XKE/XN)**TOTHRD
      BETA2M=1.-EM*EM
      SINOS=SIN(OMGASM)
      COSOS=COS(OMGASM)
      AXNM=EM*COSOS
      AYNM=EM*SINOS
      PM=AM*BETA2M
      G1=1./PM
      G2=.5*CK2*G1
      G3=G2*G1
      BETA=SQRT(BETA2M)
      G4=.25*A3COF*SINI
      G5=.25*A3COF*G1
      SNF=BETA*SINE*ZC5
      CSF=(COSE-EM)*ZC5
      FM=ATAN2(SNF,CSF)
      SNFG=SNF*COSOS+CSF*SINOS
      CSFG=CSF*COSOS-SNF*SINOS
      SN2F2G=2.*SNFG*CSFG
      CS2F2G=2.*CSFG**2-1.
      ECOSF=EM*CSF
      G10=FM-XMAM+EM*SNF
      RM=PM/(1.+ECOSF)
      AOVR=AM/RM
      G13=XN*AOVR
      G14=-G13*AOVR
      DR=G2*(UNMTH2*CS2F2G-3.*TTHMUN)-G4*SNFG
      DIWC=3.*G3*SINI*CS2F2G-G5*AYNM
      DI=DIWC*COSI
      SINI2=SIN(.5*XINC)
C
C      UPDATE FOR SHORT PERIOD PERIODICS
C
      SNI2DU=SINIO2*(
     1   G3*(.5*(1.-7.*THETA2)*SN2F2G-3.*UNM5TH*G10)-G5*SINI*CSFG*(2.+
     2         ECOSF))-.5*G5*THETA2*AXNM/COSIO2
      XLAMB=FM+OMGASM+XNODES+G3*(.5*(1.+6.*COSI-7.*THETA2)*SN2F2G-3.*
     1      (UNM5TH+2.*COSI)*G10)+G5*SINI*(COSI*AXNM/(1.+COSI)-(2.
     2      +ECOSF)*CSFG)
      Y4=SINI2*SNFG+CSFG*SNI2DU+.5*SNFG*COSIO2*DI
      Y5=SINI2*CSFG-SNFG*SNI2DU+.5*CSFG*COSIO2*DI
      R=RM+DR
      RDOT=XN*AM*EM*SNF/BETA+G14*(2.*G2*UNMTH2*SN2F2G+G4*CSFG)
      RVDOT=XN*AM**2*BETA/RM+
     1      G14*DR+AM*G13*SINI*DIWC
C
C      ORIENTATION VECTORS
C
      SNLAMB=SIN(XLAMB)
      CSLAMB=COS(XLAMB)
      TEMP=2.*(Y5*SNLAMB-Y4*CSLAMB)
      UX=Y4*TEMP+CSLAMB
      VX=Y5*TEMP-SNLAMB
      TEMP=2.*(Y5*CSLAMB+Y4*SNLAMB)
      UY=-Y4*TEMP+SNLAMB
      VY=-Y5*TEMP+CSLAMB
      TEMP=2.*SQRT(1.-Y4*Y4-Y5*Y5)
      UZ=Y4*TEMP
      VZ=Y5*TEMP
C
C      POSITION AND VELOCITY
C
      X=R*UX
      Y=R*UY
      Z=R*UZ
      XDOT=RDOT*UX+RVDOT*VX
      YDOT=RDOT*UY+RVDOT*VY
      ZDOT=RDOT*UZ+RVDOT*VZ
C
      RETURN
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SGP(IFLAG,TSINCE)
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1        XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
C  turned off. R. Noomen, March 5, 1999
C     IF  (IFLAG .EQ. 0) GO TO 19
C
C      INITIALIZATION
C
      C1= CK2*1.5
      C2= CK2/4.0
      C3= CK2/2.0
      C4= XJ3*AE**3/(4.0*CK2)
      COSIO=COS(XINCL)
      SINIO=SIN(XINCL)
      A1=(XKE/XNO)**TOTHRD
      D1=     C1/A1/A1*(3.*COSIO*COSIO-1.)/(1.-EO*EO)**1.5
      AO=A1*(1.-1./3.*D1-D1*D1-134./81.*D1*D1*D1)
      PO=AO*(1.-EO*EO)
      QO=AO*(1.-EO)
      XLO=XMO+OMEGAO+XNODEO
      D1O= C3 *SINIO*SINIO
      D2O= C2 *(7.*COSIO*COSIO-1.)
      D3O=C1*COSIO
      D4O=D3O*SINIO
      PO2NO=XNO/(PO*PO)
      OMGDT=C1*PO2NO*(5.*COSIO*COSIO-1.)
      XNODOT=-2.*D3O*PO2NO
      C5=.5*C4*SINIO*(3.+5.*COSIO)/(1.+COSIO)
      C6=C4*SINIO
      IFLAG=0
C
C      UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
C
   19 A=XNO+(2.*XNDT2O+3.*XNDD6O*TSINCE)*TSINCE
      A=AO*(XNO/A)**TOTHRD
      E=E6A
      IF(A.GT.QO) E=1.-QO/A
      P=A*(1.-E*E)
      XNODES= XNODEO+XNODOT*TSINCE
      OMGAS= OMEGAO+OMGDT*TSINCE
      XLS=MOD(XLO+(XNO+OMGDT+XNODOT+(XNDT2O+XNDD6O*TSINCE)*
     1 TSINCE)*TSINCE , TWOPI )
C
C      LONG PERIOD PERIODICS
C
      AXNSL=E*COS(OMGAS)
      AYNSL=E*SIN(OMGAS)-C6/P
      XL=MOD(XLS-C5/P*AXNSL , TWOPI )
C
C      SOLVE KEPLERS EQUATION
C
      U=MOD(XL-XNODES , TWOPI )
      ITEM3=0
      EO1=U
      TEM5=1.
   20 SINEO1=SIN(EO1)
      COSEO1=COS(EO1)
      IF(ABS(TEM5).LT.E6A) GO TO 30
      IF(ITEM3.GE.10) GO TO 30
      ITEM3=ITEM3+1
      TEM5=1.-COSEO1*AXNSL-SINEO1*AYNSL
      TEM5=(U-AYNSL*COSEO1+AXNSL*SINEO1-EO1)/TEM5
      TEM2=ABS(TEM5)
      IF(TEM2.GT.1.) TEM5=TEM2/TEM5
      EO1=EO1+TEM5
      GO TO 20
C
C      SHORT PERIOD PRELIMINARY QUANTITIES
C
   30 ECOSE=AXNSL*COSEO1+AYNSL*SINEO1
      ESINE=AXNSL*SINEO1-AYNSL*COSEO1
      EL2=AXNSL*AXNSL+AYNSL*AYNSL
      PL=A*(1.-EL2)
      PL2=PL*PL
      R=A*(1.-ECOSE)
      RDOT=XKE*SQRT(A)/R*ESINE
      RVDOT=XKE*SQRT(PL)/R
      TEMP=ESINE/(1.+SQRT(1.-EL2))
      SINU=A/R*(SINEO1-AYNSL-AXNSL*TEMP)
      COSU=A/R*(COSEO1-AXNSL+AYNSL*TEMP)
      SU=ATAN2(SINU,COSU)
C
C      UPDATE FOR SHORT PERIODICS
C
      SIN2U=(COSU+COSU)*SINU
      COS2U=1.-2.*SINU*SINU
      RK=R+D1O/PL*COS2U
      UK=SU-D2O/PL2*SIN2U
      XNODEK=XNODES+D3O*SIN2U/PL2
      XINCK =XINCL+D4O/PL2*COS2U
C
C      ORIENTATION VECTORS
C
      SINUK=SIN(UK)
      COSUK=COS(UK)
      SINNOK=SIN(XNODEK)
      COSNOK=COS(XNODEK)
      SINIK=SIN(XINCK)
      COSIK=COS(XINCK)
      XMX=-SINNOK*COSIK
      XMY=COSNOK*COSIK
      UX=XMX*SINUK+COSNOK*COSUK
      UY=XMY*SINUK+SINNOK*COSUK
      UZ=SINIK*SINUK
      VX=XMX*COSUK-COSNOK*SINUK
      VY=XMY*COSUK-SINNOK*SINUK
      VZ=SINIK*COSUK
C
C      POSITION AND VELOCITY
C
      X=RK*UX
      Y=RK*UY
      Z=RK*UZ
      XDOT=RDOT*UX
      YDOT=RDOT*UY
      ZDOT=RDOT*UZ
      XDOT=RVDOT*VX+XDOT
      YDOT=RVDOT*VY+YDOT
      ZDOT=RVDOT*VZ+ZDOT
C
      RETURN
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SGP4(IFLAG,TSINCE)
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
C  turned off. R. Noomen, March 5, 1999
C     IF  (IFLAG .EQ. 0) GO TO 100
C
C      RECOVER ORIGINAL MEAN MOTION (XNODP) AND SEMIMAJOR AXIS (AODP)
C      FROM INPUT ELEMENTS
C
*      WRITE (0,*) 'SGP4: xmo   =',XMO*180.0D0/PI
*      WRITE (0,*) 'SGP4: xnodeo [deg]=',XNODEO*180.0D0/PI
*      WRITE (0,*) 'SGP4: omegao [deg]=',OMEGAO*180.0D0/PI
*      WRITE (0,*) 'SGP4: eo    =',EO
*      WRITE (0,*) 'SGP4: xincl [deg] =',XINCL*180.0D0/PI
*      WRITE (0,*) 'SGP4: xno  [deg]  =',XNO*180.0D0/PI
*      WRITE (0,*) 'SGP4: xndt20=',XNDT2O
*      WRITE (0,*) 'SGP4: xndd60=',XNDD6O
*      WRITE (0,*) 'SGP4: iexp=',IEXP
*      WRITE (0,*) 'SGP4: bstar=',BSTAR
*      WRITE (0,*) 'SGP4: ibexp=',IBEXP
*      WRITE (0,*) 'SGP4:',X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      A1=(XKE/XNO)**TOTHRD
      COSIO=COS(XINCL)
      THETA2=COSIO*COSIO
      X3THM1=3.*THETA2-1.
      EOSQ=EO*EO
      BETAO2=1.-EOSQ
      BETAO=SQRT(BETAO2)
      DEL1=1.5*CK2*X3THM1/(A1*A1*BETAO*BETAO2)
      AO=A1*(1.-DEL1*(.5*TOTHRD+DEL1*(1.+134./81.*DEL1)))
      DELO=1.5*CK2*X3THM1/(AO*AO*BETAO*BETAO2)
      XNODP=XNO/(1.+DELO)
      AODP=AO/(1.-DELO)
C
C      INITIALIZATION
C
C      FOR PERIGEE LESS THAN 220 KILOMETERS, THE ISIMP FLAG IS SET AND
C      THE EQUATIONS ARE TRUNCATED TO LINEAR VARIATION IN SQRT A AND
C      QUADRATIC VARIATION IN MEAN ANOMALY.  ALSO, THE C3 TERM, THE
C      DELTA OMEGA TERM, AND THE DELTA M TERM ARE DROPPED.
C
      ISIMP=0
      IF((AODP*(1.-EO)/AE) .LT. (220./XKMPER+AE)) ISIMP=1
C
C      FOR PERIGEE BELOW 156 KM, THE VALUES OF
C      S AND QOMS2T ARE ALTERED
C
      S4=S
      QOMS24=QOMS2T
      PERIGE=(AODP*(1.-EO)-AE)*XKMPER
      IF(PERIGE .GE. 156.) GO TO 10
      S4=PERIGE-78.
      IF(PERIGE .GT. 98.) GO TO 9
      S4=20.
    9 QOMS24=((120.-S4)*AE/XKMPER)**4
      S4=S4/XKMPER+AE
   10 PINVSQ=1./(AODP*AODP*BETAO2*BETAO2)
      TSI=1./(AODP-S4)
      ETA=AODP*EO*TSI
      ETASQ=ETA*ETA
      EETA=EO*ETA
      PSISQ=ABS(1.-ETASQ)
      COEF=QOMS24*TSI**4
      COEF1=COEF/PSISQ**3.5
      C2=COEF1*XNODP*(AODP*(1.+1.5*ETASQ+EETA*(4.+ETASQ))+.75*
     1         CK2*TSI/PSISQ*X3THM1*(8.+3.*ETASQ*(8.+ETASQ)))
      C1=BSTAR*C2
      SINIO=SIN(XINCL)
      A3OVK2=-XJ3/CK2*AE**3
      C3=COEF*TSI*A3OVK2*XNODP*AE*SINIO/EO
      X1MTH2=1.-THETA2
      C4=2.*XNODP*COEF1*AODP*BETAO2*(ETA*
     1         (2.+.5*ETASQ)+EO*(.5+2.*ETASQ)-2.*CK2*TSI/
     2         (AODP*PSISQ)*(-3.*X3THM1*(1.-2.*EETA+ETASQ*
     3         (1.5-.5*EETA))+.75*X1MTH2*(2.*ETASQ-EETA*
     4         (1.+ETASQ))*COS(2.*OMEGAO)))
      C5=2.*COEF1*AODP*BETAO2*(1.+2.75*(ETASQ+EETA)+EETA*ETASQ)
      THETA4=THETA2*THETA2
      TEMP1=3.*CK2*PINVSQ*XNODP
      TEMP2=TEMP1*CK2*PINVSQ
      TEMP3=1.25*CK4*PINVSQ*PINVSQ*XNODP
      XMDOT=XNODP+.5*TEMP1*BETAO*X3THM1+.0625*TEMP2*BETAO*
     1         (13.-78.*THETA2+137.*THETA4)
      X1M5TH=1.-5.*THETA2
      OMGDOT=-.5*TEMP1*X1M5TH+.0625*TEMP2*(7.-114.*THETA2+
     1         395.*THETA4)+TEMP3*(3.-36.*THETA2+49.*THETA4)
      XHDOT1=-TEMP1*COSIO
      XNODOT=XHDOT1+(.5*TEMP2*(4.-19.*THETA2)+2.*TEMP3*(3.-
     1         7.*THETA2))*COSIO
      OMGCOF=BSTAR*C3*COS(OMEGAO)
      XMCOF=-TOTHRD*COEF*BSTAR*AE/EETA
      XNODCF=3.5*BETAO2*XHDOT1*C1
      T2COF=1.5*C1
      XLCOF=.125*A3OVK2*SINIO*(3.+5.*COSIO)/(1.+COSIO)
      AYCOF=.25*A3OVK2*SINIO
      DELMO=(1.+ETA*COS(XMO))**3
      SINMO=SIN(XMO)
      X7THM1=7.*THETA2-1.
      IF(ISIMP .EQ. 1) GO TO 90
      C1SQ=C1*C1
      D2=4.*AODP*TSI*C1SQ
      TEMP=D2*TSI*C1/3.
      D3=(17.*AODP+S4)*TEMP
      D4=.5*TEMP*AODP*TSI*(221.*AODP+31.*S4)*C1
      T3COF=D2+2.*C1SQ
      T4COF=.25*(3.*D3+C1*(12.*D2+10.*C1SQ))
      T5COF=.2*(3.*D4+12.*C1*D3+6.*D2*D2+15.*C1SQ*(
     1         2.*D2+C1SQ))
   90 IFLAG=0
C
C      UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
C
  100 XMDF=XMO+XMDOT*TSINCE
      OMGADF=OMEGAO+OMGDOT*TSINCE
      XNODDF=XNODEO+XNODOT*TSINCE
      OMEGA=OMGADF
      XMP=XMDF
      TSQ=TSINCE*TSINCE
      XNODE=XNODDF+XNODCF*TSQ
      TEMPA=1.-C1*TSINCE
      TEMPE=BSTAR*C4*TSINCE
      TEMPL=T2COF*TSQ
      IF(ISIMP .EQ. 1) GO TO 110
      DELOMG=OMGCOF*TSINCE
      DELM=XMCOF*((1.+ETA*COS(XMDF))**3-DELMO)
      TEMP=DELOMG+DELM
      XMP=XMDF+TEMP
      OMEGA=OMGADF-TEMP
      TCUBE=TSQ*TSINCE
      TFOUR=TSINCE*TCUBE
      TEMPA=TEMPA-D2*TSQ-D3*TCUBE-D4*TFOUR
      TEMPE=TEMPE+BSTAR*C5*(SIN(XMP)-SINMO)
      TEMPL=TEMPL+T3COF*TCUBE+
     1         TFOUR*(T4COF+TSINCE*T5COF)
  110 A=AODP*TEMPA**2
      E=EO-TEMPE
      XL=XMP+OMEGA+XNODE+XNODP*TEMPL
      BETA=SQRT(1.-E*E)
      XN=XKE/A**1.5
C
C      LONG PERIOD PERIODICS
C
      AXN=E*COS(OMEGA)
      TEMP=1./(A*BETA*BETA)
      XLL=TEMP*XLCOF*AXN
      AYNL=TEMP*AYCOF
      XLT=XL+XLL
      AYN=E*SIN(OMEGA)+AYNL
C
C      SOLVE KEPLERS EQUATION
C
      CAPU=MOD(XLT-XNODE , TWOPI )
      TEMP2=CAPU
      DO 130 I=1,10
      SINEPW=SIN(TEMP2)
      COSEPW=COS(TEMP2)
      TEMP3=AXN*SINEPW
      TEMP4=AYN*COSEPW
      TEMP5=AXN*COSEPW
      TEMP6=AYN*SINEPW
      EPW=(CAPU-TEMP4+TEMP3-TEMP2)/(1.-TEMP5-TEMP6)+TEMP2
      IF(ABS(EPW-TEMP2) .LE. E6A) GO TO 140
  130 TEMP2=EPW
C
C      SHORT PERIOD PRELIMINARY QUANTITIES
C
  140 ECOSE=TEMP5+TEMP6
      ESINE=TEMP3-TEMP4
*      WRITE (0,*) 'SGP4: ecc.anom+omega[deg]=',EPW*180.0D0/PI
*      WRITE (0,*) 'SGP4: omega[deg]=',OMEGA*180.0D0/PI
*      WRITE (0,*) 'SGP4: ecc.anom[deg]=',(EPW-OMEGA)*180.0D0/PI
      ELSQ=AXN*AXN+AYN*AYN
      TEMP=1.-ELSQ
      PL=A*TEMP
      R=A*(1.-ECOSE)
      TEMP1=1./R
      RDOT=XKE*SQRT(A)*ESINE*TEMP1
      RFDOT=XKE*SQRT(PL)*TEMP1
      TEMP2=A*TEMP1
      BETAL=SQRT(TEMP)
      TEMP3=1./(1.+BETAL)
      COSU=TEMP2*(COSEPW-AXN+AYN*ESINE*TEMP3)
      SINU=TEMP2*(SINEPW-AYN-AXN*ESINE*TEMP3)
      U=ATAN2(SINU,COSU)
*      WRITE (0,*) 'SGP4: u[deg]=',U*180.0D0/PI
      SIN2U=2.*SINU*COSU
      COS2U=2.*COSU*COSU-1.
      TEMP=1./PL
      TEMP1=CK2*TEMP
      TEMP2=TEMP1*TEMP
C
C      UPDATE FOR SHORT PERIODICS
C
      RK=R*(1.-1.5*TEMP2*BETAL*X3THM1)+.5*TEMP1*X1MTH2*COS2U
      UK=U-.25*TEMP2*X7THM1*SIN2U
      XNODEK=XNODE+1.5*TEMP2*COSIO*SIN2U
      XINCK=XINCL+1.5*TEMP2*COSIO*SINIO*COS2U
      RDOTK=RDOT-XN*TEMP1*X1MTH2*SIN2U
      RFDOTK=RFDOT+XN*TEMP1*(X1MTH2*COS2U+1.5*X3THM1)
C
C      ORIENTATION VECTORS
C
*      WRITE (0,*) 'SGP4: uk[deg]=',UK*180.0D0/PI
*      WRITE (0,*) 'SGP4: i[deg]=',XINCK*180.0D0/PI
*      WRITE (0,*) 'SGP4: OMEGA[deg]=',XNODEK*180.0D0/PI
      SINUK=SIN(UK)
      COSUK=COS(UK)
      SINIK=SIN(XINCK)
      COSIK=COS(XINCK)
      SINNOK=SIN(XNODEK)
      COSNOK=COS(XNODEK)
      XMX=-SINNOK*COSIK
      XMY=COSNOK*COSIK
      UX=XMX*SINUK+COSNOK*COSUK
      UY=XMY*SINUK+SINNOK*COSUK
      UZ=SINIK*SINUK
      VX=XMX*COSUK-COSNOK*SINUK
      VY=XMY*COSUK-SINNOK*SINUK
      VZ=SINIK*COSUK
C
C      POSITION AND VELOCITY
C
      X=RK*UX
      Y=RK*UY
      Z=RK*UZ
      XDOT=RDOTK*UX+RFDOTK*VX
      YDOT=RDOTK*UY+RFDOTK*VY
      ZDOT=RDOTK*UZ+RFDOTK*VZ
C
      RETURN
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SGP8(IFLAG,TSINCE)
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
C  turned off. R. Noomen, March 5, 1999
C     IF  (IFLAG .EQ. 0) GO TO 100
C
C      RECOVER ORIGINAL MEAN MOTION (XNODP) AND SEMIMAJOR AXIS (AODP)
C      FROM INPUT ELEMENTS --------- CALCULATE BALLISTIC COEFFICIENT
C      (B TERM) FROM INPUT B* DRAG TERM
C
      A1=(XKE/XNO)**TOTHRD
      COSI=COS(XINCL)
      THETA2=COSI*COSI
      TTHMUN=3.*THETA2-1.
      EOSQ=EO*EO
      BETAO2=1.-EOSQ
      BETAO=SQRT(BETAO2)
      DEL1=1.5*CK2*TTHMUN/(A1*A1*BETAO*BETAO2)
      AO=A1*(1.-DEL1*(.5*TOTHRD+DEL1*(1.+134./81.*DEL1)))
      DELO=1.5*CK2*TTHMUN/(AO*AO*BETAO*BETAO2)
      AODP=AO/(1.-DELO)
      XNODP=XNO/(1.+DELO)
      B=2.*BSTAR/RHO
C
C      INITIALIZATION
C
      ISIMP=0
      PO=AODP*BETAO2
      POM2=1./(PO*PO)
      SINI=SIN(XINCL)
      SING=SIN(OMEGAO)
      COSG=COS(OMEGAO)
      TEMP=.5*XINCL
      SINIO2=SIN(TEMP)
      COSIO2=COS(TEMP)
      THETA4=THETA2**2
      UNM5TH=1.-5.*THETA2
      UNMTH2=1.-THETA2
      A3COF=-XJ3/CK2*AE**3
      PARDT1=3.*CK2*POM2*XNODP
      PARDT2=PARDT1*CK2*POM2
      PARDT4=1.25*CK4*POM2*POM2*XNODP
      XMDT1=.5*PARDT1*BETAO*TTHMUN
      XGDT1=-.5*PARDT1*UNM5TH
      XHDT1=-PARDT1*COSI
      XLLDOT=XNODP+XMDT1+
     2           .0625*PARDT2*BETAO*(13.-78.*THETA2+137.*THETA4)
      OMGDT=XGDT1+
     1      .0625*PARDT2*(7.-114.*THETA2+395.*THETA4)+PARDT4*(3.-36.*
     2         THETA2+49.*THETA4)
      XNODOT=XHDT1+
     1       (.5*PARDT2*(4.-19.*THETA2)+2.*PARDT4*(3.-7.*THETA2))*COSI
      TSI=1./(PO-S)
      ETA=EO*S*TSI
      ETA2=ETA**2
      PSIM2=ABS(1./(1.-ETA2))
      ALPHA2=1.+EOSQ
      EETA=EO*ETA
      COS2G=2.*COSG**2-1.
      D5=TSI*PSIM2
      D1=D5/PO
      D2=12.+ETA2*(36.+4.5*ETA2)
      D3=ETA2*(15.+2.5*ETA2)
      D4=ETA*(5.+3.75*ETA2)
      B1=CK2*TTHMUN
      B2=-CK2*UNMTH2
      B3=A3COF*SINI
      C0=.5*B*RHO*QOMS2T*XNODP*AODP*TSI**4*PSIM2**3.5/SQRT(ALPHA2)
      C1=1.5*XNODP*ALPHA2**2*C0
      C4=D1*D3*B2
      C5=D5*D4*B3
      XNDT=C1*(
     1  (2.+ETA2*(3.+34.*EOSQ)+5.*EETA*(4.+ETA2)+8.5*EOSQ)+
     1   D1*D2*B1+   C4*COS2G+C5*SING)
      XNDTN=XNDT/XNODP
C
C      IF DRAG IS VERY SMALL, THE ISIMP FLAG IS SET AND THE
C      EQUATIONS ARE TRUNCATED TO LINEAR VARIATION IN MEAN
C      MOTION AND QUADRATIC VARIATION IN MEAN ANOMALY
C
      IF(ABS(XNDTN*XMNPDA) .LT. 2.16E-3) GO TO 50
      D6=ETA*(30.+22.5*ETA2)
      D7=ETA*(5.+12.5*ETA2)
      D8=1.+ETA2*(6.75+ETA2)
      C8=D1*D7*B2
      C9=D5*D8*B3
      EDOT=-C0*(
     1   ETA*(4.+ETA2+EOSQ*(15.5+7.*ETA2))+EO*(5.+15.*ETA2)+
     1   D1*D6*B1 +
     1   C8*COS2G+C9*SING)
      D20=.5*TOTHRD*XNDTN
      ALDTAL=EO*EDOT/ALPHA2
      TSDTTS=2.*AODP*TSI*(D20*BETAO2+EO*EDOT)
      ETDT=(EDOT+EO*TSDTTS)*TSI*S
      PSDTPS=-ETA*ETDT*PSIM2
      SIN2G=2.*SING*COSG
      C0DTC0=D20+4.*TSDTTS-ALDTAL-7.*PSDTPS
      C1DTC1=XNDTN+4.*ALDTAL+C0DTC0
      D9=ETA*(6.+68.*EOSQ)+EO*(20.+15.*ETA2)
      D10=5.*ETA*(4.+ETA2)+EO*(17.+68.*ETA2)
      D11=ETA*(72.+18.*ETA2)
      D12=ETA*(30.+10.*ETA2)
      D13=5.+11.25*ETA2
      D14=TSDTTS-2.*PSDTPS
      D15=2.*(D20+EO*EDOT/BETAO2)
      D1DT=D1*(D14+D15)
      D2DT=ETDT*D11
      D3DT=ETDT*D12
      D4DT=ETDT*D13
      D5DT=D5*D14
      C4DT=B2*(D1DT*D3+D1*D3DT)
      C5DT=B3*(D5DT*D4+D5*D4DT)
      D16=
     1     D9*ETDT+D10*EDOT +
     1     B1*(D1DT*D2+D1*D2DT) +
     1     C4DT*COS2G+C5DT*SING+XGDT1*(C5*COSG-2.*C4*SIN2G)
      XNDDT=C1DTC1*XNDT+C1*D16
      EDDOT=C0DTC0*EDOT-C0*(
     1     (4.+3.*ETA2+30.*EETA+EOSQ*(15.5+21.*ETA2))*ETDT+(5.+15.*ETA2
     '         +EETA*(31.+14.*ETA2))*EDOT +
     1     B1*(D1DT*D6+D1*ETDT*(30.+67.5*ETA2))  +
     1     B2*(D1DT*D7+D1*ETDT*(5.+37.5*ETA2))*COS2G+
     1     B3*(D5DT*D8+D5*ETDT*ETA*(13.5+4.*ETA2))*SING+XGDT1*(C9*
     '         COSG-2.*C8*SIN2G))
      D25=EDOT**2
      D17=XNDDT/XNODP-XNDTN**2
      TSDDTS=2.*TSDTTS*(TSDTTS-D20)+AODP*TSI*(TOTHRD*BETAO2*D17-4.*D20*
     '         EO*EDOT+2.*(D25+EO*EDDOT))
      ETDDT =(EDDOT+2.*EDOT*TSDTTS)*TSI*S+TSDDTS*ETA
      D18=TSDDTS-TSDTTS**2
      D19=-PSDTPS**2/ETA2-ETA*ETDDT*PSIM2-PSDTPS**2
      D23=ETDT*ETDT
      D1DDT=D1DT*(D14+D15)+D1*(D18-2.*D19+TOTHRD*D17+2.*(ALPHA2*D25
     '         /BETAO2+EO*EDDOT)/BETAO2)
      XNTRDT=XNDT*(2.*TOTHRD*D17+3.*
     1  (D25+EO*EDDOT)/ALPHA2-6.*ALDTAL**2 +
     1  4.*D18-7.*D19 )   +
     1  C1DTC1*XNDDT+C1*(C1DTC1*D16+
     1  D9*ETDDT+D10*EDDOT+D23*(6.+30.*EETA+68.*EOSQ)+
     1  ETDT*EDOT*(40.+30.*
     '  ETA2+272.*EETA)+D25*(17.+68.*ETA2) +
     1    B1*(D1DDT*D2+2.*D1DT*D2DT+D1*(ETDDT*D11+D23*(72.+54.*ETA2))) +
     1    B2*(D1DDT*D3+2.*D1DT*D3DT+D1*(ETDDT*D12+D23*(30.+30.*ETA2))) *
     1    COS2G+
     1      B3*((D5DT*D14+D5*(D18-2.*D19)) *
     1 D4+2.*D4DT*D5DT+D5*(ETDDT*D13+22.5*ETA*D23)) *SING+XGDT1*
     1         ((7.*D20+4.*EO*EDOT/BETAO2)*
     '         (C5*COSG-2.*C4*SIN2G)
     '         +((2.*C5DT*COSG-4.*C4DT*SIN2G)-XGDT1*(C5*SING+4.*
     '         C4*COS2G))))
      TMNDDT=XNDDT*1.E9
      TEMP=TMNDDT**2-XNDT*1.E18*XNTRDT
      PP=(TEMP+TMNDDT**2)/TEMP
      GAMMA=-XNTRDT/(XNDDT*(PP-2.))
      XND=XNDT/(PP*GAMMA)
      QQ=1.-EDDOT/(EDOT*GAMMA)
      ED=EDOT/(QQ*GAMMA)
      OVGPP=1./(GAMMA*(PP+1.))
      GO TO 70
   50 ISIMP=1
      EDOT=-TOTHRD*XNDTN*(1.-EO)
   70 IFLAG=0
C
C      UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
C
  100 XMAM=MOD(XMO+XLLDOT*TSINCE , TWOPI )
      OMGASM=OMEGAO+OMGDT*TSINCE
      XNODES=XNODEO+XNODOT*TSINCE
      IF(ISIMP .EQ. 1) GO TO 105
      TEMP=1.-GAMMA*TSINCE
      TEMP1=TEMP**PP
      XN=XNODP+XND*(1.-TEMP1)
      EM=EO+ED*(1.-TEMP**QQ)
      Z1=XND*(TSINCE+OVGPP*(TEMP*TEMP1-1.))
      GO TO 108
  105 XN=XNODP+XNDT*TSINCE
      EM=EO+EDOT*TSINCE
      Z1=.5*XNDT*TSINCE*TSINCE
  108 Z7=3.5*TOTHRD*Z1/XNODP
      XMAM=MOD(XMAM+Z1+Z7*XMDT1 , TWOPI )
      OMGASM=OMGASM+Z7*XGDT1
      XNODES=XNODES+Z7*XHDT1
C
C      SOLVE KEPLERS EQUATION
C
      ZC2=XMAM+EM*SIN(XMAM)*(1.+EM*COS(XMAM))
      DO 130 I=1,10
      SINE=SIN(ZC2)
      COSE=COS(ZC2)
      ZC5=1./(1.-EM*COSE)
      CAPE=(XMAM+EM*SINE-ZC2)*
     1   ZC5+ZC2
      IF(ABS(CAPE-ZC2) .LE. E6A) GO TO 140
  130 ZC2=CAPE
C
C      SHORT PERIOD PRELIMINARY QUANTITIES
C
  140 AM=(XKE/XN)**TOTHRD
      BETA2M=1.-EM*EM
      SINOS=SIN(OMGASM)
      COSOS=COS(OMGASM)
      AXNM=EM*COSOS
      AYNM=EM*SINOS
      PM=AM*BETA2M
      G1=1./PM
      G2=.5*CK2*G1
      G3=G2*G1
      BETA=SQRT(BETA2M)
      G4=.25*A3COF*SINI
      G5=.25*A3COF*G1
      SNF=BETA*SINE*ZC5
      CSF=(COSE-EM)*ZC5
      FM=ATAN2(SNF,CSF)
      SNFG=SNF*COSOS+CSF*SINOS
      CSFG=CSF*COSOS-SNF*SINOS
      SN2F2G=2.*SNFG*CSFG
      CS2F2G=2.*CSFG**2-1.
      ECOSF=EM*CSF
      G10=FM-XMAM+EM*SNF
      RM=PM/(1.+ECOSF)
      AOVR=AM/RM
      G13=XN*AOVR
      G14=-G13*AOVR
      DR=G2*(UNMTH2*CS2F2G-3.*TTHMUN)-G4*SNFG
      DIWC=3.*G3*SINI*CS2F2G-G5*AYNM
      DI=DIWC*COSI
C
C      UPDATE FOR SHORT PERIOD PERIODICS
C
      SNI2DU=SINIO2*(
     1   G3*(.5*(1.-7.*THETA2)*SN2F2G-3.*UNM5TH*G10)-G5*SINI*CSFG*(2.+
     2         ECOSF))-.5*G5*THETA2*AXNM/COSIO2
      XLAMB=FM+OMGASM+XNODES+G3*(.5*(1.+6.*COSI-7.*THETA2)*SN2F2G-3.*
     1      (UNM5TH+2.*COSI)*G10)+G5*SINI*(COSI*AXNM/(1.+COSI)-(2.
     2      +ECOSF)*CSFG)
      Y4=SINIO2*SNFG+CSFG*SNI2DU+.5*SNFG*COSIO2*DI
      Y5=SINIO2*CSFG-SNFG*SNI2DU+.5*CSFG*COSIO2*DI
      R=RM+DR
      RDOT=XN*AM*EM*SNF/BETA+G14*(2.*G2*UNMTH2*SN2F2G+G4*CSFG)
      RVDOT=XN*AM**2*BETA/RM+
     1      G14*DR+AM*G13*SINI*DIWC
C
C      ORIENTATION VECTORS
C
      SNLAMB=SIN(XLAMB)
      CSLAMB=COS(XLAMB)
      TEMP=2.*(Y5*SNLAMB-Y4*CSLAMB)
      UX=Y4*TEMP+CSLAMB
      VX=Y5*TEMP-SNLAMB
      TEMP=2.*(Y5*CSLAMB+Y4*SNLAMB)
      UY=-Y4*TEMP+SNLAMB
      VY=-Y5*TEMP+CSLAMB
      TEMP=2.*SQRT(1.-Y4*Y4-Y5*Y5)
      UZ=Y4*TEMP
      VZ=Y5*TEMP
C
C      POSITION AND VELOCITY
C
      X=R*UX
      Y=R*UY
      Z=R*UZ
      XDOT=RDOT*UX+RVDOT*VX
      YDOT=RDOT*UY+RVDOT*VY
      ZDOT=RDOT*UZ+RVDOT*VZ
C
      RETURN
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      FUNCTION THETAG(EP)
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON /E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
      YR=(EP+2.D-7)*1.D-3
      JY=YR
      YR=JY
      D=EP-YR*1.D3
      IF(JY.LT.10) JY=JY+80
      N=(JY-69)/4
      IF(JY.LT.70) N=(JY-72)/4
      DS50=7305.D0 + 365.D0*(JY-70) +N + D
      THETA=1.72944494D0 + 6.3003880987D0*DS50
      TEMP=THETA/TWOPI
      I=TEMP
      TEMP=I
      THETAG=THETA-TEMP*TWOPI
      IF(THETAG.LT.0.D0) THETAG=THETAG+TWOPI
C
      RETURN
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DPINIT(EQSQ,SINIQ,COSIQ,RTEQSQ,AO,COSQ2,SINOMO,COSOMO,
     1         BSQ,XLLDOT,OMGDT,XNODOT,XNODP)
C
C     ENTRANCE FOR DEEP SPACE INITIALIZATION
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
C     ENTRANCE FOR DEEP SPACE INITIALIZATION
C
      THGR=THETAG(EPOCH)
      EQ = EO
      XNQ = XNODP
      AQNV = 1./AO
      XQNCL = XINCL
      XMAO=XMO
      XPIDOT=OMGDT+XNODOT
      SINQ = SIN(XNODEO)
      COSQ = COS(XNODEO)
      OMEGAQ = OMEGAO
C
C     INITIALIZE LUNAR SOLAR TERMS
C
    5 DAY=DS50+18261.5D0
C     IF (DAY.EQ.PREEP)    GO TO 10
      DIFF = DAY - PREEP
      IF ( ABS ( DIFF ) .LT. EPS )    GO TO 10
      PREEP = DAY
      XNODCE=4.5236020-9.2422029E-4*DAY
      STEM= SIN (XNODCE)
      CTEM= COS (XNODCE)
      ZCOSIL=.91375164-.03568096*CTEM
      ZSINIL=SQRT (1.-ZCOSIL*ZCOSIL)
      ZSINHL= .089683511*STEM/ZSINIL
      ZCOSHL=SQRT (1.-ZSINHL*ZSINHL)
      C=4.7199672+.22997150*DAY
      GAM=5.8351514+.0019443680*DAY
      ZMOL = MOD(C-GAM , TWOPI )
      ZX= .39785416*STEM/ZSINIL
      ZY= ZCOSHL*CTEM+0.91744867*ZSINHL*STEM
      ZX=ATAN2(ZX,ZY)
      ZX=GAM+ZX-XNODCE
      ZCOSGL=COS (ZX)
      ZSINGL=SIN (ZX)
      ZMOS=6.2565837D0+.017201977D0*DAY
      ZMOS=MOD(ZMOS , TWOPI )
C
C     DO SOLAR TERMS
C
   10 LS = 0
      SAVTSN=1.D20
      ZCOSG=ZCOSGS
      ZSING=ZSINGS
      ZCOSI=ZCOSIS
      ZSINI=ZSINIS
      ZCOSH=COSQ
      ZSINH=SINQ
      CC=C1SS
      ZN=ZNS
      ZE=ZES
      ZMO=ZMOS
      XNOI=1./XNQ
C modification by R. Noomen, March 4, 1999
C     ASSIGN 30 TO LS
      LS = 30
C  eom
   20 A1=ZCOSG*ZCOSH+ZSING*ZCOSI*ZSINH
      A3=-ZSING*ZCOSH+ZCOSG*ZCOSI*ZSINH
      A7=-ZCOSG*ZSINH+ZSING*ZCOSI*ZCOSH
      A8=ZSING*ZSINI
      A9=ZSING*ZSINH+ZCOSG*ZCOSI*ZCOSH
      A10=ZCOSG*ZSINI
      A2= COSIQ*A7+ SINIQ*A8
      A4= COSIQ*A9+ SINIQ*A10
      A5=- SINIQ*A7+ COSIQ*A8
      A6=- SINIQ*A9+ COSIQ*A10
C
      X1=A1*COSOMO+A2*SINOMO
      X2=A3*COSOMO+A4*SINOMO
      X3=-A1*SINOMO+A2*COSOMO
      X4=-A3*SINOMO+A4*COSOMO
      X5=A5*SINOMO
      X6=A6*SINOMO
      X7=A5*COSOMO
      X8=A6*COSOMO
C
      Z31=12.*X1*X1-3.*X3*X3
      Z32=24.*X1*X2-6.*X3*X4
      Z33=12.*X2*X2-3.*X4*X4
      Z1=3.*(A1*A1+A2*A2)+Z31*EQSQ
      Z2=6.*(A1*A3+A2*A4)+Z32*EQSQ
      Z3=3.*(A3*A3+A4*A4)+Z33*EQSQ
      Z11=-6.*A1*A5+EQSQ *(-24.*X1*X7-6.*X3*X5)
      Z12=-6.*(A1*A6+A3*A5)+EQSQ *(-24.*(X2*X7+X1*X8)-6.*(X3*X6+X4*X5))
      Z13=-6.*A3*A6+EQSQ *(-24.*X2*X8-6.*X4*X6)
      Z21=6.*A2*A5+EQSQ *(24.*X1*X5-6.*X3*X7)
      Z22=6.*(A4*A5+A2*A6)+EQSQ *(24.*(X2*X5+X1*X6)-6.*(X4*X7+X3*X8))
      Z23=6.*A4*A6+EQSQ *(24.*X2*X6-6.*X4*X8)
      Z1=Z1+Z1+BSQ*Z31
      Z2=Z2+Z2+BSQ*Z32
      Z3=Z3+Z3+BSQ*Z33
      S3=CC*XNOI
      S2=-.5*S3/RTEQSQ
      S4=S3*RTEQSQ
      S1=-15.*EQ*S4
      S5=X1*X3+X2*X4
      S6=X2*X3+X1*X4
      S7=X2*X4-X1*X3
      SE=S1*ZN*S5
      SI=S2*ZN*(Z11+Z13)
      SL=-ZN*S3*(Z1+Z3-14.-6.*EQSQ)
      SGH=S4*ZN*(Z31+Z33-6.)
      SH=-ZN*S2*(Z21+Z23)
      IF(XQNCL.LT.5.2359877E-2) SH=0.0
      EE2=2.*S1*S6
      E3=2.*S1*S7
      XI2=2.*S2*Z12
      XI3=2.*S2*(Z13-Z11)
      XL2=-2.*S3*Z2
      XL3=-2.*S3*(Z3-Z1)
      XL4=-2.*S3*(-21.-9.*EQSQ)*ZE
      XGH2=2.*S4*Z32
      XGH3=2.*S4*(Z33-Z31)
      XGH4=-18.*S4*ZE
      XH2=-2.*S2*Z22
      XH3=-2.*S2*(Z23-Z21)
C modification by R. Noomen, March 4, 1999
C     GO TO LS
      IF ( LS .EQ. 30 ) GOTO 30
      IF ( LS .EQ. 40 ) GOTO 40
C  eom
C
C     DO LUNAR TERMS
C
   30 SSE = SE
      SSI=SI
      SSL=SL
      SSH=SH/SINIQ
      SSG=SGH-COSIQ*SSH
      SE2=EE2
      SI2=XI2
      SL2=XL2
      SGH2=XGH2
      SH2=XH2
      SE3=E3
      SI3=XI3
      SL3=XL3
      SGH3=XGH3
      SH3=XH3
      SL4=XL4
      SGH4=XGH4
      LS=1
      ZCOSG=ZCOSGL
      ZSING=ZSINGL
      ZCOSI=ZCOSIL
      ZSINI=ZSINIL
      ZCOSH=ZCOSHL*COSQ+ZSINHL*SINQ
      ZSINH=SINQ*ZCOSHL-COSQ*ZSINHL
      ZN=ZNL
      CC=C1L
      ZE=ZEL
      ZMO=ZMOL
C modification by R. Noomen, March 4, 1999
C     ASSIGN 40 TO LS
      LS = 40
C  eom
      GO TO 20
   40 SSE = SSE+SE
      SSI=SSI+SI
      SSL=SSL+SL
      SSG=SSG+SGH-COSIQ/SINIQ*SH
      SSH=SSH+SH/SINIQ
C
C     GEOPOTENTIAL RESONANCE INITIALIZATION FOR 12 HOUR ORBITS
C
      IRESFL=0
      ISYNFL=0
      IF(XNQ.LT.(.0052359877).AND.XNQ.GT.(.0034906585)) GO TO 70
      IF (XNQ.LT.(8.26E-3) .OR. XNQ.GT.(9.24E-3))    RETURN
      IF (EQ.LT.0.5)    RETURN
      IRESFL =1
      EOC=EQ*EQSQ
      G201=-.306-(EQ-.64)*.440
      IF(EQ.GT.(.65)) GO TO 45
      G211=3.616-13.247*EQ+16.290*EQSQ
      G310=-19.302+117.390*EQ-228.419*EQSQ+156.591*EOC
      G322=-18.9068+109.7927*EQ-214.6334*EQSQ+146.5816*EOC
      G410=-41.122+242.694*EQ-471.094*EQSQ+313.953*EOC
      G422=-146.407+841.880*EQ-1629.014*EQSQ+1083.435*EOC
      G520=-532.114+3017.977*EQ-5740*EQSQ+3708.276*EOC
      GO TO 55
   45 G211=-72.099+331.819*EQ-508.738*EQSQ+266.724*EOC
      G310=-346.844+1582.851*EQ-2415.925*EQSQ+1246.113*EOC
      G322=-342.585+1554.908*EQ-2366.899*EQSQ+1215.972*EOC
      G410=-1052.797+4758.686*EQ-7193.992*EQSQ+3651.957*EOC
      G422=-3581.69+16178.11*EQ-24462.77*EQSQ+12422.52*EOC
      IF(EQ.GT.(.715)) GO TO 50
      G520=1464.74-4664.75*EQ+3763.64*EQSQ
      GO TO 55
   50 G520=-5149.66+29936.92*EQ-54087.36*EQSQ+31324.56*EOC
   55 IF(EQ.GE.(.7)) GO TO 60
      G533=-919.2277+4988.61*EQ-9064.77*EQSQ+5542.21*EOC
      G521 = -822.71072+4568.6173*EQ-8491.4146*EQSQ+5337.524*EOC
      G532 = -853.666+4690.25*EQ-8624.77*EQSQ+5341.4*EOC
      GO TO 65
   60 G533=-37995.78+161616.52*EQ-229838.2*EQSQ+109377.94*EOC
      G521 = -51752.104+218913.95*EQ-309468.16*EQSQ+146349.42*EOC
      G532 = -40023.88+170470.89*EQ-242699.48*EQSQ+115605.82*EOC
   65 SINI2=SINIQ*SINIQ
      F220=.75*(1.+2.*COSIQ+COSQ2)
      F221=1.5*SINI2
      F321=1.875*SINIQ*(1.-2.*COSIQ-3.*COSQ2)
      F322=-1.875*SINIQ*(1.+2.*COSIQ-3.*COSQ2)
      F441=35.*SINI2*F220
      F442=39.3750*SINI2*SINI2
      F522=9.84375*SINIQ*(SINI2*(1.-2.*COSIQ-5.*COSQ2)
     1     +.33333333*(-2.+4.*COSIQ+6.*COSQ2))
      F523 = SINIQ*(4.92187512*SINI2*(-2.-4.*COSIQ+10.*COSQ2)
     *      +6.56250012*(1.+2.*COSIQ-3.*COSQ2))
      F542 = 29.53125*SINIQ*(2.-8.*COSIQ+COSQ2*(-12.+8.*COSIQ
     *      +10.*COSQ2))
      F543=29.53125*SINIQ*(-2.-8.*COSIQ+COSQ2*(12.+8.*COSIQ-10.*COSQ2))
      XNO2=XNQ*XNQ
      AINV2=AQNV*AQNV
      TEMP1 = 3.*XNO2*AINV2
      TEMP = TEMP1*ROOT22
      D2201 = TEMP*F220*G201
      D2211 = TEMP*F221*G211
      TEMP1 = TEMP1*AQNV
      TEMP = TEMP1*ROOT32
      D3210 = TEMP*F321*G310
      D3222 = TEMP*F322*G322
      TEMP1 = TEMP1*AQNV
      TEMP = 2.*TEMP1*ROOT44
      D4410 = TEMP*F441*G410
      D4422 = TEMP*F442*G422
      TEMP1 = TEMP1*AQNV
      TEMP = TEMP1*ROOT52
      D5220 = TEMP*F522*G520
      D5232 = TEMP*F523*G532
      TEMP = 2.*TEMP1*ROOT54
      D5421 = TEMP*F542*G521
      D5433 = TEMP*F543*G533
      XLAMO = XMAO+XNODEO+XNODEO-THGR-THGR
      BFACT = XLLDOT+XNODOT+XNODOT-THDT-THDT
      BFACT=BFACT+SSL+SSH+SSH
      GO TO 80
C
C      SYNCHRONOUS RESONANCE TERMS INITIALIZATION
C
   70 IRESFL=1
      ISYNFL=1
      G200=1.0+EQSQ*(-2.5+.8125*EQSQ)
      G310=1.0+2.0*EQSQ
      G300=1.0+EQSQ*(-6.0+6.60937*EQSQ)
      F220=.75*(1.+COSIQ)*(1.+COSIQ)
      F311=.9375*SINIQ*SINIQ*(1.+3.*COSIQ)-.75*(1.+COSIQ)
      F330=1.+COSIQ
      F330=1.875*F330*F330*F330
      DEL1=3.*XNQ*XNQ*AQNV*AQNV
      DEL2=2.*DEL1*F220*G200*Q22
      DEL3=3.*DEL1*F330*G300*Q33*AQNV
      DEL1=DEL1*F311*G310*Q31*AQNV
      FASX2=.13130908
      FASX4=2.8843198
      FASX6=.37448087
      XLAMO=XMAO+XNODEO+OMEGAO-THGR
      BFACT = XLLDOT+XPIDOT-THDT
      BFACT=BFACT+SSL+SSG+SSH
   80 XFACT=BFACT-XNQ
C
C     INITIALIZE INTEGRATOR
C
      XLI=XLAMO
      XNI=XNQ
      ATIME=0.D0
      STEPP=720.D0
      STEPN=-720.D0
      STEP2 = 259200.D0
      RETURN
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DPSEC(XLL,OMGASM,XNODES,EM,XINC,XN,T,OMGDT)
C     SUBROUTINE DPSEC(XLL,OMGASM,XNODES,EM,XINC,XN,T)
C
C     ENTRANCE FOR DEEP SPACE SECULAR EFFECTS
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
C     ENTRANCE FOR DEEP SPACE SECULAR EFFECTS
C
      XLL=XLL+SSL*T
      OMGASM=OMGASM+SSG*T
      XNODES=XNODES+SSH*T
      EM=EO+SSE*T
      XINC=XINCL+SSI*T
      IF(XINC .GE. 0.) GO TO 90
      XINC = -XINC
      XNODES = XNODES + PI
      OMGASM = OMGASM - PI
   90 IF(IRESFL .EQ. 0) RETURN
C 100 IF (ATIME.EQ.0.D0)    GO TO 170
  100 IF ( ABS ( ATIME ) .LT. EPS )    GO TO 170
      IF(T.GE.(0.D0).AND.ATIME.LT.(0.D0)) GO TO 170
      IF(T.LT.(0.D0).AND.ATIME.GE.(0.D0)) GO TO 170
  105 IF( ABS(T).GE. ABS(ATIME)) GO TO 120
      DELT=STEPP
      IF (T.GE.0.D0)    DELT = STEPN
C modification by R. Noomen, March 4, 1999
C 110 ASSIGN 100 TO IRET
  110 IRET = 100
C  eom
      GO TO 160
  120 DELT=STEPN
      IF (T.GT.0.D0)    DELT = STEPP
  125 IF ( ABS(T-ATIME).LT.STEPP)    GO TO 130
C modification by R. Noomen, March 4, 1999
C     ASSIGN 125 TO IRET
      IRET = 125
C  eom
      GO TO 160
  130 FT = T-ATIME
C modification by R. Noomen, March 4, 1999
C     ASSIGN 140 TO IRETN
      IRETN = 140
C  eom
      GO TO 150
  140 XN = XNI+XNDOT*FT+XNDDT*FT*FT*0.5
      XL = XLI+XLDOT*FT+XNDOT*FT*FT*0.5
      TEMP = -XNODES+THGR+T*THDT
      XLL = XL-OMGASM+TEMP
      IF (ISYNFL.EQ.0)    XLL = XL+TEMP+TEMP
      RETURN
C
C     DOT TERMS CALCULATED
C
  150 IF (ISYNFL.EQ.0)    GO TO 152
      XNDOT=DEL1*SIN (XLI-FASX2)+DEL2*SIN (2.*(XLI-FASX4))
     1     +DEL3*SIN (3.*(XLI-FASX6))
      XNDDT = DEL1*COS(XLI-FASX2)
     *       +2.*DEL2*COS(2.*(XLI-FASX4))
     *       +3.*DEL3*COS(3.*(XLI-FASX6))
      GO TO 154
  152 XOMI = OMEGAQ+OMGDT*ATIME
      X2OMI = XOMI+XOMI
      X2LI = XLI+XLI
      XNDOT = D2201*SIN(X2OMI+XLI-G22)
     *       +D2211*SIN(XLI-G22)
     *       +D3210*SIN(XOMI+XLI-G32)
     *       +D3222*SIN(-XOMI+XLI-G32)
     *       +D4410*SIN(X2OMI+X2LI-G44)
     *       +D4422*SIN(X2LI-G44)
     *       +D5220*SIN(XOMI+XLI-G52)
     *       +D5232*SIN(-XOMI+XLI-G52)
     *       +D5421*SIN(XOMI+X2LI-G54)
     *       +D5433*SIN(-XOMI+X2LI-G54)
      XNDDT = D2201*COS(X2OMI+XLI-G22)
     *       +D2211*COS(XLI-G22)
     *       +D3210*COS(XOMI+XLI-G32)
     *       +D3222*COS(-XOMI+XLI-G32)
     *       +D5220*COS(XOMI+XLI-G52)
     *       +D5232*COS(-XOMI+XLI-G52)
     *       +2.*(D4410*COS(X2OMI+X2LI-G44)
     *       +D4422*COS(X2LI-G44)
     *       +D5421*COS(XOMI+X2LI-G54)
     *       +D5433*COS(-XOMI+X2LI-G54))
  154 XLDOT=XNI+XFACT
      XNDDT = XNDDT*XLDOT
C modification by R. Noomen, March 4, 1999
C     GO TO IRETN
      IF ( IRETN .EQ. 140 ) GOTO 140
      IF ( IRETN .EQ. 165 ) GOTO 165
C  eom
C
C     INTEGRATOR
C
C modification by R. Noomen, March 4, 1999
C 160 ASSIGN 165 TO IRETN
  160 IRETN = 165
C  eom
      GO TO 150
  165 XLI = XLI+XLDOT*DELT+XNDOT*STEP2
      XNI = XNI+XNDOT*DELT+XNDDT*STEP2
      ATIME=ATIME+DELT
C modification by R. Noomen, March 4, 1999
C     GO TO IRET
      IF ( IRET .EQ. 100 ) GOTO 100
      IF ( IRET .EQ. 125 ) GOTO 125
C  eom
C
C     EPOCH RESTART
C
  170 IF (T.GE.0.D0)    GO TO 175
      DELT=STEPN
      GO TO 180
  175 DELT = STEPP
  180 ATIME = 0.D0
      XNI=XNQ
      XLI=XLAMO
      GO TO 125
C
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE DPPER(EM,XINC,OMGASM,XNODES,XLL)
      SUBROUTINE DPPER(EM,XINC,OMGASM,XNODES,XLL,T,COSIQ,SINIQ)
C
C     ENTRANCES FOR LUNAR-SOLAR PERIODICS
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
C
      COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,
     .          IEXP , BSTAR , IBEXP ,
     1          X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
      COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
     1           XJ3,XKE,XKMPER,XMNPDA,AE
C     COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
      COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2 , EPS
      COMMON / COMRN1 / SSE , SSG , SSH , SSI , SSL
      COMMON / COMRN2 / XLI , XNI , ATIME , STEPP , STEPN , STEP2
      COMMON / COMRN3 / DEL1 , DEL2 , DEL3
      COMMON / COMRN4 / FASX2 , FASX4 , FASX6
      COMMON / COMRN5 / D2201 , D2211 , D3210 , D3222 , D4410 , D4422 ,
     .                  D5220 , D5232 , D5421 , D5433
      COMMON / COMRN6 / PREEP
      COMMON / COMRN7 / IRESFL , ISYNFL
      COMMON / COMRN8 / THGR
      COMMON / COMRN9 / OMEGAQ
      COMMON / COMRN10 / XFACT , XNQ , XLAMO
      COMMON / COMRN11 / XL2 , XL3 , XL4
      COMMON / COMRN12 / XH2 , XH3
      COMMON / COMRN13 / XGH2 , XGH3 , XGH4
      COMMON / COMRN14 / SAVTSN , ZMOS , ZMOL
      COMMON / COMRN15 / SE2 , SE3 , SH2 , SH3 ,
     .                   SI2 , SI3 , SL2 , SL3 , SL4
      COMMON / COMRN16 / SGH2 , SGH3 , SGH4
      COMMON / COMRN17 / E3 , EE2 , XI2 , XI3
      COMMON / COMRN18 / XQNCL
      COMMON / COMRN19 / ZNS, C1SS, ZES, ZNL, C1L, ZEL,
     .                   ZCOSIS, ZSINIS, ZSINGS, ZCOSGS, ZCOSHS, ZSINHS
      COMMON / COMRN20 / Q22,Q31,Q33, G22,G32, G44,G52, G54,
     .                   ROOT22,ROOT32, ROOT44,ROOT52, ROOT54,
     .                   THDT, RHO
C
C     ENTRANCES FOR LUNAR-SOLAR PERIODICS
C
      SINIS = SIN(XINC)
      COSIS = COS(XINC)
      IF ( ABS(SAVTSN-T).LT.(30.D0))    GO TO 210
      SAVTSN=T
      ZM=ZMOS+ZNS*T
  205 ZF=ZM+2.*ZES*SIN (ZM)
      SINZF=SIN (ZF)
      F2=.5*SINZF*SINZF-.25
      F3=-.5*SINZF*COS (ZF)
      SES=SE2*F2+SE3*F3
      SIS=SI2*F2+SI3*F3
      SLS=SL2*F2+SL3*F3+SL4*SINZF
      SGHS=SGH2*F2+SGH3*F3+SGH4*SINZF
      SHS=SH2*F2+SH3*F3
      ZM=ZMOL+ZNL*T
      ZF=ZM+2.*ZEL*SIN (ZM)
      SINZF=SIN (ZF)
      F2=.5*SINZF*SINZF-.25
      F3=-.5*SINZF*COS (ZF)
      SEL=EE2*F2+E3*F3
      SIL=XI2*F2+XI3*F3
      SLL=XL2*F2+XL3*F3+XL4*SINZF
      SGHL=XGH2*F2+XGH3*F3+XGH4*SINZF
      SHL=XH2*F2+XH3*F3
      PE=SES+SEL
      PINC=SIS+SIL
      PL=SLS+SLL
  210 PGH=SGHS+SGHL
      PH=SHS+SHL
      XINC = XINC+PINC
      EM = EM+PE
      IF(XQNCL.LT.(.2)) GO TO 220
      GO TO 218
C
C     APPLY PERIODICS DIRECTLY
C
  218 PH=PH/SINIQ
      PGH=PGH-COSIQ*PH
      OMGASM=OMGASM+PGH
      XNODES=XNODES+PH
      XLL = XLL+PL
      GO TO 230
C
C     APPLY PERIODICS WITH LYDDANE MODIFICATION
C
  220 SINOK=SIN(XNODES)
      COSOK=COS(XNODES)
      ALFDP=SINIS*SINOK
      BETDP=SINIS*COSOK
      DALF=PH*COSOK+PINC*COSIS*SINOK
      DBET=-PH*SINOK+PINC*COSIS*COSOK
      ALFDP=ALFDP+DALF
      BETDP=BETDP+DBET
      XLS = XLL+OMGASM+COSIS*XNODES
      DLS=PL+PGH-PINC*XNODES*SINIS
      XLS=XLS+DLS
      XNODES=ATAN2(ALFDP,BETDP)
      XLL = XLL+PL
      OMGASM = XLS-XLL-COS(XINC)*XNODES
  230 CONTINUE
C
      RETURN
C
      END
