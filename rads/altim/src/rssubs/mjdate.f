**MJDATE -- Convert MJD to YYMMDD or v.v.
*+
      SUBROUTINE MJDATE(IND, MJD, YYMMDD, YY, MM, DD)
      INTEGER*4 IND, MJD, YYMMDD, YY, MM, DD
*
* This function converts Modified Julian Dates (MJD) to Year-Month-Day
* in the format YYMMDD, or vice versa.
* The algorithm is good for the years 1900 till 1999 only.
*
* Arguments:
*  IND      (input): IND=1, convert MJD --> YYMMDD and YY, MM, DD
*                    IND=2, convert YYMMDD --> MJD and YY, MM, DD
*                    IND=3, convert YY, MM, DD --> MJD
*                    IND=4, convert YYMMDD --> YY, MM, DD
*  YYMMDD  (in/out): Date in the form 851231.
*  YY, MM, DD (i/o): Date as separate integers, e.g. 85, 12, 31.
*  MJD     (in/out): Modified Julian Day, e.g. 46430
*-
* 27-Sep-1994 - Created from the subroutine MJDATE by Ron Noomen
*-----------------------------------------------------------------------
      REAL*8 XMJD, YEAR
      INTEGER*4 MJCENT, IH, ID, MONTH(13)
      PARAMETER (YEAR=365.25D0, MJCENT=15019, IH=100)
      DATA MONTH / 0,31,59,90,120,151,181,212,243,273,304,334,365 /
      SAVE MONTH
C
      GOTO (30, 10, 20, 10), IND
C        CONVERT YYMMDD TO YY, MM, DD
  10  MM = YYMMDD/IH
      DD = YYMMDD - IH*MM
      YY = MM/IH
      MM = MM - IH *YY
      IF (IND.EQ.4) RETURN
C
C        CONVERT YY, MM, DD TO MJD
  20  MJD = MJCENT + DD + MONTH(MM) + (1461*YY - 1)/4
      IF (YY.EQ.4*(YY/4) .AND. MM.GT.2) MJD = MJD + 1
      RETURN
C
C        CONVERT MJD TO YY, MM, DD
  30  XMJD = MJD - MJCENT
      YY = INT(XMJD/YEAR+1.0D-10)
      ID = INT(XMJD -DBLE(YY)*YEAR+1.D0+1.0D-10)
      MM = ID/30 + 1
      IF (YY.EQ.4*(YY/4)) THEN
	 IF (ID.EQ.60) MM = MM - 1
	 IF (ID.GT.60) ID = ID - 1
      ENDIF
  40  DD = ID - MONTH(MM)
      IF (DD.LE.0) THEN
         MM = MM - 1
         GOTO 40
      ENDIF
C        CONVERT YY, MM, DD TO YYMMDD
      YYMMDD = (YY*IH + MM)*IH + DD
      RETURN
C
      END
