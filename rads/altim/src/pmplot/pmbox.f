**PMBOX -- draw labeled frame around viewport
*+
      SUBROUTINE PMBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
      CHARACTER*(*) XOPT, YOPT
      REAL XTICK, YTICK
      INTEGER NXSUB, NYSUB
*
* Annotate the viewport with frame, grid, ticks, numeric labels, etc.
* PMBOX is rather similar to PGBOX, but is used especially for non-linear
* geographical projections.
* PMDEF and PMWINDOW must have been called before PMBOX.
*
* Arguments:
*  XOPT   (input)  : String of options for X (horizontal, longitude) axis of
*                    plot. Options are single letters, and may be in
*                    any order (see below).
*  XTICK  (input)  : Longitude interval (degrees) between major tick marks on
*                    X axis. If XTICK=0.0, the interval is chosen by PMBOX.
*  NXSUB  (input)  : Number of subintervals to divide the major coordinate
*                    interval into. If XTICK=0.0 or NXSUB=0, the number is
*                    chosen by PMBOX.
*  YOPT   (input)  : String of options for Y (vertical, latitude) axis of plot.
*                    Coding is the same as for XOPT.
*  YTICK  (input)  : Latitude interval (degrees) between major tick marks on Y
*                    axis. If YTICK=0.0, the interval will be the same as XTICK.
*  NYSUB  (input)  : Like NXSUB for the Y axis.
*
* Options (for parameters XOPT and YOPT) common to all projections:
*  B : Draw bottom (X) or left (Y) edge of frame.
*  C : Draw top (X) or right (Y) edge of frame.
*  G : Draw Grid of meridians (X) or parallels (Y).
*  H : (projections 32 & 34 only) Draw frame around area and remove anything
*      outside.
*  I : Invert the tick marks (plot them outside the viewport).
*  P : extend ("Project") major tick marks outside the box (if option I
*      is not set) or inside the box (if option I is set).
*  2 : Increase the number of steps used to draw curved meridians (X) or
*      parallels (Y) by a factor 2. The default number is 50 or half the
*      length of the meridians or parallels measured in degrees (whichever
*      is greater).
*  5 : Increase the number of steps used to draw curved meridians (X) or
*      parallels (Y) by a factor 5. If both options 2 and 5 are used the
*      number of steps is increased by a factor 10.
* To get a complete frame, specify BC in both XOPT and YOPT.
*
* Options used with projections other than azimuthal:
*  N : Write Numeric labels in the conventional location below the
*      viewport (X) or to the left of the viewport (Y).
*  M : Write numeric labels in the unconventional location above the
*      viewport (X) or to the right of the viewport (Y).
*  . : Write numeric labels in a decimal notation. The default notation is
*      to have arcminutes as a superscript (unless they are zero).
*  : : Write numeric labels in a notation degrees:minutes. The default
*      notation is to have arcminutes as a superscript (unless they are zero).
*  - : Write - for West or South. Default is no sign.
*  T : Draw major Tick marks at the major coordinate interval.
*      (Ignored if option G is specified. Tick marks are only drawn when B or C
*      are used.)
*  S : Draw minor tick marks (Subticks).
*  V : Orient numeric labels Vertically. This is only applicable to Y.
*      The default is to write Y-labels parallel to the axis
*
* Options used for azimuthal and polar projections only:
*  H : Draw observer's Horizon and remove all graphics outside the horizon.
*      Options B and C will have the same effect.
*  T : Omit the meridians in the vicinity of the poles. Meridians will extend
*      up to one major Y-interval from the poles. (Effective for XOPT only.)
*  S : Meridians will not be draw within one sub-interval from the poles.
*      (Used only for XOPT).
*  (Ticks or numerals can not be drawn in this mode. N, M, and V are ignored.)
*--
*  9-Jan-1991 - created [Remko Scharroo]                      
* 14-Jan-1991 - Major changes to implement conic projection.
* 16-Jan-1991 - Support azimuthal projections.
* 12-Dec-1991 - Tick interval chosen when XTICK or YTICK less then 1.E6.
* 13-Jan-1992 - Standardize PMPLOT.
* 19-Mar-1992 - Call PMRND to find interval size and make BC=H.
* 27-Apr-1993 - Use new PGNUMB routine. Include options '.' and ':'
* 11-Jun-1993 - Increase number of digits sent to PGNUMB
* 22-Jul-1994 - Include option I
* 10-Aug-1994 - Include option P. Use absolute numerals, unless - is
*               specified.
* 28-Jun-1995 - H also in effect for Azimuthal.
* 16-Jan-1996 - Bug removed. Wrong major ticks in project <> 1
* 11-Jul-1996 - Adjusted to PGPLOT 5.1
*-----------------------------------------------------------------------
      INCLUDE 'pgplot.inc'
      INCLUDE 'pmplot.inc'
      CHARACTER  CLBL*20, OPT*64
      LOGICAL  XOPTB, XOPTC, XOPTG, XOPTN, XOPTM, XOPTT, XOPTS, XOPTI
      LOGICAL  YOPTB, YOPTC, YOPTG, YOPTN, YOPTM, YOPTT, YOPTS, YOPTI
      LOGICAL  XOPTH, XOPT2, XOPT5, XOPTP, XABS
      LOGICAL  YOPTH, YOPT2, YOPT5, YOPTP, YABS
      LOGICAL  YOPTV, IN, IRANGE, MAJOR, OPTH, PGNOTO
      INTEGER XNFORM, YNFORM
      REAL     X(5), Y(5)
      REAL     LAT, LON
      INTEGER CI,FS,I,NP,NC,NV,J,PMBOX1
      REAL A,B,C
      REAL PMRND,PMX,PMY
      REAL TIKLA,TIKLB,TIKL1,TIKL2
      REAL XC,      XINT,XINT2,XVAL
      REAL YC,Y0,Y1,YINT,YINT2,YVAL
      INTEGER IX1,IX2,NXSTEP,NSUBX
      INTEGER IY1,IY2,NYSTEP,NSUBY
*
* Old common arguments
*
      REAL XBLC,XTRC,XLEN,XOFF,XSP
      REAL YBLC,YTRC,YLEN,YOFF,YSP
C
      IRANGE(A,B,C) = (A.LE.B.AND.B.LE.C) .OR. (C.LE.B.AND.B.LE.A)
C
      IF (PGNOTO('PMBOX') .OR. PMOPEN.LT.4) RETURN
      IF (PTYPE.LE.0) THEN
	 CALL PGBOX(XOPT,XTICK,NXSUB,YOPT,YTICK,NYSUB)
	 RETURN
      ENDIF
      CALL PGBBUF
C
C Decode options.
C
      CALL GRTOUP(OPT,XOPT)
      XOPTB = INDEX(OPT,'B').NE.0
      XOPTC = INDEX(OPT,'C').NE.0
      XOPTG = INDEX(OPT,'G').NE.0
      XOPTH = INDEX(OPT,'H').NE.0
      XOPTI = INDEX(OPT,'I').NE.0
      XOPTM = INDEX(OPT,'M').NE.0
      XOPTN = INDEX(OPT,'N').NE.0
      XOPTS = INDEX(OPT,'S').NE.0
      XOPTT = INDEX(OPT,'T').NE.0
      XOPT2 = INDEX(OPT,'2').NE.0
      XOPT5 = INDEX(OPT,'5').NE.0
      XOPTP = INDEX(OPT,'P').NE.0
      XNFORM = 4
      IF (INDEX(OPT,'.').NE.0) XNFORM=1
      IF (INDEX(OPT,':').NE.0) XNFORM=5
      XABS = INDEX(OPT,'-').EQ.0

      CALL GRTOUP(OPT,YOPT)
      YOPTB = INDEX(OPT,'B').NE.0
      YOPTC = INDEX(OPT,'C').NE.0
      YOPTG = INDEX(OPT,'G').NE.0
      YOPTH = INDEX(OPT,'H').NE.0
      YOPTI = INDEX(OPT,'I').NE.0
      YOPTM = INDEX(OPT,'M').NE.0
      YOPTN = INDEX(OPT,'N').NE.0
      YOPTS = INDEX(OPT,'S').NE.0
      YOPTT = INDEX(OPT,'T').NE.0
      YOPTV = INDEX(OPT,'V').NE.0
      YOPT2 = INDEX(OPT,'2').NE.0
      YOPT5 = INDEX(OPT,'5').NE.0
      YOPTP = INDEX(OPT,'P').NE.0
      YNFORM = 4
      IF (INDEX(OPT,'.').NE.0) YNFORM=1
      IF (INDEX(OPT,':').NE.0) YNFORM=5
      YABS = INDEX(OPT,'-').EQ.0

      OPTH = XOPTB.OR.XOPTC.OR.XOPTH.OR.YOPTB.OR.YOPTC.OR.YOPTH

*
* Convert a few PGPLOT 5.1 naming conventions to the old ones
*
      XBLC=PGXBLC(PGID)
      XTRC=PGXTRC(PGID)
      XLEN=PGXLEN(PGID)
      XSP =PGXSP (PGID)
      XOFF=PGXOFF(PGID)
      YBLC=PGYBLC(PGID)
      YTRC=PGYTRC(PGID)
      YLEN=PGYLEN(PGID)
      YSP =PGYSP (PGID)
      YOFF=PGYOFF(PGID)
*
* Set step sizes for drawing meridians (X) and parallels (Y)
*
      IF (XCURVE) THEN
         NXSTEP=MAX(50,NINT((LATMAX-LATMIN)/2.0))
         IF (XOPT2) NXSTEP=NXSTEP*2
         IF (XOPT5) NXSTEP=NXSTEP*5
      ELSE
         NXSTEP=1
      ENDIF
      IF (YCURVE) THEN
         NYSTEP=MAX(50,NINT((LONMAX-LONMIN)/2.0))
         IF (YOPT2) NYSTEP=NYSTEP*2
         IF (YOPT5) NYSTEP=NYSTEP*5
      ELSE
         NYSTEP=1
      ENDIF
C
C Remove window.
C
      CALL GRAREA(PGID,0.,0.,-1.,-1.)
*
* Draw observer's horizon
*
*     IF (AZMTAL.AND.(XOPTH.OR.YOPTH.OR.XOPTB.OR.XOPTC.OR.
*    .YOPTB.OR.YOPTC)) THEN
*        CALL GRMOVA(XTRC,0.)
*        XVAL=2*PI/NXSTEP
*        DO 10 I=1,NXSTEP
*  10       CALL GRLINA(XTRC*COS(I*XVAL),XTRC*SIN(I*XVAL))
*     ENDIF
C
C Draw box.
C
      IF (XOPTB) THEN
          CALL GRMOVA(XBLC, YBLC)
          CALL GRLINA(XTRC, YBLC)
      END IF
      IF (YOPTC) THEN
          CALL GRMOVA(XTRC, YBLC)
          CALL GRLINA(XTRC, YTRC)
      END IF
      IF (XOPTC) THEN
          CALL GRMOVA(XTRC, YTRC)
          CALL GRLINA(XBLC, YTRC)
      END IF
      IF (YOPTB) THEN
          CALL GRMOVA(XBLC, YTRC)
          CALL GRLINA(XBLC, YBLC)
      END IF
C
C Length of major/minor tick marks.
C
      TIKL1 = XSP*0.6*(YTRC-YBLC)/YLEN
      IF (XOPTI) TIKL1 = -TIKL1
      TIKL2 = TIKL1*0.5
C
C Choose X tick intervals. Major interval = XINT,
C minor interval = XINT2 = XINT/NSUBX.
C
      IF (ABS(XTICK).LT.1E-6) THEN
         XINT = MAX(0.05, MIN(10.0*XSP/XLEN, 0.20))*(LONMAX-LONMIN)
         XINT = PMRND(XINT,NSUBX)
      ELSE
         XINT = XTICK
         NSUBX = MAX(NXSUB,1)
      END IF
      IF (.NOT.XOPTS) NSUBX = 1
      XINT2 = XINT/NSUBX
      NP = INT(ALOG10(ABS(XINT2)))-5
      NV = NINT(XINT2/10.**NP)
      CALL PGBOX1(LONMIN,LONMAX,XINT2,IX1,IX2)
*
* Choose Y tick intervals.
*
      IF (ABS(YTICK).LT.1E-6.AND.XOPT.EQ.YOPT) THEN
         YINT = MAX(0.05, MIN(7.0*YSP/YLEN, 0.25))*(LATMAX-LATMIN)
         YINT = PMRND(YINT,NSUBY)
      ELSE IF (ABS(YTICK).LT.1E-6) THEN
         YINT = XINT
         NSUBY = NSUBX
      ELSE
         YINT = YTICK
         NSUBY = MAX(NYSUB,1)
      END IF
      IF (.NOT.YOPTS) NSUBY = 1
      YINT2 = YINT/NSUBY
      CALL PGBOX1(LATMIN,LATMAX,YINT2,IY1,IY2)
*
* Start longitude loop
*
      DO I=IX1,IX2
         MAJOR=(MOD(I,NSUBX).EQ.0)
         LON=I*XINT2
         TIKLA=0
         IF (MAJOR.AND.XOPTP) TIKLA=TIKL2
         TIKLB=TIKL2
         IF (MAJOR.AND.XOPTT) TIKLB=TIKL1
*
* Draw bottom ticks and numerals. For orthographic and perspective
* projection, cut grid one (sub-)interval from the poles.
*
         XVAL=PMX(YBLC,LON)
         XC=(XVAL-XBLC)/(XTRC-XBLC)
         IN=IRANGE(XBLC,XVAL,XTRC)
         IF ((AZMTAL.OR.PTYPE.EQ.41.OR.PTYPE.EQ.42)
     .		.AND.(YOPTT.OR.YOPTS)) THEN
            Y0=LATMIN+YINT2
            Y1=LATMAX-YINT2
         ELSE
            Y0=LATMIN
            Y1=LATMAX
         ENDIF
         IF ((XOPTT.OR.XOPTS).AND.XOPTB.AND.IN) THEN
	    Y(1)=YBLC-TIKLA
	    Y(2)=YBLC+TIKLB
	    X(1)=XVAL
	    X(2)=PMX(Y(2),LON)
	    CALL PGLINE(2,X,Y)
         ENDIF
         IF (MAJOR.AND.XOPTN.AND.IN) THEN
            CALL PGNUMB(PMBOX1(XABS,LON,NP),NP,XNFORM,CLBL,NC)           
            CALL PGMTXT('B', 1.2, XC, 0.5, CLBL(1:NC))
         ENDIF

* Draw X grid

         IF (MAJOR.AND.XOPTG) THEN
            CALL GRAREA(PGID,XOFF,YOFF,XLEN,YLEN)
            DO J=0,NYSTEP
               XVAL=LON
               YVAL=Y0+J*(Y1-Y0)/NYSTEP
               CALL PMCONV(1,XVAL,YVAL)
               IF (J.EQ.0) THEN
                  CALL GRMOVA(XVAL,YVAL)
               ELSE
                  CALL GRLINA(XVAL,YVAL)
               ENDIF
	    enddo
            CALL GRAREA(PGID,0.,0.,-1.,-1.)
         ENDIF
*
* Draw top ticks and numerals
*
         XVAL=PMX(YTRC,LON)
         XC=(XVAL-XBLC)/(XTRC-XBLC)
         IN=IRANGE(XBLC,XVAL,XTRC)
         IF ((XOPTT.OR.XOPTS).AND.XOPTC.AND.IN) THEN
	    Y(1)=YTRC+TIKLA
	    Y(2)=YTRC-TIKLB
	    X(1)=XVAL
	    X(2)=PMX(Y(2),LON)
	    CALL PGLINE(2,X,Y)
         ENDIF
         IF (MAJOR.AND.XOPTM.AND.IN) THEN
            CALL PGNUMB(PMBOX1(XABS,LON,NP),NP,XNFORM,CLBL,NC)
            CALL PGMTXT('T', 0.7, XC, 0.5, CLBL(1:NC))
         ENDIF
      enddo
*
* Length of Y tick marks.
*
      TIKL1 = XSP*0.6*(XTRC-XBLC)/XLEN
      IF (YOPTI) TIKL1 = -TIKL1
      TIKL2 = TIKL1*0.5
      NP = INT(ALOG10(ABS(YINT2)))-5
      NV = NINT(YINT2/10.**NP)
*
* Start latitude loop
*
      DO I=IY1,IY2
         MAJOR=(MOD(I,NSUBY).EQ.0)
         LAT=I*YINT2
         TIKLA=0
         IF (MAJOR.AND.YOPTP) TIKLA=TIKL2
         TIKLB=TIKL2
         IF (MAJOR.AND.YOPTT) TIKLB=TIKL1

* Draw left ticks and numerals

         YVAL=PMY(XBLC,LAT)
         YC=(YVAL-YBLC)/(YTRC-YBLC)
         IN=IRANGE(YBLC,YVAL,YTRC)
         IF ((YOPTT.OR.YOPTS).AND.YOPTB.AND.IN) THEN
	    X(1)=XBLC-TIKLA
	    X(2)=XBLC+TIKLB
	    Y(1)=YVAL
	    Y(2)=PMY(X(2),LAT)
	    CALL PGLINE(2,X,Y)
         ENDIF
         IF (MAJOR.AND.YOPTN.AND.IN) THEN
            CALL PGNUMB(PMBOX1(YABS,LAT,NP),NP,YNFORM,CLBL,NC)
            IF (YOPTV) THEN
               CALL PGMTXT('LV',0.7,YC,1.,CLBL(1:NC))
            ELSE
               CALL PGMTXT('L',0.7,YC,.5,CLBL(1:NC))
            ENDIF
         ENDIF

* Draw Y grid

         IF (MAJOR.AND.YOPTG) THEN
            CALL GRAREA(PGID,XOFF,YOFF,XLEN,YLEN)
            DO J=0,NXSTEP
               XVAL=LONMIN+J*(LONMAX-LONMIN)/NXSTEP
               YVAL=LAT
               CALL PMCONV(1,XVAL,YVAL)
               IF (J.EQ.0) THEN
                  CALL GRMOVA(XVAL,YVAL)
               ELSE
                  CALL GRLINA(XVAL,YVAL)
               ENDIF
            enddo
            CALL GRAREA(PGID,0.,0.,-1.,-1.)
         ENDIF
*
* Draw right ticks and numerals
*
         YVAL=PMY(XTRC,LAT)
         YC=(YVAL-YBLC)/(YTRC-YBLC)
         IN=IRANGE(YBLC,YVAL,YTRC)
         IF ((YOPTT.OR.YOPTS).AND.YOPTC.AND.IN) THEN
	    X(1)=XTRC+TIKLA
	    X(2)=XTRC-TIKLB
	    Y(1)=YVAL
	    Y(2)=PMY(X(2),LAT)
	    CALL PGLINE(2,X,Y)
         ENDIF
         IF (MAJOR.AND.YOPTM.AND.IN) THEN
            CALL PGNUMB(PMBOX1(YABS,LAT,NP),NP,YNFORM,CLBL,NC)
            IF (YOPTV) THEN
               CALL PGMTXT('RV',1.2,YC,0.,CLBL(1:NC))
            ELSE
               CALL PGMTXT('R',1.2,YC,.5,CLBL(1:NC))
            ENDIF
         ENDIF
      enddo
*
* For tilted rectangular projection only: draw frame and remove
*
      IF (PTYPE.EQ.34 .AND. OPTH) THEN
	 X(1)=LONUMAX
	 Y(1)=LATUMIN
	 X(2)=LONUMAX
	 Y(2)=LATUMAX
	 X(3)=LONUMIN
	 Y(3)=LATUMAX
	 X(4)=LONUMIN
	 Y(4)=LATUMIN
	 CALL PGQFS(FS)
	 CALL PGQCI(CI)
	 CALL PMCONV(4,X,Y)
	 XVAL=X(3)
	 YVAL=Y(3)
	 CALL PGSFS(1)
	 CALL PGSCI(0)
	 IF (PPARA2.GT.0) THEN
	    X(5)=XBLC
	    Y(5)=YTRC
	    CALL PGPOLY(3,X(3),Y(3))
	    X(3)=XTRC
	    Y(3)=YBLC
	    CALL PGPOLY(3,X,Y)
	 ELSE IF (PPARA2.LT.0) THEN
	    X(5)=XBLC
	    Y(5)=YBLC
	    CALL PGPOLY(3,X(3),Y(3))
	    X(3)=XTRC
	    Y(3)=YTRC
	    CALL PGPOLY(3,X,Y)
	 ENDIF
	 CALL PGSFS(2)
	 CALL PGSCI(CI)
	 X(3)=XVAL
	 Y(3)=YVAL
	 CALL PGPOLY(4,X,Y)
	 CALL PGSFS(FS)
      ENDIF
*
* For Molweide something similar
*
      IF (PTYPE.EQ.32 .AND. OPTH) THEN
	 XVAL=LONUMIN
	 YVAL=-90.
	 CALL PMCONV(1,XVAL,YVAL)
	 CALL PGMOVE(XVAL,YVAL)
	 DO I=1,90
	    XVAL=LONUMIN
	    YVAL=I*2-90.
	    CALL PMCONV(1,XVAL,YVAL)
	    CALL GRLINA(XVAL,YVAL)
	 ENDDO
	 XVAL=LONUMAX
	 YVAL=-90.
	 CALL PMCONV(1,XVAL,YVAL)
	 CALL PGMOVE(XVAL,YVAL)
	 DO I=1,90
	    XVAL=LONUMAX
	    YVAL=I*2-90.
	    CALL PMCONV(1,XVAL,YVAL)
	    CALL GRLINA(XVAL,YVAL)
	 ENDDO
      ENDIF
*
* For Polar projection: remove corners
*
      if ((ptype.eq.41 .or. ptype.eq.42) .and. opth) then
	 if (ptype.eq.41) yval=latumin
	 if (ptype.eq.42) yval=latumax
	 call pmbox2(  0+loncen, 3.,yval)
	 call pmbox2(  0+loncen,-3.,yval)
	 call pmbox2(180+loncen, 3.,yval)
	 call pmbox2(180+loncen,-3.,yval)
      endif
*
* For Azimuthal projection: remove corners
*
      if (azmtal .and. opth) then
	 call pmbox2(  0., 3.,100.)
	 call pmbox2(  0.,-3.,100.)
	 call pmbox2(180., 3.,100.)
	 call pmbox2(180.,-3.,100.)
      endif
C
C Restore window: interior of box.
C
      CALL GRAREA(PGID,XOFF,YOFF,XLEN,YLEN)
C
C Restore window: interior of box.
C
      CALL GRAREA(PGID,XOFF,YOFF,XLEN,YLEN)
      CALL PGEBUF
      END

      FUNCTION PMBOX1(ABS,X,NP)
      LOGICAL ABS
      INTEGER PMBOX1, NP
      REAL X,Y

      Y=X
      IF (ABS) THEN
         IF (Y.LT.-180) Y=Y+360
         IF (Y.GT.180) Y=Y-360
	 IF (Y.LT.0) Y=-Y
      ENDIF
      PMBOX1=NINT(Y/10.**NP)
      END

      subroutine pmbox2(angle0,dangle,lat)
      real angle0,dangle,lat
      integer i
      real pi,rad,x(0:31),y(0:31)

      pi=4*atan(1.)
      rad=pi/180
      do i=0,30
	 x(i)=angle0+i*dangle
	 y(i)=lat
      enddo
      if (abs(lat).le.90) then
         call pmconv(31,x,y)
         x(31)=x(30)
         y(31)=y(0)
      else
	 do i=0,30
	    y(i)=sin(x(i)*rad)
	    x(i)=cos(x(i)*rad)
	 enddo
         x(31)=x(0)
         y(31)=y(30)
      endif
      call pgsave
      call pgsfs(1)
      call pgsci(0)
      call pgpoly(32,x(0),y(0))
      call pgunsa
      call pgsci(1)
      call pgline(31,x(0),y(0))
      end
