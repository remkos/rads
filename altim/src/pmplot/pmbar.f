**PMBAR -- draw scale-bar
*+
      SUBROUTINE PMBAR (OPTIONS, DISP, LENGTH, NMAIN, NSUB)
      CHARACTER*(*) OPTIONS
      REAL DISP, LENGTH
      INTEGER NMAIN, NSUB
*
* Draw a horizontal scale-bar at the top or bottom of the viewport. This routine
* is only valid when the geographical projections are used. The
* scale can be defined with PMDEF; the window must have been specified
* by PMWINDOW.
*
* The bar will look like this   0        10       20        30 km
*                               II  II  II    *    IIIIIIIIII
* (* is the center of the bar)          Scale 1:40000
*
* The dimension of the various items in the bar (in units of character size):
* Kilometer scale             100%
* Height of the scale-bar      50%
* Scale                       120%
* Total height (excl. scale)  180%
* Total height (incl. scale)  310%
*
* Arguments:
*  OPTIONS  (input) : String of plot options as defined below.
*  DISP     (input) : Displacement of the center of the scale bar from
*                     the top or the bottom of the viewport, measured
*                     outwards in units of character height.
*                     Use a negative value to put the bar inside the viewport,
*                     a positive to put it outside.
*  LENGTH   (input) : Length of the scale-bar in kilometers (Example: 30.).
*                     If LENGTH=0.0, PMBAR will choose the length such
*                     that it is at least 20 character units long, but no
*                     more than 60. NMAIN and NSUB will then also be set
*                     automatically.
*          (output) : Length chosen by PMBAR.
*  NMAIN    (input) : Number of main intervals (Example: 3).
*  NSUB     (input) : Number of sub-intervals (Example: 5).
*
* Options (for parameter OPTIONS):
*  T : draw scale-bar at the top of viewport (Default).
*  B : draw scale-bar at the bottom of viewport.
*  L : align the left of the scale-bar with the left of the viewport.
*  C : align the center of the scale-bar with the center of the
*      viewport (Default).
*  R : align the right of the scale-bar with the right of the viewport.
*  S : put the scale at the bottom of the scale-bar.
*--
*  8-Jan-1991 - created [Remko Scharroo]
* 13-Dec-1991 - Select length of scale-bar automatically if LENGTH < 1 meter
* 13-Jan-1992 - Standardize PMPLOT.
*  7-Jul-1994 - All variables defined
* 18-Jul-1994 - PGPIXL used iso PGGRAY
* 18-Jan-1996 - Adjusted to PGPLOT v5.0
* 11-Jul-1996 - Adjusted to PGPLOT v5.1
*-----------------------------------------------------------------------
      CHARACTER OPT*64
      REAL HEIGHT /.5/, TOPOFF /.3/, BOTOFF /1.1/, FONTSZ /1.2/
      LOGICAL BOT, LEFT, RIGHT, SCALE, PGNOTO
      REAL X0, X1, XS0, XS1, XV0, XV1, XINT
      REAL Y0, Y1, YS0, YS1, YV0, YV1
      REAL CHARSZ, CHRHGT, BARHGT, BARLEN, PGRND
      INTEGER I, NC /0/, NXSUB, NXINT, NMAX, FILL
      PARAMETER (NMAX=10)
      INTEGER A(NMAX)
      REAL KPERIN,NEWLNG
      INCLUDE 'pgplot.inc'
      INCLUDE 'pmplot.inc'
*
      IF (PGNOTO('PMBAR') .OR. PMOPEN.LT.4) RETURN
      CALL PGBBUF
*
* Set parameters
*
      CALL GRTOUP(OPT,OPTIONS)
      BOT=INDEX(OPT,'B').NE.0
      LEFT=INDEX(OPT,'L').NE.0
      RIGHT=INDEX(OPT,'R').NE.0
      SCALE=INDEX(OPT,'S').NE.0
      CHRHGT=PGYSP(PGID)/PGYPIN(PGID)
      BARHGT=HEIGHT*CHRHGT
*
* Compute main interval if requested (KPERIN=kilometer per inch)
*
      KPERIN=25.4E-6*PSCALE
      IF (LENGTH.LE.1E-3) THEN
         XINT=PGRND(MAX(1.,4*KPERIN/PGXPIN(PGID)*PGXSP(PGID)),NXSUB)
         IF (NXSUB.EQ.2) NXSUB=4
         NXINT=5
      ELSE
         NXINT=MAX(1,MIN(NMAX,NMAIN))
         NXSUB=MAX(1,MIN(NMAX,NSUB))
         XINT=MAX(1.,ANINT(LENGTH/NXINT))
      ENDIF
      NEWLNG=NXINT*XINT
      BARLEN=NEWLNG/KPERIN
*
* Query viewport and window dimensions
*
      CALL PGQVP(1,XV0,XV1,YV0,YV1)
      CALL PGQWIN(X0,X1,Y0,Y1)
*
* Set new viewport and window
*
      IF (LEFT) THEN
         XS0=XV0
      ELSE IF (RIGHT) THEN
         XS0=XV1-BARLEN
      ELSE
         XS0=(XV0+XV1-BARLEN)/2
      ENDIF
      XS1=XS0+BARLEN
      IF (BOT) THEN
         YS0=YV0-DISP*CHRHGT-.5*BARHGT
      ELSE
         YS0=YV1+DISP*CHRHGT-.5*BARHGT
      ENDIF
      YS1=YS0+BARHGT
      CALL PGVSIZ(XS0,XS1,YS0,YS1)
      CALL PGSWIN(0.,NEWLNG,0.,1.)
*
* Fill bar with black/white blocks
*
      FILL=0
      DO I=NXSUB,1,-1
         FILL=1-FILL
         A(I)=FILL
      ENDDO
      CALL PGPIXL(A,NXSUB,1,1,NXSUB,1,1,0.,XINT,0.,1.)
      FILL=0
      DO I=1,NXINT
         FILL=1-FILL
         A(I)=FILL
      ENDDO
      CALL PGPIXL(A,NXINT,1,2,NXINT,1,1,XINT,NEWLNG,0.,1.)
*
* Draw box and numbers at main intervals
*
      CALL PGBOX('BC',1.,1,'BC',1.,1)
      DO I=0,NXINT
         CALL PGNUMB(NINT(I*XINT),0,1,OPT,NC)
         CALL PGMTXT('T',TOPOFF,I*XINT/NEWLNG,.5,OPT(1:NC))
      ENDDO
*
* Write 'km' and 'Scale 1:?????'
*
      CALL PGMTXT('T',TOPOFF,1.,-.25*NC-.15,'km')
      IF (SCALE) THEN
         OPT='Scale 1:'
         CALL PGNUMB(NINT(PSCALE),0,1,OPT(9:),NC)
         CALL PGQCH(CHARSZ)
         CALL PGSCH(CHARSZ*FONTSZ)
         CALL PGMTXT('B',BOTOFF,0.5,0.5,OPT(1:NC+8))
         CALL PGSCH(CHARSZ)
      ENDIF
*
* Restore viewport and window
*
      CALL PGVSIZ(XV0,XV1,YV0,YV1)
      CALL PGSWIN(X0,X1,Y0,Y1)
      CALL PGEBUF
      END
