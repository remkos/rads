C*PGSVPX -- set viewport (any unit)
C+
      SUBROUTINE PGSVPX (UNITS, XLEFT, XRIGHT, YBOT, YTOP)
      INTEGER UNITS
      REAL    XLEFT, XRIGHT, YBOT, YTOP
C
C Change the size and position of the viewport, specifying
C the viewport in device coordinates (several possible units specified by
C UNITS).  The viewport is the rectangle on the view surface "through"
C which one views the graph.  All the PG routines which plot lines
C etc. plot them within the viewport, and lines are truncated at
C the edge of the viewport (except for axes, labels etc drawn with
C PGBOX or PGLAB).  The region of world space (the coordinate
C space of the graph) which is visible through the viewport is
C specified by a call to PGSWIN.  It is legal to request a
C viewport larger than the view surface; only the part which
C appears on the view surface will be plotted.
C
C Arguments:
C  UNITS  (input)  : used to specify the units of the output parameters:
C                    UNITS = 0 : normalised device coordinates (as PGSVP).
C                    UNITS = 1 : inches
C                    UNITS = 2 : millimeters
C                    UNITS = 3 : pixels
C                    Other values give an error message, and are
C                    treated as 0.
C  XLEFT  (input)  : x-coordinate of left hand edge of viewport
C  XRIGHT (input)  : x-coordinate of right hand edge of viewport
C  YBOT   (input)  : y-coordinate of bottom edge of viewport
C  YTOP   (input)  : y-coordinate of top  edge of viewport
C--
* 19-Mar-1992 - Introduced to PGPLOT [Remko Scharroo]
* 11-Jul-1996 - Adjusted to PGPLOT 5.1
C-----------------------------------------------------------------------
      INCLUDE  'pgplot.inc'
      LOGICAL PGNOTO
      REAL    XS, YS
C
      IF (PGNOTO('PGSVPX'))  RETURN
      IF (XLEFT.GE.XRIGHT .OR. YBOT.GE.YTOP) THEN
	 CALL GRWARN('PGSVPX ignored: invalid arguments')
	 RETURN
      ENDIF
C
      IF (UNITS.EQ.0) THEN
         XS = PGXSZ(PGID)/PGXPIN(PGID)
         YS = PGYSZ(PGID)/PGYPIN(PGID)
         CALL PGVSIZ(XLEFT*XS, XRIGHT*XS, YBOT*YS, YTOP*YS)
      ELSE IF (UNITS.EQ.1) THEN
         CALL PGVSIZE(XLEFT, XRIGHT, YBOT, YTOP)
      ELSE IF (UNITS.EQ.2) THEN
         CALL PGVSIZE(XLEFT/25.4, XRIGHT/25.4, YBOT/25.4, YTOP/25.4)
      ELSE IF (UNITS.EQ.3) THEN
         CALL PGVSIZE(XLEFT/PGXPIN(PGID), XRIGHT/PGXPIN(PGID),
     .      YBOT/PGYPIN(PGID), YTOP/PGXPIN(PGID))
      ENDIF
      END

      SUBROUTINE PGVP (UNITS, XLEFT, XRIGHT, YBOT, YTOP)
      INTEGER UNITS
      REAL    XLEFT, XRIGHT, YBOT, YTOP
      CALL PGSVPX (UNITS, XLEFT, XRIGHT, YBOT, YTOP)
      END
