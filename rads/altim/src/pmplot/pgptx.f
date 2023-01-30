**PGPTX -- draw one or more graph markers at given size and angle
*+
      SUBROUTINE PGPTX (N, XPTS, YPTS, SIZE, ANGLE, SYMBOL)
      INTEGER N
      REAL XPTS(*), YPTS(*), SIZE(*), ANGLE(*)
      INTEGER SYMBOL
*
* Primitive routine to draw Graph Markers (polymarker). The markers
* are drawn using the current values of attributes color-index and
* line-width (character-font applies if the symbol number is >31).
* The Markers are drawn SIZE times greater than the current character size.
* The angle at which the markers are drawn is given by ANGLE.
* If the point to be marked lies outside the window, no marker is drawn.
* The "pen position" is changed to (XPTS(N),YPTS(N)) in world coordinates
* (if N > 0).
*
* Arguments:
*  N      (input)  : number of points to mark.
*  XPTS   (input)  : world x-coordinates of the points.
*  YPTS   (input)  : world y-coordinates of the points.
*  SIZE   (input)  : size of the points measured in units of current character
*                    size.
*  ANGLE  (input)  : tilt of the character, measured anti-clockwise (degrees).
*  SYMBOL (input)  : code number of the symbol to be drawn at each
*                    point:
*                    -1, -2  : a single dot (diameter = current
*                              line width).
*                    -3..-31 : a regular polygon with ABS(SYMBOL)
*                              edges (style set by current fill style).
*                    0..31   : standard marker symbols.
*                    32..127 : ASCII characters (in current font).
*                              e.g. to use letter F as a marker, let
*                              SYMBOL = ICHAR('F').
*                    > 127  :  a Hershey symbol number.
*--
* 21-Jan-1991 - created [Remko Scharroo].
*  7-Jul-1994 - adapted to PGPLOT v4.9h
* 11-Jul-1996 - adjusted to PGPLOT 5.1
*-----------------------------------------------------------------------
      INCLUDE  'pgplot.inc'
      LOGICAL  PGNOTO
*
      IF (PGNOTO('PGPTX')) RETURN
      IF (N.LT.1) RETURN
      CALL PGBBUF
*
      IF (SYMBOL.GE.0) THEN
          CALL GRMKER(SYMBOL,.FALSE.,N,XPTS,YPTS,.TRUE.,SIZE,ANGLE)
      ELSE
          CALL GRVCT0(3,.FALSE.,N,XPTS,YPTS)
      END IF
      CALL GRMOVA(XPTS(N),YPTS(N))
      CALL PGEBUF
      END
