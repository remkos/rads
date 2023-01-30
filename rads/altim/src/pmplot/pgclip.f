**PGCLIP -- Select clipping mode
*+
      SUBROUTINE PGCLIP (MODE)
      INTEGER MODE
*
* Select whether symbols or line should be clipped. Clipping is switched
* on or off with the MODE parameter (see below).
*
* Argument:
* MODE (input) : = 0 : Clip all graphics.
*                      Draw symbols completely, but only if the center is
*                      inside the viewport. (DEFAULT)
*                = 1 : Clip all graphics and all symbols.
*                = 2 : Clip all graphics.
*                      Let symbols be drawn completely, also when they are
*                      outside the viewport.
*                = 3 : Do not clip any graphics nor symbols.
*--
*  5-Nov-1993 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      include 'grclip.inc'
      if (mode.lt.0.or.mode.gt.3) then
	 call grwarn('PGCLIP: invalid mode')
      else
	 grclipmode = mode
      endif
      end
