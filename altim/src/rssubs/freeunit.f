**FREEUNIT -- Scan for free unit number for opening file
*+
      FUNCTION FREEUNIT()
      INTEGER*4 FREEUNIT
*
* This function returns a unit number that is not currently connected to
* any file. A value of 0 is returned if no free unit is found.
*
* Example:
*     INTEGER*4 UNIT, FREEUNIT
*     UNIT=FREEUNIT()
*     OPEN (UNIT=UNIT, STATUS='OLD')
*     ...
*-
* 22-May-1995 -- Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer*4 unit
      logical opened
      do unit=99,7,-1
	 inquire (unit=unit,opened=opened)
	 if (.not.opened) then
	    freeunit=unit
	    return
	 endif
      enddo
      freeunit=0
      end
