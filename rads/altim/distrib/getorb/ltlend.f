**LTLEND -- Determined byte order of integers
*+
      FUNCTION LTLEND()
      LOGICAL  LTLEND
*
* This function returns .TRUE. if your system uses the LITTLE ENDIAN
* notation of INTEGERS. For BIG ENDIAN notation, it returns .FALSE.
*
* The terms LITTLE ENDIAN and BIG ENDIAN refer to the order of the
* bytes in an INTEGER word, both for two-byte and four-byte integers.
* Within a BYTE, the BITS are always ordered most significant to
* least significant. Within a WORD the BYTES may be ordered:
* - Most significant to least significant: BIG endian
*   (IBM RS6000, Sun, SGI, Convex, MacIntosh).
* - Least significant to most significant: LITTLE endian
*   (PC, VAX, DEC).
*
* Arguments:
*  (none)
* Returned value:
*  LTLEND : .TRUE. for LITTLE endian notation of integers
*           .FALSE. for BIG endian notation of integers
*-
*  5-Jan-1996 - Created by Remko Scharroo
*  8-Aug-2000 - Changed to check integer*2 in stead of integer*4
*-----------------------------------------------------------------------
      integer*2 itest(2)
      character*2 ctest(2)
      equivalence (ctest,itest)
      data ctest /'12','21'/
      ltlend=itest(1).gt.itest(2)
      end
