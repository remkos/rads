**BUBBLE -- Make a bubble sort of an integer array
*+
      SUBROUTINE BUBBLE (INDEX, VALUE, NR)
      INTEGER*4 NR, INDEX(NR), VALUE(*)
*
* This routine bubble-sorts an integer INDEX according to the corresponding
* integer value in array VALUE. The array VALUE may be larger than NR,
* while array INDEX should contain NR values.
*
* Arguments:
*  INDEX (input) : One-dimensional array of index numbers
*       (output) : Row of index numbers sorted on corresponding value
*  VALUE (input) : One-dimensional array of integer values
*  NR    (input) : Number of indices/values
*
* Example:
* Assume you have 6 values 100, 10, 11, 21, 17, and 90, to be sorted in
* the right order, then your input can be like this:
*
*   INDEX    1    2    3    4    5    6
*   VALUE  100   10   11   21   17   90
*      NR    6
*
* After running BUBBLE, the output will look like this:
*
*   INDEX    2    3    5    4    6    1
*   VALUE  100   10   11   21   17   90
*      NR    6
*
* Note that the array VALUE hasn't changed, but that INDEX is sorted by
* the order in which VALUE should be sorted.
*-
* 24-Jan-1990 - Created
*  2-Jul-1993 - New manual
* 11-Jan-1994 - Faster version (check changes)
*  4-Aug-1998 - New manual
*-----------------------------------------------------------------------
      integer*4 i,j,k1,k2
      logical changed
      do i=2,nr
	 changed=.false.
         do j=nr,i,-1
            k2=index(j)
            k1=index(j-1)
            if (value(k2).lt.value(k1)) then
               index(j)=k1
               index(j-1)=k2
	       changed=.true.
            endif
	 enddo
	 if (.not.changed) return
      enddo
      end
