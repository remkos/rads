!*iqsort -- Make a quick sort of an INTEGER*4 array
!+
subroutine iqsort (index, value, nr)
use typesizes
integer(fourbyteint), intent(inout) :: index(nr)
integer(fourbyteint), intent(in) :: value(*), nr
!
! This routine quick-sorts an INTEGER*4 array <index> according to the
! corresponding value in array <value> of type INTEGER*4.
! The array <index> contains <nr> pointers to the values in array <value>,
! which does not have to be equally large.
!
! Arguments:
!  index (input) : One-dimensional array of index numbers (INTEGER*4)
!       (output) : Row of index numbers sorted on corresponding value
!  value (input) : One-dimensional array of values (INTEGER*4)
!  nr    (input) : Number of indices/values
!
! Example 1:
! Assume you have 6 values 100, 10, 11, 21, 17, and 90, stored in array
! <value>. These values have to be sorted in ascending order. Your input
! will be like this:
!
!   index    1    2    3    4    5    6
!   value  100   10   11   21   17   90
!      nr    6
!
! After running iqsort, the output will look like this:
!
!   index    2    3    5    4    6    1
!   value  100   10   11   21   17   90
!      nr    6
!
! Note that the array <value> has not been changed, but that INDEX is sorted by
! the order in which <value> should be sorted. You may print the values in
! the ascending order as follows:
!
!     do i=1,nr
!        write (*,*) value(index(i))
!     enddo
!
! Example 2:
! It is also possible that the indices are not in logical order. As long as
! they point to the respective locations of the values in array <value>, the
! sorting will be done accordingly. Assume you have the following input:
! (values ** are irrelevant)
!
!   index    1    2    4    5    6    8
!   value  100   10   **   21   17   90   **   1
!      nr    6
!
! Then after running iqsort, the result will be:
!
!   index    8    2    5    4    6    1
!   value  100   10   **   21   17   90   **   1
!      nr    6
!
! Printing the values in ascending order after this goes the same as in
! Example 1.
!-
! 30-Mar-1994 - Created by Marc Naeije
! 24-Sep-2010 - Adopted to Fortran 90 by Remko Scharroo
!-----------------------------------------------------------------------
integer(fourbyteint), parameter :: m = 7, nstack = 2000
integer(fourbyteint) :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
integer(fourbyteint) :: a
jstack = 0
l = 1
ir = nr
do
	if (ir-l < m) then
		jloop: do j = l+1,ir
			indxt = index(j)
			a = value(indxt)
			do i = j-1,1,-1
				if (value(index(i)) < a) then
					index(i+1) = indxt
					cycle jloop
				endif
				index(i+1) = index(i)
			enddo
			i = 0
			index(i+1)=indxt
		enddo jloop
		if (jstack == 0) return
		ir = istack(jstack)
		l = istack(jstack-1)
		jstack = jstack-2
	else
		k = (l+ir) / 2
		itemp = index(k)
		index(k) = index(l+1)
		index(l+1) = itemp
		if (value(index(l+1)) > value(index(ir))) then
			itemp = index(l+1)
			index(l+1) = index(ir)
			index(ir) = itemp
		endif
		if (value(index(l)) > value(index(ir))) then
			itemp = index(l)
			index(l) = index(ir)
			index(ir) = itemp
		endif
		if (value(index(l+1)) > value(index(l))) then
			itemp = index(l+1)
			index(l+1) = index(l)
			index(l) = itemp
		endif
		i = l+1
		j = ir
		indxt = index(l)
		a = value(indxt)
		do
			do
				i = i+1
				if (value(index(i)) >= a) exit
			enddo
			do
				j = j-1
				if (value(index(j)) <= a) exit
			enddo
			if (j < i) exit
			itemp=index(i)
			index(i)=index(j)
			index(j)=itemp
		enddo
		index(l) = index(j)
		index(j) = indxt
		jstack = jstack+2
		if (jstack > nstack) stop 'NSTACK too small in QSORT'
		if (ir-i+1 >= j-1) then
			istack(jstack) = ir
			istack(jstack-1) = i
			ir = j-1
		else
			istack(jstack) = j-1
			istack(jstack-1) = l
			l = i
		endif
	endif
enddo
end subroutine iqsort
