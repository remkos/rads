program varindex_test

use typesizes
use rads

type(rads_sat) :: S
integer(fourbyteint) :: i, j

write (*,*) 'Start rads_init'
S%debug = 3
call rads_init (S, 'e2a')
do i = 1,10000
	j = varindex (S, 'sla')
	write (*,*) j
enddo

end program varindex_test
