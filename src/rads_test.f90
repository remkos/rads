!-----------------------------------------------------------------------
! Copyright (c) 2011-2026  Remko Scharroo
! See LICENSE.TXT file for copying and redistribution conditions.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!-----------------------------------------------------------------------

program rads_test

use typesizes
use rads
use rads_misc

type(rads_sat) :: S
type(rads_pass) :: P
real(eightbytereal), allocatable :: data(:)
integer(fourbyteint) :: cycle,pass
real :: cpu(2)
integer :: clock(2),clockrate,i

call system_clock(clock(1),clockrate,i)
call cpu_time(cpu(1))

write (*,600) 'Start rads_init'
call rads_init (S)
write (*,600) 'debug =',rads_verbose
do cycle = S%cycles(1),S%cycles(2)
	do pass = S%passes(1),S%passes(2)
		call rads_open_pass (S, P, cycle, pass)
		if (P%ndata > 0) then
			write (*,600) 'cycle, pass, ndata = ',cycle, pass, P%ndata
			allocate (data(P%ndata))
			write (*,600) 'Start rads_get_var'
			do i = 1,S%nsel
				call rads_get_var (S, P, S%sel(i), data)
				if (S%sel(i)%info%boz_format) call bit_transfer (data)
				write (*,'(2a,2(1x,'//S%sel(i)%info%format//'))' ) S%sel(i)%name,'data(1), data(ndata) = ',data(1),data(P%ndata)
			enddo
			deallocate (data)
		endif
		call rads_close_pass (S, P)
	enddo
enddo
write (*,600) 'Start rads_end'
call rads_stat (S)
call rads_end (S)
write (*,600) 'End of program'

call cpu_time(cpu(2))
call system_clock(clock(2),clockrate,i)
write (*,'("CLOCK, CPU =",2f9.4)') modulo(clock(2)-clock(1),i)/dble(clockrate), cpu(2)-cpu(1)

600 format (a,3i6)

end program rads_test
