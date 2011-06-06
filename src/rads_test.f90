program rads_test

use typesizes
use rads

type(rads_sat) :: S
type(rads_pass) :: P
real(eightbytereal), allocatable :: data(:)
integer(fourbyteint) :: cycle,pass
real :: cpu(2)
integer :: clock(2),clockrate,i

call system_clock(clock(1),clockrate,i)
call cpu_time(cpu(1))

write (*,*) 'Start rads_init'
call rads_init (S)
write (*,*) 'debug =',S%debug
do cycle = S%cycles(1),S%cycles(2)
	do pass = S%passes(1),S%passes(2)
		call rads_open_pass (S, P, cycle, pass)
		if (P%ndata > 0) then
			write (*,*) 'cycle, pass, ndata = ',cycle, pass, P%ndata
			allocate (data(P%ndata))
			write (*,*) 'Start rads_get_var'
			do i = 1,S%nsel
				call rads_get_var (S, P, S%sel(i), data)
				write (*,'(a,2(1x,'//S%sel(i)%info%format//'))' ) S%sel(i)%name//'data(1), data(ndata) = ',data(1),data(P%ndata)
			enddo
			deallocate (data)
		endif
		call rads_close_pass (S, P)
	enddo
enddo
write (*,*) 'Start rads_end'
call rads_stat (S)
call rads_end (S)
write (*,*) 'End of program'

call cpu_time(cpu(2))
call system_clock(clock(2),clockrate,i)
write (*,'("CLOCK, CPU =",2f9.4)') modulo(clock(2)-clock(1),i)/dble(clockrate), cpu(2)-cpu(1)

end program rads_test
