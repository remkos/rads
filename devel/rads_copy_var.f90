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

!*rads_copy_var -- Copy variables from one set of RADS data to another
!+
! This program gets the values of variables from one flavour of RADS
! data and puts them onto files of the other flavour.
! The passes in the two flavours do not need to be equal length;
! the closest matching values within a second are copied.
!
! usage: rads_copy_var [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_copy_var

use rads
use rads_misc
use rads_devel

! Data variables

type(rads_sat) :: S(2)
type(rads_pass) :: P(2)

! Command line arguments

integer(fourbyteint) :: cyc, pass, j

! Initialise

call synopsis ('--head')
call rads_set_options (' all')
call rads_init (S)

! Check options

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('all')
!	dummy
	end select
enddo

! Read the flags.dat file and execute data flagging if required

do cyc = S(1)%cycles(1), S(1)%cycles(2), S(1)%cycles(3)
	do pass = S(1)%passes(1), S(1)%passes(2), S(1)%passes(3)
		call rads_open_pass (S(1), P(1), cyc, pass)
		if (P(1)%ndata > 0) then
			call rads_open_pass (S(2), P(2), cyc, pass, .true.)
			if (P(2)%ndata > 0) call process_pass (P(1)%ndata, P(2)%ndata)
			call rads_close_pass (S(2), P(2))
		endif
		call rads_close_pass (S(1), P(1))
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Copy variable between RADS data flavours', flag=flag)) return
call synopsis_devel ('')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  --all                     (dummy, no impact)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n1, n2)
integer(fourbyteint), intent(in) :: n1, n2
real(eightbytereal) :: var(0:n1)
integer :: i, idx(n2)

call log_pass (P(2))

! Store all data fields.

call rads_put_history (S(2), P(2))
do i = 1, S(1)%nsel
	call rads_def_var (S(2), P(2), S(1)%sel(i)%name)
enddo

! If all times match, then do a simple copy

if (n1 == n2 .and. all(P(1)%tll(:,1) == P(2)%tll(:,1))) then
	do i = 1, S(1)%nsel
		call rads_get_var (S(1), P(1), S(1)%sel(i)%name, var(1:))
		call rads_put_var (S(2), P(2), S(1)%sel(i)%name, var(1:))
	enddo
	call log_records (n2)
	return
endif

! Otherwise make an index first

idx = 0
do i = 1,n2
	j = minloc(abs(P(1)%tll(:,1)-P(2)%tll(i,1)), 1)
	if (abs(P(1)%tll(j,1)-P(2)%tll(i,1)) < 0.5d0) idx(i) = j
enddo
do i = 1, S(1)%nsel
	call rads_get_var (S(1), P(1), S(1)%sel(i)%name, var(1:))
	call rads_put_var (S(2), P(2), S(1)%sel(i)%name, var(idx))
enddo
call log_records (n2)

end subroutine process_pass

end program rads_copy_var
