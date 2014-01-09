!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2014  Remko Scharroo (Altimetrics LLC)
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

!*rads_fix_reaper -- Patch RADS altimeter files of ERS/REAPER for various anomalies
!
! This program makes numerous patches to the ERS/REAPER RADS data processed
! by rads_gen_reaper. These patches include:
!
! ptr:
! - Correct range for erroneous PTR
!
! usage: rads_fix_reaper [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_reaper

use rads
use rads_misc
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

character(len=rads_cmdl) :: path
integer(fourbyteint) :: i, cyc, pass
real(eightbytereal) :: time_ptr, range_ptr
logical :: lptr = .false.

! Scan command line for options

call synopsis
call rads_set_options (' ptr all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('ptr')
		lptr = .true.
	case ('all')
		lptr = .true.
	end select
enddo

! Load PTR data

if (lptr) then
	call parseenv ('${RADSROOT}/ext/reaper/commissioning/diff_ptrolc_ers'//S%sat(2:2)//'.dat', path)
	open (10,file=path,status='old')
endif

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

close (10)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Patch ERS/REAPER data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --ptr                     Correct range for PTR error' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n), range_ku(n)
integer :: n_changed

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename(len_trim(S%dataroot)+2:))

! Check if time range matches PTR table

do while (time_ptr + 0.5d0 < P%start_time)
	call next_ptr
enddo
if (time_ptr - 0.5d0 > P%end_time) then
	write (*,552) 0
	return
endif

! Process data records

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'range_ku', range_ku, .true.)

! Match up PTR

n_changed = 0
if (lptr) then
	do i = 1,n
		do while (time_ptr + 0.5d0 < time(i))
			call next_ptr
		enddo
		if (time_ptr - 0.5d0 > time(i)) cycle
		n_changed = n_changed + 1
		range_ku(i) = range_ku(i) - range_ptr
	enddo
endif

! If nothing changed, stop here

if (n_changed == 0) then
	write (*,552) 0
	return
endif

! Write out all the data

call rads_put_history (S, P)
if (lptr) call rads_put_var (S, P, 'range_ku', range_ku)

write (*,552) n_changed

end subroutine process_pass

!-----------------------------------------------------------------------
! Get next PTR value
!-----------------------------------------------------------------------

subroutine next_ptr
character(len=80) :: line
integer :: ios, n
real(eightbytereal) :: corr1, corr2
do
	read (10,'(a)',iostat=ios) line
	if (ios /= 0) then
		time_ptr = 1d20
		exit
	endif
	if (line(1:1) == '#') cycle
	read (line,*) time_ptr, n, corr1, corr2, range_ptr
	exit
enddo
end subroutine next_ptr

end program rads_fix_reaper
