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

!*rads_add_flags -- Add additional flags to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! based on a existing data files of passes indicating the flags to be
! changed.
!
! usage: rads_add_flags [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_flags

use rads
use rads_misc
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Command line arguments

character(rads_cmdl) :: filename, arg, cause
integer(fourbyteint) :: cyc, pass, pass0, pass1, bit, set, sel, ios
real(eightbytereal) :: val0, val1

! Formats

550 format (a)

! Initialise

call synopsis ('--head')
call rads_init (S)

! Open the data file

call parseenv ('${ALTIM}/data/tables/', filename)
filename = trim(filename) // trim(S%sat) // '_flags.dat'
open (10,file=filename, status='old', iostat=ios)
if (ios /= 0) then
	write (*,550) 'Error opening data file ' // trim(filename)
	stop
endif

! Read the flags.dat file and execute data flagging if required

do
	read (10,'(a)',iostat=ios) arg
	if (ios /= 0) exit
	if (arg(:1) == '#') cycle
	read (arg,*) bit,set,cyc,pass0,pass1,sel,val0,val1,cause
	if (cyc < S%cycles(1) .or. cyc > S%cycles(2)) cycle
	do pass = max(pass0,S%passes(1)), min(pass1,S%passes(2)), S%passes(3)
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
if (rads_version ('$Revision$', 'Add additional flags to RADS data', flag=flag)) return
call synopsis_devel ('')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  --file=name               Set the name of the flag data base to be used' / &
'                            (default is ${ALTIM}/data/tables/${SAT}_flags.dat)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: val(n), flags(n)
integer(fourbyteint) :: flag, i, changed

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename(len_trim(S%dataroot)+2:))

! Get lat, lon, flags, surface_type

call rads_get_var (S, P, 'flags', flags, .true.)
if (sel >= 0) call rads_get_var (S, P, sel, val, .true.)

! Process data records

changed = 0
do i = 1,n
	if (sel < 0 .or. (val(i) >= val0 .and. val(i) <= val1)) then
		changed = changed + 1
		flag = nint(flags(i))
		if (set == 1) then
			flag = ibset(flag, bit)
		else
			flag = ibclr(flag, bit)
		endif
		flags(i) = flag
	endif
enddo

! Store all data fields.

if (changed == 0) then
	write (*,552) 0
	return
endif

call rads_put_history (S, P)
call rads_def_var (S, P, 'flags')
call rads_put_var (S, P, 'flags', flags)
write (*,552) n
end subroutine process_pass

end program rads_add_flags
