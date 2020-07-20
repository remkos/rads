!-----------------------------------------------------------------------
! Copyright (c) 2011-2020  Remko Scharroo
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

!*rads_fix_eqtime -- Patch RADS for erroneous equator time
!
! This program makes patches the Jason and SARAL RADS data for errors
! introduced by an older version of ogdrsplit (prior to 26 Aug 2015).
! It moves the equator time one day forward when it is clearly too
! early. This would only affect OGDR files.
!
! usage: rads_fix_eqtime [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_eqtime

use rads
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

integer(fourbyteint) :: cyc, pass

! Scan command line for options

call synopsis ('--head')
call rads_init (S)

! Run process for all files

do cyc = S%cycles(1),max(0,S%cycles(2)),S%cycles(3)
	! Process passes
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Patch Jason and SARAL data for equator time bug', flag=flag)) return
call synopsis_devel ('')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n


call log_pass (P)

if (index (P%original, '_OPN_') == 0) then
	! Skip if not OGDR
	call log_records (0)
else if (P%start_time - P%equator_time < 43200d0) then
	! Skip if time is OK
	call log_records (0)
else
	! Move equator time one forward
	P%equator_time = P%equator_time + 86400d0
	call rads_put_passinfo (S, P)
	call rads_put_history (S, P)
	call log_records (n)
endif
end subroutine process_pass

end program rads_fix_eqtime
