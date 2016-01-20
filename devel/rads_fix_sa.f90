!-----------------------------------------------------------------------
! Copyright (c) 2011-2016  Remko Scharroo
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

!*rads_fix_sa -- Patch RADS altimeter files of SARAL for various anomalies
!
! This program makes numerous patches to the SARAL RADS data processed
! by rads_gen_saral. These patches include:
!
! wet:
! - Shift MWR wet prior to 2013-10-22: subtract 6.4 mm
!
! usage: rads_fix_sa [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_sa

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
! Since all files are now "L2 Library V5" (except cycle 0), we only do this
! for cycle 0

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
if (rads_version ('Patch SARAL data for several anomalies', flag=flag)) return
call synopsis_devel ('')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: wet(n)

call log_pass (P)

! Shift MWR wet tropo prior to 2013-10-22 12:15:27 (pre-patch2 data only)

call rads_get_var (S, P, 'wet_tropo_rad', wet, .true.)
wet = wet - 6.4d-3

! Write out all the data

call rads_put_history (S, P)
call rads_put_var (S, P, 'wet_tropo_rad', wet)

call log_records (n)
end subroutine process_pass

end program rads_fix_sa
