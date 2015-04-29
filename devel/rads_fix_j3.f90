!-----------------------------------------------------------------------
! Copyright (c) 2011-2015  Remko Scharroo
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

!*rads_fix_j3 -- Patch RADS altimeter files of Jason-3 for various anomalies
!
! This program makes numerous patches to the Jason-3 RADS data processed
! by rads_gen_j3. These patches include:
!
! sig0:
! - Adjust backscatter coefficient for apparent off-nadir angle
!
! wind:
! - Recompute wind speed from adjusted sigma0 based on Collard model
!
! usage: rads_fix_j3 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_j3

use rads
use rads_devel
use meteo_subs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

real(eightbytereal), parameter :: dsig0_ku = -2.40d0, dsig0_c = -0.73d0	! Ku- and C-band Sigma0 bias of Jason-1
integer(fourbyteint) :: i, cyc, pass

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('all')
	! Nothing
	end select
enddo

! Nothing to be done for the time being

stop

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
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
if (rads_version ('Patch Jason-3 data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --all                     All of the above (i.e. nothing for the moment)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
!real(eightbytereal) :: sig0_ku(n), sig0_c(n), psi2(n), swh_ku(n), u(n)
!integer(fourbyteint) :: i

call log_pass (P)

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

call log_records (n)
end subroutine process_pass

end program rads_fix_j3
