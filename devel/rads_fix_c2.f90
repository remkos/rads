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

!*rads_fix_c2 -- Patch RADS altimeter files of CryoSat for various anomalies
!
! This program makes numerous patches to the CryoSat RADS data processed
! by rads_gen_c2_op. These patches include:
!
! tbias-lrm:
! - Adjust LRM range and ssha for a pseudo time tag bias of +0.1 ms (CRYO-COP-16),
!   with value modified from 0.394 ms to 0.1 ms based on work by Marc Naeije.
!
! tbias-plrm:
! - Adjust PLRM range and ssha for a pseudo time tag bias of +1.8 ms (CRYO-COP-72/73),
!   with value confirmed by own analysis.
!
! Documentation:
! [1] CryoSat Ocean Data Quality Status Summary, ESA IDEAS+-VEG-OQC-MEM-3159,
!     Issue 4, 1 Oct 2024
!
! usage: rads_fix_c2 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_c2

use rads
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

real(eightbytereal), parameter :: tbias_lrm = 0.1d-3, tbias_plrm = 1.8d-3
integer(fourbyteint) :: i, cyc, pass
logical :: ltbias_lrm = .false., ltbias_plrm = .false.

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' tbias-lrm tbias-plrm all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('tbias-lrm')
		ltbias_lrm = .true.
	case ('tbias-plrm')
		ltbias_plrm = .true.
	case ('all')
		ltbias_lrm = .true.
		ltbias_plrm = .true.
	end select
enddo

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
if (rads_version ('Patch CryoSat-2 data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --tbias-lrm               Correct LRM range and ssh for pseudo time tag bias' / &
'  --tbias-plrm              Correct PLRM range and ssh for pseudo time tag bias' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: flag_alt_oper_mode(n), alt_rate(n), &
	ssha(n), ssha_plrm(n), range_ku(n), range_ku_plrm(n)
integer(fourbyteint) :: i
character(len=4) :: software_version

call log_pass (P)

! If nothing requested, stop here

if (.not.(ltbias_lrm .or. ltbias_plrm)) then
	call log_records (0)
	return
endif

! Determine baseline version; skip if not 'D'

i = index(P%original,'/')
software_version = P%original(i+1:i+4)
if (software_version <= '4.00') then
	call log_records(0)
	return
endif

! We always need alt_rate

call rads_get_var (S, P, 'alt_rate', alt_rate, .true.)

! Fix LRM time tag bias (CRYO-COP-16)

if (ltbias_lrm) then
	call rads_get_var (S, P, 'range_ku', range_ku, .true.)
	call rads_get_var (S, P, 'ssha', ssha, .true.)
	call rads_get_var (S, P, 'flag_alt_oper_mode', flag_alt_oper_mode, .true.)
	where (flag_alt_oper_mode == 0)
		range_ku = range_ku + tbias_lrm * alt_rate
		ssha = ssha - tbias_lrm * alt_rate
	end where
endif

! Fix PLRM time tag bias (CRYO-COP-72/73)

if (ltbias_plrm) then
	call rads_get_var (S, P, 'range_ku_plrm', range_ku_plrm, .true.)
	call rads_get_var (S, P, 'ssha_plrm', ssha_plrm, .true.)
	range_ku_plrm = range_ku_plrm + tbias_plrm * alt_rate
	ssha_plrm = ssha_plrm - tbias_plrm * alt_rate
endif

! Update history

call rads_put_history (S, P)

! Redefine the variables

if (ltbias_lrm) then
	call rads_put_var (S, P, 'range_ku', range_ku)
	call rads_put_var (S, P, 'ssha', ssha)
endif
if (ltbias_plrm) then
	call rads_put_var (S, P, 'range_ku_plrm', range_ku_plrm)
	call rads_put_var (S, P, 'ssha_plrm', ssha_plrm)
endif

call log_records (n)
end subroutine process_pass

end program rads_fix_c2
