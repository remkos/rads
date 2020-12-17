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

!*rads_fix_s6 -- Patch RADS altimeter files of Sentinel-6 for various anomalies
!
! This program makes numerous patches to the RADS data for Sentinel-6
! processed by rads_gen_s6. These patches include:
!
! range:
! - Add offset to altimeter range
!
! sig0:
! - Add offset to backscatter coefficient
!
! wind:
! - Recompute wind speed from adjusted sigma0 based on Collard model
!
! usage: rads_fix_s6 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_s6

use rads
use rads_devel
use rads_misc
use meteo_subs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

real(eightbytereal) :: drange(2) = nan, dsig0(2) = nan
integer(fourbyteint) :: i, ios, cyc, pass
logical :: lrange = .false., lsig0 = .false., lwind = .false.

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' range: sig0: wind')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('range')
		read (rads_opt(i)%arg, *, iostat=ios) drange
		if (isnan_(drange(2))) drange(2) = drange(1)
		lrange = .true.
	case ('sig0')
		read (rads_opt(i)%arg, *, iostat=ios) dsig0
		if (isnan_(dsig0(2))) dsig0(2) = dsig0(1)
		lsig0 = .true.
	case ('wind')
		lwind = .true.
	end select
enddo

! If nothing selected, stop here

if (.not.(lrange .or. lsig0 .or. lwind)) stop

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
if (rads_version ('Patch Sentinel-6 data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --range=OFFSET            Add OFFSET (m) to altimeter range' / &
'  --sig0=OFFSET             Add OFFSET (dB) to backscatter coefficient' / &
'  --wind                    Recompute wind speed from adjusted sigma0 based on Collard model')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: range_ku(n), range_ku_mle3(n), range_c(n), &
	sig0_ku(n), sig0_ku_mle3(n), sig0_c(n), swh_ku(n), swh_ku_mle3(n), &
	wind_speed_alt(n), wind_speed_alt_mle3(n)
logical :: lr
integer(fourbyteint) :: i

call log_pass (P)
lr = (index(P%original, '_LR_') > 0)

! ** Both LR and HR **

! Adjust range for offset

if (lrange) then
	call rads_get_var (S, P, 'range_ku', range_ku, .true.)
	range_ku = range_ku + drange(1)
endif

! Adjust sigma0 for offset

if (lsig0) then
	call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	sig0_ku = sig0_ku + dsig0(1)
endif

! Compute wind speed from Collard model (using unadjusted sig0_ku)

if (lwind) then
	call rads_get_var (S, P, 'swh_ku', swh_ku, .true.)
	if (.not.lsig0) call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	do i = 1,n
		wind_speed_alt(i) = wind_j1 (1, sig0_ku(i), swh_ku(i))
	enddo
endif

! ** LR only **

if (lr) then

! Adjust range for offset

	if (lrange) then
		call rads_get_var (S, P, 'range_ku_mle3', range_ku_mle3, .true.)
		range_ku_mle3 = range_ku_mle3 + drange(1)
		call rads_get_var (S, P, 'range_c', range_c, .true.)
		range_c = range_c + drange(2)
	endif

! Adjust sigma0 for offset

	if (lsig0) then
		call rads_get_var (S, P, 'sig0_ku_mle3', sig0_ku_mle3, .true.)
		sig0_ku_mle3 = sig0_ku_mle3 + dsig0(1)
		call rads_get_var (S, P, 'sig0_c', sig0_c, .true.)
		sig0_c = sig0_c + dsig0(2)
	endif

! Compute wind speed from Collard model (using unadjusted sig0_ku)

	if (lwind) then
		call rads_get_var (S, P, 'swh_ku_mle3', swh_ku_mle3, .true.)
		if (.not.lsig0) call rads_get_var (S, P, 'sig0_ku_mle3', sig0_ku_mle3, .true.)
		do i = 1,n
			wind_speed_alt_mle3(i) = wind_j1 (1, sig0_ku_mle3(i), swh_ku_mle3(i))
		enddo
	endif
endif

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Write out all the data
! ** HR and LR **

if (lrange)	call rads_put_var (S, P, 'range_ku', range_ku)
if (lsig0) call rads_put_var (S, P, 'sig0_ku', sig0_ku)
if (lwind) call rads_put_var (S, P, 'wind_speed_alt', wind_speed_alt)

! ** LR only **

if (lr) then
	if (lrange) then
		call rads_put_var (S, P, 'range_ku_mle3', range_ku_mle3)
		call rads_put_var (S, P, 'range_c', range_c)
	endif
	if (lsig0) then
		call rads_put_var (S, P, 'sig0_ku_mle3', sig0_ku_mle3)
		call rads_put_var (S, P, 'sig0_c', sig0_c)
	endif
	if (lwind) call rads_put_var (S, P, 'wind_speed_alt_mle3', wind_speed_alt_mle3)
endif

call log_records (n)
end subroutine process_pass

end program rads_fix_s6
