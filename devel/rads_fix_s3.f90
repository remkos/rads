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

!*rads_fix_s3 -- Patch RADS altimeter files of Sentinel-3 for various anomalies
!
! This program makes numerous patches to the Sentinel-3 RADS data processed
! by rads_gen_s3. These patches include:
!
!  --sig0   Adjust backscatter coefficient for apparent biases
!  --wind   Update wind speed using Envisat model
!  --tb     Adjust brightness temperatures for apparent biases
!  --mwr    Update radiometer wet parameters
!
! usage: rads_fix_s3 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_s3

use rads
use rads_misc
use rads_devel
use meteo_subs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

real(eightbytereal), parameter :: dsig0_ku = 0.0d0, dsig0_c = 3.3d0	! Ku- and C-band Sigma0 bias
real(eightbytereal), parameter :: dtb_238 = 2.0d0, dtb_365 = 3.0d0	! Rough values from MTR presentation by M. Frery.
real(eightbytereal), parameter :: drange = 59.3d-3 ! Range bias in earlier rep_002 data
integer(fourbyteint) :: i, cyc, pass
logical :: lsig0 = .false., lwind = .false., ltb = .false., lmwr = .false., lrange = .false.

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' sig0 wind tb mwr all range')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('sig0')
		lsig0 = .true.
	case ('wind')
		lwind = .true.
	case ('tb')
		ltb = .true.
	case ('mwr')
		lmwr = .true.
	case ('all')
		lsig0 = .true.
		lwind = .true.
		ltb = .true.
		lmwr = .true.
	case ('range')
		lrange = .true.
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
if (rads_version ('Patch Sentinel-3 data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --sig0                    Adjust backscatter coefficient for apparent biases and attenuation' / &
'  --wind                    Update wind speed using Envisat model' / &
'  --tb                      Adjust brightness temperatures for apparent biases' / &
'  --mwr                     Update radiometer wet parameters' / &
'  --all                     All of the above' / &
'  --range                   Subtract 59.3 mm from Ku- and C-band ranges')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: sig0_ku(n), sig0_c(n), atten_ku(n), atten_c(n), wind(n), tb_238(n), tb_365(n), &
	wet_tropo_ecmwf(n), wet_tropo_rad(n), range_ku(n), range_ku_plrm(n), range_c(n), ssha(n), ssha_plrm(n)
integer(fourbyteint) :: i

call log_pass (P)

call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
call rads_get_var (S, P, 'sig0_c', sig0_c, .true.)
call rads_get_var (S, P, 'dsig0_atmos_ku', atten_ku, .true.)
call rads_get_var (S, P, 'dsig0_atmos_c', atten_c, .true.)
call rads_get_var (S, P, 'tb_238', tb_238, .true.)
call rads_get_var (S, P, 'tb_365', tb_365, .true.)
call rads_get_var (S, P, 'wet_tropo_ecmwf', wet_tropo_ecmwf, .true.)

! Adjust backscatter for bias

if (lsig0) then
	sig0_ku = sig0_ku + dsig0_ku
	sig0_c  = sig0_c  + dsig0_c
endif

! Adjust brightness temperatures for bias

if (ltb) then
	tb_238 = tb_238 + dtb_238
	tb_365 = tb_365 + dtb_365
endif

! Adjust radiometer parameters using Envisat NN model
! Note that the original sig0 is NOT corrected for atmospheric attenuation

if (lmwr) then
	do i = 1,n
		atten_ku(i)      = nn_l2_mwr (tb_238(i), tb_365(i), sig0_ku(i), 1)
		wet_tropo_rad(i) = nn_l2_mwr (tb_238(i), tb_365(i), sig0_ku(i), 3)
	enddo
endif

! Fill in the defaulted attenuation values with functional form of ECMWF wet tropo correction
! (See notes of 2016-04-14 "Simple backup for sig0 attenuation")
! Add the attenuation correction to sigma0

if (lsig0) then
	where (isnan_(atten_ku)) atten_ku = 0.118d0 - 0.456d0 * wet_tropo_ecmwf
	where (isnan_(atten_c))  atten_c  = 0.09d0
	sig0_ku = sig0_ku + atten_ku
	sig0_c  = sig0_c  + atten_c
endif

! Adjust wind speed

if (lwind) wind = wind_ecmwf (sig0_ku)

! Update ranges

if (lrange) then
	call rads_get_var (S, P, 'range_ku', range_ku, .true.)
	call rads_get_var (S, P, 'range_ku_plrm', range_ku_plrm, .true.)
	call rads_get_var (S, P, 'range_c', range_c, .true.)
	call rads_get_var (S, P, 'ssha', ssha, .true.)
	call rads_get_var (S, P, 'ssha_plrm', ssha_plrm, .true.)
endif

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Write out all the data

if (lsig0) then
	call rads_put_var (S, P, 'sig0_ku', sig0_ku)
	call rads_put_var (S, P, 'sig0_c' , sig0_c)
	call rads_put_var (S, P, 'dsig0_atmos_ku', atten_ku)
	call rads_put_var (S, P, 'dsig0_atmos_c' , atten_c)
endif
if (ltb) then
	call rads_put_var (S, P, 'tb_238', tb_238)
	call rads_put_var (S, P, 'tb_365', tb_365)
endif
if (lwind) call rads_put_var (S, P, 'wind_speed_alt', wind)
if (lmwr) call rads_put_var (S, P, 'wet_tropo_rad', wet_tropo_rad)
if (lrange) then
	call rads_put_var (S, P, 'range_ku', range_ku - drange)
	call rads_put_var (S, P, 'range_ku_plrm', range_ku_plrm - drange)
	call rads_put_var (S, P, 'range_c', range_c - drange)
	call rads_put_var (S, P, 'ssha', ssha + drange)
	call rads_put_var (S, P, 'ssha_plrm', ssha_plrm + drange)
endif

call log_records (n)
end subroutine process_pass

end program rads_fix_s3
