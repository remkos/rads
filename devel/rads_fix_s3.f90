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
use rads_devel
use meteo_subs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

real(eightbytereal), parameter :: dsig0_ku = -30.0d0, dsig0_c = -5.5d0	! Ku- and C-band Sigma0 bias
real(eightbytereal), parameter :: dtb_238 = 2.0d0, dtb_365 = 3.0d0	! Rough values from MTR presentation by M. Frery.
integer(fourbyteint) :: i, cyc, pass
logical :: lsig0 = .false., lwind = .false., ltb = .false., lmwr = .false.

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' sig0 wind tb mwr all')
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
if (rads_version ('Patch Jason-3 data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --sig0                    Adjust backscatter coefficient for apparent biases' / &
'  --wind                    Update wind speed using Envisat model' / &
'  --tb                      Adjust brightness temperatures for apparent biases' / &
'  --mwr                     Update radiometer wet parameters' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: sig0_ku(n), sig0_c(n), atten_ku(n), wind(n), tb_238(n), tb_365(n), wet(n)
integer(fourbyteint) :: i

call log_pass (P)

! Adjust backscatter for bias

call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
call rads_get_var (S, P, 'sig0_c', sig0_c, .true.)
if (lsig0) then
	sig0_ku = sig0_ku + dsig0_ku
	sig0_c  = sig0_c  + dsig0_c
endif

! Adjust brightness temperatures for bias

call rads_get_var (S, P, 'tb_238', tb_238, .true.)
call rads_get_var (S, P, 'tb_365', tb_365, .true.)
if (ltb) then
	tb_238 = tb_238 + dtb_238
	tb_365 = tb_365 + dtb_365
endif

! Adjust radiometer parameters using Envisat NN model

if (lmwr) then
	call rads_get_var (S, P, 'dsig0_atmos_ku', atten_ku, .true.)
	sig0_ku = sig0_ku - atten_ku	! Remove applied attenuation first
	do i = 1,n
		atten_ku(i) = nn_l2_mwr (tb_238(i), tb_365(i), sig0_ku(i), 1)
		wet(i)      = nn_l2_mwr (tb_238(i), tb_365(i), sig0_ku(i), 3)
	enddo
	sig0_ku = sig0_ku + atten_ku	! Add the recomputed attenuation back
endif

! Adjust wind speed

if (lwind) wind = wind_ecmwf (sig0_ku)

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Write out all the data

if (lsig0) then
	call rads_put_var (S, P, 'sig0_ku', sig0_ku)
	call rads_put_var (S, P, 'sig0_c' , sig0_c)
endif
if (ltb) then
	call rads_put_var (S, P, 'tb_238', tb_238)
	call rads_put_var (S, P, 'tb_365', tb_365)
endif
if (lwind) call rads_put_var (S, P, 'wind_speed_alt', wind)
if (lmwr) call rads_put_var (S, P, 'wet_tropo_rad', wet)

call log_records (n)
end subroutine process_pass

end program rads_fix_s3
