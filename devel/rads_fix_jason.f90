!-----------------------------------------------------------------------
! Copyright (c) 2011-2025  Remko Scharroo
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

!*rads_fix_jason -- Patch RADS altimeter files of Jason for various anomalies
!
! This program makes numerous patches to the RADS data for Jason-1/2/3
! processed by rads_gen_jason, and originating from (O/I)GDR versions D
! or E. These patches include:
!
! sig0:
! - Adjust backscatter coefficient for apparent off-nadir angle
! - Add optional bias to Ku and C sigma0
!
! rad:
! - Add offset to radiometer wet tropo or read biases to the wet tropo
!   from a file
!
! wind:
! - Recompute wind speed from adjusted sigma0 based on Collard model
!
! usage: rads_fix_jason [data-selectors] [options]
!
! References:
! Collard, F., Algorithmes de vent et periode moyenne des vagues JASON
! a base de reseaux de neurons,
! Rapport BO-021-CLS-0407-RF, Boost Technologies, 2005.
!
! Quartly, G. D., Optimizing sigma0 information from Jason-2 altimeter,
! IEEE Geosci. Rem. Sens. Lett., 6(3), 398–402,
! doi:10.1109/LGRS.2009.2013973, 2009.
!-----------------------------------------------------------------------
program rads_fix_jason

use rads
use rads_misc
use rads_devel
use meteo_subs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

real(eightbytereal) :: dwet = 0d0, dsig0(2) = 0d0, wet_cor(254,0:999) = nan
integer(fourbyteint) :: i, ios, cyc, pass
logical :: lrad = .false., lsig0 = .false., lwind = .false.

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' sig0:: rad: wind all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('sig0')
		read (rads_opt(i)%arg, *, iostat=ios) dsig0
		lsig0 = .true.
	case ('rad')
		if (read_wet_cor(rads_opt(i)%arg)) then
			read (rads_opt(i)%arg, *, iostat=ios) dwet
			wet_cor = dwet * 1d-3
		endif
		lrad = .true.
	case ('wind')
		lwind = .true.
	case ('all')
		if (S%sat(:2) == 'j3') then
			call parseenv ("${RADSROOT}/ext/j3/JASON_3_PD_CORRECTION_20230925.txt", rads_opt(i)%arg)
			if (read_wet_cor(rads_opt(i)%arg)) call rads_message ('Error loading radiometer correction file')
			wet_cor(:,349:364) = 1.8d-3
			wet_cor(:,365:999) = 0.8d-3
			lrad = .true.
		else
			lsig0 = .true.
		endif
	end select
enddo

! If nothing selected, stop here

if (.not.(lsig0 .or. lwind .or. lrad)) stop

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
if (rads_version ('Patch Jason data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --sig0[=BIAS_KU,BIAS_C]   Adjust backscatter coefficient for apparent off-nadir angle, and' / &
'                            optionally add biases to the Ku and C band values' / &
'  --rad=OFFSET              Add OFFSET mm to radiometer wet tropo' / &
'  --rad=FILENAME            Correct radiometer wet tropo according to correction file' / &
'  --all                     JA1/JA2: --sig0 (without applying a bias)' / &
'                            JA3 cycle   0-348: --rad=${RADSROOT}/ext/j3/JASON_3_PD_CORRECTION_20230925.txt' / &
'                                cycle 349-364: --rad=1.8' / &
'                                cycle 365-   : --rad=0.8' / &
'  --wind                    Recompute wind speed from adjusted sigma0 based on Collard model')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: sig0_ku(n), sig0_c(n), psi2(n), swh_ku(n), u(n), wet(n)
integer(fourbyteint) :: i

call log_pass (P)

! Adjust radiometer wet tropo by 2 mm because of uncorrected drift.
! This should be only used for OGDRs or IGDRs.

if (lrad) then
	call rads_get_var (S, P, 'wet_tropo_rad', wet, .true.)
	wet = wet + wet_cor (pass, cyc)
endif

! Adjust backscatter for correlation with off-nadir angle (See Quartly [2009])

if (lsig0) then
	call rads_get_var (S, P, 'off_nadir_angle2_wf_ku', psi2, .true.)
	call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	call rads_get_var (S, P, 'sig0_c', sig0_c, .true.)
	sig0_ku = sig0_ku - 11.34d0 * psi2 + dsig0(1)
	sig0_c  = sig0_c  -  2.01d0 * psi2 + dsig0(2)
endif

! Compute wind speed from Collard model (using unadjusted sig0_ku)

if (lwind) then
	call rads_get_var (S, P, 'swh_ku', swh_ku, .true.)
	if (.not.lsig0) call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	do i = 1,n
		u(i) = wind_j1 (1, sig0_ku(i) - dsig0(1), swh_ku(i))
	enddo
endif

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Redefine the variables

if (lrad) call rads_def_var (S, P, 'wet_tropo_rad')
if (lsig0) then
	call rads_def_var (S, P, 'sig0_ku')
	call rads_def_var (S, P, 'sig0_c' )
endif
if (lwind) call rads_def_var (S, P, 'wind_speed_alt')

! Write out all the data

if (lrad) call rads_put_var (S, P, 'wet_tropo_rad', wet)
if (lsig0) then
	call rads_put_var (S, P, 'sig0_ku', sig0_ku)
	call rads_put_var (S, P, 'sig0_c' , sig0_c)
endif
if (lwind) call rads_put_var (S, P, 'wind_speed_alt', u)

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Read the radiometer wet tropo correction file
!-----------------------------------------------------------------------

logical function read_wet_cor (filenm)
character(len=*), intent(in) :: filenm
integer :: ios, pass, cyc
character(len=17) :: date
character(len=80) :: line
real(eightbytereal) :: dwet

read_wet_cor = .true.
open (10, file=filenm, form='formatted', status='old', iostat=ios)
if (ios /= 0) return
read_wet_cor = .false.
! Read the data while skipping headers
do
	read (10, '(a)', iostat = ios) line
	if (ios /= 0) exit
	if (line(:3) == 'HDR') cycle
	read (line, *) cyc, pass, date, dwet
	wet_cor (pass, cyc) = dwet * 1d-2
enddo
! Fill missing data
dwet = 0d0
do cyc = 0,348
	do pass = 1,254
		if (isnan_(wet_cor (pass,cyc))) wet_cor (pass,cyc) = dwet
		dwet = wet_cor (pass, cyc)
	enddo
enddo
wet_cor(:,349:999) = 0d0
close (10)
end function read_wet_cor

end program rads_fix_jason
