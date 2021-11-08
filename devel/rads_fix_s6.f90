!-----------------------------------------------------------------------
! Copyright (c) 2011-2021  Remko Scharroo
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
! processed by rads_gen_s6.
!
! usage: rads_fix_s6 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_s6

use rads
use rads_devel
use rads_devel_netcdf
use rads_misc
use rads_grid

! Other local variables

character(len=rads_cmdl) :: aux_wind = '', aux_ssbk = '', aux_ssbc = '', aux_rain = '', aux_cal1 = ''
integer(fourbyteint) :: i, cyc, pass, ios
type(grid) :: info_wind, info_ssbk, info_ssbc
logical :: lcal1 = .false., lsideB_range = .false., lsideB_sig0 = .false., lrange = .false., lsig0 = .false., &
	lwind = .false., lssb = .false., lrain = .false.
integer, parameter :: sig0_nx = 500
real(eightbytereal) :: exp_ku_sigma0(sig0_nx), rms_exp_ku_sigma0(sig0_nx)
real(eightbytereal) :: bias_range(2) = 0d0, bias_sig0(2) = 0d0
real(eightbytereal), parameter :: sig0_dx = 0.1d0, gate_width = 0.3795d0, sign_error = 2 * 0.528d0
! For loading CAL1 file
integer :: ncal
real(eightbytereal) :: cal1_interval(2) = (/ -1d0, 1d0 /), &
	drange_sideB(2) = (/ -2d-3, 0d-3 /), dsig0_sideB(2) = (/ -0.07d0, 0.49d0 /)
real(eightbytereal), allocatable :: cal1_time(:), cal1(:,:), cal1_flags(:)
logical, allocatable :: cal1_mask(:)

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' cal1:: sideB_range:: sideB_sig0:: range sig0 wind ssb rain all bias_range: bias_sig0:')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('cal1')
		read (rads_opt(i)%arg, *, iostat=ios) cal1_interval
		lcal1 = .true.
	case ('sideB_range')
		read (rads_opt(i)%arg, *, iostat=ios) drange_sideB
		lsideB_range = .true.
	case ('sideB_sig0')
		read (rads_opt(i)%arg, *, iostat=ios) dsig0_sideB
		lsideB_sig0 = .true.
	case ('range')
		lrange = .true.
	case ('sig0')
		lsig0 = .true.
	case ('wind')
		lwind = .true.
	case ('ssb')
		lssb = .true.
	case ('rain')
		lrain = .true.
	case ('all')
		lcal1 = .true.
		lsideB_range = .true.
		lsideB_sig0 = .true.
		lrange = .true.
		lsig0 = .true.
		lwind = .true.
		lssb = .true.
		lrain = .true.
	case ('bias_range')
		read (rads_opt(i)%arg, *, iostat=ios) bias_range
	case ('bias_sig0')
		read (rads_opt(i)%arg, *, iostat=ios) bias_sig0
	end select
enddo
cal1_interval = cal1_interval * 86400d0

! If nothing selected, stop here

if (.not.(lcal1 .or. lrange .or. lsig0 .or. lwind .or. lssb .or. lrain)) stop

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
'  --cal1[=T0,T1]            Replace CAL1 range and power by values from LTM file, averaged over time' / &
'                            interval T0,T1 days around measurement time (default: -1,1) (L2 CONF < 011)' / &
'  --sideB_range[=KU,C]      Add biases to range (Ku, C, in m) for Side B and CHDR < 005 (def: -0.002,0.000)' / &
'  --sideB_sig0[=KU,C]       Add biases to sig0 (Ku, C, in dB) for Side B and CHDR < 005 (def: -0.07,+0.49)' / &
'  --range                   Fix range biases known for PDAP v3.0 and v3.1' / &
'                            Also fix result of temporary error in radar data base (RMC only, 34 gates)' / &
'  --sig0                    Add -7.41 dB to HR sigma0' / &
'  --rain                    Add biases to sigma0 before calling rain model (L2 CONF < 008 or with --sig0)' / &
'  --wind                    Add biases to sigma0 before calling wind model (L2 CONF < 008 or 009 or with --sig0)' / &
'  --ssb                     Update SSB (with --wind)' / &
'  --all                     All of the above' / &
'  --bias_range=KU,C         Add additional bias to range (Ku, C, in m)' / &
'  --bias_sig0=KU,C          Add additional bias to sig0 (Ku, C, in dB)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: latency(n), range_ku(n), range_ku_mle3(n), range_c(n), &
	sig0_ku(n), sig0_ku_mle3(n), sig0_c(n), dsig0_atmos_ku(n), dsig0_atmos_c(n), dsig0_atten(n), &
	swh_ku(n), swh_ku_mle3(n), wind_speed_alt(n), wind_speed_alt_mle3(n), qual_alt_rain_ice(n), flags(n), &
	ssb_cls(n), ssb_cls_mle3(n), ssb_cls_c(n), dum(n)
real(eightbytereal) :: drange(2), dsig0(2), dwind(2), drain(2), cal1_old(4), cal1_new(4)
logical :: lr, redundant, val, do_range, do_sig0, do_wind, do_ssb, do_rain
character(len=3) :: chd_ver, cnf_ver

! Initialise

call log_pass (P)
drange = 0d0
dsig0 = 0d0
dwind = 0d0
drain = 0d0

! Determine if LR/HR, OPE/VAL, CHD and CONF versions

lr = (index(P%original, '_LR_') > 0)
val = (index(P%original, '_VAL_') > 0)
i = index(P%original, 'CHD')
chd_ver = P%original(i+5:i+7)
redundant = (P%original(i+3:i+3) == 'R')
i = index(P%original, 'CONF')
cnf_ver = P%original(i+5:i+7)

! Get latency

call rads_get_var (S, P, 'latency', latency, .true.)

! cal1: Replace the values with those from the LTM (prior to Baseline F04)

if (lcal1 .and. cnf_ver < '011') then
	if (val) then
		if (need_file ('VAL/S6A_P4_1__C1HR_AX.nc', aux_cal1)) call load_cal1
	else
		if (need_file ('OPE/S6A_P4_1__C1HR_AX.nc', aux_cal1)) call load_cal1
	endif
	call cal1_average (P%equator_time, cal1_new, redundant)
	call rads_get_var (S, P, 'cal1_range_ku', dum)
	cal1_old(1:1) = dum(1)
	call rads_get_var (S, P, 'cal1_power_ku', dum)
	cal1_old(2:2) = dum(1)
	if (lr) then
		call rads_get_var (S, P, 'cal1_range_c',  dum)
		cal1_old(3:3) = dum(1)
		call rads_get_var (S, P, 'cal1_power_c',  dum)
		cal1_old(4:4) = dum(1)
	endif
endif

! sideB: Add biases to range and sigma0 (prior to CHDR version 005)

if (redundant .and. chd_ver < '005') then
	if (lsideB_range) drange = drange + drange_sideB
	if (lsideB_sig0) dsig0 = dsig0 + dsig0_sideB
endif

! range: Replace CAL1 and/or add additional bias
! range: Determine offsets to solve known biases
! range: 2 * 0.528 m: sign error in COG offset, present in PDAP v3.0 and still present in STC/NTC
!                     through the platform file
!          -0.0435 m: error in characterisation file
!           24 gates: error in radar data base

if (lrange) then
	if (lcal1 .and. all(isan_(cal1_old))) drange = drange + (cal1_new(1:3:2) - cal1_old(1:3:2))
	drange = drange + bias_range
	if (latency(1) == rads_nrt) then
		if (chd_ver < '003') drange = drange + sign_error - 0.0435d0
	else if (latency(1) == rads_stc) then
		if (cnf_ver < '005') drange = drange + sign_error
	else if (val) then
		if (P%cycle * 1000 + P%pass < 9002) drange = drange + sign_error
	else if (lr) then
		if (P%cycle * 1000 + P%pass < 10004) drange = drange + sign_error
	else
		if (P%cycle * 1000 + P%pass < 9113) drange = drange + sign_error
	endif
	if (.not.lr .and. P%cycle == 8 .and. (P%pass >= 12 .and. P%pass <= 59)) drange = drange + 24 * gate_width
endif
do_range = any(drange /= 0d0)

! sigma0: Replace CAL1 and/or add additional bias

if (lsig0) then
	if (lcal1 .and. all(isan_(cal1_old))) dsig0 = dsig0 + (cal1_old(2:4:2) - cal1_new(2:4:2))
	dsig0 = dsig0 + bias_sig0
endif

! wind: apply biases before calling wind model
! L2 CONF < 008: did not use proper wind model and/or sig0 bias
! L2 CONF = 009: did not use proper sig0 bias in HR

if (lwind) then
	if (cnf_ver < '008' .or. cnf_ver == '009' .or. dsig0(1) /= 0) dwind = (/ 1.29d0, 1.37d0 /)
endif
do_wind = any(dwind /= 0d0)

! ssb: do when requested and wind has changed

do_ssb = (lssb .and. do_wind)

! rain: apply biases before calling rain model
! L2 CONF < 008: did not use proper sig0 bias

if (lrain .and. lr) then
	if (cnf_ver < '008' .or. any(dsig0 /= 0d0)) drain = (/ 1.23d0, 1.64d0 /)
endif
do_rain = any(drain /= 0d0)

! sigma0: alignment of HR with LR

if (lsig0) then
	if (.not.lr) then
		if (cnf_ver < '009') then
			dsig0(1) = dsig0(1) - 7.41d0
		else
			dsig0(1) = dsig0(1) - 5.67d0	! Changed by 1.74 dB since L2 CONF 009
		endif
	endif
endif
do_sig0 = any(dsig0 /= 0d0)

! If nothing to change, skip

if (.not.do_range .and. .not.do_sig0 .and. .not.do_wind .and. .not.do_ssb) then
	call log_records(0)
	return
endif

! Adjust range for offset

if (do_range) then
	call rads_get_var (S, P, 'range_ku', range_ku, .true.)
	range_ku = range_ku + drange(1)
	if (lr) then
		call rads_get_var (S, P, 'range_ku_mle3', range_ku_mle3, .true.)
		range_ku_mle3 = range_ku_mle3 + drange(1)
		call rads_get_var (S, P, 'range_c', range_c, .true.)
		range_c = range_c + drange(2)
	endif
endif

! Adjust sigma0 for offset

if (do_sig0) then
	call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	sig0_ku = sig0_ku + dsig0(1)
	if (lr) then
		call rads_get_var (S, P, 'sig0_ku_mle3', sig0_ku_mle3, .true.)
		sig0_ku_mle3 = sig0_ku_mle3 + dsig0(1)
		call rads_get_var (S, P, 'sig0_c', sig0_c, .true.)
		sig0_c = sig0_c + dsig0(2)
	endif
endif

! Compute wind speed from 2D wind model after adding biases
! Load wind model if required

if (do_wind) then
	if (lr) then
		if (need_file('AUX_WNDL_S6A_002.nc', aux_wind)) then
			if (grid_load(aux_wind,info_wind) /= 0) call rads_exit ('Error loading '//trim(aux_wind))
		endif
		call rads_get_var (S, P, 'swh_ku_mle3', swh_ku_mle3, .true.)
		if (.not.do_sig0) call rads_get_var (S, P, 'sig0_ku_mle3', sig0_ku_mle3, .true.)
		call grid_inter (info_wind, n, sig0_ku_mle3 + dwind(2), swh_ku_mle3, wind_speed_alt_mle3)
	else if (need_file('AUX_WNDH_S6A_002.nc', aux_wind)) then
		if (grid_load(aux_wind,info_wind) /= 0) call rads_exit ('Error loading '//trim(aux_wind))
	endif
	call rads_get_var (S, P, 'swh_ku', swh_ku, .true.)
	if (.not.do_sig0) call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	call grid_inter (info_wind, n, sig0_ku + dwind(1), swh_ku, wind_speed_alt)
endif

! Compute SSB from 2D model after updating wind
! Load SSB model if required

if (do_ssb) then
	if (lr) then
		if (need_file('AUX_SBLK_S6A_002.nc', aux_ssbk)) then
			if (grid_load(aux_ssbk,info_ssbk) /= 0) call rads_exit ('Error loading '//trim(aux_ssbk))
		endif
		if (need_file('AUX_SBLC_S6A_002.nc', aux_ssbc)) then
			if (grid_load(aux_ssbc,info_ssbc) /= 0) call rads_exit ('Error loading '//trim(aux_ssbc))
		endif
	else
		if (need_file('AUX_SBHK_S6A_002.nc', aux_ssbk)) then
			if (grid_load(aux_ssbk,info_ssbk) /= 0) call rads_exit ('Error loading '//trim(aux_ssbk))
		endif
	endif
	call grid_inter (info_ssbk, n, wind_speed_alt, swh_ku, ssb_cls)
	if (lr) then
		call grid_inter (info_ssbk, n, wind_speed_alt_mle3, swh_ku_mle3, ssb_cls_mle3)
		call grid_inter (info_ssbc, n, wind_speed_alt, swh_ku, ssb_cls_c)
	endif
endif

! Determine rain flag after adding biases

if (do_rain) then
	if (need_file('AUX_SIGL_S6A_001.nc', aux_rain)) call rain_table
	if (dsig0(1) == 0d0) call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	if (dsig0(2) == 0d0) call rads_get_var (S, P, 'sig0_c', sig0_c, .true.)
	call rads_get_var (S, P, 'dsig0_atmos_ku', dsig0_atmos_ku, .true.)
	call rads_get_var (S, P, 'dsig0_atmos_c', dsig0_atmos_c, .true.)
	call rads_get_var (S, P, 'qual_alt_rain_ice', qual_alt_rain_ice, .true.)
	call rads_get_var (S, P, 'flags', flags, .true.)
	call compute_rain_flag (n, sig0_ku - dsig0_atmos_ku + drain(1), sig0_c - dsig0_atmos_c + drain(2), &
		dsig0_atten, qual_alt_rain_ice, flags)
endif

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Write out all the data

if (lcal1) then
	call rads_put_var (S, P, 'cal1_range_ku', cal1_new(1))
	call rads_put_var (S, P, 'cal1_power_ku', cal1_new(2))
	if (lr) then
		call rads_put_var (S, P, 'cal1_range_c', cal1_new(3))
		call rads_put_var (S, P, 'cal1_power_c', cal1_new(4))
	endif
endif

if (do_range) then
	call rads_put_var (S, P, 'range_ku', range_ku)
	if (lr) then
		call rads_put_var (S, P, 'range_ku_mle3', range_ku_mle3)
		call rads_put_var (S, P, 'range_c', range_c)
	endif
endif

if (do_sig0) then
	call rads_put_var (S, P, 'sig0_ku', sig0_ku)
	if (lr) then
		call rads_put_var (S, P, 'sig0_ku_mle3', sig0_ku_mle3)
		call rads_put_var (S, P, 'sig0_c', sig0_c)
	endif
endif

if (do_wind) then
	call rads_put_var (S, P, 'wind_speed_alt', wind_speed_alt)
	if (lr) call rads_put_var (S, P, 'wind_speed_alt_mle3', wind_speed_alt_mle3)
endif

if (do_ssb) then
	call rads_put_var (S, P, 'ssb_cls', ssb_cls)
	if (lr) then
		call rads_put_var (S, P, 'ssb_cls_mle3', ssb_cls_mle3)
		call rads_put_var (S, P, 'ssb_cls_c', ssb_cls_c)
	endif
endif

if (do_rain) then
	call rads_put_var (S, P, 'qual_alt_rain_ice', qual_alt_rain_ice)
	call rads_put_var (S, P, 'dsig0_atten', dsig0_atten)
	call rads_put_var (S, P, 'flags', flags)
endif

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Determine if file needs to be loaded
!-----------------------------------------------------------------------

logical function need_file (filenm, pathnm)
character(len=*), intent(in) :: filenm
character(len=*), intent(inout) :: pathnm

need_file = (index(pathnm, filenm) == 0)
if (need_file) then
	if (filenm(:3) == 'AUX') then
		call parseenv ('${ALTIM}/data/models/' // trim(filenm), pathnm)
	else
		call parseenv ('${RADSROOT}/ext/6a/' // trim(filenm), pathnm)
	endif
	write (*,600) trim(filenm)
endif
600 format ('(Loading ',a,') ... ', $)
end function need_file

!-----------------------------------------------------------------------
! Load CAL1 values
!-----------------------------------------------------------------------

subroutine load_cal1
use netcdf
use rads_netcdf
integer(fourbyteint) :: ncid, ncidk, ncidc, dimid

call nfs(nf90_open(aux_cal1, nf90_nowrite, ncid))
call nfs(nf90_inq_ncid(ncid, 'ku', ncidk))
call nfs(nf90_inq_dimid(ncidk, 'time', dimid))
call nfs(nf90_inquire_dimension(ncidk, dimid, len=ncal))
allocate (cal1_time(ncal), cal1(ncal,4), cal1_flags(ncal), cal1_mask(ncal))
call get_var(ncidk, 'time 473299200 ADD', cal1_time)
call get_var(ncidk, 'p4_instrument_configuration_flags', cal1_flags)
call get_var(ncidk, 'cog_delay -149896229 MUL', cal1(:,1))
call get_var(ncidk, 'total_power', cal1(:,2))
call nfs(nf90_inq_ncid(ncid, 'c', ncidc))
call get_var(ncidc, 'cog_delay -149896229 MUL', cal1(:,3))
call get_var(ncidc, 'total_power', cal1(:,4))
call nfs(nf90_close(ncid))
end subroutine load_cal1

!-----------------------------------------------------------------------
! Average CAL1 values from the CAL1 LTM file and round to 0.1 mm and 0.01 dB, resp.
!-----------------------------------------------------------------------

subroutine cal1_average (time, cal1_avg, redundant)
real(eightbytereal), intent(in) :: time
real(eightbytereal), intent(out) :: cal1_avg(:)
logical, intent(in) :: redundant
integer(fourbyteint) :: navg, flags
real(eightbytereal), parameter :: cal1_scale(4) = (/ 1d-4, 1d-2, 1d-4, 1d-2 /)

! Take only the data from the appropriate side of the altimeter and within the required interval
if (redundant) then
	flags = 33
else
	flags = 1
endif
cal1_mask = (cal1_time >= time + cal1_interval(1) .and. cal1_time <= time + cal1_interval(2) .and. cal1_flags == flags)
navg = count(cal1_mask)

! If there are no CAL1 values matching the criterion, take the closest in time
if (navg == 0) then
	i = minloc (abs(cal1_time(:) - time), 1, (cal1_flags(:) == flags))
	cal1_avg = cal1(i,:)
	return
endif

! Otherwise, average and round
do i = 1,4
	cal1_avg(i) = nint(sum(cal1(:,i), cal1_mask) / navg / cal1_scale(i)) * cal1_scale(i)
enddo
end subroutine cal1_average

! Interpolate grid

subroutine grid_inter (info, n, xval, yval, zval)
type(grid) :: info
integer(fourbyteint), intent(in) :: n
real(eightbytereal), intent(in) :: xval(n), yval(n)
real(eightbytereal), intent(out) :: zval(n)
integer(fourbyteint) :: i
real(eightbytereal) :: x, y

do i = 1,n
	x = xval(i)
	if (x < info%xmin) x = info%xmin
	if (x > info%xmax) x = info%xmax
	y = yval(i)
	if (y < info%ymin) y = info%ymin
	if (y > info%ymax) y = info%ymax
	zval(i) = grid_lininter (info, x, y)
enddo
end subroutine grid_inter

! Load the rain tables

subroutine rain_table
use rads_netcdf
use netcdf
integer :: ncid, varid
call nfs(nf90_open(aux_rain, nf90_nowrite, ncid))
call nfs(nf90_inq_varid(ncid, 'exp_ku_sigma0', varid))
call nfs(nf90_get_var(ncid, varid, exp_ku_sigma0))
call nfs(nf90_inq_varid(ncid, 'RMS_exp_ku_sigma0', varid))
call nfs(nf90_get_var(ncid, varid, rms_exp_ku_sigma0))
call nfs(nf90_close(ncid))
end subroutine rain_table

! Compute the rain flag

subroutine compute_rain_flag (n, sig0_ku, sig0_c, rain_attenuation, rain_flag, flags)
integer(fourbyteint), intent(in) :: n
real(eightbytereal), intent(in) :: sig0_ku(:), sig0_c(:)
real(eightbytereal), intent(out) :: rain_attenuation(:), rain_flag(:)
real(eightbytereal), intent(inout) :: flags(:)
real(eightbytereal) :: surface_class(n), climato_use_flag(n), rad_cloud_liquid_water(n)
integer, parameter :: no_rain = 0, rain = 1, high_rain_probability_from_altimeter = 2, &
	high_probability_of_no_rain_from_altimeter = 3, ambiguous_situation_possibility_of_ice = 4, &
	evaluation_not_possible = 5
integer, parameter :: open_ocean = 0
real(eightbytereal), parameter :: lat_thres = 50d0, rain_cloud_liquid = 0.2d0, rain_flag_coef = 1.8d0, &
	delta_sigma0_diff_threshold = 0.5d0
integer :: i, j
real(eightbytereal) :: rain_threshold, x

call rads_get_var (S, P, 'surface_class', surface_class, .true.)
call rads_get_var (S, P, 'qual_dsig0_atmos', climato_use_flag, .true.)
call rads_get_var (S, P, 'liquid_water_rad', rad_cloud_liquid_water, .true.)

rain_attenuation = nan
rain_flag = evaluation_not_possible
do i = 1, n
	if (.not.(sig0_ku(i) > 0d0 .and. sig0_c(i) > 0d0)) cycle
	if (.not.(surface_class(i) == open_ocean)) cycle !.or. surface_class(i) == continental_water)) cycle

	x = sig0_c(i) / sig0_dx
	j = max(1, min(int(x), sig0_nx - 1))
	x = x - j
	rain_attenuation(i) = (exp_ku_sigma0(j) * (1-x) + exp_ku_sigma0(j+1) * x) - sig0_ku(i)
	rain_threshold = min(delta_sigma0_diff_threshold, rain_flag_coef * (rms_exp_ku_sigma0(j) * (1-x) + rms_exp_ku_sigma0(j+1) * x))

	if (climato_use_flag(i) == 1d0) then
		if (rain_attenuation(i) > rain_threshold) then
			rain_flag(i) = high_rain_probability_from_altimeter
		else
			rain_flag(i) = high_probability_of_no_rain_from_altimeter
		endif
	else if (rain_attenuation(i) > rain_threshold .and. rad_cloud_liquid_water(i) > rain_cloud_liquid) then
		if (abs(P%tll(i,2)) > lat_thres) then
			rain_flag(i) = ambiguous_situation_possibility_of_ice
		else
			rain_flag(i) = rain
		endif
	else
		rain_flag(i) = no_rain
	endif
enddo
do i = 1,n
	j = nint(flags(i))
	if (rain_flag(i) == rain .or. rain_flag(i) == high_probability_of_no_rain_from_altimeter .or. &
		rain_flag(i) == ambiguous_situation_possibility_of_ice) then
		j = ibset(j,7)
	else
		j = ibclr(j,7)
	endif
	flags(i) = j
enddo
end subroutine compute_rain_flag

end program rads_fix_s6
