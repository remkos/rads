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

character(len=rads_cmdl) :: aux_wind = '', aux_ssbk = '', aux_ssbc = '', aux_rain = ''
integer(fourbyteint) :: i, cyc, pass, ios, lrain = 0, lwind = 0
type(grid) :: info_wind, info_ssbk, info_ssbc
logical :: lsig0 = .false., lssb = .false., nr_only = .false., lflag = .false.
integer, parameter :: sig0_nx = 500
real(eightbytereal) :: exp_ku_sigma0(sig0_nx), rms_exp_ku_sigma0(sig0_nx)
real(eightbytereal) :: bias_range(2) = 0d0, bias_sig0(2) = 0d0, &
	dwind(2) = 0.57d0, drain(2) = (/ 0.51d0, 0.72d0 /)
real(eightbytereal), parameter :: sig0_dx = 0.1d0, gate_width = 0.3795d0, sign_error = 2 * 0.528d0

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' range sig0 wind:: ssb rain:: all bias-range: bias-sig0: nr-only p2p' // &
	' flag-bit0')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('sig0')
		lsig0 = .true.
	case ('wind')
		read (rads_opt(i)%arg, *, iostat=ios) dwind
		lwind = 1
		if (ios == 0) lwind = 2
	case ('ssb')
		lssb = .true.
	case ('rain')
		read (rads_opt(i)%arg, *, iostat=ios) drain
		lrain = 1
		if (ios == 0) lrain = 2
	case ('all')
		lsig0 = .true.
	case ('bias-range')
		read (rads_opt(i)%arg, *, iostat=ios) bias_range
	case ('bias-sig0')
		read (rads_opt(i)%arg, *, iostat=ios) bias_sig0
	case ('nr-only')
		nr_only = .true.
	case ('flag-bit0')
		lflag = .true.
	case ('p2p')
		bias_sig0(1) = 10d0 * log10(4d0)	! Impact of reducing the waveform accumulation by factor 4
		lrain = 1
		lwind = 1
		lssb = .true.
	end select
enddo

! If nothing selected, stop here

if (.not.(lsig0 .or. lwind > 0 .or. lssb .or. lrain > 0 .or. lflag .or. any(bias_sig0 /= 0d0))) stop

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
'  --sig0                    Add -5.67 dB (Baseline < F09) or +0.47 dB (Baseline F09) to HR sigma0' / &
'  --all                     All of the above' / &
'  --p2p                     Counter effects of reduced 2.2 kHz LR waveform accumulation,' / &
'                            implies --rain --wind --ssb)' / &
'  --rain[=KU,C]             Add biases to sigma0 (KU, C, in dB) before calling rain model' / &
'                            (with --sig0 use default 0.51,0.72)' / &
'  --wind[=MLE4,MLE3]        Add biases to sigma0 (MLE4, MLE3, in dB) before calling wind model' / &
'                            (with --sig0 use default 0.57,0.57)' / &
'  --ssb                     Update SSB (with --wind)' / &
'  --bias-range=KU,C         Add additional bias to range (Ku, C, in m)' / &
'  --bias-sig0=KU,C          Add additional bias to sig0 (Ku, C, in dB)' / &
'  --nr-only                 Only update numerical retracker values' / &
'  --flag-bit0               Clear (VAL=0) or set (VAL=1) flag bit 0 and update attributes')
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
	flags_mle3(n), flags_nr(n), ssb_cls(n), ssb_cls_mle3(n), ssb_cls_c(n)
real(eightbytereal) :: drange(2), dsig0(2)
logical :: lr, has_nr, has_mle3, has_c, redundant, val, do_range = .false., do_sig0 = .false., do_wind = .false., &
	do_ssb = .false., do_rain = .false., do_flag = .false.
character(len=3) :: chd_ver, cnf_ver, baseline

! Initialise

call log_pass (P)
drange = bias_range
dsig0 = bias_sig0
drain = 0d0

! Determine if LR/HR, OPE/VAL, CHD and CONF versions

lr = (index(P%original, '_LR_') > 0)
val = (index(P%original, '_VAL_') > 0)
i = index(P%original, 'Baseline')
baseline = P%original(i+9:i+11)
i = index(P%original, 'CHD')
chd_ver = P%original(i+5:i+7)
redundant = (P%original(i+3:i+3) == 'R')
i = index(P%original, 'CONF')
cnf_ver = P%original(i+5:i+7)
has_nr = (lr .and. baseline > 'F07') .or. baseline > 'F08'
has_c = (lr .and. .not.nr_only)
has_mle3 = (lr .and. .not.nr_only .and. baseline < 'G01')

! Get latency

call rads_get_var (S, P, 'latency', latency, .true.)

do_range = any(drange /= 0d0)

! sigma0: alignment of HR with LR

if (lsig0 .and. .not.lr) then
	if (baseline < 'F09') then
		dsig0(1) = dsig0(1) - 5.67d0
	else if (baseline < 'G01') then
		dsig0(1) = dsig0(1) + 0.47d0	! Changed by +6.14 dB in Baseline F09, removed in G01
	endif
endif

do_sig0 = any(dsig0 /= 0d0)

! wind: apply biases before calling wind model

if (lwind > 0) then
	do_wind = (dsig0(1) /= 0 .or. lwind == 2)
endif

! ssb: do when requested and wind has changed

do_ssb = (lssb .and. do_wind)

! rain: apply biases before calling rain model

if (lrain > 0 .and. lr .and. .not.nr_only) then
	do_rain = (do_sig0 .or. lrain == 2)
endif

! flag bit 0: not needed for HR Side B

do_flag = lflag .and. (lr .or. .not.redundant)

! If nothing to change, skip

if (.not.(do_range .or. do_sig0 .or. do_wind .or. do_ssb .or. do_rain .or. do_flag)) then
	call log_records(0)
	return
endif

! Adjust range for offset

if (do_range) then
	if (nr_only) then
		call rads_get_var (S, P, 'range_ku_nr', range_ku, .true.)
	else
		call rads_get_var (S, P, 'range_ku', range_ku, .true.)
		if (has_mle3) then
			call rads_get_var (S, P, 'range_ku_mle3', range_ku_mle3, .true.)
			range_ku_mle3 = range_ku_mle3 + drange(1)
		endif
		if (has_c) then
			call rads_get_var (S, P, 'range_c', range_c, .true.)
			range_c = range_c + drange(2)
		endif
	endif
	range_ku = range_ku + drange(1)
endif

! Adjust sigma0 for offset

if (do_sig0) then
	if (nr_only) then
		call rads_get_var (S, P, 'sig0_ku_nr', sig0_ku, .true.)
	else
		call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
		if (has_mle3) then
			call rads_get_var (S, P, 'sig0_ku_mle3', sig0_ku_mle3, .true.)
			sig0_ku_mle3 = sig0_ku_mle3 + dsig0(1)
		endif
		if (has_c) then
			call rads_get_var (S, P, 'sig0_c', sig0_c, .true.)
			sig0_c = sig0_c + dsig0(2)
		endif
	endif
	sig0_ku = sig0_ku + dsig0(1)
endif

! Compute wind speed from 2D wind model after adding biases
! Load wind model if required

if (do_wind) then
	if (lr) then
		if (need_file('AUX_WNDL_S6A_002.nc', aux_wind)) then
			if (grid_load(aux_wind,info_wind) /= 0) call rads_exit ('Error loading '//trim(aux_wind))
		endif
		if (has_mle3) then
			call rads_get_var (S, P, 'swh_ku_mle3', swh_ku_mle3, .true.)
			if (.not.do_sig0) call rads_get_var (S, P, 'sig0_ku_mle3', sig0_ku_mle3, .true.)
			call grid_inter (info_wind, n, sig0_ku_mle3 + dwind(2), swh_ku_mle3, wind_speed_alt_mle3)
		endif
	else if (need_file('AUX_WNDH_S6A_002.nc', aux_wind)) then
		if (grid_load(aux_wind,info_wind) /= 0) call rads_exit ('Error loading '//trim(aux_wind))
	endif
	if (nr_only) then
		call rads_get_var (S, P, 'swh_ku_nr', swh_ku, .true.)
		if (.not.do_sig0) call rads_get_var (S, P, 'sig0_ku_nr', sig0_ku, .true.)
	else
		call rads_get_var (S, P, 'swh_ku', swh_ku, .true.)
		if (.not.do_sig0) call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	endif
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
	if (has_mle3) call grid_inter (info_ssbk, n, wind_speed_alt_mle3, swh_ku_mle3, ssb_cls_mle3)
	if (has_c) call grid_inter (info_ssbc, n, wind_speed_alt, swh_ku, ssb_cls_c)
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

! Set the flag bit 0 appropriately

if (do_flag) then
	if (.not.do_rain) call rads_get_var (S, P, 'flags', flags, .true.)
	call set_flag (n, flags, redundant)
	if (has_mle3) then
		call rads_get_var (S, P, 'flags_mle3', flags_mle3, .true.)
		call set_flag (n, flags_mle3, redundant)
	endif
	if (has_nr) then
		call rads_get_var (S, P, 'flags_nr', flags_nr, .true.)
		call set_flag (n, flags_nr, redundant)
	endif
endif

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! (Re)define the variables

if (do_range) then
	if (nr_only) then
		call rads_def_var (S, P, 'range_ku_nr')
	else
		call rads_def_var (S, P, 'range_ku')
		if (has_mle3) call rads_def_var (S, P, 'range_ku_mle3')
		if (has_c) call rads_def_var (S, P, 'range_c')
	endif
endif

if (do_sig0) then
	if (nr_only) then
		call rads_def_var (S, P, 'sig0_ku_nr')
	else
		call rads_def_var (S, P, 'sig0_ku')
		if (has_mle3) call rads_def_var (S, P, 'sig0_ku_mle3')
		if (has_c) call rads_def_var (S, P, 'sig0_c')
	endif
endif

if (do_wind) then
	if (nr_only) then
		call rads_def_var (S, P, 'wind_speed_alt_nr')
	else
		call rads_def_var (S, P, 'wind_speed_alt')
		if (has_mle3) call rads_def_var (S, P, 'wind_speed_alt_mle3')
	endif
endif

if (do_ssb) then
	if (nr_only) then
		call rads_def_var (S, P, 'ssb_cls_nr')
		if (lr) call rads_def_var (S, P, 'ssb_cls_c_nr')
	else
		call rads_def_var (S, P, 'ssb_cls')
		if (has_mle3) call rads_def_var (S, P, 'ssb_cls_mle3')
		if (has_c) call rads_def_var (S, P, 'ssb_cls_c')
	endif
endif

if (do_rain) then
	call rads_def_var (S, P, 'qual_alt_rain_ice')
	call rads_def_var (S, P, 'dsig0_atten')
	call rads_def_var (S, P, 'flags')
else if (do_flag) then
	call rads_def_var (S, P, 'flags')
endif

if (do_flag) then
	if (has_mle3) call rads_def_var (S, P, 'flags_mle3')
	if (has_nr) call rads_def_var (S, P, 'flags_nr')
endif

! Write out all the data

if (do_range) then
	if (nr_only) then
		call rads_put_var (S, P, 'range_ku_nr', range_ku)
	else
		call rads_put_var (S, P, 'range_ku', range_ku)
		if (has_mle3) call rads_put_var (S, P, 'range_ku_mle3', range_ku_mle3)
		if (has_c) call rads_put_var (S, P, 'range_c', range_c)
	endif
endif

if (do_sig0) then
	if (nr_only) then
		call rads_put_var (S, P, 'sig0_ku_nr', sig0_ku)
	else
		call rads_put_var (S, P, 'sig0_ku', sig0_ku)
		if (has_mle3) call rads_put_var (S, P, 'sig0_ku_mle3', sig0_ku_mle3)
		if (has_c) call rads_put_var (S, P, 'sig0_c', sig0_c)
	endif
endif

if (do_wind) then
	if (nr_only) then
		call rads_put_var (S, P, 'wind_speed_alt_nr', wind_speed_alt)
	else
		call rads_put_var (S, P, 'wind_speed_alt', wind_speed_alt)
		if (has_mle3) call rads_put_var (S, P, 'wind_speed_alt_mle3', wind_speed_alt_mle3)
	endif
endif

if (do_ssb) then
	if (nr_only) then
		call rads_put_var (S, P, 'ssb_cls_nr', ssb_cls)
		if (lr) call rads_put_var (S, P, 'ssb_cls_c_nr', ssb_cls_c)
	else
		call rads_put_var (S, P, 'ssb_cls', ssb_cls)
		if (has_mle3) call rads_put_var (S, P, 'ssb_cls_mle3', ssb_cls_mle3)
		if (has_c) call rads_put_var (S, P, 'ssb_cls_c', ssb_cls_c)
	endif
endif

if (do_rain) then
	call rads_put_var (S, P, 'qual_alt_rain_ice', qual_alt_rain_ice)
	call rads_put_var (S, P, 'dsig0_atten', dsig0_atten)
	call rads_put_var (S, P, 'flags', flags)
else if (do_flag) then
	call rads_put_var (S, P, 'flags', flags)
endif

if (do_flag) then
	if (has_mle3) call rads_put_var (S, P, 'flags_mle3', flags_mle3)
	if (has_nr) call rads_put_var (S, P, 'flags_nr', flags_nr)
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

! Set flag bit 0

subroutine set_flag (n, flags, set)
integer(fourbyteint), intent(in) :: n
real(eightbytereal), intent(inout) :: flags(:)
logical, intent(in) :: set
integer(fourbyteint) :: i, j
if (set) then
	do i = 1,n
		j = nint(flags(i))
		flags(i) = ibset(j,0)
	enddo
else
	do i = 1,n
		j = nint(flags(i))
		flags(i) = ibclr(j,0)
	enddo
endif
end subroutine set_flag

end program rads_fix_s6
