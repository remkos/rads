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
logical :: lsig0 = .false., lssb = .false., liono = .false., lflag = .false., lp2p = .false., lg02 = .false., lamr = .false.
integer, parameter :: sig0_nx = 500
real(eightbytereal) :: exp_ku_sigma0(sig0_nx), rms_exp_ku_sigma0(sig0_nx), f
real(eightbytereal) :: bias_range(3) = 0d0, bias_swh(3) = 0d0, bias_sig0(3) = 0d0, bias_off_nadir(2) = 0d0, &
	dwind(2) = 0d0, drain(3) = (/ 0.51d0, 0.52d0, 0.72d0 /)
real(eightbytereal), parameter :: sig0_dx = 0.1d0, gate_width = 0.3795d0, sign_error = 2 * 0.528d0

! AMR variables to copy

real(eightbytereal), parameter :: sec2000 = 473299200d0
integer, parameter :: n_amr_var = 7
type :: amr_type
	character(len=rads_varl) :: amr_name, alt_name
end type
type(amr_type) :: amr_var(n_amr_var)
! Assuming no change in flags
amr_var(1) = amr_type ('rad_wind_speed', 'wind_speed_rad')
amr_var(2) = amr_type ('rad_wet_tropo_cor', 'wet_tropo_rad')
amr_var(3) = amr_type ('rad_tb_187', 'tb_187')
amr_var(4) = amr_type ('rad_tb_238', 'tb_238')
amr_var(5) = amr_type ('rad_tb_340', 'tb_340')
amr_var(6) = amr_type ('rad_cloud_liquid_water', 'liquid_water_rad')
amr_var(7) = amr_type ('rad_water_vapor', 'water_vapor_rad')

! Scan command line

call synopsis ('--head')
call rads_set_options ( &
	' range sig0 wind:: ssb rain:: all iono bias-range: bias-swh: bias-sig0: bias-off-nadir: p2p g02 amr flag-bit0')
call rads_init (S)

! Set defaults

select case (S%branch(1)(:5))
case ('6a', '6a.lr')
	dwind = (/ 0.57d0, 0.58d0 /)
case ('6a.hr')
	dwind = (/ 1.89d0, 2.02d0 /)
case ('6b', '6b.lr')
	dwind = (/ 0.65d0, 0.66d0 /)
case ('6b.hr')
	dwind = (/ 1.92d0, 2.07d0 /)
case default
	call rads_exit ('Unknown case '//S%branch(1)(:5))
end select

! Determine conversion factor from range difference to ionospheric correction

f = 1d0/(1d0-(S%frequency(1)/S%frequency(2))**2)

! Scan the command line for command-specific options

do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('wind')
		read (rads_opt(i)%arg, *, iostat=ios) dwind
		if (ios == 0) lwind = 2
	case ('ssb')
		lssb = .true.
	case ('rain')
		read (rads_opt(i)%arg, *, iostat=ios) drain
		lrain = 1
		if (ios == 0) lrain = 2
	case ('all')
		lrain = 1
		lwind = 1
		lssb = .true.
		liono = .true.
	case ('iono')
		liono = .true.
	case ('bias-range')
		read (rads_opt(i)%arg, *, iostat=ios) bias_range
	case ('bias-swh')
		read (rads_opt(i)%arg, *, iostat=ios) bias_swh
	case ('bias-sig0')
		read (rads_opt(i)%arg, *, iostat=ios) bias_sig0
	case ('bias-off-nadir')
		read (rads_opt(i)%arg, *, iostat=ios) bias_off_nadir
	case ('flag-bit0')
		lflag = .true.
	case ('p2p')
		lp2p = .true.
		if (S%sat /= '6b') call rads_exit ('--p2p can be used only with Sentinel-6B')
		lrain = 1
		lwind = 1
		lssb = .true.
	case ('g02')
		if (S%sat /= '6b') call rads_exit ('--g02 can be used only with Sentinel-6B')
		lg02 = .true.
		lrain = 1
		lwind = 1
		lssb = .true.
		liono = .true.
	case ('amr')
		if (S%sat /= '6b') call rads_exit ('--amr can be used only with Sentinel-6B')
		lamr = .true.
	end select
enddo

! If nothing selected, stop here

if (.not.(lsig0 .or. lwind > 0 .or. lssb .or. lrain > 0 .or. lflag .or. lp2p .or. lg02 .or. lamr .or. &
	any(bias_range /= 0d0) .or. any(bias_swh /= 0d0) .or. any(bias_sig0 /= 0d0) .or. any(bias_off_nadir /= 0d0))) stop

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
'  --rain[=KU,NR,C]          Add biases to sigma0 (Ku conv, Ku NR, C, in dB) before calling rain model' / &
'                            (with --rain use default 0.51, 0.52, 0.72 for LR only)' / &
'  --wind[=KU,NR]            Add bias to sigma0 (Ku conv, Ku NR, in dB) before calling wind model' / &
'                            (with --wind use defaults:' / &
'                            S6A LR: 0.57, 0.58  S6A HR: 1.89, 2.02' / &
'                            S6B LR: 0.65, 0.66  S6B HR: 1.92, 2.07)' / &
'  --ssb                     Update SSB (with --wind)' / &
'  --iono                    Correct ionospheric corrections also for range biases (LR only)' / &
'  --all                     Same as --rain --wind --ssb --iono' / &
'  --bias-range=KU,NR,C      Add additional bias to range (Ku conv, Ku NR, C, in m)' / &
'  --bias-swh=KU,NR,C        Add additional bias to SWH (Ku conv, Ku NR, C, in m)' / &
'  --bias-sig0=KU,NR,C       Add additional bias to sig0 (Ku conv, Ku NR, C, in dB)' / &
'  --bias-off-nadir=KU,NR    Add additional bias to off nadir angle (Ku conv, Ku NR, in deg**2)' / &
'  --p2p                     Counter effects of reduced 2.2 kHz LR waveform accumulation,' / &
'                            only for S6B, implies --rain --wind --ssb)' / &
'  --g02                     Upgrade S6B OPE Side A data to G02 standards, implies --all' / &
'  --amr                     Update S6B AMR-C data' / &
'  --flag-bit0               Clear (VAL=0) or set (VAL=1) flag bit 0 and update attributes')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
use netcdf
use rads_netcdf
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n), latency(n), range_ku(n), range_ku_nr(n), range_c(n), &
	sig0_ku(n), sig0_ku_nr(n), sig0_c(n), dsig0_atmos_ku(n), dsig0_atmos_c(n), dsig0_atten(n), &
	off_nadir_angle2_wf_ku(n), off_nadir_angle2_wf_ku_nr(n), swh_ku(n), swh_ku_nr(n), swh_c(n), &
	wind_speed_alt(n), wind_speed_alt_nr(n), qual_alt_rain_ice(n), flags(n), flags_nr(n), &
	ssb_cls(n), ssb_cls_nr(n), ssb_cls_c(n), ssb_cls_c_nr(n), &
	iono_alt(n), iono_alt_smooth(n), iono_alt_nr(n), iono_alt_smooth_nr(n), tmp(n)
real(eightbytereal) :: drange(3), dswh(3), dsig0(3), diono(3)
real(eightbytereal), parameter :: dsig0_p2p = 10d0 * log10(4d0)	! Impact of reducing the waveform accumulation by factor 4
real(eightbytereal), parameter :: dsig0_rmc = 0.016d0	! Change of residual_onboard_proc_rmc_chd
real(eightbytereal), parameter :: p2p_start = 1295016747d0, p2p_end = 1295872205d0 ! Start and end time of P2P period (S6B Cycle 7)

logical :: lr, redundant, val, do_p2p, do_range, do_swh, do_sig0, do_off_nadir, do_wind, do_ssb, do_rain, do_flag, do_iono
character(len=3) :: chd_ver, cnf_ver, baseline
character(len=rads_cmdl) :: amr_fname
integer :: ncid, i, j, l, idx(n), n_amr_val
real(eightbytereal), allocatable :: amr_val(:)

! Initialise

call log_pass (P)
drange = bias_range
dsig0 = bias_sig0
dswh = bias_swh
drain = 0d0
diono = f * (drange(3) - drange)

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

! Get latency

call rads_get_var (S, P, 'latency', latency, .true.)

! Update of S6B Side A OPE data to G02

if (lg02) then
	drange = 0d0
	dsig0 = 0d0
	dswh = 0d0
	if (val) call rads_exit ('--g02 only supported for S6B OPE data')

	! Side A updates
	if (.not.redundant) then
		! AR4322: 3.4 mm for Side A with CHAN 001, with 5.4 mm subtracted to make Side B unchanged with the fix
		drange(2) = drange(2) - 2.05d-3 ! All NR ranges
		! CHAN 002: combined effect on range and sigma0
		drange = drange + 0.15d-3 ! On all ranges
		dsig0 = dsig0 - 0.033d0
		diono = drange(3) - drange
		! CHDN 002: add -0.8 mm to HR range, -6/8 mm to HR SWH and -0.10 dB to Sigma0 for Side A with CHDN 001
		if (lr) then
			dsig0(1:2) = dsig0(1:2) - 0.10d0 ! All LR Ku-band sigma0
			bias_off_nadir = (/ -0.0072d0, -0.0075d0 /) ! All LR Ku-band off-nadir angles
		else
			drange(1:2) = drange(1:2) - 0.8d-3 ! All HR Ku-band ranges
			dswh(1:2) = dswh(1:2) - (/ 8d-3, 6d-3 /) ! HR Ku-band SWH: -8 mm for SAMOSA, -6 mm for NR
			dsig0(1:2) = dsig0(1:2) - 0.014d0 ! All HR Ku-band sigma0 (RMC to be adjusted later)
		endif
		! L2 CONF 026: change of sig0 biases for wind speed
		if (lr) then
			dwind = (/ 0.593d0, 0.608d0 /)	! Wind_System_Bias_Sigma0* for LR
		else
			dwind = (/ 1.866d0, 2.009d0 /)	! Wind_System_Bias_Sigma0* for HR
		endif
	endif
endif

! p2p: only for LR Cycle 7

do_p2p = (lp2p .and. lr .and. P%cycle == 7)

! bias_range: apply bias to ranges

do_range = any(drange /= 0d0)

! iono: correct ionospheric corrections for range_bias

do_iono = (liono .and. any(diono /= 0))

! bias_swh: apply bias to sigma0

do_swh = any(dswh /= 0d0)

! bias_sig0: apply bias to sigma0

do_sig0 = any(dsig0 /= 0d0)

! bias_off_nadir: apply bias to off-nadir angle

do_off_nadir = (lr .and. any(bias_off_nadir /= 0d0))

! wind: apply biases before calling wind model

do_wind = ((lwind == 1 .and. do_sig0) .or. lwind == 2)

! ssb: do when requested and wind or SWH has changed

do_ssb = (lssb .and. (do_wind .or. do_swh))

! rain: apply biases before calling rain model

do_rain = lr .and. ((lrain == 1 .and. do_sig0) .or. lrain == 2)

! flag bit 0: not needed for HR Side B

do_flag = lflag .and. (lr .or. .not.redundant)

! If nothing to change, skip

if (.not.(do_range .or. do_iono .or. do_swh .or. do_sig0 .or. do_off_nadir .or. &
	do_wind .or. do_ssb .or. do_rain .or. lamr .or. do_flag)) then
	call log_records(0)
	return
endif

! Adjust range for offset

if (do_range) then
	call rads_get_var (S, P, 'range_ku', range_ku, .true.)
	call rads_get_var (S, P, 'range_ku_nr', range_ku_nr, .true.)
	if (lr) then
		call rads_get_var (S, P, 'range_c', range_c, .true.)
		range_c = range_c + drange(3)
	endif
	range_ku = range_ku + drange(1)
	range_ku_nr = range_ku_nr + drange(2)
endif

! Adjust SWH for offset

if (do_swh) then
	call rads_get_var (S, P, 'swh_ku', swh_ku, .true.)
	call rads_get_var (S, P, 'swh_ku_nr', swh_ku_nr, .true.)
	if (lr) then
		call rads_get_var (S, P, 'swh_c', swh_c, .true.)
		swh_c = swh_c + dswh(3)
	endif
	swh_ku = swh_ku + dswh(1)
	swh_ku_nr = swh_ku_nr + dswh(2)

! Other cases for which the retreive wave height

else if (do_wind .or. do_ssb) then
	call rads_get_var (S, P, 'swh_ku', swh_ku, .true.)
	call rads_get_var (S, P, 'swh_ku_nr', swh_ku_nr, .true.)
endif

! Insert reprocessed AMR data

if (lamr) then
	! Open the corresponding (reprocessed) AMR data file
	l = len_trim(P%fileinfo(1)%name)
	call parseenv ('${RADSROOT}/ext/6b/amr_rep/' // P%fileinfo(1)%name(l-18:l), amr_fname)
	l = index(amr_fname, 'amr_rep/')
	write (*, 600) trim(amr_fname(l:))
	call nfs(nf90_open (amr_fname, nf90_nowrite, ncid))
	call nfs(nf90_inquire_dimension (ncid, 1, len=n_amr_val))
	allocate (amr_val(0:n_amr_val))
	amr_val = nan
	! Match RADS and AMR passes based on time
	call get_var (ncid, 'time', amr_val(1:))
	amr_val = amr_val + sec2000
	call rads_get_var (S, P, 'time', time, .true.)
	idx = 0
	do i = 1,n
		j = minloc(abs(amr_val(1:)-time(i)), 1)
		if (abs(amr_val(j)-time(i)) < 0.5d0) idx(i) = j
	enddo
	! Get original sigma0 and atmospheric correction
	call rads_get_var (S, P, 'sig0_ku', sig0_ku)
	call rads_get_var (S, P, 'sig0_ku_nr', sig0_ku_nr)
	call rads_get_var (S, P, 'dsig0_atmos_ku', tmp)
	! Correct sigma0 for the new atmospheric correction
	call get_var (ncid, 'rad_atm_cor_sig0_ku', amr_val(1:))
	dsig0_atmos_ku = amr_val(idx)
	sig0_ku = sig0_ku - tmp + dsig0_atmos_ku
	sig0_ku_nr = sig0_ku_nr - tmp + dsig0_atmos_ku
	if (lr) then
		call rads_get_var (S, P, 'sig0_c', sig0_c)
		call rads_get_var (S, P, 'dsig0_atmos_c', tmp)
		call get_var (ncid, 'rad_atm_cor_sig0_c', amr_val(1:))
		dsig0_atmos_c = amr_val(idx)
		sig0_c = sig0_c - tmp + dsig0_atmos_c
	endif
	deallocate (amr_val)
600 format ('(Reading ',a,') ... ', $)

! Other cases for which to retrieve sigma0

else if (do_sig0 .or. do_wind .or. do_rain .or. do_p2p) then
	call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	call rads_get_var (S, P, 'sig0_ku_nr', sig0_ku_nr, .true.)
	call rads_get_var (S, P, 'dsig0_atmos_ku', dsig0_atmos_ku, .true.)
	if (lr) then
		call rads_get_var (S, P, 'sig0_c', sig0_c, .true.)
		call rads_get_var (S, P, 'dsig0_atmos_c', dsig0_atmos_c, .true.)
	endif
endif

! Correct sigma0 for P2P effect

if (do_p2p) then
	call rads_get_var (S, P, 'time', time, .true.)
	where (time >= p2p_start .and. time <= p2p_end)
		sig0_ku = sig0_ku + dsig0_p2p
		sig0_ku_nr = sig0_ku_nr + dsig0_p2p
	endwhere
endif

! Adjust RMC sigma0

if (lg02 .and. .not.lr .and. .not.redundant) then
	call rads_get_var (S, P, 'flag_alt_oper_mode', tmp)
	where (tmp == 2)
		sig0_ku = sig0_ku + dsig0_rmc
		sig0_ku_nr = sig0_ku_nr + dsig0_rmc
	endwhere
endif

! Adjust sigma0 for offset

if (do_sig0) then
	sig0_ku = sig0_ku + dsig0(1)
	sig0_ku_nr = sig0_ku_nr + dsig0(2)
	if (lr) sig0_c = sig0_c + dsig0(3)
endif
if (lamr .or. do_p2p) do_sig0 = .true. ! To write out corrected sigma0

! Adjust off nadir angle for offset

if (do_off_nadir) then
	call rads_get_var (S, P, 'off_nadir_angle2_wf_ku', off_nadir_angle2_wf_ku)
	call rads_get_var (S, P, 'off_nadir_angle2_wf_ku_nr', off_nadir_angle2_wf_ku_nr)
	off_nadir_angle2_wf_ku = off_nadir_angle2_wf_ku + bias_off_nadir(1)
	off_nadir_angle2_wf_ku_nr = off_nadir_angle2_wf_ku_nr + bias_off_nadir(2)
endif

! Compute wind speed from 2D wind model after adding biases
! Load wind model if required

if (do_wind) then
	if (lr) then
		if (need_file('AUX_WNDL_S6A_002.nc', aux_wind)) then
			if (grid_load(aux_wind,info_wind) /= 0) call rads_exit ('Error loading '//trim(aux_wind))
		endif
	else if (need_file('AUX_WNDH_S6A_002.nc', aux_wind)) then
		if (grid_load(aux_wind,info_wind) /= 0) call rads_exit ('Error loading '//trim(aux_wind))
	endif
	call grid_inter (info_wind, n, sig0_ku + dwind(1), swh_ku, wind_speed_alt)
	call grid_inter (info_wind, n, sig0_ku_nr + dwind(2), swh_ku_nr, wind_speed_alt_nr)

! Other cases for which to retireve wind speed

else if (do_ssb) then
	call rads_get_var (S, P, 'wind_speed_alt', wind_speed_alt, .true.)
	call rads_get_var (S, P, 'wind_speed_alt_nr', wind_speed_alt_nr, .true.)
endif

! Compute SSB from 2D model after updating wind
! Load SSB model if required

if (do_ssb) then
	if (lr) then
		if (need_file('AUX_SBLC_S6A_002.nc', aux_ssbc)) then
			if (grid_load(aux_ssbc,info_ssbc) /= 0) call rads_exit ('Error loading '//trim(aux_ssbc))
		endif
		call grid_inter (info_ssbc, n, wind_speed_alt, swh_ku, ssb_cls_c)
		call grid_inter (info_ssbc, n, wind_speed_alt_nr, swh_ku_nr, ssb_cls_c_nr)
		if (need_file('AUX_SBLK_S6A_002.nc', aux_ssbk)) then
			if (grid_load(aux_ssbk,info_ssbk) /= 0) call rads_exit ('Error loading '//trim(aux_ssbk))
		endif
	else
		if (need_file('AUX_SBHK_S6A_002.nc', aux_ssbk)) then
			if (grid_load(aux_ssbk,info_ssbk) /= 0) call rads_exit ('Error loading '//trim(aux_ssbk))
		endif
	endif
	call grid_inter (info_ssbk, n, wind_speed_alt, swh_ku, ssb_cls)
	call grid_inter (info_ssbk, n, wind_speed_alt_nr, swh_ku_nr, ssb_cls_nr)
	if (lr) then ! Compute average contributions to ionospheric correction
		call rads_get_var (S, P, 'ssb_cls_c', tmp)
		diono(1) = diono(1) + mean(ssb_cls_c - tmp)
		call rads_get_var (S, P, 'ssb_cls_c_nr', tmp)
		diono(2) = diono(2) + mean(ssb_cls_c_nr - tmp)
		call rads_get_var (S, P, 'ssb_cls', tmp)
		diono(1) = diono(1) - mean(ssb_cls - tmp)
		call rads_get_var (S, P, 'ssb_cls_nr', tmp)
		diono(2) = diono(2) - mean(ssb_cls_nr - tmp)
	endif
endif

! Adjust ionospheric corrections for range bias

if (do_iono) then
	call rads_get_var (S, P, 'iono_alt', iono_alt, .true.)
	call rads_get_var (S, P, 'iono_alt_smooth', iono_alt_smooth, .true.)
	call rads_get_var (S, P, 'iono_alt_nr', iono_alt_nr, .true.)
	call rads_get_var (S, P, 'iono_alt_smooth_nr', iono_alt_smooth_nr, .true.)
	diono = f * diono
	iono_alt = iono_alt + diono(1)
	iono_alt_smooth = iono_alt_smooth + diono(1)
	iono_alt_nr = iono_alt_nr + diono(2)
	iono_alt_smooth_nr = iono_alt_smooth_nr + diono(2)
endif

! Determine rain flag after adding biases

if (do_rain) then
	if (need_file('AUX_SIGL_S6A_001.nc', aux_rain)) call rain_table
	call rads_get_var (S, P, 'qual_alt_rain_ice', qual_alt_rain_ice, .true.)
	call rads_get_var (S, P, 'flags', flags, .true.)
	call rads_get_var (S, P, 'flags_nr', flags_nr, .true.)
	call compute_rain_flag (n, sig0_ku - dsig0_atmos_ku + drain(1), sig0_c - dsig0_atmos_c + drain(3), &
		dsig0_atten, qual_alt_rain_ice, flags)
	call compute_rain_flag (n, sig0_ku_nr - dsig0_atmos_ku + drain(2), sig0_c - dsig0_atmos_c + drain(3), &
		dsig0_atten, qual_alt_rain_ice, flags_nr)
endif

! Set the flag bit 0 appropriately

if (do_flag) then
	if (.not.do_rain) then
		call rads_get_var (S, P, 'flags', flags, .true.)
		call rads_get_var (S, P, 'flags_nr', flags_nr, .true.)
	endif
	call set_flag (n, flags, redundant)
	call set_flag (n, flags_nr, redundant)
endif

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! (Re)define the variables

if (do_range) then
	call rads_def_var (S, P, 'range_ku')
	call rads_def_var (S, P, 'range_ku_nr')
	if (lr) call rads_def_var (S, P, 'range_c')
endif

if (do_iono) then
	call rads_def_var (S, P, 'iono_alt')
	call rads_def_var (S, P, 'iono_alt_smooth')
	call rads_def_var (S, P, 'iono_alt_nr')
	call rads_def_var (S, P, 'iono_alt_smooth_nr')
endif

if (do_swh) then
	call rads_def_var (S, P, 'swh_ku')
	call rads_def_var (S, P, 'swh_ku_nr')
	if (lr) call rads_def_var (S, P, 'swh_c')
endif

if (do_sig0) then
	call rads_def_var (S, P, 'sig0_ku')
	call rads_def_var (S, P, 'sig0_ku_nr')
	if (lr) call rads_def_var (S, P, 'sig0_c')
endif

if (do_off_nadir) then
	call rads_def_var (S, P, 'off_nadir_angle2_wf_ku')
	call rads_def_var (S, P, 'off_nadir_angle2_wf_ku_nr')
endif

if (do_wind) then
	call rads_def_var (S, P, 'wind_speed_alt')
	call rads_def_var (S, P, 'wind_speed_alt_nr')
endif

if (do_ssb) then
	call rads_def_var (S, P, 'ssb_cls')
	call rads_def_var (S, P, 'ssb_cls_nr')
	if (lr) then
		call rads_def_var (S, P, 'ssb_cls_c')
		call rads_def_var (S, P, 'ssb_cls_c_nr')
	endif
endif

if (do_rain) then
	call rads_def_var (S, P, 'qual_alt_rain_ice')
	call rads_def_var (S, P, 'dsig0_atten')
endif

if (do_rain .or. do_flag) then
	call rads_def_var (S, P, 'flags')
	call rads_def_var (S, P, 'flags_nr')
endif

if (lamr) then
	do i = 1,n_amr_var
		call rads_def_var (S, P, amr_var(i)%alt_name)
	enddo
	call rads_def_var (S, P, 'dsig0_atmos_ku')
	if (lr) call rads_def_var (S, P, 'dsig0_atmos_c')
endif

! Write out all the data

if (do_range) then
	call rads_put_var (S, P, 'range_ku', range_ku)
	call rads_put_var (S, P, 'range_ku_nr', range_ku_nr)
	if (lr) call rads_put_var (S, P, 'range_c', range_c)
endif

if (do_iono) then
	call rads_put_var (S, P, 'iono_alt', iono_alt)
	call rads_put_var (S, P, 'iono_alt_smooth', iono_alt_smooth)
	call rads_put_var (S, P, 'iono_alt_nr', iono_alt_nr)
	call rads_put_var (S, P, 'iono_alt_smooth_nr', iono_alt_smooth_nr)
endif

if (do_swh) then
	call rads_put_var (S, P, 'swh_ku', swh_ku)
	call rads_put_var (S, P, 'swh_ku_nr', swh_ku_nr)
	if (lr) call rads_put_var (S, P, 'swh_c', swh_c)
endif

if (do_sig0) then
	call rads_put_var (S, P, 'sig0_ku', sig0_ku)
	call rads_put_var (S, P, 'sig0_ku_nr', sig0_ku_nr)
	if (lr) call rads_put_var (S, P, 'sig0_c', sig0_c)
endif

if (do_off_nadir) then
	call rads_put_var (S, P, 'off_nadir_angle2_wf_ku', off_nadir_angle2_wf_ku)
	call rads_put_var (S, P, 'off_nadir_angle2_wf_ku_nr', off_nadir_angle2_wf_ku_nr)
endif

if (do_wind) then
	call rads_put_var (S, P, 'wind_speed_alt', wind_speed_alt)
	call rads_put_var (S, P, 'wind_speed_alt_nr', wind_speed_alt_nr)
endif

if (do_ssb) then
	call rads_put_var (S, P, 'ssb_cls', ssb_cls)
	call rads_put_var (S, P, 'ssb_cls_nr', ssb_cls_nr)
	if (lr) then
		call rads_put_var (S, P, 'ssb_cls_c', ssb_cls_c)
		call rads_put_var (S, P, 'ssb_cls_c_nr', ssb_cls_c_nr)
	endif
endif

if (do_rain) then
	call rads_put_var (S, P, 'qual_alt_rain_ice', qual_alt_rain_ice)
	call rads_put_var (S, P, 'dsig0_atten', dsig0_atten)
endif

if (do_rain .or. do_flag) then
	call rads_put_var (S, P, 'flags', flags)
	call rads_put_var (S, P, 'flags_nr', flags_nr)
endif

if (lamr) then
	allocate (amr_val(0:n_amr_val))
	do i = 1,n_amr_var
		call get_var (ncid, amr_var(i)%amr_name, amr_val(1:))
		call rads_put_var (S, P, amr_var(i)%alt_name, amr_val(idx))
	enddo
	deallocate (amr_val)
	call rads_put_var (S, P, 'dsig0_atmos_ku', dsig0_atmos_ku)
	if (lr) call rads_put_var (S, P, 'dsig0_atmos_c', dsig0_atmos_c)
	call nfs(nf90_close(ncid))
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

! Computre mean

pure function mean (x)
real(eightbytereal), intent(in) :: x(:)
real(eightbytereal) :: mean, y
integer(fourbyteint) :: n
call mean_variance (x, mean, y, n)
end function mean

end program rads_fix_s6
