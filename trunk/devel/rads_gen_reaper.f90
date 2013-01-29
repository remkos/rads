!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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

!*rads_gen_reaper -- Converts REAPER ERS-1/2 data to RADS
!+
program rads_gen_reaper

! This program reads REAPER files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/eE/F.r/eEpPPPPcCCC.nc.
!     E = 1 or 2
!     F = mission phase
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: reaperraw [options] < list_of_REAPER_file_names
!
! This program handles only the REAPER ALT_2M files in netCDF format.
!-----------------------------------------------------------------------
!
!  Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_reaper - Orbit altitude
! alt_rate - Orbit altitude rate
! range_ku - Ocean range (retracked)
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_rad - Radiometer wet tropo correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! iono_gim - GIM ionosphetic correction
! iono_nic09 - NIC09 ionospheric correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! tide_solid - Solid earth tide
! tide_ocean_fes04 - FES2008 ocean tide
! tide_ocean_got47 - GOT4.7 ocean tide
! tide_load_fes04 - FES2008 load tide
! tide_load_got47 - GOT4.7 load tide
! tide_pole - Pole tide
! ssb_bm3 - SSB
! mss_cls01 - CLS01 MSS
! geoid_egm2008 - EGM2008 geoid
! mss_ucl04 - UCL04 MSS
! swh_ku - Significant wave height
! sig0_ku - Sigma0
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! range_rms_ku - Std dev of range
! range_numval_ku - Nr of averaged range measurements
! topo_macess - MACESS topography
! tb_238 - Brightness temperature (23.8 GHz)
! tb_365 - Brightness temperature (36.5 GHz)
! peakiness_ku - Peakiness
! flags - Engineering flags
! swh_rms_ku - Std dev of SWH
! sig0_rms_ku - Std dev of sigma0
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! liquid_water - Liquid water content
! water_vapor_content - Water vapor content
! tide_equil - Long-period equilibrium tide
! tide_non_equil - Long-period non-equilibrium tide
!
! On ALT_2_ and ALT_2S only:
! wind_speed_alt - Altimeter wind speed
! drange_cal - Internal calibration correction to range (appied)
! drange_fm - Doppler correction (applied)
! dsig0_atmos_ku - Sigma0 attenuation
! mqe - Mean quadratic error of waveform fit
!-----------------------------------------------------------------------
use typesizes
use netcdf
use rads
use rads_misc
use rads_time
use rads_netcdf
use rads_devel

! Command line arguments

integer(fourbyteint) :: verbose=0,c0=0,c1=999,ios
real(eightbytereal) :: t0, t1
character(160) :: arg
character(20) :: sat='',newsat='',oldsat='',optopt,optarg
character(80), parameter :: optlist='v debug: sat: cycle: t: mjd: sec: ymd: doy:'

! Header variables

character(2) :: mission
character(1) :: phasenm(2)
character(80) :: l2_proc_time, l2_version
logical :: alt_2m
real(eightbytereal) :: tnode(2),lnode(2)
integer(fourbyteint) :: datanr=0,nrec,ncid,varid,orbitnr(2),cyclenr(2),passnr(2),ers

! Data variables

integer(fourbyteint), parameter :: mrec = 10000
real(eightbytereal) :: time(mrec), lat(mrec), lon(mrec), alt_reaper(mrec), alt_rate(mrec), range_ku(mrec), &
	dry_tropo_ecmwf(mrec), wet_tropo_rad(mrec), wet_tropo_ecmwf(mrec), iono_gim(mrec), iono_nic09(mrec), &
	inv_bar_static(mrec), inv_bar_mog2d(mrec), tide_solid(mrec), tide_ocean_fes04(mrec), tide_ocean_got47(mrec), &
	tide_load_fes04(mrec), tide_load_got47(mrec), tide_pole(mrec), ssb_bm3(mrec), &
	mss_cls01(mrec), geoid_egm2008(mrec), mss_ucl04(mrec), swh_ku(mrec), sig0_ku(mrec), wind_speed_alt(mrec), &
	wind_speed_ecmwf_u(mrec), wind_speed_ecmwf_v(mrec), range_rms_ku(mrec), range_numval_ku(mrec), &
	topo_macess(mrec), tb_238(mrec), tb_365(mrec), peakiness_ku(mrec), drange_cal(mrec), drange_fm(mrec), &
	swh_rms_ku(mrec), sig0_rms_ku(mrec), off_nadir_angle2_wf_ku(mrec), dsig0_atmos_ku(mrec), &
	liquid_water(mrec), water_vapor_content(mrec), mqe(mrec), tide_equil(mrec), tide_non_equil(mrec)
integer(twobyteint) :: flags(mrec)
type(rads_sat) :: S
type(rads_pass) :: P
character(640) :: original = ''

! Other local variables

real(eightbytereal), parameter :: sec1990=157766400d0	! UTC seconds from 1 Jan 1985 to 1 Jan 1990
real(eightbytereal), parameter :: picosec_to_mm=0.5d-12*299792458d3	! picoseconds of 2-way range to mm 1-way
real(eightbytereal) :: nan, sum_d_applied
real(eightbytereal), allocatable :: a(:),b(:),c(:),d(:,:),sum_c_applied(:)
integer(fourbyteint) :: i,k,i0,i1,flag
logical :: new

! Initialise

nan = 0d0
nan = nan / nan
t0 = nan
t1 = nan
550 format (a)
551 format (a,' ...')

! Scan command line for options

call synopsis ()
do
	call getopt (optlist, optopt, optarg)
	select case (optopt)
	case ('!')
		exit
	case ('v')
		verbose = 1
	case ('debug')
		read (optarg,*) verbose
	case ('sat')
		sat = optarg
	case ('cycle')
		c1 = -1
		read (optarg,*,iostat=ios) c0,c1
		if (c1 < c0) c1 = c0
	case default
		if (.not.dateopt (optopt, optarg, t0, t1)) then
			call synopsis ()
			stop
		endif
	end select
enddo

!----------------------------------------------------------------------
! Read all file names for standard input
!----------------------------------------------------------------------

files: do
	read (*,550,iostat=ios) arg
	if (ios /= 0) exit files

! Check input file name

	i = index(arg,'ERS_ALT_2')
	if (i <= 0) then
		write (*,550) 'Error: Wrong input file'
		cycle files
	endif
	alt_2m = index(arg,'ERS_ALT_2M') > 0

! Open input file

	write (*,551) trim(arg)
	if (nf90_open(arg,nf90_nowrite,ncid) /= nf90_noerr) then
		write (*,550) 'Error opening file'
		cycle files
	endif

! Reduce file name to basename only

	i = index(arg,'/',.true.)
	arg = arg(i+1:)

! Check for ERS-1 or -2
! Do not trust 'mission' attribute. It is always 'E1'.

	if (arg(:2) == 'E1') then
		ers = 1
		mission = 'e1'
	else if (arg(:2) == 'E2') then
		ers = 2
		mission = 'e2'
	else
		write (*,550) 'Error: Unknown file type'
		cycle files
	endif

! Read header records

	call nfs(nf90_inq_dimid(ncid,'Record',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
	if (nrec > mrec) then
		write (*,'("Error: Too many measurements:",i5)') nrec
		cycle files
	endif
	call nfs(nf90_get_att(ncid,nf90_global,'l2_proc_time',l2_proc_time))
	call nfs(nf90_get_att(ncid,nf90_global,'l2_software_ver',l2_version))

	i0 = datanr+1
	i1 = datanr+nrec
	if (i1 > mrec) then
		write (*,'("Error: Too many accumulated measurements:",i5)') i1
		cycle files
	endif

! Allocate arrays

	allocate (a(nrec),b(nrec),c(nrec),d(20,nrec),sum_c_applied(nrec))

! Time and orbit: Low rate

	call get_var_1d ('time_day_1hz',a)
	call get_var_1d ('time_milsec_1hz',b)
	call get_var_1d ('time_micsec_1hz',c)
	time(i0:i1) = a * 86400d0 + b * 1d-3 + c * 1d-6 + sec1990
	call get_var_1d ('latitude_1hz', lat(i0:i1))
	call get_var_1d ('longitude_1hz', lon(i0:i1))
	call get_var_1d ('altitude_1hz', alt_reaper(i0:i1))
	call get_var_1d ('altitude_rate_1hz', alt_rate(i0:i1))
	call get_var_1d ('wf_attitude_1hz', a)
	off_nadir_angle2_wf_ku(i0:i1) = a * a

! Range data: Low rate

	call get_var_1d ('ocean_range_1hz', range_ku(i0:i1))
	call get_var_1d ('ocean_stdev_1hz', range_rms_ku(i0:i1))
	call get_var_1d ('ocean_valid_num_1hz', range_numval_ku(i0:i1))
	call invalidate (range_numval_ku(i0:i1) == 0, range_ku(i0:i1))
	call invalidate (range_numval_ku(i0:i1) <= 1, range_rms_ku(i0:i1))
	if (alt_2m) then ! Different name in ALT_2M
		call get_var_1d ('ocean_valid_bitmap_1hz', a)
	else
		call get_var_1d ('f_ocean_valid_bitmap_1hz', a)
	endif
	! For some reason ALT_2M has 1-Hz peakiness
	! where ALT_2S and ALT_2_ have ocean_wind_1hz
	if (alt_2m) then
		call get_var_1d ('wf_pk_1hz', peakiness_ku(i0:i1))
	else
		call get_var_1d ('ocean_wind_1hz', wind_speed_alt(i0:i1))
	endif

! Range data: High rate

	if (.not.alt_2m) then
		call get_var_2d ('ocean_mean_quadratic_error', d)
		call mean_1hz (d,mqe(i0:i1),a)
		call get_var_2d ('wf_pk', d)
		call mean_1hz (d,peakiness_ku(i0:i1),a)
		call get_var_2d ('dop_c+delta_dop_c', d)
		call mean_1hz (d,drange_fm(i0:i1),a)
		call get_var_2d ('sptr_jumps_c', d)
		drange_cal(i0:i1) = d(1,:) * picosec_to_mm
		call get_var_2d ('sum_c_applied', d)
		sum_c_applied = d(1,:)
		do i = 1,nrec
			do k = 2,19
				if (d(k,i) /= d(1,i)) write (*,*) 'Error: sum_c_applied changed:',i,k,d(1,i),d(k,i)
			enddo
		enddo
	endif

! Sigma zero: Low rate

	call get_var_1d ('ocean_sig0_1hz', sig0_ku(i0:i1))
	call get_var_1d ('ocean_sig0_stdev_1hz', swh_rms_ku(i0:i1))
	call get_var_1d ('ocean_sig0_valid_num_1hz', a)
	call invalidate (a(1:nrec) == 0, sig0_ku(i0:i1))
	call invalidate (a(1:nrec) == 0, wind_speed_alt(i0:i1))
	call invalidate (a(1:nrec) == 0, peakiness_ku(i0:i1))
	call invalidate (a(1:nrec) <= 1, sig0_rms_ku(i0:i1))

! SWH: Low rate

	call get_var_1d ('swh_signed_1hz', swh_ku(i0:i1))
	call get_var_1d ('swh_stdev_1hz', swh_rms_ku(i0:i1))
	call get_var_1d ('swh_valid_num_1hz', a)
	call invalidate (a(1:nrec) == 0, swh_ku(i0:i1))
	call invalidate (a(1:nrec) <= 1, swh_rms_ku(i0:i1))
	call get_var_1d ('em_bias_1hz', ssb_bm3(i0:i1))

! MWR Flags: Low rate

	call get_var_1d ('f_sea_ice_flag_1hz', a)
	flags(i0:i1) = 0
	call flag_set (a(1:nrec) == 1, flags(i0:i1), 8)

! MWR: Low rate

	call get_var_1d ('tb_23_8_1hz', tb_238(i0:i1))
	call get_var_1d ('tb_36_5_1hz', tb_365(i0:i1))
	call get_var_1d ('f_mwr_srf_typ_1hz', a)
	call flag_set (a(1:nrec) == 1, flags(i0:i1), 5)
	call get_var_1d ('f_mwr_interp_qual_1hz', a)
	call invalidate (a(1:nrec) == 3, tb_238(i0:i1))
	call invalidate (a(1:nrec) == 3, tb_365(i0:i1))
	if (alt_2m) then ! Different name in ALT_2M
		call get_var_1d ('f_mwr_valid_1hz', a)
	else
		call get_var_1d ('f_MWR_valid_1hz', a)
	endif
	call invalidate (a(1:nrec) == 1, tb_238(i0:i1))
	call invalidate (a(1:nrec) == 1, tb_365(i0:i1))

! Atmospheric and geophysical: Low rate

	call get_var_1d ('dry_c_1hz', dry_tropo_ecmwf(i0:i1))
	call get_var_1d ('ib_c_1hz', inv_bar_static(i0:i1))
	call get_var_1d ('mog2d_c_1hz', inv_bar_mog2d(i0:i1))
	call get_var_1d ('wet_c_mod_1hz', wet_tropo_ecmwf(i0:i1))
	call get_var_1d ('wet_c_mwr_1hz', wet_tropo_rad(i0:i1))
	call get_var_1d ('water_vapor_content_1hz', water_vapor_content(i0:i1))
	call get_var_1d ('liquid_water_content_1hz', liquid_water(i0:i1))
	call get_var_1d ('u_wind_1hz', wind_speed_ecmwf_u(i0:i1))
	call get_var_1d ('v_wind_1hz', wind_speed_ecmwf_v(i0:i1))
	call get_var_1d ('iono_c_mod_1hz', iono_nic09(i0:i1))
	call get_var_1d ('iono_c_gps_1hz', iono_gim(i0:i1))
	call get_var_1d ('h_mss_cls01_1hz', mss_cls01(i0:i1))
	call get_var_1d ('h_mss_ucl04_1hz', mss_ucl04(i0:i1))
	call get_var_1d ('h_geo_1hz', geoid_egm2008(i0:i1))
	call get_var_1d ('h_ot_1hz', tide_ocean_got47(i0:i1))
	call get_var_1d ('h_ot2_1hz', tide_ocean_fes04(i0:i1))
	call get_var_1d ('h_olt_1hz', tide_load_got47(i0:i1))
	call get_var_1d ('h_olt2_1hz', tide_load_fes04(i0:i1))
	call get_var_1d ('h_lpt_1hz', tide_equil(i0:i1))
	call get_var_1d ('h_lptne_1hz', tide_non_equil(i0:i1))
	call get_var_1d ('h_set_1hz', tide_solid(i0:i1))
	call get_var_1d ('h_pol_1hz', tide_pole(i0:i1))

	call get_var_1d ('f_srf_typ_1hz', a)
	call flag_set (a(1:nrec) >= 3, flags(i0:i1), 2)
	call flag_set (a(1:nrec) >= 2, flags(i0:i1), 4)
	call flag_set (a(1:nrec) >= 1, flags(i0:i1), 5)

	call get_var_1d ('h_odle_1hz', topo_macess(i0:i1))
	if (.not.alt_2m) call get_var_1d ('sig0_attn_c_1hz', dsig0_atmos_ku(i0:i1))

	call get_var_1d ('f_corr_error_1hz', a)
	do i = 1,nrec
		k = datanr + i
		flag = nint(a(i))
		if (btest(flag, 1)) dry_tropo_ecmwf(k) = nan
		if (btest(flag, 2)) inv_bar_static(k) = nan
		if (btest(flag, 3)) inv_bar_mog2d(k) = nan
		if (btest(flag, 5)) wet_tropo_ecmwf(k) = nan
		if (btest(flag, 6)) then
			wet_tropo_rad(k) = nan
			liquid_water(k) = nan
			water_vapor_content(k) = nan
		endif
		if (btest(flag, 7)) wind_speed_ecmwf_u(k) = nan
		if (btest(flag, 8)) wind_speed_ecmwf_v(k) = nan
		if (btest(flag, 9)) iono_nic09(k) = nan
		if (btest(flag,10)) iono_gim(k) = nan
		if (btest(flag,11)) mss_cls01(k) = nan
		if (btest(flag,12)) geoid_egm2008(k) = nan
		if (btest(flag,13)) tide_ocean_got47(k) = nan
		if (btest(flag,14)) tide_ocean_fes04(k) = nan
		if (btest(flag,15)) tide_load_got47(k) = nan
		if (btest(flag,16)) tide_load_fes04(k) = nan
		if (btest(flag,17)) tide_equil(k) = nan
		if (btest(flag,18)) tide_non_equil(k) = nan
		if (btest(flag,19)) tide_solid(k) = nan
		if (btest(flag,20)) tide_pole(k) = nan
		if (btest(flag,22)) topo_macess(k) = nan
		if (btest(flag,23)) ssb_bm3(k) = nan
		if (btest(flag,27)) mss_ucl04(k) = nan
		if (btest(flag,28)) dsig0_atmos_ku(k) = nan
	enddo

! Remove applied corrections from range

	call get_var_1d ('f_corr_applied_1hz', a)
	do i = 1,nrec
		k = datanr + i
		flag = nint(a(i))
		sum_d_applied = 0d0
		if (btest(flag, 1)) sum_d_applied = sum_d_applied + dry_tropo_ecmwf(k)
		if (btest(flag, 2)) sum_d_applied = sum_d_applied + inv_bar_static(k)
		if (btest(flag, 3)) sum_d_applied = sum_d_applied + inv_bar_mog2d(k)
		if (btest(flag, 5)) sum_d_applied = sum_d_applied + wet_tropo_ecmwf(k)
		if (btest(flag, 6)) sum_d_applied = sum_d_applied + wet_tropo_rad(k)
		if (btest(flag, 7)) sum_d_applied = sum_d_applied + wind_speed_ecmwf_u(k)
		if (btest(flag, 8)) sum_d_applied = sum_d_applied + wind_speed_ecmwf_v(k)
		if (btest(flag, 9)) sum_d_applied = sum_d_applied + iono_nic09(k)
		if (btest(flag,10)) sum_d_applied = sum_d_applied + iono_gim(k)
		if (btest(flag,11)) sum_d_applied = sum_d_applied + mss_cls01(k)
		if (btest(flag,12)) sum_d_applied = sum_d_applied + geoid_egm2008(k)
		if (btest(flag,13)) sum_d_applied = sum_d_applied + tide_ocean_got47(k)
		if (btest(flag,14)) sum_d_applied = sum_d_applied + tide_ocean_fes04(k)
		if (btest(flag,15)) sum_d_applied = sum_d_applied + tide_load_got47(k)
		if (btest(flag,16)) sum_d_applied = sum_d_applied + tide_load_fes04(k)
		if (btest(flag,17)) sum_d_applied = sum_d_applied + tide_equil(k)
		if (btest(flag,18)) sum_d_applied = sum_d_applied + tide_non_equil(k)
		if (btest(flag,19)) sum_d_applied = sum_d_applied + tide_solid(k)
		if (btest(flag,20)) sum_d_applied = sum_d_applied + tide_pole(k)
		if (btest(flag,22)) sum_d_applied = sum_d_applied + topo_macess(k)
		if (btest(flag,23)) sum_d_applied = sum_d_applied + ssb_bm3(k)
		if (btest(flag,27)) sum_d_applied = sum_d_applied + mss_ucl04(k)
		range_ku(k) = range_ku(k) - sum_d_applied
		if (.not.alt_2m .and. .not.(sum_d_applied == sum_c_applied(i))) &
			write (*,*) "Error: sum_c_applied wrong: ",i,sum_d_applied,sum_c_applied(i)
	enddo

! Close this input file

	deallocate (a,b,c,d,sum_c_applied)

	call nfs(nf90_close(ncid))

	datanr = i1

! Add input file name to "original" string

	original = trim(original) // rads_linefeed // arg

! Dump whatever data we can

	do
		new = erspass (ers, time(1), orbitnr(1), phasenm(1), cyclenr(1), passnr(1), tnode(1), lnode(1))
		do i = 2,datanr
			if (erspass (ers, time(i), orbitnr(2), phasenm(2), cyclenr(2), passnr(2), tnode(2), lnode(2))) exit
		enddo
		if (i > datanr) exit
		nrec = i - 1
		call write_data (nrec)
		datanr = datanr - nrec
		call move_data (nrec, datanr)
	enddo

! If we have data left, keep last filename in "original" string, otherwise clear it

	if (datanr == 0) then
		original = ''
	else
		original = rads_linefeed // arg
	endif

enddo files ! Each file

! Dump whatever remains

new = erspass (ers, time(1), orbitnr(1), phasenm(1), cyclenr(1), passnr(1), tnode(1), lnode(1))
call write_data (datanr)

call rads_end (S)

contains

subroutine synopsis ()
if (rads_version ('$Revision$', 'Write REAPER data to RADS')) return
write (*,1310)
1310 format (/ &
'syntax: reaperraw [options] < list_of_REAPER_file_names'// &
'This program converts REAPER ALT_2M files to RADS data'/ &
'files with the name $RADSDATAROOT/data/eE/F.r/pPPPP/eEpPPPPcCCC.nc.'/ &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

subroutine move_data (bottom, n)
integer(fourbyteint), intent(in) :: bottom, n
integer(fourbyteint) :: i0, i1
if (n == 0) return
i0 = bottom + 1
i1 = bottom + n
time(1:n) = time(i0:i1)
lat(1:n) = lat(i0:i1)
lon(1:n) = lon(i0:i1)
alt_reaper(1:n) = alt_reaper(i0:i1)
alt_rate(1:n) = alt_rate(i0:i1)
range_ku(1:n) = range_ku(i0:i1)
dry_tropo_ecmwf(1:n) = dry_tropo_ecmwf(i0:i1)
wet_tropo_rad(1:n) = wet_tropo_rad(i0:i1)
wet_tropo_ecmwf(1:n) = wet_tropo_ecmwf(i0:i1)
iono_gim(1:n) = iono_gim(i0:i1)
iono_nic09(1:n) = iono_nic09(i0:i1)
inv_bar_static(1:n) = inv_bar_static(i0:i1)
inv_bar_mog2d(1:n) = inv_bar_mog2d(i0:i1)
tide_solid(1:n) = tide_solid(i0:i1)
tide_ocean_fes04(1:n) = tide_ocean_fes04(i0:i1)
tide_ocean_got47(1:n) = tide_ocean_got47(i0:i1)
tide_load_fes04(1:n) = tide_load_fes04(i0:i1)
tide_load_got47(1:n) = tide_load_got47(i0:i1)
tide_pole(1:n) = tide_pole(i0:i1)
ssb_bm3(1:n) = ssb_bm3(i0:i1)
mss_cls01(1:n) = mss_cls01(i0:i1)
geoid_egm2008(1:n) = geoid_egm2008(i0:i1)
mss_ucl04(1:n) = mss_ucl04(i0:i1)
swh_ku(1:n) = swh_ku(i0:i1)
sig0_ku(1:n) = sig0_ku(i0:i1)
wind_speed_ecmwf_u(1:n) = wind_speed_ecmwf_u(i0:i1)
wind_speed_ecmwf_v(1:n) = wind_speed_ecmwf_v(i0:i1)
range_rms_ku(1:n) = range_rms_ku(i0:i1)
range_numval_ku(1:n) = range_numval_ku(i0:i1)
topo_macess(1:n) = topo_macess(i0:i1)
tb_238(1:n) = tb_238(i0:i1)
tb_365(1:n) = tb_365(i0:i1)
flags(1:n) = flags(i0:i1)
drange_cal(1:n) = drange_cal(i0:i1)
swh_rms_ku(1:n) = swh_rms_ku(i0:i1)
sig0_rms_ku(1:n) = sig0_rms_ku(i0:i1)
off_nadir_angle2_wf_ku(1:n) = off_nadir_angle2_wf_ku(i0:i1)
liquid_water(1:n) = liquid_water(i0:i1)
water_vapor_content(1:n) = water_vapor_content(i0:i1)
tide_equil(1:n) = tide_equil(i0:i1)
tide_non_equil(1:n) = tide_non_equil(i0:i1)
peakiness_ku(1:n) = peakiness_ku(i0:i1)
wind_speed_alt(1:n) = wind_speed_alt(i0:i1)
drange_fm(1:n) = drange_fm(i0:i1)
dsig0_atmos_ku(1:n) = dsig0_atmos_ku(i0:i1)
mqe(1:n) = mqe(i0:i1)
end subroutine move_data

subroutine write_data (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i
real(eightbytereal) :: dhellips, dh

550 format (a)

if (datanr == 0) return	! Skip empty data sets
if (cyclenr(1) < c0 .or. cyclenr(1) > c1) return	! Skip chunks that are not of the selected cycle
if (tnode(1) < t0 .or. tnode(1) > t1) return	! Skip equator times that are not of selected range

! Determine mission and phase
! Use sat from command line if given, otherwise autogenerate

if (sat /= '') then
	newsat = sat
else
	newsat = mission(1:2)//'/'//strtolower(phasenm(1))//'.r'
endif

! Initialise settings for current mission.
! Set the list of variables to be considered.

if (newsat /= oldsat) then
	oldsat = newsat
	call rads_init (S, newsat, verbose)
	call rads_parse_varlist (S, &
		'time,lat,lon,alt_reaper,alt_rate,range_ku,dry_tropo_ecmwf,wet_tropo_rad,wet_tropo_ecmwf,' // &
		'iono_gim,iono_nic09,inv_bar_static,inv_bar_mog2d,tide_solid,tide_ocean_fes04,tide_ocean_got47,' // &
		'tide_load_fes04,tide_load_got47,tide_pole,ssb_bm3,mss_cls01,geoid_egm2008,mss_ucl04,swh_ku,sig0_ku,' // &
		'wind_speed_ecmwf_u,wind_speed_ecmwf_v,range_rms_ku,range_numval_ku,topo_macess,tb_238,tb_365,' // &
		'peakiness_ku,flags,swh_rms_ku,sig0_rms_ku,off_nadir_angle2_wf_ku,liquid_water,water_vapor_content,' // &
		'tide_equil,tide_non_equil')
	if (.not.alt_2m) call rads_parse_varlist (S, 'wind_speed_alt,drange_cal,drange_fm,dsig0_atmos_ku,mqe')
endif

! Store relevant info

call rads_init_pass_struct (S, P)

P%cycle = cyclenr(1)
P%pass = passnr(1)
P%start_time = time(1)
P%end_time = time(n)
P%equator_time = tnode(1)
P%equator_lon = lnode(1)
P%original = 'REAPER '//trim(l2_version)//' data of '//trim(l2_proc_time)//trim(original)

! Open output file

call rads_create_pass (S, P, n)

! Change reference of altitude, MSS and geoid from WGS84 to TOPEX ellipsoid

do i = 1,n
	dh = dhellips(1,lat(i)*1d-6)*1d3
	alt_reaper(i) = alt_reaper(i) + dh
	mss_cls01(i) = mss_cls01(i) + dh
	geoid_egm2008(i) = geoid_egm2008(i) + dh
	mss_ucl04(i) = mss_ucl04(i) + dh
enddo

! Switch off (for output) iono_gim, and MWR variables when all NaN

S%sel%noedit = .false.
if (all(isnan(iono_gim(1:n)))) then
	write (*,550,advance='no') 'No iono_gim '
	call switch_off (S, 'iono_gim')
endif
if (all(isnan(wet_tropo_rad(1:n)))) then
	write (*,550,advance='no') 'No wet_tropo_rad '
	call switch_off (S, 'wet_tropo_rad')
	call switch_off (S, 'tb_238')
	call switch_off (S, 'tb_365')
	call switch_off (S, 'water_vapor_content')
	call switch_off (S, 'liquid_water')
endif

! Define the output variables

call rads_def_var (S, P, S%sel)

! Fill all the data fields in order as specified above.

call put_var (time(1:n) * 1d0 ,.true.)
call put_var (lat(1:n) * 1d-6)
call put_var (lon(1:n) * 1d-6)
call put_var (alt_reaper(1:n) * 1d-3)
call put_var (alt_rate(1:n) * 1d-3)
call put_var (range_ku(1:n) * 1d-3)
call put_var (dry_tropo_ecmwf(1:n) * 1d-3)
call put_var (wet_tropo_rad(1:n) * 1d-3)
call put_var (wet_tropo_ecmwf(1:n) * 1d-3)
call put_var (iono_gim(1:n) * 1d-3)
call put_var (iono_nic09(1:n) * 1d-3)
call put_var (inv_bar_static(1:n) * 1d-3)
call put_var (inv_bar_mog2d(1:n) * 1d-3)
call put_var (tide_solid(1:n) * 1d-3)
call put_var (tide_ocean_fes04(1:n) * 1d-3)
call put_var (tide_ocean_got47(1:n) * 1d-3)
call put_var (tide_load_fes04(1:n) * 1d-3)
call put_var (tide_load_got47(1:n) * 1d-3)
call put_var (tide_pole(1:n) * 1d-3)
call put_var (ssb_bm3(1:n) * 1d-3)
call put_var (mss_cls01(1:n) * 1d-3)
call put_var (geoid_egm2008(1:n) * 1d-3)
call put_var (mss_ucl04(1:n) * 1d-3)
call put_var (swh_ku(1:n) * 1d-3)
call put_var (sig0_ku(1:n) * 1d-2)
call put_var (wind_speed_ecmwf_u(1:n) * 1d-3)
call put_var (wind_speed_ecmwf_v(1:n) * 1d-3)
call put_var (range_rms_ku(1:n) * 1d-3)
call put_var (range_numval_ku(1:n) * 1d0 )
call put_var (topo_macess(1:n) * 1d-3)
call put_var (tb_238(1:n) * 1d-2)
call put_var (tb_365(1:n) * 1d-2)
call put_var (peakiness_ku(1:n) * 1d-3)
call put_var (flags(1:n) * 1d0 )
call put_var (swh_rms_ku(1:n) * 1d-3)
call put_var (sig0_rms_ku(1:n) * 1d-2)
call put_var (off_nadir_angle2_wf_ku(1:n) * 1d-4)
call put_var (liquid_water(1:n) * 1d-2)
call put_var (water_vapor_content(1:n) * 1d-2)
call put_var (tide_equil(1:n) * 1d-3)
call put_var (tide_non_equil(1:n) * 1d-3)
if (.not.alt_2m) then
	call put_var (wind_speed_alt(1:n) * 1d-3)
	call put_var (drange_cal(1:n) * 1d-3)
	call put_var (drange_fm(1:n) * 1d-3)
	call put_var (dsig0_atmos_ku(1:n) * 1d-2)
	call put_var (mqe(1:n) * 1d-4)
endif

! Write to the data file

552 format ('...',i5,' records written to ',a)
write (*,552) n,trim(P%filename(len_trim(S%dataroot)+2:))
call rads_close_pass (S, P)

end subroutine write_data

subroutine switch_off (S, name)
! Set variable name to 'noedit'
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: name
integer :: i
do i = 1, S%nsel
	if (S%sel(i)%name == name) then
		S%sel(i)%noedit = .true.
		exit
	endif
enddo
end subroutine switch_off

subroutine put_var (data, first)
! Write variables one after the other to the output file
real(eightbytereal), intent(in) :: data(:)
logical, optional, intent(in) :: first
integer :: start(1) = (/ 1 /), i
save i
if (present(first)) then
	if (first) i = 0
endif
i = i + 1
call rads_put_var (S, P, S%sel(i), data, start)
end subroutine put_var

subroutine get_var_1d (varnm, array)
character(*), intent(in) :: varnm
real(eightbytereal), intent(out) :: array(:)
real(eightbytereal) :: array2(mrec)
integer(fourbyteint) :: i0,i1,l,varid,constant
i1 = 0
l = len_trim(varnm)
do
	if (i1 > l) exit
	i0 = i1
	i1 = scan(varnm(i0+1:), '+-') + i0
	if (i1 == i0) i1 = l + 1
	if (nf90_inq_varid(ncid,varnm(i0+1:i1-1),varid) /= nf90_noerr) then
		write (*,'("No such variable: ",a)') varnm(i0+1:i1-1)
		return
	endif
	if (i0 == 0) then
		call nfs(nf90_get_var(ncid,varid,array(1:nrec)))
	else
		call nfs(nf90_get_var(ncid,varid,array2(1:nrec)))
		constant = 0
		if (varnm(i0:i0) == '-') constant = -1
		if (varnm(i0:i0) == '+') constant = 1
		array(1:nrec) = array(1:nrec) + constant * array2(1:nrec)
	endif
enddo
end subroutine get_var_1d

subroutine get_var_2d (varnm, array)
character(*), intent(in) :: varnm
real(eightbytereal), intent(out) :: array(:,:)
real(eightbytereal) :: array2(20,mrec)
integer(fourbyteint) :: i0,i1,l,varid,constant
i1 = 0
l = len_trim(varnm)
do
	if (i1 > l) exit
	i0 = i1
	i1 = scan(varnm(i0+1:), '+-') + i0
	if (i1 == i0) i1 = l + 1
	if (nf90_inq_varid(ncid,varnm(i0+1:i1-1),varid) /= nf90_noerr) then
		write (*,'("No such variable: ",a)') varnm(i0+1:i1-1)
		return
	endif
	if (i0 == 0) then
		call nfs(nf90_get_var(ncid,varid,array(:,1:nrec)))
	else
		call nfs(nf90_get_var(ncid,varid,array2(:,1:nrec)))
		constant = 0
		if (varnm(i0:i0) == '-') constant = -1
		if (varnm(i0:i0) == '+') constant = 1
		array(:,1:nrec) = array(:,1:nrec) + constant * array2(:,1:nrec)
	endif
enddo
end subroutine get_var_2d

subroutine mean_1hz (y, mean, rms)
real(eightbytereal), intent(in) :: y(:,:)
real(eightbytereal), intent(out) :: mean(:), rms(:)
integer(fourbyteint) :: i, j, n
do j = 1,nrec
	mean(j) = 0d0
	rms(j) = 0d0
	n = 20
	do i = 1,n
		mean(j) = mean(j) + y(i,j)
		rms(j) = rms(j) + y(i,j)**2
	enddo
	mean(j) = mean(j) / n
	rms(j) = sqrt ((rms(j) - n * mean(j)**2) / (n - 1))
enddo
end subroutine mean_1hz

subroutine flag_set (a, flags, bit)
logical, intent(in) :: a(:)
integer(twobyteint), intent(inout) :: flags(:)
integer(fourbyteint), intent(in) :: bit
integer(fourbyteint) :: i
integer(twobyteint) :: j
if (size(a) /= size(flags)) stop "Error in flag_set"
j = int(bit,twobyteint)
do i = 1,size(a)
	if (a(i)) flags(i) = ibset(flags(i),j)
enddo
end subroutine flag_set

subroutine invalidate (a, b)
logical, intent(in) :: a(:)
real(eightbytereal), intent(inout) :: b(:)
if (size(a) /= size(b)) stop "Error in invalidate"
where (a) b = nan
end subroutine invalidate

end program rads_gen_reaper
