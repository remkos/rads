!-----------------------------------------------------------------------
! Copyright (c) 2011-2023  Remko Scharroo and Eric Leuliette
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

!*rads_gen_tp_gdrf -- Converts TOPEX/Poseidon GDR-F data to RADS
!+
program rads_gen_tp_gdrf

! This program reads Jason-1/2/3 (O/I)GDR files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/SS/F/SSpPPPPcCCC.nc.
!    SS = satellite abbreviation
!     F = phase (a)
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_tp_gdrf [options] < list_of_file_names
!
! where [options] include:
!  --min-rec <min_rec> : Specify minimum number of records per pass to process.
!
! This program handles TOPEX/Poseidon GDR-F Level 2 products in NetCDF format.
! The format is described in:
!
! [1] TOPEX/POSEIDON GDR-F Products Handbook, SALP-MU-MAO-OP-17607-CN (CNES)
!     Version 1.0, 16 Jun 2023
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdrf - Orbital altitude
! alt_rate - Orbital altitude rate
! range_* - Ocean range (retracked)
! range_rms_* - Std dev of range
! range_numval_* - Nr of averaged range measurements
! qual_range - Quality of range measurement
! swh_* - Significant wave height
! swh_rms_* - Std dev of SWH
! qual_swh - Quality of SWH measurement
! sig0_* - Sigma0
! sig0_rms_* - Std dev of sigma0
! qual_sig0 - Quality of sigma0 measurement
! dsig0_atmos_* - Atmospheric attenuation of sigma0
! wind_speed_alt - Altimeter wind speed
! wind_speed_alt_mle3 - Altimeter wind speed (MLE3)
! wind_speed_rad - Radiometer wind speed
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! qual_alt_rain_ice - Altimeter rain flag
! qual_rad_rain_ice - Radiometer rain flag
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! off_nadir_angle2_wf_rms_ku - RMS of mispointing from waveform squared
! qaul_attitude - Quality of attitude measurement
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! wet_tropo_rad - Radiometer wet tropo correction
! iono_alt - Dual-frequency ionospheric correction
! iono_alt_smooth - Filtered dual-frequency ionospheric correction
! iono_alt_mle3 - Dual-frequency ionospheric correction (MLE3)
! iono_alt_smooth_mle3 - Filtered dual-frequency ionospheric correction (MLE3)
! qual_iono_alt - Quality of dual-frequency ionosphere correction
! iono_gim - GIM ionosphetic correction
! ssb_cls_* - SSB
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! tide_ocean/load_got410 - GOT4.10c ocean and load tide
! tide_ocean/load_fes14 - FES2014 ocean and load tide
! tide_non_equal - Long-period non-equilibrium tide
! tide_solid - Solid earth tide
! tide_pole - Pole tide
! geoid_egm2008 - EGM2008 geoid
! cnescls15 - CNES/CLS15 mean sea surface
! mss_dtu18 - DTU13 mean sea surface
! topo_ace2 - ACE2 topography
! surface_class - Surgace classification
! surface_type_rad - Radiometer surface type
! dist_coast - Distance to the coast
! angle_coast - Angle to the coast
! rads_dist_coast - Radiometer distance to the coast
! tb_187 - Brightness temperature (18.7 GHz)
! tb_238 - Brightness temperature (23.8 GHz)
! tb_365 - Brightness temperature (36.5 GHz)
! qual_rad_tb - Quality of brightness temperatures
! rad_liquid_water - Liquid water content
! rad_water_vapor - Water vapor content
! mean_wave_period - Mean wave period (t02)
! mean_wave_direction - Mean wave direction
! ssha - Sea surface height anomaly
! ssha_mle3 - Sea surface height anomaly (MLE3)
! latency - Latency (NRT, STC, NTC)
! flags - Engineering flags
!
! Extensions _* are:
! _ku:      Ku-band
! _ku_mle3  Ku-band MLE3
! _c:       C-band
!-----------------------------------------------------------------------
use rads
use rads_devel
use rads_devel_netcdf
use rads_gen
use rads_misc
use rads_netcdf
use rads_time
use rads_geo
use netcdf

! Command line arguments

integer(fourbyteint) :: ios
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: ncid, cyclenr, passnr, varid
logical :: tx = .true.
real(eightbytereal) :: equator_time

! Data variables

integer(twobyteint), allocatable :: flags_mle3(:), flags_save(:)
character(len=2) :: mss_cnescls_ver = '15', mss_dtu_ver = '18', tide_fes_ver = '14'
character(len=3) :: tide_got_ver = '410'

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000

! Initialise

call synopsis
call rads_gen_getopt ('', ' min-rec:')
call synopsis ('--head')
if (sat == '') call rads_exit ('Need to specify -S')
call rads_init (S, sat)

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

do
	read (*,'(a)',iostat=ios) infile
	if (ios /= 0) exit

! Open input file

	call log_string (basename(infile))
	if (nft(nf90_open(infile,nf90_nowrite,ncid))) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Check if input is a Jason GDR-F Level 2 data set

	if (index(infile,'TP_GPN_') .eq. 0) then
		call log_string ('Error: this is not TOPEX/Poseidon GDR-F', .true.)
		call nfs(nf90_close(ncid))
		cycle
	endif

	call nfs(nf90_get_att(ncid,nf90_global,'title',arg))
	if (arg /= 'TOPEX/POSEIDON GDR-F') then
		call log_string ('Error: this is not TOPEX/Poseidon GDR-F', .true.)
		call nfs(nf90_close(ncid))
		cycle
	endif

! Get the mission name and initialise RADS (if not done before)

	call nfs(nf90_get_att(ncid,nf90_global,'altimeter_sensor_name',arg))
	if (sat(:2) == 'tx' .and. arg(:5) == 'TOPEX') then
		tx = .true.
	else if (sat(:2) == 'pn' .and. arg == 'POSEIDON') then
		tx = .false.
	else
		call nfs(nf90_close(ncid))
		call log_string ('Skipped', .true.)
		cycle
	endif

! Read global attributes

	call nfs(nf90_inq_dimid(ncid, 'time', varid))
	call nfs(nf90_inquire_dimension(ncid, varid, len=nrec))
	if (nrec == 0) then
		call log_string ('Error: file skipped: no measurements', .true.)
		cycle
	else if (nrec < min_rec) then
		call log_string ('Warning: file skipped: too few measurements', .true.)
		cycle
	else if (nrec > mrec) then
		call log_string ('Error: file skipped: too many measurements', .true.)
		cycle
	endif

! Read cycle, pass, equator time

	call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cyclenr))
	call nfs(nf90_get_att(ncid,nf90_global,'pass_number',passnr))
	call nfs(nf90_get_att(ncid,nf90_global,'equator_time',arg))
	equator_time = strp1985f (arg)

! Skip passes of which the cycle number or equator crossing time is outside the specified interval

	if (equator_time < times(1) .or. equator_time > times(2) .or. cyclenr < cycles(1) .or. cyclenr > cycles(2)) then
		call nfs(nf90_close(ncid))
		call log_string ('Skipped', .true.)
		cycle
	endif

! Set mission phase based on equator_time

	call rads_set_phase (S, equator_time)

! Store relevant info

	call rads_init_pass_struct (S, P)
	P%cycle = cyclenr
	P%pass = passnr
	P%equator_time = equator_time
	call nfs(nf90_get_att(ncid,nf90_global,'equator_longitude',P%equator_lon))
	call nfs(nf90_get_att(ncid,nf90_global,'time_coverage_start',arg))
	P%start_time = strp1985f(arg)
	call nfs(nf90_get_att(ncid,nf90_global,'time_coverage_end',arg))
	P%end_time = strp1985f(arg)

! Store input file name

	P%original = trim(basename(infile))

! Allocate variables

	allocate (a(nrec),flags(nrec),flags_mle3(nrec),flags_save(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
	call nc2f (ncid, 'alt_state_flag_oper', 0, eq=1)			! bit  0: Altimeter Side A/B
	call nc2f (ncid, 'off_nadir_angle_wf_ku_qual', 1)			! bit  1: Quality off-nadir pointing
	call nc2f (ncid, 'surface_classification_flag', 2, eq=4)	! bit  2: Continental ice
	if (tx) call nc2f (ncid, 'range_c_qual', 3)					! bit  3: Quality dual-frequency iono
	call nc2f (ncid, 'surface_classification_flag', 4, eq=1)
	call nc2f (ncid, 'surface_classification_flag', 4, ge=3)	! bit  4: Water/land
	call nc2f (ncid, 'surface_classification_flag', 5, ge=1)	! bit  5: Ocean/other
	call nc2f (ncid, 'rad_surface_type_flag', 6, ge=2)			! bit  6: Radiometer land flag
	if (tx) call nc2f (ncid, 'rain_flag', 7)
	call nc2f (ncid, 'ice_flag', 7)								! bit  7: Altimeter rain or ice flag
	call nc2f (ncid, 'rad_rain_flag', 8)
	call nc2f (ncid, 'rad_sea_ice_flag', 8)						! bit  8: Radiometer rain or ice flag
	call nc2f (ncid, 'rad_tb_18_qual',  9)
	call nc2f (ncid, 'rad_tb_21_qual',  9)						! bit  9: Quality 18.7 or 23.8 GHz channel
	call nc2f (ncid, 'rad_tb_37_qual', 10)						! bit 10: Quality 34.0 GHz channel

! Now do specifics for MLE3

	if (tx) then
		flags_save = flags	! Keep flags for later
		call nc2f (ncid, 'range_ku_mle3_qual', 11)					! bit 11: Quality range
		call nc2f (ncid, 'swh_ku_mle3_qual', 12)					! bit 12: Quality SWH
		call nc2f (ncid, 'sig0_ku_mle3_qual', 13)			! bit 13: Quality Sigma0
		flags_mle3 = flags	! Copy result for MLE3

! Redo the last ones for standard retracker

		flags = flags_save	! Reset to general flags
		call nc2f (ncid, 'range_ku_qual', 11)						! bit 11: Quality range
		call nc2f (ncid, 'swh_ku_qual', 12)							! bit 12: Quality SWH
		call nc2f (ncid, 'sig0_ku_qual', 13)				! bit 13: Quality Sigma0
endif

! Time and location

	call get_var (ncid, 'time', a)
	call new_var ('time', a + sec2000)
	call cpy_var (ncid, 'latitude', 'lat')
	call cpy_var (ncid, 'longitude','lon')
	call cpy_var (ncid, 'altitude delta_ellipsoid_tp_wgs84 SUB', 'alt_std1808')
	call cpy_var (ncid, 'altitude_cnes delta_ellipsoid_tp_wgs84 SUB', 'alt_gdrf')
	call cpy_var (ncid, 'altitude_rate_mean_sea_surface', 'alt_rate')

! Range

	if (tx) then
		call cpy_var (ncid, 'range_ku', 'range_ku')
		call cpy_var (ncid, 'range_rms_ku', 'range_rms_ku')
		call cpy_var (ncid, 'range_numval_ku', 'range_numval')
		call cpy_var (ncid, 'range_ku_qual', 'qual_range')
		call cpy_var (ncid, 'range_ku_mle3', 'range_ku_mle3')
		call cpy_var (ncid, 'range_rms_ku_mle3', 'range_rms_ku_mle3')
		call cpy_var (ncid, 'range_numval_ku_mle3', 'range_numval_ku_mle3')
		call cpy_var (ncid, 'range_ku_mle3_qual', 'qual_range_mle3')
		call cpy_var (ncid, 'range_c', 'range_c')
		call cpy_var (ncid, 'range_rms_c', 'range_rms_c')
		call cpy_var (ncid, 'range_numval_c', 'range_numval_c')
!		call cpy_var (ncid, 'range_c_qual', 'qual_range_c')
		call cpy_var (ncid, 'net_instr_cor_range_ku', 'drange_ku')
		call cpy_var (ncid, 'net_instr_cor_range_ku_mle3', 'drange_ku_mle3')
		call cpy_var (ncid, 'net_instr_cor_range_c', 'drange_c')
		call cpy_var (ncid, 'cg_to_altimeter_timevarying_offset', 'drange_cg')
		call cpy_var (ncid, 'osc_drift_cor', 'drange_uso')
		call cpy_var (ncid, 'range_cor_doppler_ku', 'drange_fm')
	else
		call cpy_var (ncid, 'range_ku_mgdr', 'range_ku')
		call cpy_var (ncid, 'range_rms_ku_mgdr', 'range_rms_ku')
		call cpy_var (ncid, 'range_numval_ku_mgdr', 'range_numval')
		call cpy_var (ncid, 'net_instr_cor_range_ku_mgdr', 'drange_ku')
	endif

! SWH

	if (tx) then
		call cpy_var (ncid, 'swh_ku', 'swh_ku')
		call cpy_var (ncid, 'swh_rms_ku', 'swh_rms_ku')
		call cpy_var (ncid, 'swh_ku_qual', 'qual_swh')
		call cpy_var (ncid, 'swh_ku_mle3', 'swh_ku_mle3')
		call cpy_var (ncid, 'swh_rms_ku_mle3', 'swh_rms_ku_mle3')
		call cpy_var (ncid, 'swh_c', 'swh_c')
		call cpy_var (ncid, 'swh_rms_c', 'swh_rms_c')
	else
		call cpy_var (ncid, 'swh_ku_mgdr', 'swh_ku')
	endif

! Backscatter

	if (tx) then
		call cpy_var (ncid, 'sig0_ku', 'sig0_ku')
		call cpy_var (ncid, 'sig0_rms_ku', 'sig0_rms_ku')
		call cpy_var (ncid, 'sig0_ku_qual', 'qual_sig0')
		call cpy_var (ncid, 'sig0_ku_mle3', 'sig0_ku_mle3')
		call cpy_var (ncid, 'sig0_rms_ku_mle3', 'sig0_rms_ku_mle3')
		call cpy_var (ncid, 'rad_atm_cor_sig0_ku', 'dsig0_atmos_ku')
		call cpy_var (ncid, 'net_instr_cor_sig0_ku', 'dsig0_ku')
		call cpy_var (ncid, 'net_instr_cor_sig0_ku_mle3', 'dsig0_ku_mle3')
		call cpy_var (ncid, 'sig0_c', 'sig0_c')
		call cpy_var (ncid, 'sig0_rms_c', 'sig0_rms_c')
		call cpy_var (ncid, 'rad_atm_cor_sig0_c', 'dsig0_atmos_c')
		call cpy_var (ncid, 'net_instr_cor_sig0_c', 'dsig0_c')
	else
		call cpy_var (ncid, 'sig0_ku_mgdr', 'sig0_ku')
		call cpy_var (ncid, 'sig0_rms_ku_mgdr', 'sig0_rms_ku')
		call cpy_var (ncid, 'rad_atm_cor_sig0_ku', 'dsig0_atmos_ku')
		call cpy_var (ncid, 'net_instr_cor_sig0_ku_mgdr', 'dsig0_ku')
	endif

! Wind speed

	if (tx) then
		call cpy_var (ncid, 'wind_speed_alt', 'wind_speed_alt')
		call cpy_var (ncid, 'wind_speed_alt_mle3', 'wind_speed_alt_mle3')
	else
		call cpy_var (ncid, 'wind_speed_alt_mgdr', 'wind_speed_alt')
	endif
	call cpy_var (ncid, 'rad_wind_speed', 'wind_speed_rad')
	call cpy_var (ncid, 'wind_speed_mod_u', 'wind_speed_ecmwf_u')
	call cpy_var (ncid, 'wind_speed_mod_v', 'wind_speed_ecmwf_v')

! Rain or ice

	if (tx) then
		call cpy_var (ncid, 'rain_flag ice_flag IOR', 'qual_alt_rain_ice')
	else
		call cpy_var (ncid, 'ice_flag', 'qual_alt_rain_ice')
	endif
	call cpy_var (ncid, 'rad_rain_flag rad_sea_ice_flag IOR', 'qual_rad_rain_ice')

! Off-nadir angle

	if (tx) then
		call cpy_var (ncid, 'off_nadir_angle_wf_ku_smoothed', 'off_nadir_angle2_wf_ku')
!		call cpy_var (ncid, 'off_nadir_angle_wf_ku', 'off_nadir_angle2_wf_ku')
		call cpy_var (ncid, 'off_nadir_angle_wf_rms_ku', 'off_nadir_angle2_wf_rms_ku')
		call cpy_var (ncid, 'off_nadir_angle_wf_ku_qual', 'qual_attitude')
	else
		call cpy_var (ncid, 'off_nadir_angle_wf_ku_mgdr', 'off_nadir_angle2_wf_ku')
	endif

! Path delay

	call cpy_var (ncid, 'model_dry_tropo_cor_measurement_altitude', 'dry_tropo_ecmwf')
	call cpy_var (ncid, 'model_wet_tropo_cor_zero_altitude', 'wet_tropo_ecmwf')
	call cpy_var (ncid, 'rad_wet_tropo_cor', 'wet_tropo_rad')
	call cpy_var (ncid, 'composite_wet_tropo_gpd', 'gpd_wet_tropo_cor')
	if (tx) then
		call cpy_var (ncid, 'iono_cor_alt_ku', 'iono_alt')
		call cpy_var (ncid, 'iono_cor_alt_ku_mle3', 'iono_alt_mle3')
		call cpy_var (ncid, 'range_c_qual', 'qual_iono_alt')
	else
		call cpy_var (ncid, 'iono_cor_doris', 'iono_doris')
	endif
	call cpy_var (ncid, 'iono_cor_gim_ku', 'iono_gim')

! SSB

	if (tx) then
		call cpy_var (ncid, 'sea_state_bias_ku', 'ssb_tx')
		call cpy_var (ncid, 'sea_state_bias_ku_mle3', 'ssb_tx_mle3')
		call cpy_var (ncid, 'sea_state_bias_c', 'ssb_tx_c')
	else
		call cpy_var (ncid, 'sea_state_bias_ku_mgdr', 'ssb_bm3')
	endif

! IB

	call cpy_var (ncid, 'inv_bar_cor', 'inv_bar_static')
	call cpy_var (ncid, 'dac', 'inv_bar_mog2d')

! Tides

	call cpy_var (ncid, 'ocean_tide_got load_tide_got SUB', 'tide_ocean_got' // tide_got_ver)
	call cpy_var (ncid, 'ocean_tide_fes load_tide_fes SUB ocean_tide_non_eq ADD', 'tide_ocean_fes' // tide_fes_ver)
	call cpy_var (ncid, 'load_tide_got', 'tide_load_got' // tide_got_ver)
	call cpy_var (ncid, 'load_tide_fes', 'tide_load_fes' // tide_fes_ver)
	call cpy_var (ncid, 'ocean_tide_eq', 'tide_equil')
	call cpy_var (ncid, 'ocean_tide_non_eq', 'tide_non_equil')
	call cpy_var (ncid, 'internal_tide_hret', 'tide_internal')
	call cpy_var (ncid, 'solid_earth_tide', 'tide_solid')
	call cpy_var (ncid, 'pole_tide', 'tide_pole')

! Geoid and MSS

	call cpy_var (ncid, 'geoid delta_ellipsoid_tp_wgs84 SUB', 'geoid_egm2008')
	call cpy_var (ncid, 'mean_sea_surface_cnescls delta_ellipsoid_tp_wgs84 SUB', 'mss_cnescls' // mss_cnescls_ver)
	call cpy_var (ncid, 'mean_sea_surface_dtu delta_ellipsoid_tp_wgs84 SUB', 'mss_dtu' // mss_dtu_ver)

! Surface type and coastal proximity

	call cpy_var (ncid, 'depth_or_elevation', 'topo_ace2')
	call cpy_var (ncid, 'surface_classification_flag', 'surface_class')
	call get_var (ncid, 'rad_surface_type_flag', a)
	where (a > 0) a = a + 1
	call new_var ('surface_type_rad', a)
	call cpy_var (ncid, 'distance_to_coast 1e-3 MUL', 'dist_coast') ! Convert m to km
	call cpy_var (ncid, 'rad_distance_to_land 1e-3 MUL', 'rad_dist_coast') ! Convert m to km

! Other flags

	call cpy_var (ncid, 'alt_state_flag_oper', 'flag_alt_oper_mode')

! Bit flags

	call new_var ('flags', dble(flags))
	if (tx) call new_var ('flags_mle3', dble(flags_mle3))

! Other radiometer measurements
! Selected smoothed TBs

	call cpy_var (ncid, 'rad_tb_18', 'tb_180')
	call cpy_var (ncid, 'rad_tb_21', 'tb_210')
	call cpy_var (ncid, 'rad_tb_37', 'tb_370')
	call cpy_var (ncid, 'rad_tb_18_qual 2 MUL rad_tb_21_qual ADD 2 MUL rad_tb_37_qual ADD', 'qual_rad_tb')
	call cpy_var (ncid, 'rad_cloud_liquid_water', 'liquid_water_rad')
	call cpy_var (ncid, 'rad_water_vapor', 'water_vapor_rad')

! Wave model

	call cpy_var (ncid, 'swh_model', 'swh_mfwam')
	call cpy_var (ncid, 'mean_wave_period_t02', 'mean_wave_period')

! SSHA

	if (tx) then
		call cpy_var (ncid, 'ssha', 'ssha')
		call cpy_var (ncid, 'ssha_mle3', 'ssha_mle3')
	else
		call cpy_var (ncid, 'ssha_mgdr', 'ssha')
	endif

! Dump the data

	call nfs(nf90_close(ncid))
	call put_rads
	deallocate (a, flags, flags_mle3, flags_save)

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Write TOPEX/Poseidon GDR-F data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_file_names')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --min-rec=MIN_REC         Specify minimum number of records per pass to process' // &
'This program converts TOPEX/Poseidon GDR-F Level 2 products to RADS data' / &
'files with the name $RADSDATAROOT/data/SS/F/pPPPP/SSpPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_tp_gdrf
