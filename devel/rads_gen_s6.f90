!-----------------------------------------------------------------------
! Copyright (c) 2011-2021  Remko Scharroo and Eric Leuliette
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

!*rads_gen_s6 -- Converts Sentinel-6 data to RADS
!+
program rads_gen_s6

! This program reads Sentinel-6 NRT/STC/NTC files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/SS/F/SSpPPPPcCCC.nc.
!    SS = satellite (6a or 6b)
!     F = phase (a)
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_s6 [options] < list_of_Sentinel6_file_names
!
! where [options] include:
!  --min-rec <min_rec> : Specify minimum number of records per pass to process.
!
! This program handles Sentinel-6 Level 2 products in NetCDF format.
! The format is described in:
!!
! [1] Sentinel-6/Jason-CS Level 2 Product Format Specification (L2 ALT PFS)
!     EUM/LEO-JASCS/SPE/17/901187, v4B
!
! [2] Sentinel-6/Jason-CS Generic Product Format Specification (GPFS)
!     EUM/LEO-JASCS/SPE/17/897975, v4B
!
! [3] Sentinel-6/Jason-CS Level 2 NetCDF Dump
!     EUM/LEO-JASCS/SPE/17/957846, v4B
!
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
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! wet_tropo_rad - Radiometer wet tropo correction
! iono_alt - Dual-frequency ionospheric correction
! iono_alt_smooth - Filtered dual-frequency ionospheric correction
! iono_gim - GIM ionosphetic correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! ssb_cls_* - SSB
! swh_* - Significant wave height
! swh_rms_* - Std dev of SWH
! sig0_* - Sigma0
! sig0_rms_* - Std dev of sigma0
! dsig0_atmos_* - Atmospheric attenuation of sigma0
! wind_speed_alt - Altimeter wind speed
! wind_speed_rad - Radiometer wind speed
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! tide_ocean/load_got410 - GOT4.10c ocean and load tide
! tide_ocean/load_fes14 - FES2014 ocean and load tide
! tide_pole - Pole tide
! tide_solid - Solid earth tide
! topo_ace2 - ACE2 topography
! geoid_egm2008 - EGM2008 geoid
! cnescls15 - CNES/CLS15 mean sea surface
! mss_dtu18 - DTU18 mean sea surface
! tb_238 - Brightness temperature (23.8 GHz)
! tb_365 - Brightness temperature (36.5 GHz)
! flags - Engineering flags
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! rad_liquid_water - Liquid water content
! rad_water_vapor - Water vapor content
! ssha - Sea surface height anomaly
!
! Extensions _* are:
! _ku:      Ku-band
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

integer(fourbyteint) :: ios, i
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: ncid, ncid1, ncidk, ncidc, ncid20, ncid20k, ncid20c, cyclenr, passnr, varid, nrec20, &
	cyclenr_orig, passnr_orig
real(eightbytereal) :: equator_time

! Data variables

integer(twobyteint), allocatable :: flags_mle3(:), flags_save(:)
character(len=2) :: mss_cnescls_ver = '15', mss_dtu_ver = '18', tide_fes_ver = '14'
character(len=3) :: tide_got_ver = '410', env
character(len=8) :: chd_ver, cha_ver, cnf_ver
integer :: latency = rads_nrt, nsat
logical :: lr, lcal1 = .false.

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000
real(eightbytereal), allocatable :: b(:), telemetry_type_flag(:)

! Initialise

call synopsis
call rads_gen_getopt ('', ' cal1 min-rec:')
call synopsis ('--head')
call rads_set_options (' cal1 min-rec:')
if (sat /= '') call rads_init (S, sat)
call rads_load_options (nsat)

do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('cal1')
		lcal1 = .true.
	end select
enddo

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

do
	read (*,'(a)',iostat=ios) infile
	if (ios /= 0) exit

! Open input file

	call log_string (infile)
	if (nft(nf90_open(infile,nf90_nowrite,ncid))) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Check if input is a Sentinel-6 Level 2 data set

	if (nft(nf90_get_att(ncid,nf90_global,'product_name',arg)) .or. &
		arg(5:10) /= 'P4_2__') then
		call log_string ('Error: this is not a Sentinel-6 Level 2 data set', .true.)
		cycle
	endif

! Determine latency (NRT, STC, NTC) and resolution (LR, HR)

	if (index(arg, '_NR_') > 0) then
		latency = rads_nrt
	else if (index(arg, '_ST_') > 0) then
		latency = rads_stc
	else if (index(arg, '_NT_') > 0) then
		latency = rads_ntc
	else
		call log_string ('Error: file skipped: unknown latency', .true.)
		cycle
	endif
	lr = (index(arg, '_LR_') > 0)

! Get the mission name and initialise RADS (if not done before)

	select case (arg(:3))
	case ('S6A')
		arg = '6a'
	case ('S6B')
		arg = '6b'
	case default
		call log_string ('Error: wrong misson: must be S6A or S6B', .true.)
		cycle
	end select
	if (sat == '') then
		call rads_init (S, arg)
	else if (S%sat /= arg(:2)) then
		call log_string ('Error: wrong misson: must be S'//S%sat, .true.)
		cycle
	endif

! Get NetCDF ID for 1-Hz data

	call nfs(nf90_inq_ncid(ncid, 'data_01', ncid1))
	call nfs(nf90_inq_ncid(ncid1, 'ku', ncidk))
	if (lr) call nfs(nf90_inq_ncid(ncid1, 'c', ncidc))

! Read global attributes

	call nfs(nf90_inq_dimid(ncid1, 'time', varid))
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

	call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cyclenr_orig))
	call nfs(nf90_get_att(ncid,nf90_global,'pass_number',passnr_orig))
	call nfs(nf90_get_att(ncid,nf90_global,'equator_time',arg))
	equator_time = nint(strp1985f (arg) * 1d3, eightbyteint) * 1d-3 ! Round nicely to nearest millisecond

! Recompute the cycle and pass number since they cannot fully be relied on [AR 1751]

	call rads_time_to_cycle_pass (S, equator_time, cyclenr, passnr)
	if (cyclenr /= cyclenr_orig .or. passnr /= passnr_orig) call log_string('Warning: cycle/pass nr adjusted')

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

! Determine L2 processing baseline, characterisation and configuration versions

	call nfs(nf90_get_att(ncid,nf90_global,'xref_altimeter_characterization',arg))
	chd_ver = arg(11:14) // ' ' // arg(96:98)
	call nfs(nf90_get_att(ncid,nf90_global,'xref_altimeter_characterization_array',arg))
	cha_ver = arg(11:14) // ' ' // arg(96:98)
	call nfs(nf90_get_att(ncid,nf90_global,'xref_processor_configuration',arg))
	cnf_ver = arg(11:14) // ' ' // arg(96:98)

	call nfs(nf90_get_att(ncid,nf90_global,'source',arg))

! Store input file name

	i = index(infile, '/', .true.) + 1
	P%original = trim(infile(i:)) // ' (' // trim(arg) // ', ' // chd_ver // ', ' // cha_ver // ', ' // cnf_ver // ')'

! Allocate variables

	allocate (a(nrec),b(nrec),flags(nrec),flags_mle3(nrec),flags_save(nrec))
	nvar = 0

! Get NetCDF ID for 20-Hz data (if available)

	if (nft(nf90_inq_ncid(ncid, 'data_20', ncid20))) then
		ncid20 = 0
		ncid20k = 0
		nrec20 = 0
	else
		call nfs(nf90_inq_ncid(ncid20, 'ku', ncid20k))
		call nfs(nf90_inq_dimid(ncid20k, 'time', varid))
		call nfs(nf90_inquire_dimension(ncid20k, varid, len=nrec20))
		allocate (telemetry_type_flag(nrec20))
	endif

! Compile flag bits

	flags = 0
	if (.not.lr) flags = 1										! bit  0: Altimeter mode (LR = 0, HR = 1)
	if (lr) call nc2f (ncidk, 'off_nadir_angle_wf_ocean_qual', 1)		! bit  1: Quality off-nadir pointing
	call nc2f (ncid1, 'surface_classification_flag', 2, eq=4)	! bit  2: Continental ice
	if (lr) call nc2f (ncidc, 'range_ocean_qual', 3)			! bit  3: Quality dual-frequency iono
	call nc2f (ncid1, 'surface_classification_flag', 4, eq=1)
	call nc2f (ncid1, 'surface_classification_flag', 4, ge=3)	! bit  4: Water/land
	call nc2f (ncid1, 'surface_classification_flag', 5, ge=1)	! bit  5: Ocean/other
	call nc2f (ncid1, 'rad_surface_type_flag', 6, ge=2)			! bit  6: Radiometer land flag
	call nc2f (ncid1, 'rain_flag', 7, eq=1)
	call nc2f (ncid1, 'rain_flag', 7, eq=2)
	call nc2f (ncid1, 'rain_flag', 7, eq=4)						! bit  7: Altimeter rain or ice flag
	call nc2f (ncid1, 'rad_rain_flag', 8)
	call nc2f (ncid1, 'rad_sea_ice_flag', 8)					! bit  8: Radiometer rain or ice flag
	call nc2f (ncid1, 'rad_tmb_187_qual',  9)
	call nc2f (ncid1, 'rad_tmb_238_qual',  9)					! bit  9: Quality 18.7 or 23.8 GHz channel
	call nc2f (ncid1, 'rad_tmb_340_qual', 10)					! bit 10: Quality 34.0 GHz channel
	call nc2f (ncid1, 'orbit_type_flag', 15, le=2)				! bit 15: Quality of orbit

! Now do specifics for MLE3

	if (lr) then
		flags_save = flags	! Keep flags for later
		call nc2f (ncidk, 'range_ocean_mle3_qual', 11)			! bit 11: Quality range
		call nc2f (ncidk, 'swh_ocean_mle3_qual', 12)			! bit 12: Quality SWH
		call nc2f (ncidk, 'sig0_ocean_mle3_qual', 13)			! bit 13: Quality Sigma0
		flags_mle3 = flags	! Copy result for MLE3
		flags = flags_save	! Continue with nominal ocean retracking flags
	endif

! Redo the last ones for standard retracker

	call nc2f (ncidk, 'range_ocean_qual', 11)					! bit 11: Quality range
	call nc2f (ncidk, 'swh_ocean_qual', 12)						! bit 12: Quality SWH
	call nc2f (ncidk, 'sig0_ocean_qual', 13)					! bit 13: Quality Sigma0

! Time and location

	call get_var (ncid1, 'time', a)
	a = a + sec2000
	P%start_time = a(1)  ! Using 1-Hz time stamps, not the ones from the header which are
	P%end_time = a(nrec) ! the time stamps of the 20-Hz measurements
	call new_var ('time', a)
	call cpy_var (ncid1, 'latitude', 'lat')
	call cpy_var (ncid1, 'longitude','lon')
	call cpy_var (ncid1, 'altitude delta_ellipsoid_tp_wgs84 SUB', 'alt_gdrf')
	call cpy_var (ncid1, 'altitude_rate_mean_sea_surface', 'alt_rate')
	call cpy_var (ncid1, 'orbit_type_flag', 'qual_orbit')

! Telemetry type (copied from 20-Hz to 1-Hz)

	if (ncid20k /= 0) then
		call get_var (ncidk, 'index_first_20hz_measurement', a)
		! Note: index_first_20hz_measurement starts with 0, so add 1 in Fortran
		call get_var (ncidk, 'numtotal_20hz_measurement', b)
		call get_var (ncid20k, 'telemetry_type_flag', telemetry_type_flag)
		call new_var ('flag_alt_oper_mode', telemetry_type_flag(nint(a+b/2)) - 1)
		! Note: 1 subtracted because in PGF: 0 = null_record.
	endif

! Range

	call cpy_var (ncidk, 'range_ocean', 'range_ku')
	call cpy_var (ncidk, 'range_ocean_rms', 'range_rms_ku')
	call cpy_var (ncidk, 'range_ocean_numval', 'range_numval')
	call cpy_var (ncidk, 'range_ocean_qual', 'qual_range')
	call cpy_var (ncidk, 'net_instr_cor_range_ocean', 'drange_ku')
	if (lr) then
		call cpy_var (ncidk, 'range_ocean_mle3', 'range_ku_mle3')
		call cpy_var (ncidk, 'range_ocean_mle3_rms', 'range_rms_ku_mle3')
		call cpy_var (ncidk, 'range_ocean_mle3_numval', 'range_numval_ku_mle3')
		call cpy_var (ncidk, 'net_instr_cor_range_ocean_mle3', 'drange_ku_mle3')
		call cpy_var (ncidc, 'range_ocean', 'range_c')
		call cpy_var (ncidc, 'range_ocean_rms', 'range_rms_c')
		call cpy_var (ncidc, 'range_ocean_numval', 'range_numval_c')
		call cpy_var (ncidc, 'net_instr_cor_range_ocean', 'drange_c')
	endif

! SWH

	call cpy_var (ncidk, 'swh_ocean', 'swh_ku')
	call cpy_var (ncidk, 'swh_ocean_rms', 'swh_rms_ku')
	call cpy_var (ncidk, 'swh_ocean_qual', 'qual_swh')
	call cpy_var (ncidk, 'net_instr_cor_swh_ocean', 'dswh_ku')
	if (lr) then
		call cpy_var (ncidk, 'swh_ocean_mle3', 'swh_ku_mle3')
		call cpy_var (ncidk, 'swh_ocean_mle3_rms', 'swh_rms_ku_mle3')
		call cpy_var (ncidk, 'net_instr_cor_swh_ocean_mle3', 'dswh_ku_mle3')
		call cpy_var (ncidc, 'swh_ocean', 'swh_c')
		call cpy_var (ncidc, 'swh_ocean_rms', 'swh_rms_c')
		call cpy_var (ncidc, 'net_instr_cor_swh_ocean', 'dswh_c')
	endif

! Backscatter

	call cpy_var (ncidk, 'sig0_ocean', 'sig0_ku')
	call cpy_var (ncidk, 'sig0_ocean_rms', 'sig0_rms_ku')
	call cpy_var (ncidk, 'sig0_ocean_qual', 'qual_sig0')
	call cpy_var (ncidk, 'net_instr_cor_sig0_ocean', 'dsig0_ku')
	call cpy_var (ncidk, 'atm_cor_sig0', 'dsig0_atmos_ku')
	if (lr) then
		call cpy_var (ncidk, 'sig0_ocean_mle3', 'sig0_ku_mle3')
		call cpy_var (ncidk, 'sig0_ocean_mle3_rms', 'sig0_rms_ku_mle3')
		call cpy_var (ncidk, 'net_instr_cor_sig0_ocean_mle3', 'dsig0_ku_mle3')
		call cpy_var (ncidc, 'sig0_ocean', 'sig0_c')
		call cpy_var (ncidc, 'sig0_ocean_rms', 'sig0_rms_c')
		call cpy_var (ncidc, 'atm_cor_sig0', 'dsig0_atmos_c')
		call cpy_var (ncidc, 'net_instr_cor_sig0_ocean', 'dsig0_c')
	endif
	call cpy_var (ncid1, 'climato_use_flag', 'qual_dsig0_atmos')

! Wind speed

	call cpy_var (ncid1, 'wind_speed_alt', 'wind_speed_alt')
	if (lr) call cpy_var (ncid1, 'wind_speed_alt_mle3', 'wind_speed_alt_mle3')
	call cpy_var (ncid1, 'rad_wind_speed', 'wind_speed_rad')
	call cpy_var (ncid1, 'wind_speed_mod_u', 'wind_speed_ecmwf_u')
	call cpy_var (ncid1, 'wind_speed_mod_v', 'wind_speed_ecmwf_v')

! Rain or ice

	call cpy_var (ncid1, 'rain_flag', 'qual_alt_rain_ice')
	call cpy_var (ncid1, 'rad_rain_flag rad_sea_ice_flag IOR', 'qual_rad_rain_ice')
	call cpy_var (ncid1, 'rain_attenuation', 'dsig0_atten')

! Off-nadir angle

	if (lr) then
		call cpy_var (ncidk, 'off_nadir_angle_wf_ocean', 'off_nadir_angle2_wf_ku')
		call cpy_var (ncidk, 'off_nadir_angle_wf_ocean_rms', 'off_nadir_angle2_wf_rms_ku')
		call cpy_var (ncidk, 'off_nadir_angle_wf_ocean_qual', 'qual_attitude')
	endif
	call cpy_var (ncid1, 'off_nadir_pitch_angle_pf', 'attitude_pitch')
	call cpy_var (ncid1, 'off_nadir_roll_angle_pf', 'attitude_roll')
	call cpy_var (ncid1, 'off_nadir_yaw_angle_pf', 'attitude_yaw')

! Path delay

	call cpy_var (ncid1, 'model_dry_tropo_cor_measurement_altitude', 'dry_tropo_ecmwf')
	call cpy_var (ncid1, 'model_wet_tropo_cor_measurement_altitude', 'wet_tropo_ecmwf')
	call cpy_var (ncid1, 'rad_wet_tropo_cor', 'wet_tropo_rad')
	call cpy_var (ncid1, 'iono_cor_alt', 'iono_alt')
	call cpy_var (ncid1, 'iono_cor_alt_filtered', 'iono_alt_smooth')
	if (lr) then
		call cpy_var (ncid1, 'iono_cor_alt_mle3', 'iono_alt_mle3')
		call cpy_var (ncid1, 'iono_cor_alt_filtered_mle3', 'iono_alt_smooth_mle3')
		call cpy_var (ncidc, 'range_ocean_qual', 'qual_iono_alt')
	endif
	call cpy_var (ncidk, 'iono_cor_gim', 'iono_gim')

! SSB

	call cpy_var (ncidk, 'sea_state_bias', 'ssb_cls')
	if (lr) then
		call cpy_var (ncidk, 'sea_state_bias_mle3', 'ssb_cls_mle3')
		call cpy_var (ncidc, 'sea_state_bias', 'ssb_cls_c')
	endif

! IB

	call cpy_var (ncid1, 'inv_bar_cor', 'inv_bar_static')
	if (latency == rads_nrt) then
		call cpy_var (ncid1, 'inv_bar_cor', 'inv_bar_mog2d')
	else
		call cpy_var (ncid1, 'dac', 'inv_bar_mog2d')
	endif

! Tides

	call cpy_var (ncid1, 'ocean_tide_sol1 load_tide_sol1 SUB', 'tide_ocean_got' // tide_got_ver)
	call cpy_var (ncid1, 'ocean_tide_sol2 load_tide_sol2 SUB ocean_tide_non_eq ADD', 'tide_ocean_fes' // tide_fes_ver)
	call cpy_var (ncid1, 'load_tide_sol1', 'tide_load_got' // tide_got_ver)
	call cpy_var (ncid1, 'load_tide_sol2', 'tide_load_fes' // tide_fes_ver)
	call cpy_var (ncid1, 'ocean_tide_eq', 'tide_equil')
	call cpy_var (ncid1, 'ocean_tide_non_eq', 'tide_non_equil')
	call cpy_var (ncid1, 'internal_tide', 'tide_internal')
	call cpy_var (ncid1, 'solid_earth_tide', 'tide_solid')
	call cpy_var (ncid1, 'pole_tide', 'tide_pole')

! Geoid and MSS

	call cpy_var (ncid1, 'geoid delta_ellipsoid_tp_wgs84 SUB', 'geoid_egm2008')
	call cpy_var (ncid1, 'mean_sea_surface_sol1 delta_ellipsoid_tp_wgs84 SUB', 'mss_cnescls' // mss_cnescls_ver)
	call cpy_var (ncid1, 'mean_sea_surface_sol2 delta_ellipsoid_tp_wgs84 SUB', 'mss_dtu' // mss_dtu_ver)

! Surface type and coastal proximity

	call cpy_var (ncid1, 'depth_or_elevation', 'topo_ace2')
	call cpy_var (ncid1, 'surface_classification_flag', 'surface_class')
	call cpy_var (ncid1, 'rad_surface_type_flag', 'surface_type_rad')
	call cpy_var (ncid1, 'distance_to_coast 1e-3 MUL', 'dist_coast') ! Convert m to km
	call cpy_var (ncid1, 'angle_of_approach_to_coast', 'angle_coast')
	call cpy_var (ncid1, 'rad_distance_to_land 1e-3 MUL', 'rad_dist_coast') ! Convert m to km

! Bit flags

	call new_var ('flags', dble(flags))
	if (lr) call new_var ('flags_mle3', dble(flags_mle3))

! Other radiometer measurements

	call cpy_var (ncid1, 'rad_tmb_187', 'tb_187')
	call cpy_var (ncid1, 'rad_tmb_238', 'tb_238')
	call cpy_var (ncid1, 'rad_tmb_340', 'tb_340')
	call cpy_var (ncid1, 'rad_tmb_340_qual 2 MUL rad_tmb_238_qual ADD 2 MUL rad_tmb_187_qual ADD', 'qual_rad_tb')
	call cpy_var (ncid1, 'rad_cloud_liquid_water', 'liquid_water_rad')
	call cpy_var (ncid1, 'rad_water_vapor', 'water_vapor_rad')

! SSHA

	call cpy_var (ncidk, 'ssha', 'ssha')
	if (lr) call cpy_var (ncidk, 'ssha_mle3', 'ssha_mle3')

! Misc

	a = latency
	call new_var ('latency', a)

! Get the L1B file name so we can look for the calibration values

	call nfs(nf90_get_att(ncid,nf90_global,'xref_altimeter_level1b',arg))
	call nfs(nf90_close(ncid))

! When requested, load the calibration values from the L1B file (Baseline < F04)

	if (lcal1 .and. cnf_ver(6:8) < '011') then
		call log_string(arg)
		if (index(infile, 'REP/') > 0) then
			env = 'REP'
		else if (index(infile, 'RMC/') > 0) then
			env = 'RMC'
		else if (index(infile, 'DEV/') > 0) then
			env = 'DEV'
		else
			env = arg(89:91)
		endif
		call parseenv('${RADSROOT}/ext/' // S%sat // '/' // env // '/' // arg(11:12) // '/' // &
			arg(93:94) // '/L1B/c' // arg(72:74) // '/' // trim(arg) // '/measurement.nc', infile)
		if (nft(nf90_open(infile,nf90_nowrite,ncid))) then
			call log_string ('Error: failed to open L1B file')
		else
			call nfs(nf90_inq_ncid(ncid, 'data_20', ncid20))
			call nfs(nf90_inq_ncid(ncid20, 'ku', ncid20k))
			call get_var (ncid20k, 'range_cor_internal_delay_cal', a(1:1))
			call new_var ('cal1_range_ku', a(1:1))
			call get_var (ncid20k, 'cal1_power', a(1:1))
			call new_var ('cal1_power_ku', a(1:1))
			if (lr) then
				call nfs(nf90_inq_ncid(ncid20, 'c', ncid20c))
				call get_var(ncid20c, 'range_cor_internal_delay_cal', a(1:1))
				call new_var ('cal1_range_c', a(1:1))
				call get_var(ncid20c, 'cal1_power', a(1:1))
				call new_var ('cal1_power_c', a(1:1))
			endif
			call nfs(nf90_close(ncid))
		endif
	endif

! Dump the data

	call put_rads
	deallocate (a, b, flags, flags_mle3, flags_save)
	if (nrec20 /= 0) deallocate (telemetry_type_flag)

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional, intent(in) :: flag
if (rads_version ('Write Sentinel-6 data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_Sentinel6_file_names')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --min-rec=MIN_REC         Specify minimum number of records per pass to process' / &
'  --cal1                    Read the CAL1 values from the appropriate L1B file (L2 CONF < 011)' // &
'This program converts Sentinel-6 NRT/STC/NTC Level 2 products to RADS data' / &
'files with the name $RADSDATAROOT/data/SS/F/pPPPP/s6pPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_s6
