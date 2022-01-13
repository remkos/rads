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

!*rads_gen_jason_gdrf -- Converts Jason GDR-F data to RADS
!+
program rads_gen_jason_gdrf

! This program reads Jason-1/2/3 (O/I)GDR files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/jJ/F/jJpPPPPcCCC.nc.
!    jJ = Jason satellite abbreviation
!     F = phase (a)
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_jason_gdrf [options] < list_of_JASON3_file_names
!
! where [options] include:
!  --min-rec <min_rec> : Specify minimum number of records per pass to process.
!
! This program handles Jason-1, -2 and -3 OGDR, IGDR and GDR-F Level 2 products in NetCDF format.
! The format is described in:
!
! [1] Jason-3 Products Handbook, SALP-MU-M-OP-16118-CN (CNES)
!     Version 2.6, 21 Sept 2020
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

integer(fourbyteint) :: ios, i = 0
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: ncid, ncid1, ncidk, ncidc, ncid20, ncid20k, cyclenr, passnr, varid, nrec20
logical :: mle3 = .true.
real(eightbytereal) :: equator_time

! Data variables

integer(twobyteint), allocatable :: flags_mle3(:), flags_adaptive(:), flags_save(:)
character(len=2) :: mss_cnescls_ver = '15', mss_dtu_ver = '18', tide_fes_ver = '14'
character(len=3) :: tide_got_ver = '410'
integer :: latency = rads_nrt

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000
real(eightbytereal), allocatable :: dh(:)

! Initialise

call synopsis
call rads_gen_getopt ('', ' min-rec:')
call synopsis ('--head')
if (sat /= '') call rads_init (S, sat)

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

	if (index(infile,'_2Pf') .eq. 0) then
		call log_string ('Error: this is not GDR-F', .true.)
		cycle
	endif

! Determine latency (OGDR, IGDR, GDR)
	call nfs(nf90_get_att(ncid,nf90_global,'title',arg))
	if (arg(:4) == 'OGDR') then
		latency = rads_nrt
	else if (arg(:4) == 'IGDR') then
		latency = rads_stc
	else if (arg(:3) == 'GDR') then
		latency = rads_ntc
	else
		call log_string ('Error: file skipped: unknown latency', .true.)
		cycle
	endif

! Get the mission name and initialise RADS (if not done before)

	call nfs(nf90_get_att(ncid,nf90_global,'mission_name',arg))
	mle3 = (arg /= 'Jason-1')
	select case (arg)
	case ('Jason-1')
		arg = 'j1'
	case ('OSTM/Jason-2')
		arg = 'j2'
	case ('Jason-3')
		arg = 'j3'
	case default
		call log_string ('Error: wrong mission_name found in header', .true.)
		cycle
	end select
	if (sat == '') then
		call rads_init (S, arg)
	else if (S%sat /= arg(:2)) then
		call log_string ('Error: wrong mission_name found in header', .true.)
		cycle
	endif

! Get NetCDF ID for 1-Hz data

	call nfs(nf90_inq_ncid(ncid, 'data_01', ncid1))
	call nfs(nf90_inq_ncid(ncid1, 'ku', ncidk))
	call nfs(nf90_inq_ncid(ncid1, 'c', ncidc))

! Read global attributes

	call nfs(nf90_inq_dimid(ncid1, 'time', varid))
	call nfs(nf90_inquire_dimension(ncid1, varid, len=nrec))
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

! If pass_number is 0, get cycle and pass number from the file name

	if (passnr == 0) then
		read (infile(i+17:i+23),'(i3,1x,i3)') cyclenr, passnr
	endif

! Skip passes of which the cycle number or equator crossing time is outside the specified interval

	if (equator_time < times(1) .or. equator_time > times(2) .or. cyclenr < cycles(1) .or. cyclenr > cycles(2)) then
		call nfs(nf90_close(ncid))
		call log_string ('Skipped', .true.)
		cycle
	endif

! Set mission phase based on equator_time

!	call rads_set_phase (S, equator_time)

! Store relevant info

	call rads_init_pass_struct (S, P)
	P%cycle = cyclenr
	P%pass = passnr
	P%equator_time = equator_time
	call nfs(nf90_get_att(ncid,nf90_global,'equator_longitude',P%equator_lon))
	call nfs(nf90_get_att(ncid,nf90_global,'first_meas_time',arg))
	P%start_time = strp1985f(arg)
	call nfs(nf90_get_att(ncid,nf90_global,'last_meas_time',arg))
	P%end_time = strp1985f(arg)

! Determine L2 processing version (currently not used)

	call nfs(nf90_get_att(ncid,nf90_global,'source',arg))

! Store input file name

	P%original = trim(basename(infile)) // ' (' // trim(arg) // ')'

! Allocate variables

	allocate (a(nrec),dh(nrec),flags(nrec),flags_mle3(nrec),flags_adaptive(nrec),flags_save(nrec))
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
	endif

! Compile flag bits

	flags = 0
	call nc2f (ncid1, 'alt_state_oper_flag',0)				! bit  0: Altimeter Side A/B
	call nc2f (ncidk, 'off_nadir_angle_wf_ocean_compression_qual', 1)		! bit  1: Quality off-nadir pointing
	call nc2f (ncid1, 'surface_classification_flag', 2, eq=4)	! bit  2: Continental ice
	call nc2f (ncidc, 'range_ocean_compression_qual', 3)			! bit  3: Quality dual-frequency iono
	call nc2f (ncid1, 'surface_classification_flag', 4, eq=1)
	call nc2f (ncid1, 'surface_classification_flag', 4, ge=3)	! bit  4: Water/land
	call nc2f (ncid1, 'surface_classification_flag', 5, ge=1)	! bit  5: Ocean/other
	call nc2f (ncid1, 'rad_surface_type_flag', 6, ge=2)			! bit  6: Radiometer land flag
	call nc2f (ncid1, 'rain_flag', 7, eq=1)
	call nc2f (ncid1, 'rain_flag', 7, eq=2)
	call nc2f (ncid1, 'rain_flag', 7, eq=4)						! bit  7: Altimeter rain or ice flag
	call nc2f (ncid1, 'rad_rain_flag', 8)
	call nc2f (ncid1, 'rad_sea_ice_flag', 8)					! bit  8: Radiometer rain or ice flag
	call nc2f (ncid1, 'rad_tb_187_qual',  9)
	call nc2f (ncid1, 'rad_tb_238_qual',  9)					! bit  9: Quality 18.7 or 23.8 GHz channel
	call nc2f (ncid1, 'rad_tb_340_qual', 10)					! bit 10: Quality 34.0 GHz channel

! Now do specifics for MLE3

	flags_save = flags	! Keep flags for later
	call nc2f (ncidk, 'range_ocean_mle3_compression_qual', 11)		! bit 11: Quality range
	call nc2f (ncidk, 'swh_ocean_mle3_compression_qual', 12)		! bit 12: Quality SWH
	call nc2f (ncidk, 'sig0_ocean_mle3_compression_qual', 13)		! bit 13: Quality Sigma0
	flags_mle3 = flags	! Copy result for MLE3

! Now do specifics for adaptive retracking

	if (latency == rads_ntc) then
		flags = flags_save	! Reset to general flags
		call nc2f (ncidk, 'range_adaptive_compression_qual', 11)	! bit 11: Quality range
		call nc2f (ncidk, 'swh_adaptive_compression_qual', 12)		! bit 12: Quality SWH
		call nc2f (ncidk, 'sig0_adaptive_compression_qual', 13)		! bit 13: Quality Sigma0
		flags_adaptive = flags	! Copy result for adaptive retracking
	endif

! Redo the last ones for standard retracker

	flags = flags_save	! Reset to general flags
	call nc2f (ncidk, 'range_ocean_compression_qual', 11)			! bit 11: Quality range
	call nc2f (ncidk, 'swh_ocean_compression_qual', 12)				! bit 12: Quality SWH
	call nc2f (ncidk, 'sig0_ocean_compression_qual', 13)			! bit 13: Quality Sigma0

! Time and location

	call get_var (ncid1, 'time', a)
	call new_var ('time', a + sec2000)
	call cpy_var (ncid1, 'latitude', 'lat')
	call cpy_var (ncid1, 'longitude','lon')
	! Compute ellipsoid corrections
	do i = 1,nrec
		dh(i) = dhellips(1,a(i))
	enddo
	call get_var (ncid1, 'altitude', a)
	call new_var ('alt_gdrf', a + dh)
	call cpy_var (ncid1, 'altitude_rate', 'alt_rate')

! Range

	call cpy_var (ncidk, 'range_ocean', 'range_ku')
	call cpy_var (ncidk, 'range_ocean_rms', 'range_rms_ku')
	call cpy_var (ncidk, 'range_ocean_numval', 'range_numval')
	call cpy_var (ncidk, 'range_ocean_compression_qual', 'qual_range')
	call cpy_var (ncidk, 'range_ocean_mle3', 'range_ku_mle3')
	call cpy_var (ncidk, 'range_ocean_mle3_rms', 'range_rms_ku_mle3')
	call cpy_var (ncidk, 'range_ocean_mle3_numval', 'range_numval_ku_mle3')
	call cpy_var (ncidc, 'range_ocean', 'range_c')
	call cpy_var (ncidc, 'range_ocean_rms', 'range_rms_c')
	call cpy_var (ncidc, 'range_ocean_numval', 'range_numval_c')
	call cpy_var (ncidk, 'range_cor_ocean_net_instr', 'drange_ku')
	call cpy_var (ncidk, 'range_cor_ocean_mle3_net_instr', 'drange_ku_mle3')
	call cpy_var (ncidc, 'range_cor_ocean_net_instr', 'drange_c')
	if (latency == rads_ntc) then
		call cpy_var (ncidk, 'range_adaptive', 'range_ku_adaptive')
		call cpy_var (ncidk, 'range_adaptive_rms', 'range_rms_ku_adaptive')
		call cpy_var (ncidk, 'range_adaptive_numval', 'range_numval_ku_adaptive')
		call cpy_var (ncidk, 'range_cor_adaptive_net_instr', 'drange_ku_adaptive')
	endif

! SWH

	call cpy_var (ncidk, 'swh_ocean', 'swh_ku')
	call cpy_var (ncidk, 'swh_ocean_rms', 'swh_rms_ku')
	call cpy_var (ncidk, 'swh_ocean_compression_qual', 'qual_swh')
	call cpy_var (ncidk, 'swh_ocean_mle3', 'swh_ku_mle3')
	call cpy_var (ncidk, 'swh_ocean_mle3_rms', 'swh_rms_ku_mle3')
	call cpy_var (ncidc, 'swh_ocean', 'swh_c')
	call cpy_var (ncidc, 'swh_ocean_rms', 'swh_rms_c')
	call cpy_var (ncidk, 'swh_cor_ocean_net_instr', 'dswh_ku')
	call cpy_var (ncidk, 'swh_cor_ocean_mle3_net_instr', 'dswh_ku_mle3')
	call cpy_var (ncidc, 'swh_cor_ocean_net_instr', 'dswh_c')
	if (latency == rads_ntc) then
		call cpy_var (ncidk, 'swh_adaptive', 'swh_ku_adaptive')
		call cpy_var (ncidk, 'swh_adaptive_rms', 'swh_rms_ku_adaptive')
		call cpy_var (ncidk, 'swh_cor_adaptive_net_instr', 'dswh_ku_adaptive')
	endif

! Backscatter

	call cpy_var (ncidk, 'sig0_ocean', 'sig0_ku')
	call cpy_var (ncidk, 'sig0_ocean_rms', 'sig0_rms_ku')
	call cpy_var (ncidk, 'sig0_ocean_compression_qual', 'qual_sig0')
	call cpy_var (ncidk, 'sig0_cor_atm', 'dsig0_atmos_ku')
	call cpy_var (ncidk, 'sig0_ocean_mle3', 'sig0_ku_mle3')
	call cpy_var (ncidk, 'sig0_ocean_mle3_rms', 'sig0_rms_ku_mle3')
	call cpy_var (ncidc, 'sig0_ocean', 'sig0_c')
	call cpy_var (ncidc, 'sig0_ocean_rms', 'sig0_rms_c')
	call cpy_var (ncidc, 'sig0_cor_atm', 'dsig0_atmos_c')
	call cpy_var (ncidk, 'sig0_cor_ocean_net_instr', 'dswh_ku')
	call cpy_var (ncidk, 'sig0_cor_ocean_mle3_net_instr', 'dsig0_ku_mle3')
	call cpy_var (ncidc, 'sig0_cor_ocean_net_instr', 'dsig0_c')
	if (latency == rads_ntc) then
		call cpy_var (ncidk, 'sig0_adaptive', 'sig0_ku_adaptive')
		call cpy_var (ncidk, 'sig0_adaptive_rms', 'sig0_rms_ku_adaptive')
		call cpy_var (ncidk, 'sig0_cor_adaptive_net_instr', 'dsig0_ku_adaptive')
	endif

! Wind speed

	call cpy_var (ncid1, 'wind_speed_alt', 'wind_speed_alt')
	call cpy_var (ncid1, 'wind_speed_alt_mle3', 'wind_speed_alt_mle3')
	call cpy_var (ncid1, 'rad_wind_speed', 'wind_speed_rad')
	call cpy_var (ncid1, 'wind_speed_mod_u', 'wind_speed_ecmwf_u')
	call cpy_var (ncid1, 'wind_speed_mod_v', 'wind_speed_ecmwf_v')
	if (latency == rads_ntc) call cpy_var (ncid1  , 'wind_speed_alt_adaptive', 'wind_speed_alt_adaptive')

! Rain or ice

	call cpy_var (ncid1, 'rain_flag', 'qual_alt_rain_ice')
	call cpy_var (ncid1, 'rad_rain_flag rad_sea_ice_flag IOR', 'qual_rad_rain_ice')

! Off-nadir angle

	call cpy_var (ncidk, 'off_nadir_angle_wf_ocean', 'off_nadir_angle2_wf_ku')
	call cpy_var (ncidk, 'off_nadir_angle_wf_ocean_rms', 'off_nadir_angle2_wf_rms_ku')
	call cpy_var (ncidk, 'off_nadir_angle_wf_ocean_compression_qual', 'qual_attitude')

! Path delay

	call cpy_var (ncid1, 'model_dry_tropo_cor_measurement_altitude', 'dry_tropo_ecmwf')
	call cpy_var (ncid1, 'model_wet_tropo_cor_measurement_altitude', 'wet_tropo_ecmwf')
	call cpy_var (ncid1, 'rad_wet_tropo_cor', 'wet_tropo_rad')
	call cpy_var (ncidk, 'iono_cor_alt', 'iono_alt')
	call cpy_var (ncidk, 'iono_cor_alt_filtered', 'iono_alt_smooth')
	call cpy_var (ncidk, 'iono_cor_alt_mle3', 'iono_alt_mle3')
	call cpy_var (ncidk, 'iono_cor_alt_filtered_mle3', 'iono_alt_smooth_mle3')
	call cpy_var (ncidc, 'range_ocean_compression_qual', 'qual_iono_alt')
	call cpy_var (ncidk, 'iono_cor_gim', 'iono_gim')
	if (latency == rads_ntc) then
		call cpy_var (ncidk, 'iono_cor_alt_adaptive', 'iono_alt_adaptive')
		call cpy_var (ncidk, 'iono_cor_alt_filtered_adaptive', 'iono_alt_smooth_adaptive')
	endif

! SSB

	call cpy_var (ncidk, 'sea_state_bias', 'ssb_cls')
	call cpy_var (ncidk, 'sea_state_bias_mle3', 'ssb_cls_mle3')
	call cpy_var (ncidk, 'sea_state_bias_3d_mp2', 'ssb_cls_3d')
	call cpy_var (ncidc, 'sea_state_bias', 'ssb_cls_c')
	if (latency == rads_ntc) then
		call cpy_var (ncidk, 'sea_state_bias_adaptive', 'ssb_cls_adaptive')
		call cpy_var (ncidk, 'sea_state_bias_adaptive_3d_mp2', 'ssb_cls_3d_adaptive')
		call cpy_var (ncidc, 'sea_state_bias_adaptive', 'ssb_cls_c_adaptive')
	endif

! IB

	call cpy_var (ncid1, 'inv_bar_cor', 'inv_bar_static')
	call cpy_var (ncid1, 'dac', 'inv_bar_mog2d')

! Tides

	call cpy_var (ncid1, 'ocean_tide_got load_tide_got SUB', 'tide_ocean_got' // tide_got_ver)
	call cpy_var (ncid1, 'ocean_tide_fes load_tide_fes SUB ocean_tide_non_eq ADD', 'tide_ocean_fes' // tide_fes_ver)
	call cpy_var (ncid1, 'load_tide_got', 'tide_load_got' // tide_got_ver)
	call cpy_var (ncid1, 'load_tide_fes', 'tide_load_fes' // tide_fes_ver)
	call cpy_var (ncid1, 'ocean_tide_eq', 'tide_equil')
	call cpy_var (ncid1, 'ocean_tide_non_eq', 'tide_non_equil')
	call cpy_var (ncid1, 'internal_tide', 'tide_internal')
	call cpy_var (ncid1, 'solid_earth_tide', 'tide_solid')
	call cpy_var (ncid1, 'pole_tide', 'tide_pole')

! Geoid and MSS

	call get_var (ncid1, 'geoid', a)
	call new_var ('geoid_egm2008', a + dh)
	call get_var (ncid1, 'mean_sea_surface_cnescls', a)
	call new_var ('mss_cnescls' // mss_cnescls_ver, a + dh)
	call get_var (ncid1, 'mean_sea_surface_dtu', a)
	call new_var ('mss_dtu' // mss_dtu_ver, a + dh)

! Surface type and coastal proximity

	call cpy_var (ncid1, 'depth_or_elevation', 'topo_ace2')
	call cpy_var (ncid1, 'surface_classification_flag', 'surface_class')
	call cpy_var (ncid1, 'rad_surface_type_flag', 'surface_type_rad')
	call cpy_var (ncid1, 'distance_to_coast 1e-3 MUL', 'dist_coast') ! Convert m to km
	call cpy_var (ncid1, 'angle_of_approach_to_coast', 'angle_coast')
	call cpy_var (ncid1, 'rad_distance_to_land 1e-3 MUL', 'rad_dist_coast') ! Convert m to km

! Bit flags

	call new_var ('flags', dble(flags))
	call new_var ('flags_mle3', dble(flags_mle3))
	if (latency == rads_ntc) call new_var ('flags_adaptive', dble(flags_adaptive))

! Other radiometer measurements
! Selected smoothed TBs

	call cpy_var (ncid1, 'rad_tb_187', 'tb_187')
	call cpy_var (ncid1, 'rad_tb_238', 'tb_238')
	call cpy_var (ncid1, 'rad_tb_340', 'tb_340')
	call cpy_var (ncid1, 'rad_tb_340_qual 2 MUL rad_tb_238_qual ADD 2 MUL rad_tb_187_qual ADD', 'qual_rad_tb')
	call cpy_var (ncid1, 'rad_cloud_liquid_water', 'liquid_water_rad')
	call cpy_var (ncid1, 'rad_water_vapor', 'water_vapor_rad')

! Wave model

	call cpy_var (ncid1, 'mean_wave_period_t02', 'mean_wave_period')
	call cpy_var (ncid1, 'mean_wave_direction', 'mean_wave_direction')

! SSHA

	call cpy_var (ncidk, 'ssha', 'ssha')
	call cpy_var (ncidk, 'ssha_mle3', 'ssha_mle3')
	if (latency == rads_ntc) call cpy_var (ncidk, 'ssha_adaptive', 'ssha_adaptive')

! Misc

	a = latency
	call new_var ('latency', a)

! Dump the data

	call nfs(nf90_close(ncid))
	call put_rads
	deallocate (a, dh, flags, flags_mle3, flags_adaptive, flags_save)

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Write Jason-1/2/3 GDR-F data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_jason_file_names')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --min-rec=MIN_REC         Specify minimum number of records per pass to process' // &
'This program converts Jason-1/2/3 OGDR/IGDR/GDR GDR-F Level 2 products to RADS data' / &
'files with the name $RADSDATAROOT/data/jJ/F/pPPPP/jJpPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_jason_gdrf
