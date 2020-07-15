!-----------------------------------------------------------------------
! Copyright (c) 2011-2019  Remko Scharroo and Eric Leuliette
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

!*rads_gen_s3 -- Converts Sentinel-3 data to RADS
!+
program rads_gen_s3

! This program reads Sentinel-3 NRT/STC/NTC files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/SS/F/SSpPPPPcCCC.nc.
!    SS = satellite (3a or 3b)
!     F = phase (a)
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_s3 [options] < list_of_Sentinel3_file_names
!
! where [options] include:
!  --min-rec <min_rec> : Specify minimum number of records per pass to process.
!
! This program handles Sentinel-3 standard_measurement files in NetCDF format.
! The format is described in:
!
! [1] SENTINEL-3 Products and Algorithms Definition (S3PAD)
!     Surface Topography Mission (STM) SRAL/MWR L2 Products Specifications
!     CLS-DOS-NT-09-033, 9 rev 0, 12 July 2011
!     S3PAD-RS-CLS-SD02-00013
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
! iono_alt_* - Dual-frequency ionospheric correction (not _c)
! iono_gim - GIM ionosphetic correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! ssb_cls_* - SSB
! swh_* - Significant wave height
! swh_rms_* - Std dev of SWH
! sig0_* - Sigma0 (not corrected for attenuation ... will be done in rads_fix_s3)
! sig0_rms_* - Std dev of sigma0
! dsig0_atmos_* - Atmospheric attenuation of sigma0 (not _plrm)
! wind_speed_alt_* - Altimeter wind speed (not _c)
! wind_speed_rad - Radiometer wind speed
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! tide_ocean/load_got410 - GOT4.10c ocean and load tide
! tide_ocean/load_fes04/fes14 - FES2004 or FES2014 ocean and load tide
! tide_pole - Pole tide
! tide_solid - Solid earth tide
! topo_ace2 - ACE2 topography
! geoid_egm2008 - EGM2008 geoid
! mss_cnescls11/cnescls15 - CNES/CLS11 or CNES/CLS15 mean sea surface
! mss_dtu13 - DTU13 mean sea surface
! tb_238 - Brightness temperature (23.8 GHz)
! tb_365 - Brightness temperature (36.5 GHz)
! flags, flags_plrm - Engineering flags
! off_nadir_angle2_wf_* - Mispointing from waveform squared (not _c)
! rad_liquid_water - Liquid water content
! rad_water_vapor - Water vapor content
! ssha_* - Sea surface height anomaly (not _c)
!
! Extensions _* are:
! _ku:      Ku-band retracked from SAR
! _ku_plrm: Ku-band retracked with PLRM
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

integer(fourbyteint) :: ncid, cyclenr, passnr, varid, orbit_type_offset
real(eightbytereal) :: equator_time

! Data variables

integer(twobyteint), allocatable :: flags_plrm(:), flags_save(:)
character(len=2) :: mss_cnescls_ver, mss_dtu_ver, tide_fes_ver
integer :: latency = rads_nrt
logical :: iono_alt_smooth = .false.

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

	call log_string (infile)
	if (nf90_open(infile,nf90_nowrite,ncid) /= nf90_noerr) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Check if input is a Sentinel-3 Level 2 data set

	if (nf90_get_att(ncid,nf90_global,'title',arg) /= nf90_noerr .or. &
		arg /= 'IPF SRAL/MWR Level 2 Measurement') then
		call log_string ('Error: this is not a Sentinel-3 SRAL Level 2 data set', .true.)
		cycle
	endif

! Get the mission name and initialise RADS (if not done before)

	call nfs(nf90_get_att(ncid,nf90_global,'mission_name',arg))
	select case (arg)
	case ('Sentinel 3A')
		arg = '3a'
	case ('Sentinel 3B')
		arg = '3b'
	case default
		call log_string ('Error: wrong misson_name found in header', .true.)
		cycle
	end select
	if (sat == '') then
		call rads_init (S, arg)
	else if (S%sat /= arg(:2)) then
		call log_string ('Error: wrong misson_name found in header', .true.)
		cycle
	endif

! Read global attributes

	call nfs(nf90_inq_dimid(ncid,'time_01',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
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

! Determine latency (NRT, STC, NTC)

	call nfs(nf90_get_att(ncid,nf90_global,'product_name',arg))
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

! Detect REP_006 (S3A) or REP_007 (S3B) or REP_008

	if (index(arg, '_R_NT_003') > 0) then
		if (S%sat == '3a') then
			latency = rads_ntc + 6 ! REP_006
		else
			latency = rads_ntc + 7 ! REP_007
		endif
	else if (index(arg, '_R_NT_004') > 0) then
		latency = rads_ntc + 8 ! REP_008
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
	call nfs(nf90_get_att(ncid,nf90_global,'first_meas_time',arg))
	P%start_time = strp1985f(arg)
	call nfs(nf90_get_att(ncid,nf90_global,'last_meas_time',arg))
	P%end_time = strp1985f(arg)

! Check for orbit_type_01, otherwise take old orbit_data_type

	orbit_type_offset = -1
	if (nft(nf90_inq_varid(ncid,'orbit_type_01',varid))) then
		i = nf90_inq_varid(ncid,'orbit_data_type',varid)
		i = nf90_get_att(ncid,varid,'flag_meanings',arg)
		if (arg(1:4) == 'scen') then
			orbit_type_offset = 0
		else
			orbit_type_offset = 1
		endif
	endif

! Determine L2 processing version

	call nfs(nf90_get_att(ncid,nf90_global,'source',arg))

! Update the versions for MSS CNES-CLS and FES tides for PB 2.19 (IPF-SM-2 06.08) and later

	if (arg(10:14) < '06.08') then
		mss_cnescls_ver = '11'
		mss_dtu_ver = '13'
		tide_fes_ver = '04'
	else
		mss_cnescls_ver = '15'
		mss_dtu_ver = '15'
		tide_fes_ver = '14'
	endif

! Update the version of MSS DTU and switch on smoothed iono for PB 2.49 (IPF-SM-2 06.16) and later

	if (arg(10:17) >= '06.16   ') then
		mss_dtu_ver = '18'
		iono_alt_smooth = .true.
	endif

! Store input file name

	i = index(infile, '/', .true.) + 1
	P%original = trim(infile(i:)) // ' (' // trim(arg) // ')'

! Allocate variables

	allocate (a(nrec),dh(nrec),flags(nrec),flags_plrm(nrec),flags_save(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
	call nc2f (ncid, 'instr_op_mode_01',0,ge=1)			! bit  0: Altimeter mode
	call nc2f (ncid, 'val_alt_off_nadir_angle_wf_ocean_01_plrm_ku',1)	! bit  1: Quality off-nadir pointing
	call nc2f (ncid, 'surf_class_01',2,eq=4)				! bit  2: Continental ice
	call nc2f (ncid, 'range_ocean_qual_01_c',3)			! bit  3: Quality dual-frequency iono
	call nc2f (ncid, 'surf_class_01',4,eq=1)
	call nc2f (ncid, 'surf_class_01',4,ge=3)				! bit  4: Water/land
	call nc2f (ncid, 'surf_class_01',5,ge=1)				! bit  5: Ocean/other
	call nc2f (ncid, 'rad_surf_type_01',6,ge=2)			! bit  6: Radiometer land flag
	call nc2f (ncid, 'rain_flag_01_ku',7)
	call nc2f (ncid, 'open_sea_ice_flag_01_ku',7)			! bit  7: Altimeter rain or ice flag
	call nc2f (ncid, 'rain_flag_01_ku',8)					! bit  8: Altimeter rain flag
	call nc2f (ncid, 'tb_238_quality_flag_01',9)			! bit  9: Quality 23.8 GHz channel
	call nc2f (ncid, 'tb_365_quality_flag_01',10)			! bit 10: Quality 36.5 GHz channel
	if (orbit_type_offset < 0) then					! bit 15: Quality of orbit
		call nc2f (ncid, 'orbit_type_01',15,le=2-orbit_type_offset)
	else
		call nc2f (ncid, 'orbit_data_type',15,le=2-orbit_type_offset)
	endif

! Now do specifics for PLRM

	flags_save = flags	! Keep flags for later
	call nc2f (ncid, 'range_ocean_qual_01_plrm_ku',11)	! bit 11: Quality range
	call nc2f (ncid, 'swh_ocean_qual_01_plrm_ku',12)		! bit 12: Quality SWH
	call nc2f (ncid, 'sig0_ocean_qual_01_plrm_ku',13)		! bit 13: Quality Sigma0
	flags_plrm = flags	! Copy result for PLRM
	flags = flags_save	! Continue with SAR flags

! Redo the last ones for SAR

	call nc2f (ncid, 'range_ocean_qual_01_ku',11)			! bit 11: Quality range
	call nc2f (ncid, 'swh_ocean_qual_01_ku',12)			! bit 12: Quality SWH
	call nc2f (ncid, 'sig0_ocean_qual_01_ku',13)			! bit 13: Quality Sigma0

! Convert all the necessary fields to RADS
	call get_var (ncid, 'time_01', a)
	call new_var ('time', a + sec2000)
	call cpy_var (ncid, 'lat_01', 'lat')
	! Compute ellipsoid corrections
	do i = 1,nrec
		dh(i) = dhellips(1,a(i))
	enddo
	call cpy_var (ncid, 'lon_01','lon')
	call get_var (ncid, 'alt_01', a)
	call new_var ('alt_gdrf', a + dh)
	! Note that GDR-E orbit standards were used for the earlier part of the mission until reprocessing.
	! However this field is updated by rads_add_orbit hereafter.
	call cpy_var (ncid, 'orb_alt_rate_01', 'alt_rate')
	call cpy_var (ncid, 'range_ocean_01_ku','range_ku')
	call cpy_var (ncid, 'range_ocean_01_plrm_ku','range_ku_plrm')
	call cpy_var (ncid, 'range_ocean_01_c','range_c')
! Add zero or meas altitude tropo measurements?
	call cpy_var (ncid, 'mod_dry_tropo_cor_meas_altitude_01', 'dry_tropo_ecmwf')
	call cpy_var (ncid, 'rad_wet_tropo_cor_01_ku', 'wet_tropo_rad')
	call cpy_var (ncid, 'rad_wet_tropo_cor_01_plrm_ku', 'wet_tropo_rad_plrm')
	call cpy_var (ncid, 'rad_wet_tropo_cor_sst_gam_01_ku', 'wet_tropo_rad_sst_gam')
	call cpy_var (ncid, 'mod_wet_tropo_cor_meas_altitude_01', 'wet_tropo_ecmwf')
	call cpy_var (ncid, 'comp_wet_tropo_cor_01_ku', 'wet_tropo_comp')
	call cpy_var (ncid, 'iono_cor_alt_01_ku', 'iono_alt')
	call cpy_var (ncid, 'iono_cor_alt_01_plrm_ku', 'iono_alt_plrm')
	if (iono_alt_smooth) then
		call cpy_var (ncid, 'iono_cor_alt_filtered_01_ku', 'iono_alt_smooth')
		call cpy_var (ncid, 'iono_cor_alt_filtered_01_plrm_ku', 'iono_alt_smooth_plrm')
	endif
	call cpy_var (ncid, 'iono_cor_gim_01_ku', 'iono_gim')
	call cpy_var (ncid, 'inv_bar_cor_01', 'inv_bar_static')
	if (latency == rads_nrt) then
		call cpy_var (ncid, 'inv_bar_cor_01', 'inv_bar_mog2d')
	else
		call cpy_var (ncid, 'inv_bar_cor_01 hf_fluct_cor_01 ADD', 'inv_bar_mog2d')
	endif
	call cpy_var (ncid, 'solid_earth_tide_01', 'tide_solid')
	call cpy_var (ncid, 'ocean_tide_sol1_01 load_tide_sol1_01 SUB', 'tide_ocean_got410')
	call cpy_var (ncid, 'ocean_tide_sol2_01 load_tide_sol2_01 SUB ocean_tide_non_eq_01 ADD', 'tide_ocean_fes' // tide_fes_ver)
	call cpy_var (ncid, 'load_tide_sol1_01', 'tide_load_got410')
	call cpy_var (ncid, 'load_tide_sol2_01', 'tide_load_fes' // tide_fes_ver)
	call cpy_var (ncid, 'ocean_tide_eq_01', 'tide_equil')
	call cpy_var (ncid, 'ocean_tide_non_eq_01', 'tide_non_equil')
	call cpy_var (ncid, 'pole_tide_01', 'tide_pole')
	call cpy_var (ncid, 'sea_state_bias_01_ku', 'ssb_cls')
	call cpy_var (ncid, 'sea_state_bias_01_plrm_ku', 'ssb_cls_plrm')
	call cpy_var (ncid, 'sea_state_bias_01_c', 'ssb_cls_c')
	call get_var (ncid, 'geoid_01', a)
	call new_var ('geoid_egm2008', a + dh)
	call get_var (ncid, 'mean_sea_surf_sol1_01', a)
	call new_var ('mss_cnescls' // mss_cnescls_ver, a + dh)
	call get_var (ncid, 'mean_sea_surf_sol2_01', a)
	call new_var ('mss_dtu' // mss_dtu_ver, a + dh)
	call cpy_var (ncid, 'swh_ocean_01_ku', 'swh_ku')
	call cpy_var (ncid, 'swh_ocean_01_plrm_ku', 'swh_ku_plrm')
	call cpy_var (ncid, 'swh_ocean_01_c', 'swh_c')
	call cpy_var (ncid, 'sig0_ocean_01_ku', 'sig0_ku')
	call cpy_var (ncid, 'sig0_ocean_01_plrm_ku', 'sig0_ku_plrm')
	call cpy_var (ncid, 'sig0_ocean_01_c', 'sig0_c')
	call cpy_var (ncid, 'wind_speed_alt_01_ku', 'wind_speed_alt')
	call cpy_var (ncid, 'wind_speed_alt_01_plrm_ku', 'wind_speed_alt_plrm')
	call cpy_var (ncid, 'wind_speed_mod_u_01', 'wind_speed_ecmwf_u')
	call cpy_var (ncid, 'wind_speed_mod_v_01', 'wind_speed_ecmwf_v')
	call cpy_var (ncid, 'range_ocean_rms_01_ku', 'range_rms_ku')
	call cpy_var (ncid, 'range_ocean_rms_01_plrm_ku', 'range_rms_ku_plrm')
	call cpy_var (ncid, 'range_ocean_rms_01_c', 'range_rms_c')
	call cpy_var (ncid, 'range_ocean_numval_01_ku', 'range_numval_ku')
	call cpy_var (ncid, 'range_ocean_numval_01_plrm_ku', 'range_numval_ku_plrm')
	call cpy_var (ncid, 'range_ocean_numval_01_c', 'range_numval_c')
	call cpy_var (ncid, 'odle_01', 'topo_ace2')
	call cpy_var (ncid, 'tb_238_01','tb_238')
	call cpy_var (ncid, 'tb_365_01','tb_365')
	call new_var ('flags', dble(flags))
	call new_var ('flags_plrm', dble(flags_plrm))
	call cpy_var (ncid, 'swh_ocean_rms_01_ku', 'swh_rms_ku')
	call cpy_var (ncid, 'swh_ocean_rms_01_plrm_ku', 'swh_rms_ku_plrm')
	call cpy_var (ncid, 'swh_ocean_rms_01_c', 'swh_rms_c')
	call cpy_var (ncid, 'sig0_ocean_rms_01_ku', 'sig0_rms_ku')
	call cpy_var (ncid, 'sig0_ocean_rms_01_plrm_ku', 'sig0_rms_ku_plrm')
	call cpy_var (ncid, 'sig0_ocean_rms_01_c', 'sig0_rms_c')
	! In case of SAR, there is no estimated attitude, so we rely on the PLRM values.
	! It may not be good to use these for editing (?)
	call cpy_var (ncid, 'corrected_off_nadir_angle_wf_ocean_01_ku corrected_off_nadir_angle_wf_ocean_01_plrm_ku AND', &
		'off_nadir_angle2_wf_ku')
	call cpy_var (ncid, 'off_nadir_angle_rms_01_ku off_nadir_angle_rms_01_plrm_ku AND', 'off_nadir_angle2_wf_rms_ku')
	call cpy_var (ncid, 'off_nadir_pitch_angle_pf_01', 'attitude_pitch')
	call cpy_var (ncid, 'off_nadir_roll_angle_pf_01', 'attitude_roll')
	call cpy_var (ncid, 'off_nadir_yaw_angle_pf_01', 'attitude_yaw')
	call cpy_var (ncid, 'atm_cor_sig0_01_ku', 'dsig0_atmos_ku')
	call cpy_var (ncid, 'atm_cor_sig0_01_c', 'dsig0_atmos_c')
	call cpy_var (ncid, 'rad_liquid_water_01_ku', 'liquid_water_rad')
	call cpy_var (ncid, 'rad_water_vapor_01_ku', 'water_vapor_rad')
	call cpy_var (ncid, 'ssha_01_ku', 'ssha')
	call cpy_var (ncid, 'ssha_01_plrm_ku', 'ssha_plrm')
	if (orbit_type_offset < 0) then	! Straight copy of orbit_type_01
		call cpy_var (ncid, 'orbit_type_01', 'orbit_type')
	else ! Convert orbit_data_type (old school, before PB 2.19)
		call get_var (ncid, 'orbit_data_type', a)
		a = a + orbit_type_offset
		where (a > 4) a = a * 2 - 4
		call new_var ('orbit_type', a)
	endif
	call cpy_var (ncid, 'surf_class_01', 'surface_class')
	a = latency
	call new_var ('latency', a)

! Dump the data

	call nfs(nf90_close(ncid))
	call put_rads
	deallocate (a, dh, flags, flags_plrm, flags_save)

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Write Sentinel-3 data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_Sentinel3_file_names')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --min-rec=MIN_REC         Specify minimum number of records per pass to process' // &
'This program converts Sentinel-3 NRT/STC/NTC files to RADS data' / &
'files with the name $RADSDATAROOT/data/SS/F/pPPPP/s3pPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_s3
