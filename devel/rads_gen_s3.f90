!-----------------------------------------------------------------------
! Copyright (c) 2011-2016  Remko Scharroo and Eric Leuliette
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
! written into files $RADSDATAROOT/data/s3/a/s3pPPPPcCCC.nc.
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_s3 [options] < list_of_Sentinel3_file_names
!
! This program handles Sentinel-3 standard_measurement files in netCDF format.
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
! alt_gdrd - Orbit altitude
! alt_rate - Orbit altitude rate
! range_* - Ocean range (retracked)
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! wet_tropo_rad - Radiometer wet tropo correction
! iono_alt_* - Dual-frequency ionospheric correction (not _c)
! iono_gim - GIM ionosphetic correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! ssb_cls_* - SSB
! swh_* - Significant wave height
! sig0_* - Sigma0
! wind_speed_alt_* - Altimeter wind speed (not _c)
! wind_speed_rad - Radiometer wind speed
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! range_rms_* - Std dev of range
! range_numval_* - Nr of averaged range measurements
! topo_dtm2000 - Bathymetry
! tb_238 - Brightness temperature (23.8 GHz)
! tb_340 - Brightness temperature (34.0 GHz)
! flags, flags_plrm - Engineering flags
! swh_rms_* - Std dev of SWH
! sig0_rms_* - Std dev of sigma0
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! rad_liquid_water - Liquid water content
! rad_water_vapor - Water vapor content
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
use netcdf

! Command line arguments

integer(fourbyteint) :: ios, i
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: cyclenr, passnr, varid
real(eightbytereal) :: equator_time
logical :: nrt

! Data variables

integer(twobyteint), allocatable :: flags_plrm(:), flags_save(:)

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000

! Initialise

call synopsis
call rads_gen_getopt ('3a')
call synopsis ('--head')
call rads_init (S, sat)

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

! Check if input is standard_measurement.nc

	if (index(infile,'standard_measurement') <= 0) then
		call log_string ('Error: this is not standard_measurement.nc', .true.)
		cycle
	endif

! Read global attributes

	call nfs(nf90_inq_dimid(ncid,'time_01',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
	if (nrec == 0) then
		cycle
	else if (nrec > mrec) then
		call log_string ('Error: too many measurements', .true.)
		cycle
	endif
	call nfs(nf90_get_att(ncid,nf90_global,'mission_name',arg))
	if (arg /= 'Sentinel 3A') then
		call log_string ('Error: wrong misson-name found in header', .true.)
		cycle
	endif

	call nfs(nf90_get_att(ncid,nf90_global,'title',arg))
	nrt = (arg(:4) == 'nrt')
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

! Determine L2 processing version

	call nfs(nf90_get_att(ncid,nf90_global,'references',arg))
	i = index(infile, '/', .true.) + 1
	P%original = trim(infile(i:)) // ' (' // arg // ')'

! Allocate variables

	allocate (a(nrec),flags(nrec),flags_plrm(nrec),flags_save(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
	call nc2f ('instr_op_mode_01',0,lim=2)			! bit  0: Altimeter mode
	call nc2f ('corrected_off_nadir_angle_wf_ocean_01_ku',1)	! bit  1: Quality off-nadir pointing
	call nc2f ('surf_type_01',2,val=2)				! bit  2: Continental ice
	call nc2f ('range_ocean_qual_01_c',3)			! bit  3: Quality dual-frequency iono
	call nc2f ('surf_type_01',4,lim=2)				! bit  4: Water/land
	call nc2f ('surf_type_01',5,lim=1)				! bit  5: Ocean/other
	call nc2f ('rad_surf_type_01',6,lim=2)			! bit  6: Radiometer land flag
	call nc2f ('rain_flag_01_ku',7)
	call nc2f ('open_sea_ice_flag_01_ku',7)			! bit  7: Altimeter rain or ice flag
	call nc2f ('rain_flag_01_ku',8)					! bit  8: Altimeter rain flag
	call nc2f ('tb_238_quality_flag_01',9)			! bit  9: Quality 23.8 GHz channel
	call nc2f ('tb_365_quality_flag_01',10)			! bit 10: Quality 36.5 GHz channel
	call nc2f ('range_ocean_qual_01_ku',11)			! bit 11: Quality range
	call nc2f ('swh_ocean_qual_01_ku',12)			! bit 12: Quality SWH
	call nc2f ('sig0_ocean_qual_01_ku',13)			! bit 13: Quality Sigma0

! Need to add PLRM flags

! Convert all the necessary fields to RADS
	call get_var (ncid, 'time_01', a)
	call new_var ('time', a + sec2000)
	call cpy_var ('lat_01','lat')
	call cpy_var ('lon_01','lon')
!Which alt field?
	call cpy_var ('alt_01', 'alt')
	call cpy_var ('orb_alt_rate_01', 'alt_rate')
	call cpy_var ('range_ocean_01_ku','range_ku')
	call cpy_var ('range_ocean_01_c','range_c')
!Add PLRM ranges?
!Add zero or meas altitude tropo measurements?
	call cpy_var ('mod_dry_tropo_cor_meas_altitude_01', 'dry_tropo_ecmwf')
	call cpy_var ('rad_wet_tropo_cor_01_ku', 'wet_tropo_rad')
	call cpy_var ('mod_wet_tropo_cor_meas_altitude_01', 'wet_tropo_ecmwf')
	call cpy_var ('iono_cor_alt_01_ku', 'iono_alt')
	if (.not.nrt) call cpy_var ('iono_cor_gim_01_ku', 'iono_gim')
	call cpy_var ('inv_bar_cor_01', 'inv_bar_static')
	if (nrt) then
		call cpy_var ('inv_bar_cor_01', 'inv_bar_mog2d')
	else
		call cpy_var ('inv_bar_cor_01 hf_fluct_cor_01 ADD', 'inv_bar_mog2d')
	endif
	call cpy_var ('solid_earth_tide_01', 'tide_solid')
	call cpy_var ('ocean_tide_sol1_01 load_tide_sol1_01 SUB', 'tide_ocean_got48')
	call cpy_var ('ocean_tide_sol2_01 load_tide_sol2_01 SUB', 'tide_ocean_fes04')
	call cpy_var ('load_tide_sol1_01', 'tide_load_got48')
	call cpy_var ('load_tide_sol2_01', 'tide_load_fes04')
	call cpy_var ('pole_tide_01', 'tide_pole')
	call cpy_var ('sea_state_bias_01_ku', 'ssb_cls')
	call cpy_var ('sea_state_bias_01_c', 'ssb_cls_c')
	call cpy_var ('geoid_01', 'geoid_egm2008')
	call cpy_var ('mean_sea_surf_sol1_01', 'mss_cnescls11')
	call cpy_var ('mean_sea_surf_sol2_01', 'mss_dtu10')
	call cpy_var ('swh_ocean_01_ku','swh_ku')
	call cpy_var ('swh_ocean_01_c','swh_c')
	call cpy_var ('sig0_ocean_01_ku','sig0_ku')
	call cpy_var ('sig0_ocean_01_c','sig0_c')
	call cpy_var ('wind_speed_alt_01_ku','wind_speed_alt')
	call cpy_var ('wind_speed_mod_u_01', 'wind_speed_ecmwf_u')
	call cpy_var ('wind_speed_mod_v_01', 'wind_speed_ecmwf_v')
	call cpy_var ('range_ocean_rms_01_ku','range_rms_ku')
	call cpy_var ('range_ocean_rms_01_c','range_rms_c')
	call cpy_var ('range_ocean_numval_01_ku','range_numval_ku')
	call cpy_var ('range_ocean_numval_01_c','range_numval_c')
	call cpy_var ('odle_01', 'topo_ace2')
	call cpy_var ('tb_238_01','tb_238')
	call cpy_var ('tb_365_01','tb_365')
	a = flags
	call new_var ('flags', a)
	a = flags_plrm
	call new_var ('flags_plrm', a)
	call cpy_var ('swh_ocean_rms_01_ku','swh_rms_ku')
	call cpy_var ('swh_ocean_rms_01_plrm_ku','swh_rms_ku_plrm')
	call cpy_var ('swh_ocean_rms_01_c','swh_rms_c')
	call cpy_var ('sig0_ocean_rms_01_ku','sig0_rms_ku')
	call cpy_var ('sig0_ocean_rms_01_plrm_ku','sig0_rms_ku_plrm')
	call cpy_var ('sig0_ocean_rms_01_c','sig0_rms_c')
	call cpy_var ('corrected_off_nadir_angle_wf_ocean_01_ku', 'off_nadir_angle2_wf_ku')
	call cpy_var ('atm_cor_sig0_01_ku', 'dsig0_atmos_ku')
	call cpy_var ('atm_cor_sig0_01_c', 'dsig0_atmos_c')
	call cpy_var ('rad_liquid_water_01_ku', 'liquid_water_rad')
	call cpy_var ('rad_wet_tropo_cor_01_ku', 'water_vapor_rad')

! Dump the data

	call nfs(nf90_close(ncid))
	call put_rads
	deallocate (a, flags, flags_plrm, flags_save)

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
'This program converts Sentinel-3 NRT/STC/NTC files to RADS data' / &
'files with the name $RADSDATAROOT/data/s3/a/pPPPP/s3pPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_s3
