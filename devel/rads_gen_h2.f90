!-----------------------------------------------------------------------
! Copyright (c) 2011-2026  Remko Scharroo and Eric Leuliette
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

!*rads_gen_h2 -- Converts HY-2A data to RADS
!+
program rads_gen_h2

! This program reads HY-2A IGDR files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/2a/a/2apPPPPcCCC.nc.
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_h2 [options] < list_of_Sentinel3_file_names
!
! This program handles HY-2A IGDR files in NetCDF format.
! The format is described in:
!
!
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdre - Orbital altitude
! orb_alt_rate - Orbit altitude rate
! range_* - Ocean range (retracked)
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! wet_tropo_rad - Radiometer wet tropo correction
! iono_alt_ku - Dual-frequency ionospheric correction
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
! topo_etopo1 - Bathymetry
! tb_187 - Brightness temperature (18.7 GHz)
! tb_238 - Brightness temperature (23.8 GHz)
! tb_370 - Brightness temperature (37.0 GHz)
! flags - Engineering flags
! swh_rms_* - Std dev of SWH
! sig0_rms_* - Std dev of sigma0
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! off_nadir_angle2_pf - Mispointing from platform squared
! rad_liquid_water - Liquid water content
! rad_water_vapor - Water vapor content
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
use netcdf

! Command line arguments

integer(fourbyteint) :: ios
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: ncid, cyclenr, passnr, varid
real(eightbytereal) :: equator_time, x0, x1

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000

! Initialise

call synopsis
call rads_gen_getopt ('2a')
call synopsis ('--head')
call rads_init (S, sat)

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

do
	read (*,'(a)',iostat=ios) infile
	if (ios /= 0) exit
	call log_string (basename(infile))

! Check if input is a HY-2A IGDR data set

	if (index(infile, 'H2A_RA1_IDR_2PT_') == 0) then
		call log_string ('Error: this is not a HY-2A IGDR file', .true.)
		cycle
	endif

! Open input file

	if (nf90_open(infile,nf90_nowrite,ncid) /= nf90_noerr) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Read global attributes

	call nfs(nf90_inq_dimid(ncid,'time',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
	if (nrec == 0) then
		cycle
	else if (nrec > mrec) then
		call log_string ('Error: too many measurements', .true.)
		cycle
	endif
	call nfs(nf90_get_att(ncid,nf90_global,'Mission_Name',arg))
	if (index(arg, 'HY-2') == 0) then
		call log_string ('Error: wrong misson-name found in header', .true.)
		cycle
	endif

	call nfs(nf90_get_att(ncid,nf90_global,'Cycle_Number',cyclenr))
	call nfs(nf90_get_att(ncid,nf90_global,'Pass_Number',passnr))

! Equator_Time is wrong, so we need to approximate

	call nfs(nf90_get_att(ncid,nf90_global,'First_Measurement_Time',arg))
	x0 = strp1985f (arg)
	call nfs(nf90_get_att(ncid,nf90_global,'Last_Measurement_Time',arg))
	x1 = strp1985f (arg)
	equator_time = 0.5d0 * (x0 + x1)

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
	P%start_time = x0
	P%end_time = x1

! Equator_Longitude is wrong, so we need to approximate

	call nfs(nf90_get_att(ncid,nf90_global,'First_Measurement_Longitude',x0))
	call nfs(nf90_get_att(ncid,nf90_global,'Last_Measurement_Longitude',x1))
	if (x0 < x1) x0 = x0 + 360d0
	P%equator_lon = 0.5d0 * (x0 + x1)

! Store input file name

	P%original = trim(basename(infile)) // ' (' // arg // ')'

! Allocate variables

	allocate (a(nrec), flags(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
	call nc2f (ncid, 'alt_state_flag',0)			! bit  0: Altimeter mode
!	call nc2f (ncid, 'val_alt_off_nadir_angle_wf_ocean_01_plrm_ku',1)	! bit  1: Quality off-nadir pointing
!	call nc2f (ncid, 'surf_type_01',2,eq=2)				! bit  2: Continental ice
!	call nc2f (ncid, 'range_ocean_qual_01_c',3)			! bit  3: Quality dual-frequency iono
	call nc2f (ncid, 'surface_type',4,ge=2)				! bit  4: Water/land
	call nc2f (ncid, 'surface_type',5,ge=1)				! bit  5: Ocean/other
	call nc2f (ncid, 'rad_surf_type',6,ge=2)				! bit  6: Radiometer land flag
	call nc2f (ncid, 'rain_flag',7)
	call nc2f (ncid, 'ice_flag',7)						! bit  7: Altimeter rain or ice flag
	call nc2f (ncid, 'rain_flag',8)						! bit  8: Altimeter rain flag
!	call nc2f (ncid, 'tb_238_quality_flag_01',9)			! bit  9: Quality 23.8 GHz channel
!	call nc2f (ncid, 'tb_365_quality_flag_01',10)			! bit 10: Quality 36.5 GHz channel
	call nc2f (ncid, 'alt_echo_type',11)					! bit 11: Quality range
!	call nc2f (ncid, 'swh_ocean_qual_01_plrm_ku',12)		! bit 12: Quality SWH
!	call nc2f (ncid, 'sig0_ocean_qual_01_plrm_ku',13)		! bit 13: Quality Sigma0
!	call nc2f (ncid, 'orb_state_flag',15,ge=1)			! bit 15: Quality orbit
	a = flags
	call new_var ('flags', a)

! Convert all the necessary fields to RADS
	call get_var (ncid, 'time', a)
	call new_var ('time', a + sec2000)
	call cpy_var (ncid, 'latitude', 'lat')
	call cpy_var (ncid, 'longitude','lon')
	call cpy_var (ncid, 'altitude', 'alt_gdrd')
	call get_var (ncid, 'orb_alt_rate', a)
	call new_var ('alt_rate', a * 1d-2)	! Convert cm/s to m/s
	call cpy_var (ncid, 'range_ku')
	call cpy_var (ncid, 'range_c')
	call cpy_var (ncid, 'range_rms_ku')
	call cpy_var (ncid, 'range_rms_c')
	call cpy_var (ncid, 'range_numval_ku')
	call cpy_var (ncid, 'range_numval_c')
	call cpy_var (ncid, 'net_instr_corr_ku', 'drange_ku')
	call cpy_var (ncid, 'net_instr_corr_c', 'drange_c')

	call cpy_var (ncid, 'model_dry_tropo_corr', 'dry_tropo_ecmwf')
	call cpy_var (ncid, 'model_wet_tropo_corr', 'wet_tropo_ecmwf')
	call cpy_var (ncid, 'rad_wet_tropo_corr', 'wet_tropo_rad')
	call cpy_var (ncid, 'iono_corr_alt_ku', 'iono_alt')
	call cpy_var (ncid, 'sea_state_bias_ku', 'ssb_cls')
	call cpy_var (ncid, 'sea_state_bias_c', 'ssb_cls_c')

	call cpy_var (ncid, 'swh_ku')
	call cpy_var (ncid, 'swh_c')
	call cpy_var (ncid, 'swh_rms_ku')
	call cpy_var (ncid, 'swh_rms_c')
!	call cpy_var (ncid, 'swh_numval_ku')
!	call cpy_var (ncid, 'swh_numval_c')
	call cpy_var (ncid, 'net_instr_corr_swh_ku', 'dswh_ku')
	call cpy_var (ncid, 'net_instr_corr_swh_c', 'dswh_c')

	call cpy_var (ncid, 'sig0_ku')
	call cpy_var (ncid, 'sig0_c')
	call cpy_var (ncid, 'sig0_rms_ku')
	call cpy_var (ncid, 'sig0_rms_c')
!	call cpy_var (ncid, 'sig0_numval_ku')
!	call cpy_var (ncid, 'sig0_numval_c')
	call cpy_var (ncid, 'net_instr_sig0_corr_ku', 'dsig0_ku')
	call cpy_var (ncid, 'net_instr_sig0_corr_c', 'dsig0_c')
	call cpy_var (ncid, 'atmos_sig0_corr_ku', 'dsig0_atmos_ku')
	call cpy_var (ncid, 'atmos_sig0_corr_c', 'dsig0_atmos_c')

	call cpy_var (ncid, 'off_nadir_angle_ku_wvf', 'off_nadir_angle2_wf_ku')
	call cpy_var (ncid, 'off_nadir_angle_ptf', 'off_nadir_angle2_pf')

	call cpy_var (ncid, 'tb_187')
	call cpy_var (ncid, 'tb_238')
	call cpy_var (ncid, 'tb_370')

	call cpy_var (ncid, 'mean_sea_surface', 'mss_cnescls11')
	call cpy_var (ncid, 'geoid', 'geoid_egm2008')
	call cpy_var (ncid, 'bathymetry', 'topo_etopo1')

	call cpy_var (ncid, 'inv_bar_corr', 'inv_bar_static')
	call cpy_var (ncid, 'inv_bar_corr hf_fluctuations_corr ADD', 'inv_bar_mog2d')

	call cpy_var (ncid, 'ocean_tide_sol1 load_tide_sol1 SUB', 'tide_ocean_got00')
	call cpy_var (ncid, 'ocean_tide_sol2 load_tide_sol2 SUB', 'tide_ocean_fes04')
	call cpy_var (ncid, 'ocean_tide_eq_lp', 'tide_equil')
	call cpy_var (ncid, 'ocean_tide_neq_lp', 'tide_non_equil')
	call cpy_var (ncid, 'load_tide_sol1', 'tide_load_got00')
	call cpy_var (ncid, 'load_tide_sol2', 'tide_load_fes04')
	call cpy_var (ncid, 'solid_earth_tide', 'tide_solid')
	call cpy_var (ncid, 'pole_tide', 'tide_pole')

	call cpy_var (ncid, 'wind_speed_model_u', 'wind_speed_ecmwf_u')
	call cpy_var (ncid, 'wind_speed_model_v', 'wind_speed_ecmwf_v')
	call cpy_var (ncid, 'wind_speed_alt')
	call cpy_var (ncid, 'wind_speed_rad')
	call cpy_var (ncid, 'rad_water_vapor', 'water_vapor_rad')
	call cpy_var (ncid, 'rad_liquid_water', 'liquid_water_rad')

	call cpy_var (ncid, 'ssha')

! Dump the data

	call nfs(nf90_close(ncid))
	call put_rads
	deallocate (a, flags)

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Write HY-2A IGDR data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_Sentinel3_file_names')
write (*,1310)
1310 format (/ &
'This program converts Sentinel-3 NRT/STC/NTC files to RADS data' / &
'files with the name $RADSDATAROOT/data/2a/a/pPPPP/2apPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_h2
