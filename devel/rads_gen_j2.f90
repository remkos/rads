!-----------------------------------------------------------------------
! Copyright (c) 2011-2016  Remko Scharroo
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

!*rads_gen_j2 -- Converts Jason-2 data to RADS
!+
program rads_gen_j2

! This program reads Jason-2 (O/I)GDR files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/j2/a/j2pPPPPcCCC.nc.
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_j2 [options] < list_of_JASON2_file_names
!
! This program handles Jason-2 OGDR, IGDR and GDR files in netCDF format,
! version D (also known as GDR-D).
! The format is described in:
!
! [1] OSTM/Jason-2 Products Handbook, SALP-MU-M-OP-15815-CN (CNES),
!     EUM/OPS-JAS/MAN/08/0041 (EUMETSAT), OSTM-29-1237 (JPL),
!     Version 1.8, 1 Dec 2011.
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdrd/alt_gdre - Orbit altitude
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
! tb_187 - Brightness temperature (18.7 GHz)
! tb_238 - Brightness temperature (23.8 GHz)
! tb_340 - Brightness temperature (34.0 GHz)
! flags, flags_mle3 - Engineering flags
! swh_rms_* - Std dev of SWH
! sig0_rms_* - Std dev of sigma0
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! rad_liquid_water - Liquid water content
! rad_water_vapor - Water vapor content
!
! Extionsions _* are:
! _ku:      Ku-band retracked with MLE4
! _ku_mle3: Ku-band retracked with MLE3
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

integer(fourbyteint) :: i, j, ios
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: cyclenr, passnr, varid
real(eightbytereal) :: equator_time
logical :: ogdr

! Data variables

integer(twobyteint), allocatable :: flags_mle3(:), flags_save(:)

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000

! Initialise

call synopsis
call rads_gen_getopt ('j2')
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

! Check if input is GDR-D

	if (index(infile,'_2Pd') <= 0) then
		call log_string ('Error: this is not GDR-D', .true.)
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
	call nfs(nf90_get_att(ncid,nf90_global,'mission_name',arg))
	if (arg /= 'OSTM/Jason-2') then
		call log_string ('Error: wrong misson-name found in header', .true.)
		cycle
	endif

	call nfs(nf90_get_att(ncid,nf90_global,'title',arg))
	ogdr = (arg(:4) == 'OGDR')
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
	j = index(arg, 'L2 library=')
	P%original = trim(infile(i:)) // ' (' // arg(j+11:)
	j = index(P%original, ',')
	P%original(j:) = ')'

! Allocate variables

	allocate (a(nrec),flags(nrec),flags_mle3(nrec),flags_save(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
	call nc2f ('alt_state_flag_oper',0)				! bit  0: Altimeter Side A/B
	call nc2f ('qual_alt_1hz_off_nadir_angle_wf_ku',1)	! bit  1: Quality off-nadir pointing
	call nc2f ('surface_type',2,val=2)				! bit  2: Continental ice
	call nc2f ('qual_alt_1hz_range_c',3)			! bit  3: Quality dual-frequency iono
	call nc2f ('surface_type',4,lim=2)				! bit  4: Water/land
	call nc2f ('surface_type',5,lim=1)				! bit  5: Ocean/other
	call nc2f ('rad_surf_type',6,lim=2)				! bit  6: Radiometer land flag
	call nc2f ('rain_flag',7)
	call nc2f ('ice_flag',7)						! bit  7: Altimeter rain or ice flag
	call nc2f ('rad_rain_flag',8)
	call nc2f ('rad_sea_ice_flag',8)				! bit  8: Radiometer rain or ice flag
	call nc2f ('qual_rad_1hz_tb187',9)
	call nc2f ('qual_rad_1hz_tb238',9)				! bit  9: Quality 18.7 and 23.8 GHz channel
	call nc2f ('qual_rad_1hz_tb340',10)				! bit 10: Quality 34.0 GHz channel

! Now do specifics for MLE3 (not used for the time being)

	flags_save = flags	! Keep flags for later
	call nc2f ('qual_alt_1hz_range_ku_mle3',3)		! bit  3: Quality dual-frequency iono
	call nc2f ('qual_alt_1hz_range_ku_mle3',11)		! bit 11: Quality range
	call nc2f ('qual_alt_1hz_swh_ku_mle3',12)		! bit 12: Quality SWH
	call nc2f ('qual_alt_1hz_sig0_ku_mle3',13)		! bit 13: Quality Sigma0
	flags_mle3 = flags	! Copy result for MLE3
	flags = flags_save	! Continue with MLE4 flags

! Redo the last ones for MLE4

	call nc2f ('qual_alt_1hz_range_ku',3)			! bit  3: Quality dual-frequency iono
	call nc2f ('qual_alt_1hz_range_ku',11)			! bit 11: Quality range
	call nc2f ('qual_alt_1hz_swh_ku',12)			! bit 12: Quality SWH
	call nc2f ('qual_alt_1hz_sig0_ku',13)			! bit 13: Quality Sigma0

! Convert all the necessary fields to RADS

	call get_var (ncid, 'time', a)
	call new_var ('time', a + sec2000)
	call cpy_var ('lat')
	call cpy_var ('lon')
	if (cyclenr < 254) then ! Before Cycle 254 was GDR-D orbit
		call cpy_var ('alt', 'alt_gdrd')
	else
		call cpy_var ('alt', 'alt_gdre')
	endif
	call cpy_var ('orb_alt_rate', 'alt_rate')
	call cpy_var ('range_ku')
	call cpy_var ('range_ku_mle3')
	call cpy_var ('range_c')
	call cpy_var ('model_dry_tropo_corr', 'dry_tropo_ecmwf')
	call cpy_var ('rad_wet_tropo_corr', 'wet_tropo_rad')
	call cpy_var ('model_wet_tropo_corr', 'wet_tropo_ecmwf')
	call cpy_var ('iono_corr_alt_ku', 'iono_alt')
	call cpy_var ('iono_corr_alt_ku_mle3', 'iono_alt_mle3')
	if (.not.ogdr) call cpy_var ('iono_corr_gim_ku', 'iono_gim')
	call cpy_var ('inv_bar_corr', 'inv_bar_static')
	if (ogdr) then
		call cpy_var ('inv_bar_corr', 'inv_bar_mog2d')
	else
		call cpy_var ('inv_bar_corr hf_fluctuations_corr ADD', 'inv_bar_mog2d')
	endif
	call cpy_var ('solid_earth_tide', 'tide_solid')
	call cpy_var ('ocean_tide_sol1 load_tide_sol1 SUB', 'tide_ocean_got48')
	call cpy_var ('ocean_tide_sol2 load_tide_sol2 SUB', 'tide_ocean_fes04')
	call cpy_var ('load_tide_sol1', 'tide_load_got48')
	call cpy_var ('load_tide_sol2', 'tide_load_fes04')
	call cpy_var ('pole_tide', 'tide_pole')
	call cpy_var ('sea_state_bias_ku', 'ssb_cls')
	call cpy_var ('sea_state_bias_ku_mle3', 'ssb_cls_mle3')
	call cpy_var ('sea_state_bias_c', 'ssb_cls_c')
	call cpy_var ('geoid', 'geoid_egm96')
	call cpy_var ('mean_sea_surface', 'mss_cnescls11')
	call cpy_var ('swh_ku')
	call cpy_var ('swh_ku_mle3')
	call cpy_var ('swh_c')
	call cpy_var ('sig0_ku')
	call cpy_var ('sig0_ku_mle3')
	call cpy_var ('sig0_c')
	call cpy_var ('wind_speed_alt')
	call cpy_var ('wind_speed_alt_mle3')
	call cpy_var ('wind_speed_rad')
	call cpy_var ('wind_speed_model_u', 'wind_speed_ecmwf_u')
	call cpy_var ('wind_speed_model_v', 'wind_speed_ecmwf_v')
	call cpy_var ('range_rms_ku')
	call cpy_var ('range_rms_ku_mle3')
	call cpy_var ('range_rms_c')
	call cpy_var ('range_numval_ku')
	call cpy_var ('range_numval_ku_mle3')
	call cpy_var ('range_numval_c')
	call cpy_var ('bathymetry', 'topo_dtm2000')
	call cpy_var ('tb_187')
	call cpy_var ('tb_238')
	call cpy_var ('tb_340')
	a = flags
	call new_var ('flags', a)
	a = flags_mle3
	call new_var ('flags_mle3', a)
	call cpy_var ('swh_rms_ku')
	call cpy_var ('swh_rms_ku_mle3')
	call cpy_var ('swh_rms_c')
	call cpy_var ('sig0_rms_ku')
	call cpy_var ('sig0_rms_ku_mle3')
	call cpy_var ('sig0_rms_c')
	call cpy_var ('off_nadir_angle_wf_ku', 'off_nadir_angle2_wf_ku')
	call cpy_var ('atmos_corr_sig0_ku', 'dsig0_atmos_ku')
	call cpy_var ('atmos_corr_sig0_c', 'dsig0_atmos_c')
	call cpy_var ('rad_liquid_water', 'liquid_water_rad')
	call cpy_var ('rad_water_vapor', 'water_vapor_rad')
	call cpy_var ('ssha', 'ssha')
	call cpy_var ('ssha_mle3', 'ssha_mle3')

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
if (rads_version ('Write Jason-2 data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_JASON2_file_names')
write (*,1310)
1310 format (/ &
'This program converts Jason-2 OGDR/IGDR/GDR files to RADS data' / &
'files with the name $RADSDATAROOT/data/j2/a/pPPPP/j2pPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_j2
