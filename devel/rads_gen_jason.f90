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

!*rads_gen_jason -- Converts Jason-1/2/3 data to RADS
!+
program rads_gen_jason

! This program reads Jason-1/2/3 (O/I)GDR files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/jJ/a/jJpPPPPcCCC.nc.
!    jJ = Jason satellite abbreviation
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_jason [options] < list_of_JASON3_file_names
!
! This program handles Jason-1, -2 and -3 OGDR, IGDR and GDR files in NetCDF format.
! The format is described in:
!
! [1] IGDR and GDR Jason Products, AVISO and PODAAC User Handbook,
!     SMM-MU-M5-OP-13184-CN (AVISO), JPL D-21352 (PODAAC),
!     Edition 4.0, June 2008
!
! [2] OSTM/Jason-2 Products Handbook, SALP-MU-M-OP-15815-CN (CNES),
!     EUM/OPS-JAS/MAN/08/0041 (EUMETSAT), OSTM-29-1237 (JPL),
!     Version 1.8, 1 Dec 2011.
!
! [3] Jason-3 Products Handbook, SALP-MU-M-OP-16118-CN (CNES)
!     Version 1.0, 9 July 2012
!
! [4] SALP Products Specification - Volume 30: Jason-3 Products
!     SALP-ST-M-EA-16122-CN, Version 1.2, 9 December 2013
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdrd, alt_gdre - Orbit altitude
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
! Extensions _* are:
! _ku:      Ku-band retracked with MLE4
! _ku_mle3: Ku-band retracked with MLE3 (not for Jason-1)
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

integer(fourbyteint) :: ios, j, q, r
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: ncid, cyclenr, passnr, varid
real(eightbytereal) :: equator_time
logical :: gdre = .false., mle3 = .true.
integer :: latency = rads_nrt

! Data variables

integer(twobyteint), allocatable :: flags_mle3(:), flags_save(:)

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000

! Initialise

call synopsis
call rads_gen_getopt ('')
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
	if (nf90_open(infile,nf90_nowrite,ncid) /= nf90_noerr) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Check if input is GDR-T/D or GDR-E

	if (index(infile,'_2PT') > 0 .or. index(infile,'_2Pd') > 0) then
		gdre = .false.
	else if (index(infile,'_2Pe') > 0) then
		gdre = .true.
	else
		call log_string ('Error: this is not GDR-T, GDR-d, or GDR-e', .true.)
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

! Read rest of global attributes

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

! Update mission phase if required

	call rads_set_phase (S, equator_time)

! Use special conversion of cycle numbers for geodetic mission

	if (S%sat == 'j1' .and. S%phase%name(1:1) == 'c') then
		! Redetermine subcycle numbering from the GDR cycle and pass numbers
		r = (cyclenr - 500) * 280 + (passnr - 1) ! GDRs of Phase C start with Cycle 500 and have 280 passes per cycle
		q = (r/10438) * 43
		r = modulo(r,10438)
		q = q + (r/4608) * 19
		r = modulo(r,4608)
		q = q + (r/1222) * 5
		r = modulo(r,1222)
		q = q + (r/280)
		r = modulo(r,280)
		cyclenr = q + 382
		passnr = r + 1
	else if (S%sat == 'j2' .and. S%phase%name(1:1) == 'c') then
		! Redetermine subcycle numbering from the GDR cycle and pass numbers
		r = (cyclenr - 500) * 254 + (passnr - 1) ! GDRs of Phase C start with Cycle 500 and have 254 passes per cycle
		q = (r/9472) * 23
		r = modulo(r,9472)
		q = q + (r/3702) * 9
		r = modulo (r,3702)
		q = q + (r/2068) * 5
		r = modulo(r,2068)
		q = q + (r/434)
		r = modulo(r,434)
		cyclenr = q + 332
		passnr = r + 1
	else if (S%sat == 'j2' .and. S%phase%name(1:1) == 'd') then
		! Redetermine subcycle numbering from the GDR cycle and pass numbers
		r = (cyclenr - 600) * 254 + (passnr - 1) ! GDRs of Phase D start with Cycle 600 and have 254 passes per cycle
		r = r + 142 ! This is to shift the pass numbering in order to align with Phase C
		q = (r/9472) * 23
		r = modulo(r,9472)
		q = q + (r/3702) * 9
		r = modulo (r,3702)
		q = q + (r/2068) * 5
		r = modulo(r,2068)
		q = q + (r/434)
		r = modulo(r,434)
		cyclenr = q + 356
		passnr = r + 1
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
	j = index(arg, 'L2 library=')
	P%original = trim(basename(infile)) // ' (' // arg(j+11:)
	j = index(P%original, ',')
	P%original(j:) = ')'

! Allocate variables

	allocate (a(nrec),flags(nrec),flags_mle3(nrec),flags_save(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
	call nc2f (ncid, 'alt_state_flag_oper',0)				! bit  0: Altimeter Side A/B
	call nc2f (ncid, 'qual_alt_1hz_off_nadir_angle_wf_ku',1)	! bit  1: Quality off-nadir pointing
	call nc2f (ncid, 'surface_type',2,eq=2)				! bit  2: Continental ice
	call nc2f (ncid, 'qual_alt_1hz_range_c',3)			! bit  3: Quality dual-frequency iono
	call nc2f (ncid, 'surface_type',4,ge=2)				! bit  4: Water/land
	call nc2f (ncid, 'surface_type',5,ge=1)				! bit  5: Ocean/other
	call nc2f (ncid, 'rad_surf_type',6,ge=2)				! bit  6: Radiometer land flag
	call nc2f (ncid, 'rain_flag',7)
	call nc2f (ncid, 'ice_flag',7)						! bit  7: Altimeter rain or ice flag
	call nc2f (ncid, 'rad_rain_flag',8)
	call nc2f (ncid, 'rad_sea_ice_flag',8)				! bit  8: Radiometer rain or ice flag
	call nc2f (ncid, 'qual_rad_1hz_tb187',9)
	call nc2f (ncid, 'qual_rad_1hz_tb238',9)				! bit  9: Quality 18.7 and 23.8 GHz channel
	call nc2f (ncid, 'qual_rad_1hz_tb340',10)				! bit 10: Quality 34.0 GHz channel

! Now do specifics for MLE3 (not used for the time being)

	if (mle3) then
		flags_save = flags	! Keep flags for later
		call nc2f (ncid, 'qual_alt_1hz_range_ku_mle3',3)		! bit  3: Quality dual-frequency iono
		call nc2f (ncid, 'qual_alt_1hz_range_ku_mle3',11)		! bit 11: Quality range
		call nc2f (ncid, 'qual_alt_1hz_swh_ku_mle3',12)		! bit 12: Quality SWH
		call nc2f (ncid, 'qual_alt_1hz_sig0_ku_mle3',13)		! bit 13: Quality Sigma0
		flags_mle3 = flags	! Copy result for MLE3
		flags = flags_save	! Continue with MLE4 flags
	endif

! Redo the last ones for MLE4

	call nc2f (ncid, 'qual_alt_1hz_range_ku',3)			! bit  3: Quality dual-frequency iono
	call nc2f (ncid, 'qual_alt_1hz_range_ku',11)			! bit 11: Quality range
	call nc2f (ncid, 'qual_alt_1hz_swh_ku',12)			! bit 12: Quality SWH
	call nc2f (ncid, 'qual_alt_1hz_sig0_ku',13)			! bit 13: Quality Sigma0

! Convert all the necessary fields to RADS

	call get_var (ncid, 'time', a)
	call new_var ('time', a + sec2000)
	call cpy_var (ncid, 'lat')
	call cpy_var (ncid, 'lon')
	if (S%sat == 'j2' .and. cyclenr < 254) then ! Jason-2 before Cycle 254 was GDR-D orbit
		call cpy_var (ncid, 'alt', 'alt_gdrd')
	else ! All Jason-1 and Jason-3 data, as well as more recent Jason-2 data
		call cpy_var (ncid, 'alt', 'alt_gdre')
	endif
	call cpy_var (ncid, 'orb_alt_rate', 'alt_rate')
	call cpy_var (ncid, 'range_ku')
	call cpy_var (ncid, 'range_ku_mle3', mle3)
	call cpy_var (ncid, 'range_c')
	call cpy_var (ncid, 'model_dry_tropo_corr', 'dry_tropo_ecmwf')
	call cpy_var (ncid, 'rad_wet_tropo_corr', 'wet_tropo_rad')
	call cpy_var (ncid, 'model_wet_tropo_corr', 'wet_tropo_ecmwf')
	call cpy_var (ncid, 'iono_corr_alt_ku', 'iono_alt')
	call cpy_var (ncid, 'iono_corr_alt_ku_mle3', 'iono_alt_mle3', mle3)
	call cpy_var (ncid, 'iono_corr_gim_ku', 'iono_gim', latency /= rads_nrt)
	call cpy_var (ncid, 'inv_bar_corr', 'inv_bar_static')
	if (latency == rads_nrt) then
		call cpy_var (ncid, 'inv_bar_corr', 'inv_bar_mog2d')
	else
		call cpy_var (ncid, 'inv_bar_corr hf_fluctuations_corr ADD', 'inv_bar_mog2d')
	endif
	call cpy_var (ncid, 'solid_earth_tide', 'tide_solid')
	if (gdre) then
		call cpy_var (ncid, 'ocean_tide_sol1 load_tide_sol1 SUB', 'tide_ocean_got410')
		call cpy_var (ncid, 'ocean_tide_sol2 load_tide_sol2 SUB ocean_tide_non_equil ADD', 'tide_ocean_fes14')
		call cpy_var (ncid, 'load_tide_sol1', 'tide_load_got410')
		call cpy_var (ncid, 'load_tide_sol2', 'tide_load_fes14')
!	else
!		call cpy_var (ncid, 'ocean_tide_sol1 load_tide_sol1 SUB', 'tide_ocean_got48')
!		call cpy_var (ncid, 'ocean_tide_sol2 load_tide_sol2 SUB ocean_tide_non_equil ADD', 'tide_ocean_fes04')
!		call cpy_var (ncid, 'load_tide_sol1', 'tide_load_got48')
!		call cpy_var (ncid, 'load_tide_sol2', 'tide_load_fes04')
	endif
	call cpy_var (ncid, 'ocean_tide_equil', 'tide_equil')
	call cpy_var (ncid, 'ocean_tide_non_equil', 'tide_non_equil')
	call cpy_var (ncid, 'pole_tide', 'tide_pole')
	call cpy_var (ncid, 'sea_state_bias_ku', 'ssb_cls')
	call cpy_var (ncid, 'sea_state_bias_ku_mle3', 'ssb_cls_mle3', mle3)
	call cpy_var (ncid, 'sea_state_bias_c', 'ssb_cls_c')
	if (gdre) then
		call cpy_var (ncid, 'geoid', 'geoid_egm2008')
!	else
!		call cpy_var (ncid, 'geoid', 'geoid_egm96')
	endif
	if (S%sat == 'j2' .and. S%phase%name >= 'c') then
		call cpy_var (ncid, 'mean_sea_surface', 'mss_cnescls15')
!	else
!		call cpy_var (ncid, 'mean_sea_surface', 'mss_cnescls11')
	endif
	call cpy_var (ncid, 'swh_ku')
	call cpy_var (ncid, 'swh_ku_mle3', mle3)
	call cpy_var (ncid, 'swh_c')
	call cpy_var (ncid, 'sig0_ku')
	call cpy_var (ncid, 'sig0_ku_mle3', mle3)
	call cpy_var (ncid, 'sig0_c')
	call cpy_var (ncid, 'wind_speed_alt')
	call cpy_var (ncid, 'wind_speed_alt_mle3', mle3)
	call cpy_var (ncid, 'wind_speed_rad')
	call cpy_var (ncid, 'wind_speed_model_u', 'wind_speed_ecmwf_u')
	call cpy_var (ncid, 'wind_speed_model_v', 'wind_speed_ecmwf_v')
	call cpy_var (ncid, 'range_rms_ku')
	call cpy_var (ncid, 'range_rms_ku_mle3', mle3)
	call cpy_var (ncid, 'range_rms_c')
	call cpy_var (ncid, 'range_numval_ku')
	call cpy_var (ncid, 'range_numval_ku_mle3', mle3)
	call cpy_var (ncid, 'range_numval_c')
!	call cpy_var (ncid, 'bathymetry', 'topo_dtm2000')
	call cpy_var (ncid, 'tb_187')
	call cpy_var (ncid, 'tb_238')
	call cpy_var (ncid, 'tb_340')
	a = flags
	call new_var ('flags', a)
	if (mle3) then
		a = flags_mle3
		call new_var ('flags_mle3', a)
	endif
	call cpy_var (ncid, 'swh_rms_ku')
	call cpy_var (ncid, 'swh_rms_ku_mle3', mle3)
	call cpy_var (ncid, 'swh_rms_c')
	call cpy_var (ncid, 'sig0_rms_ku')
	call cpy_var (ncid, 'sig0_rms_ku_mle3', mle3)
	call cpy_var (ncid, 'sig0_rms_c')
	call cpy_var (ncid, 'off_nadir_angle_wf_ku', 'off_nadir_angle2_wf_ku')
	call cpy_var (ncid, 'atmos_corr_sig0_ku', 'dsig0_atmos_ku')
	call cpy_var (ncid, 'atmos_corr_sig0_c', 'dsig0_atmos_c')
	call cpy_var (ncid, 'rad_liquid_water', 'liquid_water_rad')
	call cpy_var (ncid, 'rad_water_vapor', 'water_vapor_rad')
	call cpy_var (ncid, 'ssha')
	call cpy_var (ncid, 'ssha_mle3', mle3)
	a = latency
	call new_var ('latency', a)

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
if (rads_version ('Write Jason-1/2/3 data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_jason_file_names')
write (*,1310)
1310 format (/ &
'This program converts Jason-1/2/3 OGDR/IGDR/GDR files to RADS data' / &
'files with the name $RADSDATAROOT/data/jJ/a/pPPPP/jJpPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_jason
