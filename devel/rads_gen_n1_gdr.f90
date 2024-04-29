!-----------------------------------------------------------------------
! Copyright (c) 2011-2024  Remko Scharroo and Eric Leuliette
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

!*rads_gen_n1_gdr -- Converts CryoSat Ocean Products data to RADS
!+
program rads_gen_n1_gdr

! This program reads CryoSat NOP/IOP/GOP files and converts them to the RADS
! format, written into files $RADSDATAROOT/data/c2.[nig]op/a/c2pPPPPcCCC.nc.
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_n1_gdr [options] < list_of_Envisat_GDR_file_names
!
! where [options] include:
!  --min-rec <min_rec> : Specify minimum number of records per pass to process.
!
! This program handles CryoSat Ocean Products standard_measurement files in
! netCDF format. The format is described in:
!
! [1] ENVISAT-1 PRODUCTS SPECIFICATIONS
! VOLUME 14: RA2 PRODUCTS SPECIFICATIONS LEVEL 2
! ESA Doc Ref PO-RS-MDA-GS-2009 Issue: 5i 13 April 2018
!
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdre - Orbital altitude
! alt_rate - Orbital altitude rate
! range_* - Ocean range (retracked)
! range_rms_* - Std dev of range
! range_numval_* - Nr of averaged range measurements
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! wet_tropo_comp - Composite wet tropo correction
! iono_alt_* - Dual-frequency ionospheric correction (not _c)
! iono_gim - GIM ionosphetic correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! peakiness_ku_* - Peakiness
! ssb_cls_* - SSB
! swh_* - Significant wave height
! swh_rms_* - Std dev of SWH
! sig0_* - Sigma0 (not corrected for attenuation ... will be done in rads_fix_s3)
! sig0_rms_* - Std dev of sigma0
! wind_speed_alt_* - Altimeter wind speed (not _c)
! wind_speed_rad - Radiometer wind speed
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! tide_ocean/load_got410 - GOT4.10c ocean and load tide
! tide_ocean/load_fes04/fes14 - FES2014 ocean and load tide
! tide_pole - Pole tide
! tide_solid - Solid earth tide
! topo_ace2 - ACE2 topography
! geoid_egm2008 - EGM2008 geoid
! mss_cnescls11/cnescls15 - CNES/CLS11 or CNES/CLS15 mean sea surface
! mss_dtu13 - DTU13 mean sea surface
! flags - Engineering flags
! off_nadir_angle2_wf_* - Mispointing from waveform squared (not _c)
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
use rads_devel_misc
use rads_gen
use rads_misc
use rads_netcdf
use rads_time
use rads_geo
use netcdf

! Command line arguments

integer(fourbyteint) :: ios, i
character(len=rads_cmdl) :: filename, arg

! Header variables

integer(fourbyteint) :: pass_number, cycle_number, ncid, varid
real(eightbytereal) :: equator_time
logical :: s_ok

! Data variables

integer :: latency = rads_ntc + 2

! Struct for orbit info

integer(fourbyteint) :: ipass
integer(fourbyteint), parameter :: mpass = 120000 ! Enough for 12 years
type(orfinfo) :: orf(mpass)


! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000
real(eightbytereal), allocatable :: dh(:)

! Initialise

call synopsis
call rads_gen_getopt ('', ' min-rec:')
call synopsis ('--head')
call rads_init (S, sat)

! Load the ORF file

call read_orf ('EN1', orf)
ipass = 1

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

do
	read (*,'(a)',iostat=ios) filename
	if (ios /= 0) exit

! Open input file

	call log_string (basename(filename))
	if (nf90_open(filename,nf90_nowrite,ncid) /= nf90_noerr) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Check if input is a Cryosat Ocean Products Level 2 data set

	if (nf90_get_att(ncid,nf90_global,'title',arg) /= nf90_noerr .or. &
		arg(:23) /= 'Envisat RA2/MWR Level 2') then
		call log_string ('Error: this is not an Envisat RA2/MWR Level 2 data set', .true.)
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

! Allocate variables

	allocate (a(nrec),dh(nrec),flags(nrec))
	nvar = 0

! Cycle and pass numbers

	call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cycle_number))
	call nfs(nf90_get_att(ncid,nf90_global,'pass_number',pass_number))
	call get_var (ncid, 'time_01', a)
	call which_pass ((a(1) + a(nrec)) / 2)

	if (orf(ipass)%cycle /= cycle_number) write (*,630) 'cycle', orf(ipass)%cycle, cycle_number, trim(filename)
	if (orf(ipass)%pass /= pass_number) write (*,630) ' pass', orf(ipass)%pass, pass_number, trim(filename)
630 format ('Warning: ',a,' number; output:',i4,'; input:',i4,' from ',a)

	cycle_number = orf(ipass)%cycle
	pass_number = orf(ipass)%pass
	equator_time = orf(ipass)%eqtime + sec2000

! Skip passes of which the cycle number or equator crossing time is outside the specified interval

	if (equator_time < times(1) .or. equator_time > times(2) .or. cycle_number < cycles(1) .or. cycle_number > cycles(2)) then
		call nfs(nf90_close(ncid))
		call log_string ('Skipped', .true.)
		cycle
	endif

! Set mission phase based on equator_time

	call rads_set_phase (S, equator_time)

! Store relevant pass info

	call rads_init_pass_struct (S, P)
	P%cycle = cycle_number
	P%pass = pass_number
	P%equator_time = equator_time
	P%equator_lon = orf(ipass)%eqlon
	P%start_time = a(1) + sec2000
	P%end_time = a(nrec) + sec2000
	P%original = trim(basename(filename))

! Switch s_ok off for the whole Side-B period

    i = cycle_number*10000 + pass_number
    if (i >= 470794 .and. i < 470937) then
    	! Switch to side B on 15 May 2006 14:21:50
    	s_ok = .false.
    else if (i >= 470937 .and. i <= 480847) then
        ! S-band failure in Side B on 20 May 2006 14:08:06
        s_ok = .false.
    else if (i >= 650290) then
        ! S-band failure in Side A on 17 Jan 2008 23:23:40
        s_ok = .false.
    else
        s_ok = .true.
    endif

! Compile flag bits

	flags = 0
	call nc2f (ncid, 'off_nadir_angle_wf_ocean_qual_01_ku', 1)	! bit  1: Quality off-nadir pointing
	call nc2f (ncid, 'surf_class_01', 2, eq=4)					! bit  2: Continental ice
	call nc2f (ncid, 'range_ocean_qual_01_ku', 3)
	call nc2f (ncid, 'range_ocean_qual_01_s',  3)				! bit  3: Quality dual-freq iono
	call nc2f (ncid, 'surf_class_01', 4, eq=1)
	call nc2f (ncid, 'surf_class_01', 4, ge=3)					! bit  4: Water/land
	call nc2f (ncid, 'surf_class_01', 5, ge=1)					! bit  5: Ocean/other
	call nc2f (ncid, 'rad_surf_type_01', 6)						! bit  6: Radiometer land flag
	call nc2f (ncid, 'rain_flag_01_ku', 7, eq=1)
	call nc2f (ncid, 'rain_flag_01_ku', 7, eq=2)
	call nc2f (ncid, 'rain_flag_01_ku', 7, eq=4)				! bit  7: Altimeter rain or ice flag
	call nc2f (ncid, 'mcd_flag_tb_238_01', 8)					! bit  8: Quality TB23.8
	call nc2f (ncid, 'mcd_flag_tb_365_01', 9)					! bit  9: Quality TB36.5
	call nc2f (ncid, 'range_ocean_qual_01_ku', 11)				! bit 11: Quality range
	call nc2f (ncid, 'mcd_flag_retracking_ocean_01_ku', 12)		! bit 12: Quality SWH
	call nc2f (ncid, 'mcd_flag_retracking_ocean_01_ku', 13)		! bit 13: Quality Sigma0
	call nc2f (ncid, 'nominal_tracking_01', 14)					! bit 14: Tracking mode
	call nc2f (ncid, 'mcd_flag_orb_state_rest_01', 15, ge=4)	! bit 15: Quality orbit

! Convert all the necessary fields to RADS

! Tine and location

	call get_var (ncid, 'time_01', a)
	call new_var ('time', a + sec2000)
	call cpy_var (ncid, 'lat_01', 'lat')
! Compute ellipsoid corrections
	do i = 1,nrec
		dh(i) = dhellips(1,a(i))
	enddo
	call cpy_var (ncid, 'lon_01', 'lon')
	call get_var (ncid, 'alt_01', a)
	call new_var ('alt_gdrf', a + dh)
	call cpy_var (ncid, 'orb_alt_rate_01', 'alt_rate')

! Range

	call cpy_var (ncid, 'range_ocean_01_ku', 'range_ku')
	call cpy_var (ncid, 'range_ocean_rms_01_ku', 'range_rms_ku')
	call cpy_var (ncid, 'range_ocean_numval_01_ku', 'range_numval')
	call cpy_var (ncid, 'range_ocean_qual_01_ku', 'qual_range')

	if (s_ok) then
		call cpy_var (ncid, 'range_ocean_01_s', 'range_s')
		call cpy_var (ncid, 'range_ocean_rms_01_s', 'range_rms_s')
		call cpy_var (ncid, 'range_ocean_numval_01_s', 'range_numval_s')
		call cpy_var (ncid, 'range_ocean_qual_01_s', 'qual_range_s')
	endif

! Sigma0

	call cpy_var (ncid, 'sig0_ocean_01_ku', 'sig0_ku')
	call cpy_var (ncid, 'sig0_ocean_rms_01_ku', 'sig0_rms_ku')
	call cpy_var (ncid, 'sig0_ocean_qual_01_ku', 'qual_sig0')

	if (s_ok) then
		call cpy_var (ncid, 'sig0_ocean_01_s', 'sig0_s')
		call cpy_var (ncid, 'sig0_ocean_rms_01_s', 'sig0_rms_s')
		call cpy_var (ncid, 'sig0_ocean_qual_01_s', 'qual_sig0_s')
	endif

! SWH

	call get_var (ncid, 'square_swh_ocean_01_ku', a)
	where (a < 0)
		a = -sqrt(-a)
	elsewhere
		a = sqrt(a)
	endwhere
	call new_var ('swh_ku', a)
	call cpy_var (ncid, 'swh_ocean_rms_01_ku', 'swh_rms_ku')
	call cpy_var (ncid, 'swh_ocean_qual_01_ku', 'qual_swh')

	if (s_ok) then
		call get_var (ncid, 'square_swh_ocean_01_s', a)
		where (a < 0)
			a = -sqrt(-a)
		elsewhere
			a = sqrt(a)
		endwhere
		call new_var ('swh_s', a)
		call cpy_var (ncid, 'swh_ocean_rms_01_s', 'swh_rms_s')
		call cpy_var (ncid, 'swh_ocean_qual_01_s', 'qual_swh_s')
	endif

! SSHA

	call cpy_var (ncid, 'ssha_01_ku', 'ssha')

! Instrument corrections

	call cpy_var (ncid, 'uso_cor_01', 'drange_uso')

! Environmental corrections

	call cpy_var (ncid, 'mod_dry_tropo_cor_01', 'dry_tropo_ecmwf')
	call cpy_var (ncid, 'mod_dry_tropo_cor_reanalysis_01', 'dry_tropo_era')

	call cpy_var (ncid, 'mod_wet_tropo_cor_01', 'wet_tropo_ecmwf')
	call cpy_var (ncid, 'mod_wet_tropo_cor_reanalysis_01', 'wet_tropo_era')
	call cpy_var (ncid, 'rad_wet_tropo_cor_01', 'wet_tropo_rad')
	call cpy_var (ncid, 'rad_wet_tropo_cor_sst_gam_01', 'wet_tropo_rad_sst_gam')

	if (s_ok) then
		call cpy_var (ncid, 'iono_cor_alt_01_ku', 'iono_alt')
		call cpy_var (ncid, 'filtered_iono_cor_alt_01_ku', 'iono_alt_smooth')
	endif
	call cpy_var (ncid, 'iono_cor_gim_01_ku', 'iono_gim')

! SSB

	call cpy_var (ncid, 'sea_state_bias_01_ku', 'ssb_cls')
	if (s_ok) call cpy_var (ncid, 'sea_state_bias_01_s', 'ssb_cls_s')

! Atmospheric attenuation

	call cpy_var (ncid, 'atm_cor_sig0_01_ku', 'dsig0_atmos_ku')
	if (s_ok) call cpy_var (ncid, 'atm_cor_sig0_01_s', 'dsig0_atmos_s')

! MSS, MDT, Geoid, elevation

	call get_var (ncid, 'mean_sea_surf_sol1_01', a)
	call new_var ('mss_cnescls15', a + dh)
	call get_var (ncid, 'mean_sea_surf_sol2_01', a)
	call new_var ('mss_dtu15', a + dh)
	!call cpy_var (ncid, 'mean_dyn_topo_01', 'mdt_cnescls13')
	call cpy_var (ncid, 'geoid_01', 'geoid_egm2008')
	call cpy_var (ncid, 'odle_01', 'topo_ace2')

! DAC

	call cpy_var (ncid, 'inv_bar_cor_01', 'inv_bar_static')
	call cpy_var (ncid, 'inv_bar_cor_reanalysis_01', 'inv_bar_static_era')
	call cpy_var (ncid, 'inv_bar_cor_01 hf_fluct_cor_01 ADD', 'inv_bar_mog2d')
	call cpy_var (ncid, 'inv_bar_cor_reanalysis_01 hf_fluct_cor_reanalysis_01 ADD', 'inv_bar_mog2d_era')

! Tides

	call cpy_var (ncid, 'ocean_tide_sol1_01 load_tide_sol1_01 SUB', 'tide_ocean_got410')
	call cpy_var (ncid, 'ocean_tide_sol2_01 load_tide_sol2_01 SUB ocean_tide_non_eq_01 ADD', 'tide_ocean_fes14')
	call cpy_var (ncid, 'ocean_tide_eq_01', 'tide_equil')
	call cpy_var (ncid, 'ocean_tide_non_eq_01', 'tide_non_equil')

	call cpy_var (ncid, 'load_tide_sol1_01', 'tide_load_got410')
	call cpy_var (ncid, 'load_tide_sol2_01', 'tide_load_fes14')

	call cpy_var (ncid, 'solid_earth_tide_01', 'tide_solid')
	call cpy_var (ncid, 'pole_tide_01', 'tide_pole')

! Wind speed

	call cpy_var (ncid, 'wind_speed_mod_u_01', 'wind_speed_ecmwf_u')
!	call cpy_var (ncid, 'wind_speed_mod_u_reanalysis_01', 'wind_speed_era_u')
	call cpy_var (ncid, 'wind_speed_mod_v_01', 'wind_speed_ecmwf_v')
!	call cpy_var (ncid, 'wind_speed_mod_v_reanalysis_01', 'wind_speed_era_v')

	call cpy_var (ncid, 'wind_speed_alt_01_ku', 'wind_speed_alt')

! Radiometer

	call cpy_var (ncid, 'rad_liquid_water_01', 'liquid_water_rad')

! Mispointing

	call cpy_var (ncid, 'off_nadir_angle_wf_ocean_01_ku', 'off_nadir_angle2_wf_ku')
	call cpy_var (ncid, 'off_nadir_angle_wf_ocean_qual_01_ku', 'qual_attitude')

! Brightness temperatures

	call cpy_var (ncid, 'tb_238_01', 'tb_238')
	call cpy_var (ncid, 'tb_365_01', 'tb_365')

! Flags

	if (s_ok) then
		call cpy_var (ncid, 'rain_flag_01_ku', 'qual_alt_rain_ice')
		call cpy_var (ncid, 'open_sea_ice_flag_01_ku', 'seaice_class')
	endif
	call new_var ('flags', dble(flags))

! Peakiness

	call cpy_var (ncid, 'peakiness_01_ku', 'peakiness_ku')
	if (s_ok) call cpy_var (ncid, 'peakiness_01_s', 'peakiness_s')

	call cpy_var (ncid, 'gpd_wet_tropo_cor_01', 'gpd_wet_tropo_cor')

	a = latency
	call new_var ('latency', a)

! Dump the data

	call nfs(nf90_close(ncid))
	call put_rads
	deallocate (a, dh, flags)

enddo

call rads_end (S)

contains

!***********************************************************************
! Determine the corresponding record from ORF

subroutine which_pass (time)
real(eightbytereal), intent(in) :: time
do while (time < orf(ipass)%starttime)
	ipass = ipass - 1
	if (ipass < 1) call rads_exit ('Time is before the start of the ORF file')
enddo
do while (time > orf(ipass+1)%starttime)
	ipass = ipass + 1
	if (orf(ipass)%cycle < 0) call rads_exit ('Time is after the end of the ORF file')
enddo
end subroutine which_pass

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Write Envisat GDR v3 Products data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_Envisat_file_names')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --min-rec=MIN_REC         Specify minimum number of records per pass to process' // &
'This program converts Envisat GDR v3 files to RADS data' / &
'files with the name $RADSDATAROOT/data/SS/a/pPPPP/SSpPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_n1_gdr
