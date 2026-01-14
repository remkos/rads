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

!*rads_gen_c2_op -- Converts CryoSat Ocean Products data to RADS
!+
program rads_gen_c2_op

! This program reads CryoSat NOP/IOP/GOP files from Baselines C and D
! and converts them to the RADS format, written into files
! $RADSDATAROOT/data/c2.[nig]op/a/c2pPPPPcCCC.nc.
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_c2_op [options] < list_of_CryoSat_OP_file_names
!
! where [options] include:
!  --min-rec <min_rec> : Specify minimum number of records per pass to process.
!
! This program handles CryoSat Ocean Products standard_measurement files in
! netCDF format. The format is described in:
!
!
! [1] Baseline-C CryoSat Ocean Processor: Ocean Product Handbook
! Version 4.1, ESA/ESRIN, 5 Dec 2019
!
! [2] Baseline-D CryoSat Ocean Processor: Ocean Product Handbook
! Version 5.3, ESA/ESRIN, 26 Sep 2024
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
! tide_ocean/load_fes14 - FES2014 ocean and load tide
! tide_pole - Pole tide
! tide_solid - Solid earth tide
! topo_ace2 - ACE2 topography
! geoid_egm2008 - EGM2008 geoid
! mss_cnescls15/mss_cnescls22 - CNES/CLS15 or CNES/CLS22 mean sea surface
! mss_dtu15/mss_dtu21 - DTU15 or DTU21 mean sea surface
! flags, flags_plrm - Engineering flags
! off_nadir_angle2_wf_* - Mispointing from waveform squared (not _c)
! ssha_* - Sea surface height anomaly (not _c)
!
! Extensions _* are:
! _ku:      Ku-band retracked from SAR
! _ku_plrm: Ku-band retracked with PLRM
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

! Struct for orbit info

integer(fourbyteint), parameter :: mpass = 300000 ! Enough for 30 years
type(orfinfo) :: orf(mpass)
integer :: ipass = 1

! Command line arguments

integer(fourbyteint) :: ios, i
character(len=rads_cmdl) :: filename, arg, product_name

! Header variables

integer(fourbyteint) :: cycle_number, ncid, varid
real(eightbytereal) :: equator_time, tai_utc_difference

! Data variables

integer(twobyteint), allocatable :: flags_plrm(:), flags_save(:)
character(len=4) :: software_version
integer :: latency = rads_nrt

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000
real(eightbytereal), allocatable :: dh(:)
type(quicksort_pair), allocatable :: pair(:)

! Initialise

call synopsis
call rads_gen_getopt ('', ' min-rec:')
call synopsis ('--head')
call rads_init (S, sat)

call read_orf ('CS2', orf)

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

	if (nf90_get_att(ncid,nf90_global,'mission',arg) /= nf90_noerr .or. &
		arg /= 'Cryosat') then
		call log_string ('Error: this is not a Cryosat Ocean Products Level 2 data set', .true.)
		cycle
	endif

! Read global attributes

	call nfs(nf90_inq_dimid(ncid, 'time_01',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
	call nfs(nf90_inq_varid(ncid, 'time_01',varid))
	call nfs(nf90_get_att(ncid,varid,'tai_utc_difference',tai_utc_difference))
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
	call nfs(nf90_get_att(ncid,nf90_global,'mission',arg))
	if (arg /= 'Cryosat') then
		call log_string ('Error: file skipped: wrong mission-name found in header', .true.)
		cycle
	endif

	! Determine latency (NOP, IOP, GOP) and software version

	call nfs(nf90_get_att(ncid,nf90_global,'product_name',product_name))
	if (index(product_name, '_NOP') > 0) then
		latency = rads_nrt
	else if (index(product_name, '_IOP') > 0) then
		latency = rads_stc
	else if (index(product_name, '_GOP') > 0) then
		latency = rads_ntc
	else
		call log_string ('Error: file skipped: unknown latency', .true.)
		cycle
	endif

! Try to read the equator crossing time set by rads_pre_sort_passes.
! If not, take the middle of the sensing start and stop times from the ESA products

	if (nff(nf90_get_att(ncid,nf90_global,'equator_time',arg))) then
		equator_time = strp1985f(arg)
	else
		equator_time = (strp1985f(product_name(20:34)) + strp1985f(product_name(36:50))) / 2
	endif

! Get the proper cycle and pass numbers and equator crossing information from the ORF file.

	call which_pass (equator_time - sec2000)
	cycle_number = orf(ipass)%cycle
	equator_time = orf(ipass)%eqtime + sec2000

! Skip passes of which the cycle number or equator crossing time is outside the specified interval

	if (equator_time < times(1) .or. equator_time > times(2) .or. cycle_number < cycles(1) .or. cycle_number > cycles(2)) then
		call nfs(nf90_close(ncid))
		call log_string ('Skipped', .true.)
		cycle
	endif

! Initialise pass struct and store pass info

	call rads_set_phase (S, equator_time)
	call rads_init_pass_struct (S, P)
	P%cycle = cycle_number
	P%pass = orf(ipass)%pass
	P%equator_time = equator_time
	P%equator_lon = orf(ipass)%eqlon

! Read start and end times. First try the format set by pre_sort_passes.

	if (nff(nf90_get_att(ncid,nf90_global,'first_meas_time',arg))) then
		P%start_time = strp1985f(arg)
		call nfs(nf90_get_att(ncid,nf90_global,'last_meas_time',arg))
		P%end_time = strp1985f(arg)
	else
		! Start time and end time in the ESA files are in TAI; need to be converted to UTC.
		call nfs(nf90_get_att(ncid,nf90_global,'first_record_time',arg))
		P%start_time = strp1985f(arg(5:)) + tai_utc_difference
		call nfs(nf90_get_att(ncid,nf90_global,'last_record_time',arg))
		P%end_time = strp1985f(arg(5:)) + tai_utc_difference
	endif

! Get sortware version

	call nfs(nf90_get_att(ncid,nf90_global,'software_version',arg))
	i = index(arg,'/')
	software_version = arg(i+1:i+4)

! Store input file name

	P%original = trim(basename(filename)) // ' (' // trim(arg) // ')'

! Allocate variables

	allocate (a(nrec),dh(nrec),flags(nrec),flags_plrm(nrec),flags_save(nrec),pair(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
!	call nc2f (ncid, 'flag_instr_op_mode_01',0,ge=1)			! bit  0: Altimeter mode
!	call nc2f ('surf_class_01',2,eq=4)				! bit  2: Continental ice
!	call nc2f ('surf_class_01',4,eq=1)
!	call nc2f ('surf_class_01',4,ge=3)				! bit  4: Water/land
!	call nc2f ('surf_class_01',5,ge=1)				! bit  5: Ocean/other

! Now do specifics for PLRM

	flags_save = flags	! Keep flags for later
	call nc2f (ncid, 'qual_ssha_01_plrm_ku',11)	! bit 11: Quality range
!	call nc2f (ncid, 'range_ocean_qual_01_plrm_ku',11)	! bit 11: Quality range
!	call nc2f (ncid, 'swh_ocean_qual_01_plrm_ku',12)		! bit 12: Quality SWH
!	call nc2f (ncid, 'sig0_ocean_qual_01_plrm_ku',13)		! bit 13: Quality Sigma0
!	flags_plrm = flags	! Copy result for PLRM
	flags = flags_save	! Continue with SAR flags

! Redo the last ones for SAR

	call nc2f (ncid, 'qual_ssha_01_ku',11)			! bit 11: Quality range
!	call nc2f (ncid, 'swh_ocean_qual_01_ku',12)			! bit 12: Quality SWH
!	call nc2f (ncid, 'sig0_ocean_qual_01_ku',13)			! bit 13: Quality Sigma0
! Redo the last ones for SAR

! Convert all the necessary fields to RADS
	call get_var (ncid, 'time_01', a)
	pair%order = (/(i, i=1, nrec)/)

! Check for time reversal
	if (any(a(2:nrec) - a(1:nrec-1) < 0)) then
		call log_string ('Error: Time reversal detected', .true.)
		forall (i = 1:nrec)
			pair(i) = quicksort_pair (i, a(i))
		end forall
		call quicksort (pair)
	endif
	call new_var ('time', a(pair%order) + sec2000)
	call cpy_var (ncid, 'lat_01', 'lat')

! Compute ellipsoid corrections
	do i = 1,nrec
		dh(i) = dhellips(1,a(pair(i)%order))
	enddo

	call get_var (ncid, 'lon_01', a)
	call new_var ('lon' , a(pair%order))
	call get_var (ncid, 'alt_01', a)
	call new_var ('alt_gdre', a(pair%order) + dh)
	call get_var (ncid, 'orb_alt_rate_01', a)
	call new_var ('alt_rate' , a(pair%order))
	call get_var (ncid, 'range_ocean_01_ku', a)
	call new_var ('range_ku' , a(pair%order))
	call get_var (ncid, 'range_ocean_01_plrm_ku', a)
	call new_var ('range_ku_plrm' , a(pair%order))
	! Add zero or meas altitude tropo measurements?
	call get_var (ncid, 'mod_dry_tropo_cor_01', a)
	call new_var ('dry_tropo_ecmwf' , a(pair%order))
	call get_var (ncid, 'mod_wet_tropo_cor_01', a)
	call new_var ('wet_tropo_ecmwf' , a(pair%order))
	call get_var (ncid, 'gpd_wet_tropo_cor_01', a)
	call new_var ('gpd_wet_tropo_cor' , a(pair%order))
	call get_var (ncid, 'gpd_wet_tropo_cor_qual_01', a)
	call new_var ('gpd_source_flag' , a(pair%order))
	call get_var (ncid, 'iono_cor_gim_01', a)
	call new_var ('iono_gim' , a(pair%order))
	call get_var (ncid, 'inv_bar_cor_01', a)
	call new_var ('inv_bar_static' , a(pair%order))
! Note that ESA Ocean Products define hf_fluct_cor_01 as including inv_bar_cor_01
	if (latency == rads_nrt) then
		call get_var (ncid, 'inv_bar_cor_01', a)
		call new_var ('inv_bar_mog2d' , a(pair%order))
	else
		call get_var (ncid, 'hf_fluct_cor_01', a)
		call new_var ('inv_bar_mog2d' , a(pair%order))
	endif
	call get_var (ncid, 'solid_earth_tide_01', a)
	call new_var ('tide_solid' , a(pair%order))
	call get_var (ncid, 'ocean_tide_sol1_01 load_tide_sol1_01 SUB', a)
	call new_var ('tide_ocean_got410' , a(pair%order))
	call get_var (ncid, 'ocean_tide_sol2_01 load_tide_sol2_01 SUB ocean_tide_non_eq_01 ADD', a)
	call new_var ('tide_ocean_fes14', a(pair%order))
	call get_var (ncid, 'load_tide_sol1_01', a)
	call new_var ('tide_load_got410' , a(pair%order))
	call get_var (ncid, 'load_tide_sol2_01', a)
	call new_var ('tide_load_fes14', a(pair%order))
	call get_var (ncid, 'ocean_tide_eq_01', a)
	call new_var ('tide_equil' , a(pair%order))
	call get_var (ncid, 'ocean_tide_non_eq_01', a)
	call new_var ('tide_non_equil' , a(pair%order))
	call get_var (ncid, 'pole_tide_01', a)
	call new_var ('tide_pole' , a(pair%order))
	call get_var (ncid, 'sea_state_bias_01_ku', a)
	call new_var ('ssb_cls' , a(pair%order))
	call get_var (ncid, 'sea_state_bias_01_plrm_ku', a)
	call new_var ('ssb_cls_plrm' , a(pair%order))
	call get_var (ncid, 'geoid_01', a)
	call new_var ('geoid_egm2008', a(pair%order) + dh)

	if (software_version < '4.00') then
		call get_var (ncid, 'mean_sea_surf_sol1_01', a)
		call new_var ('mss_cnescls15', a(pair%order) + dh)
		call get_var (ncid, 'mean_sea_surf_sol2_01', a)
		call new_var ('mss_dtu15', a(pair%order) + dh)
	else
		! We do not copy the MSS CNES-CLS22, marked as obsolete
		call get_var (ncid, 'mean_sea_surf_sol2_01', a)
		call new_var ('mss_dtu21', a(pair%order) + dh)
	endif
	call get_var (ncid, 'swh_ocean_01_ku', a)
	call new_var ('swh_ku', a(pair%order))
	call get_var (ncid, 'swh_ocean_01_plrm_ku', a)
	call new_var ('swh_ku_plrm', a(pair%order))
	call get_var (ncid, 'sig0_ocean_01_ku', a)
	call new_var ('sig0_ku', a(pair%order))
	call get_var (ncid, 'sig0_ocean_01_plrm_ku', a)
	call new_var ('sig0_ku_plrm', a(pair%order))
	call get_var (ncid, 'wind_speed_alt_01_ku', a)
	call new_var ('wind_speed_alt', a(pair%order))
	call get_var (ncid, 'wind_speed_alt_01_plrm_ku', a)
	call new_var ('wind_speed_alt_plrm', a(pair%order))
	call get_var (ncid, 'wind_speed_mod_u_01', a)
	call new_var ('wind_speed_ecmwf_u', a(pair%order))
	call get_var (ncid, 'wind_speed_mod_v_01', a)
	call new_var ('wind_speed_ecmwf_v', a(pair%order))
	! In case of SAR, there is no estimated attitude, so we rely on the PLRM values.
	! It may not be good to use these for editing (?)
	call get_var (ncid, 'off_nadir_angle_wf_ocean_01_ku off_nadir_angle_wf_ocean_01_plrm_ku AND', a)
	call new_var ('off_nadir_angle2_wf_ku', a(pair%order))
	call get_var (ncid, 'range_ocean_rms_01_ku', a)
	call new_var ('range_rms_ku', a(pair%order))
	call get_var (ncid, 'range_ocean_rms_01_plrm_ku', a)
	call new_var ('range_rms_ku_plrm', a(pair%order))
	call get_var (ncid, 'range_ocean_numval_01_ku', a)
	call new_var ('range_numval_ku', a(pair%order))
	call get_var (ncid, 'range_ocean_numval_01_plrm_ku', a)
	call new_var ('range_numval_ku_plrm', a(pair%order))
	call get_var (ncid, 'peakiness_01_ku', a)
	call new_var ('peakiness_ku' , a(pair%order)/10.0)
	call get_var (ncid, 'odle_01', a)
	call new_var ('topo_ace2' , a(pair%order))
	call new_var ('flags', dble(flags(pair%order)))
	call new_var ('flags_plrm', dble(flags_plrm(pair%order)))
	call get_var (ncid, 'swh_ocean_rms_01_ku', a)
	call new_var ('swh_rms_ku', a(pair%order))
	call get_var (ncid, 'swh_ocean_rms_01_plrm_ku', a)
	call new_var ('swh_rms_ku_plrm', a(pair%order))
	call get_var (ncid, 'sig0_ocean_rms_01_ku', a)
	call new_var ('sig0_rms_ku', a(pair%order))
	call get_var (ncid, 'sig0_ocean_rms_01_plrm_ku', a)
	call new_var ('sig0_rms_ku_plrm', a(pair%order))
	call get_var (ncid, 'ssha_01_ku', a)
	call new_var ('ssha', a(pair%order))
	call get_var (ncid, 'ssha_01_plrm_ku', a)
	call new_var ('ssha_plrm', a(pair%order))
	call get_var (ncid, 'surf_type_01', a)
	call new_var ('surface_class', a(pair%order))
	call get_var (ncid, 'flag_instr_op_mode_01', a)
	call new_var ('flag_alt_oper_mode', a(pair%order) - 1)
	a = latency
	call new_var ('latency', a)

	! Dump the data

	call nfs(nf90_close(ncid))
	call put_rads
	deallocate (a, dh, flags, flags_plrm, flags_save, pair)

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
if (rads_version ('Write CryoSat Ocean Products data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_CryoSat_OP_file_names')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --min-rec=MIN_REC         Specify minimum number of records per pass to process' // &
'This program converts CryoSat NOP/IOP/GOP files to RADS data' / &
'files with the name $RADSDATAROOT/data/c2.[nig]op/a/pPPPP/s3pPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_c2_op
