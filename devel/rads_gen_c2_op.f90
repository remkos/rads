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
! alt_{gdre,gdrf} - Orbital altitude
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
! tide_{ocean,load}_got410 - GOT4.10c ocean and load tide
! tide_{ocean,load}_{fes14,fes22} - FES2014b or FES2022b ocean and load tide
! tide_pole - Pole tide
! tide_solid - Solid earth tide
! topo_ace2 - ACE2 topography
! geoid_egm2008 - EGM2008 geoid
! mss_{cnescls15,cnescls22} - CNES/CLS15 or CNES/CLS22 mean sea surface
! mss_{dtu15,dtu21} - DTU15 or DTU21 mean sea surface
! surface_slope_cor - surface slope correction
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

real(eightbytereal), parameter :: sec2000 = 473299200d0 ! UTC seconds from 1 Jan 1985 to 1 Jan 2000
real(eightbytereal), parameter :: sec20190417 = 1082073600d0 ! UTC seconds from 1 Jan 1985 to 17 Apr 2019
real(eightbytereal), allocatable :: dh(:)

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

	allocate (a(nrec),dh(nrec),flags(nrec),flags_plrm(nrec),flags_save(nrec))
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
	call new_var ('time', a + sec2000)

	call cpy_var (ncid, 'lat_01', 'lat')
! Compute ellipsoid corrections
	do i = 1,nrec
		dh(i) = dhellips(1,a(i))
	enddo
	call cpy_var (ncid, 'lon_01', 'lon')

	call get_var (ncid, 'alt_01', a)
	! POE-E until 2019-04-27, POE-F thereafter
	if (equator_time < sec20190417) then
		call new_var ('alt_gdre', a + dh)
	else
		call new_var ('alt_gdrf', a + dh)
	endif
	call cpy_var (ncid, 'orb_alt_rate_01', 'alt_rate')

	call cpy_var (ncid, 'range_ocean_01_ku', 'range_ku')
	call cpy_var (ncid, 'range_ocean_01_plrm_ku', 'range_ku_plrm')
	call cpy_var (ncid, 'range_ocean_rms_01_ku', 'range_rms_ku')
	call cpy_var (ncid, 'range_ocean_rms_01_plrm_ku', 'range_rms_ku_plrm')
	call cpy_var (ncid, 'range_ocean_numval_01_ku', 'range_numval_ku')
	call cpy_var (ncid, 'range_ocean_numval_01_plrm_ku', 'range_numval_ku_plrm')
	call cpy_var (ncid, 'dop_cor_01_ku', 'drange_fm')

	! Add zero or meas altitude tropo measurements?
	call cpy_var (ncid, 'mod_dry_tropo_cor_01', 'dry_tropo_ecmwf')
	call cpy_var (ncid, 'mod_wet_tropo_cor_01', 'wet_tropo_ecmwf')
	call cpy_var (ncid, 'gpd_wet_tropo_cor_01', 'gpd_wet_tropo_cor')
	call cpy_var (ncid, 'gpd_wet_tropo_cor_qual_01', 'gpd_source_flag')
	call cpy_var (ncid, 'iono_cor_gim_01', 'iono_gim')
	call cpy_var (ncid, 'inv_bar_cor_01', 'inv_bar_static')
! Note that ESA Ocean Products define hf_fluct_cor_01 as including inv_bar_cor_01
	if (latency == rads_nrt) then
		call cpy_var (ncid, 'inv_bar_cor_01', 'inv_bar_mog2d')
	else
		call cpy_var (ncid, 'hf_fluct_cor_01', 'inv_bar_mog2d')
	endif

	call cpy_var (ncid, 'solid_earth_tide_01', 'tide_solid')
	call cpy_var (ncid, 'ocean_tide_eq_01', 'tide_equil')
	call cpy_var (ncid, 'ocean_tide_non_eq_01', 'tide_non_equil' )
	call cpy_var (ncid, 'ocean_tide_sol1_01 load_tide_sol1_01 SUB', 'tide_ocean_got410')
	call cpy_var (ncid, 'ocean_tide_sol2_01 load_tide_sol2_01 SUB ocean_tide_non_eq_01 ADD', 'tide_ocean_fes14')
	call cpy_var (ncid, 'load_tide_sol1_01', 'tide_load_got410')
	call cpy_var (ncid, 'load_tide_sol2_01', 'tide_load_fes14')
	call cpy_var (ncid, 'ocean_tide_eq_01', 'tide_equil')
	call cpy_var (ncid, 'ocean_tide_non_eq_01', 'tide_non_equil')
	call cpy_var (ncid, 'pole_tide_01', 'tide_pole')
	if (software_version >= '4.00') call cpy_var (ncid, 'internal_tide_01', 'tide_internal')

	call cpy_var (ncid, 'sea_state_bias_01_ku', 'ssb_cls')
	call cpy_var (ncid, 'sea_state_bias_01_plrm_ku', 'ssb_cls_plrm')

	call get_var (ncid, 'geoid_01', a)
	call new_var ('geoid_egm2008', a + dh)
	call cpy_var (ncid, 'mean_dyn_topo_sol1_01', 'mdt_cnescls22')
	if (software_version < '4.00') then
		call get_var (ncid, 'mean_sea_surf_sol1_01', a)
		call new_var ('mss_cnescls15', a + dh)
		call get_var (ncid, 'mean_sea_surf_sol2_01', a)
		call new_var ('mss_dtu15', a + dh)
	else
		! We do not copy the MSS CNES-CLS22, marked as obsolete
		call get_var (ncid, 'mean_sea_surf_sol2_01', a)
		call new_var ('mss_dtu21', a + dh)
	endif

	call cpy_var (ncid, 'swh_ocean_01_ku', 'swh_ku')
	call cpy_var (ncid, 'swh_ocean_01_plrm_ku', 'swh_ku_plrm')
	call cpy_var (ncid, 'swh_ocean_rms_01_ku', 'swh_rms_ku')
	call cpy_var (ncid, 'swh_ocean_rms_01_plrm_ku', 'swh_rms_ku_plrm')

	call cpy_var (ncid, 'sig0_ocean_01_ku', 'sig0_ku')
	call cpy_var (ncid, 'sig0_ocean_01_plrm_ku', 'sig0_ku_plrm')
	call cpy_var (ncid, 'sig0_ocean_rms_01_ku', 'sig0_rms_ku')
	call cpy_var (ncid, 'sig0_ocean_rms_01_plrm_ku', 'sig0_rms_ku_plrm')
	call cpy_var (ncid, 'atm_cor_sig0_01', 'dsig0_atmos_ku')

	call cpy_var (ncid, 'wind_speed_alt_01_ku', 'wind_speed_alt')
	call cpy_var (ncid, 'wind_speed_alt_01_plrm_ku', 'wind_speed_alt_plrm')
	call cpy_var (ncid, 'wind_speed_mod_u_01', 'wind_speed_ecmwf_u')
	call cpy_var (ncid, 'wind_speed_mod_v_01', 'wind_speed_ecmwf_v')

	! In case of SAR, there is no estimated attitude, so we rely on the PLRM values.
	! It may not be good to use these for editing (?)
	call cpy_var (ncid, 'off_nadir_angle_wf_ocean_01_ku off_nadir_angle_wf_ocean_01_plrm_ku AND', 'off_nadir_angle2_wf_ku')
	call cpy_var (ncid, 'peakiness_01_ku 10 DIV', 'peakiness_ku')

	call cpy_var (ncid, 'odle_01', 'topo_ace2')
	call new_var ('flags', dble(flags))
	call new_var ('flags_plrm', dble(flags_plrm))

	call cpy_var (ncid, 'ssha_01_ku', 'ssha')
	call cpy_var (ncid, 'ssha_01_plrm_ku', 'ssha_plrm')
	call cpy_var (ncid, 'surf_type_01', 'surface_class')
	if (software_version >= '4.00') call cpy_var (ncid, 'surface_slope_cor_01', 'surface_slope_cor')
	call cpy_var (ncid, 'flag_instr_op_mode_01 1 SUB', 'flag_alt_oper_mode')
	a = latency
	call new_var ('latency', a)

	! Dump the data

	call put_rads_with_index (0.5d0 * S%dt1hz, 'range_ocean_numval_01_ku')
	call nfs(nf90_close(ncid))
	deallocate (a, dh, flags, flags_plrm, flags_save)

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads_with_index (time_slop, refvar)
real(eightbytereal), intent(in) :: time_slop
character(len=*), intent(in) :: refvar
real(eightbytereal), allocatable :: time(:), val(:)
integer(fourbyteint), allocatable :: idx(:)
integer :: i, j, ndup

! Work out time consistency since some measurements are duplicated in the GOP products.
allocate (time(nrec), val(nrec))
call get_var (ncid, 'time_01', time)
call get_var (ncid, refvar, val)

! Check for time duplicated (or time reversal)
ndup = count (time(2:nrec) - time(1:nrec-1) < time_slop)
if (ndup > 0) write (rads_log_unit,550,advance='no') ndup
550	format (i3,' duplicate times detected ... ')

allocate (idx(nrec-ndup))
j = 1
idx(1) = 1
do i = 2,nrec
	if (time(i) - time(idx(j)) < time_slop) then
		! If we have a duplicate, keep the one of which the number of elementary
		! range measurements is the largest. If the same, keep the latest.
		if (val(i) >= val(idx(j))) idx(j) = i
	else
		j = j + 1
		idx(j) = i
	endif
enddo
if (j /= nrec - ndup) call log_string ('Error: index length does not match')

! Check which variables are empty or all zero
forall (i = 1:nvar)
	var(i)%empty = all(isnan_(var(i)%d(idx)))
	var(i)%zero = all(var(i)%d(idx) == 0d0)
end forall

! Write out the empty variables to be kept
if (any(var(1:nvar)%empty)) then
	write (rads_log_unit,551,advance='no') 'Empty:'
	do i = 1,nvar
		if (var(i)%empty) write (rads_log_unit,551,advance='no') trim(var(i)%v%name)
	enddo
	write (rads_log_unit,551,advance='no') ' ...'
endif
551 format (a,1x)

! Do the same for records that are all zero
if (any(var(1:nvar)%zero)) then
	write (rads_log_unit,551,advance='no') 'All zero:'
	do i = 1,nvar
		if (var(i)%zero) write (rads_log_unit,551,advance='no') trim(var(i)%v%name)
	enddo
	write (rads_log_unit,551,advance='no') '...'
endif

! Create the output file
call rads_create_pass (S, P, nrec-ndup)

! Define all variables
do i = 1,nvar
	call rads_def_var (S, P, var(i)%v)
enddo

! Fill all the data fields
do i = 1,nvar
	call rads_put_var (S, P, var(i)%v, var(i)%d(idx))
enddo
deallocate (time, val, idx)

! Close the data file
call log_records (nrec-ndup, P)
call rads_close_pass (S, P)

end subroutine put_rads_with_index

!-----------------------------------------------------------------------
! Determine the corresponding record from ORF
!-----------------------------------------------------------------------

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
