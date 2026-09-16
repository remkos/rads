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

!*rads_gen_reaper -- Converts REAPER ERS-1/2 data to RADS
!+
program rads_gen_reaper

! This program reads REAPER files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/eE.VVVV/F/eEpPPPPcCCC.nc.
!     E = 1 or 2
!  VVVV = version ID
!     F = mission phase
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_reaper [options] < list_of_REAPER_file_names
!
! where [options] include:
!  --min-rec <min_rec> : Specify minimum number of records per pass to process.
!
! This program handles only the REAPER ERS_ALT_2 files in NetCDF format.
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_reaper - Orbit altitude
! alt_rate - Orbit altitude rate
! range_ku - Ocean range (retracked)
! dry_tropo_era - ERA dry tropospheric correction
! wet_tropo_rad - Radiometer wet tropo correction
! wet_tropo_era - ERA wet tropo correction
! iono_gim - GIM ionosphetic correction
! iono_nic09 - NIC09 ionospheric correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! tide_solid - Solid earth tide
! tide_pole - Pole tide
! ssb_hyb - Hybrid SSB solution for REAPER
! geoid_egm2008 - EGM2008 geoid
! swh_ku - Significant wave height
! sig0_ku - Sigma0
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! range_rms_ku - Std dev of range
! range_numval_ku - Nr of averaged range measurements
! tb_238 - Brightness temperature (23.8 GHz)
! tb_365 - Brightness temperature (36.5 GHz)
! peakiness_ku - Peakiness
! flags - Engineering flags
! swh_rms_ku - Std dev of SWH
! sig0_rms_ku - Std dev of sigma0
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! rad_liquid_water - Liquid water content
! rad_water_vapor - Water vapor content
! tide_equil - Long-period equilibrium tide
! tide_non_equil - Long-period non-equilibrium tide
!
! On (S)GDR only:
! wind_speed_alt - Altimeter wind speed
! drange_cal - Internal calibration correction to range (appied)
! drange_fm - Doppler correction (applied)
! dsig0_atmos_ku - Sigma0 attenuation
! mqe - Mean quadratic error of waveform fit
!-----------------------------------------------------------------------
use rads
use rads_devel
use rads_devel_misc
use rads_gen
use rads_netcdf
use rads_misc
use rads_time
use rads_geo
use netcdf

! Command line arguments

integer(fourbyteint) :: ios
character(len=rads_cmdl) :: infile, filenm, old_filenm = ''

! Header variables

character(len=2) :: mission
character(len=rads_varl) :: l2_proc_time, l2_version
character(len=4) :: mle
logical :: alt_2m
integer(fourbyteint) :: varid

! Data variables

integer(fourbyteint), parameter :: mrec=20000, mvar=50
integer(fourbyteint) :: nvar, nrec_buf=0, nrec_in=0, nrec_out=0, ncid, ers=0
real(eightbytereal) :: start_time, end_time, last_time = 0d0
real(eightbytereal), allocatable :: tmp(:)
type(rads_sat) :: S
type(rads_pass) :: P
type :: var_
	type(rads_var), pointer :: v ! Pointer to rads_var struct
	real(eightbytereal) :: d(mrec) ! Data array
	logical :: empty, zero ! .true. if all NaN or all zero
endtype
type(var_) :: var(mvar)

! Struct for orbit info

integer(fourbyteint) :: ipass
integer(fourbyteint), parameter :: mpass = 170000 ! Enough for 17 years
type(orfinfo) :: orf(mpass)

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000
real(eightbytereal), parameter :: sec1990=157766400d0	! UTC seconds from 1 Jan 1985 to 1 Jan 1990
real(eightbytereal), parameter :: sec_to_m=0.5d0*299792458d0	! seconds of 2-way range to meters 1-way
real(eightbytereal), parameter :: uso_wap=15000000.05d0	! Nominal USO frequency used in WAP products
integer :: i

! Initialise

call synopsis
call rads_gen_getopt ('', ' min-rec:')

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

! Start reading at least first file

read (*,'(a)',iostat=ios) infile
if (ios /= 0) then
	call synopsis ('--help')
else
	call synopsis ('--head')
endif
call rads_init(S, sat)
call get_reaper

! Load the ORF file and change time to sec1985

call read_orf ('ER' // mission(2:2), orf)
orf(:)%starttime = orf(:)%starttime + sec2000
orf(:)%eqtime = orf(:)%eqtime + sec2000
ipass = 1

do
	! Read the next file as long as buffer is empty or less than one orbit in memory

	do while (nrec_buf == 0 .or. var(1)%d(nrec_buf) - var(1)%d(1) < 6100d0)
		read (*,'(a)',iostat=ios) infile
		if (ios /= 0) exit
		call get_reaper
	enddo

	! To which passes do the data belong?
	call which_pass (var(1)%d(1))

	! Look where to split this chunk of data
	do i = 2,nrec_buf
		if (var(1)%d(i) >= orf(ipass+1)%starttime) exit
	enddo
	! It is OK to exit this loop with i = nrec_buf + 1. This means we dump all of the memory.

	! Write out the data
	nrec_out = i - 1 ! Number of measurements to be written out
	call put_rads

	! Number of measurements remaining
	nrec_buf = nrec_buf - nrec_out
	if (ios /= 0 .and. nrec_buf == 0) exit ! We are out of data

	! Move the data to be beginning
	forall (i = 1:nvar)
		var(i)%d(1:nrec_buf) = var(i)%d(nrec_out+1:nrec_out+nrec_buf)
	end forall
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
if (rads_version ('Write REAPER data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_REAPER_file_names')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --min-rec=MIN_REC         Specify minimum number of records per pass to process' // &
'This program converts REAPER ERS_ALT_2 files to RADS data' / &
'files with the name $RADSDATAROOT/data/eE.VVVV/F/pPPPP/eEpPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Store the content of a REAPER file into memory
!-----------------------------------------------------------------------

subroutine get_reaper
integer(fourbyteint) :: i

552 format (i4,' records ...')

! Reduce file name to basename only

i = index(infile,'/',.true.) + 1
old_filenm = filenm
filenm = infile(i:)

! Check input file name

call log_string (basename(infile))
if (filenm(8:17) /= '_ERS_ALT_2') then
	call log_string ('Error: wrong input file type', .true.)
	return
endif
alt_2m = (filenm(18:18) == 'M')

mission = filenm(1:2)
if (mission == 'E1') then
	ers = 1
else if (mission == 'E2') then
	ers = 2
else
	call log_string ('Error: Wrong file type: '//mission, .true.)
	return
endif

! Set start and end time from file name
! Add a little bit of slop at the end because sometimes the times overrun the end time

start_time = strp1985f (filenm(20:34))
end_time   = strp1985f (filenm(36:50)) + 2d0

! Open input file

i = nf90_open(infile,nf90_nowrite,ncid)
if (i /= nf90_noerr) then
	call log_string (nf90_strerror(i), .true.)
	stop
	return
endif

! Read header records

call nfs(nf90_inq_dimid(ncid,'time',varid))
call nfs(nf90_inquire_dimension(ncid,varid,len=nrec_in))
write (rads_log_unit,552) nrec_in
if (nrec_in == 0) then	! Skip empty input files
	call nfs(nf90_close(ncid))
	return
else if (nrec_buf+nrec_in > mrec) then
	call log_string ('Error: too many input measurements', .true.)
	stop
endif
call nfs(nf90_get_att(ncid,nf90_global,'l2_proc_time',l2_proc_time))
call nfs(nf90_get_att(ncid,nf90_global,'l2_software_ver',l2_version))
call nfs(nf90_get_att(ncid,nf90_global,'ocean_retracker_version_for_ocean',mle))
!call nfs(nf90_get_att(ncid,nf90_global,'mission',mission)) ! Actually still is not correct -> bypass

call get_reaper_data (nrec_in)
call nfs(nf90_close(ncid))
end subroutine get_reaper

subroutine get_reaper_data (nrec)
integer(fourbyteint), intent(inout) :: nrec
real(eightbytereal) :: a(nrec), b(nrec), c(nrec), d(20,nrec), dh(nrec)
integer(twobyteint) :: flags(nrec)
logical :: valid(20,nrec)
integer(fourbyteint) :: k, kerr(4)
real(eightbytereal) :: t(3)
character(len=34) :: strerr(4) = (/ &
'measurements out of time range   ', &
'measurements out of time sequence', &
'measurements with time overlap   ', &
'measurements with time reversal  ' /)

553 format ('Warning: Removed ',a,':',i6)

! Time and orbit: Low rate

call get_var (ncid, 'time', a)
a = a + sec1990

! Because first and last time record can be wrong, we use the ones from the file name, which
! seem to be reliable
t(1) = last_time
t(2) = start_time
t(3) = end_time

! There are significant overlaps between files
! Here we remove all new files that fall entirely before the end of the previous file, and
if (end_time < last_time + 1) then
	write (rads_log_unit,553) 'file because of time reversal     ', nrec
	return
endif

! Initialize

flags = 0
nvar = 0
allocate (tmp(nrec))

! Time and orbit: Low rate (cont'd)

call new_var ('time', a)
call cpy_var (ncid, 'lat', 'lat')
! Compute ellipsoid corrections
do i = 1,nrec
	dh(i) = dhellips(1,tmp(i))
enddo
call cpy_var (ncid, 'lon', 'lon')
call get_var (ncid, 'alt', a)
call new_var ('alt_reaper', a+dh)
call cpy_var (ncid, 'orb_alt_rate', 'alt_rate')

! Range data: Low rate

call get_var (ncid, 'ocean_range_numval', c)
call cpy_var (ncid, 'ocean_range', 'range_ku')
call cpy_var (ncid, 'ocean_range_rms', 'range_rms_ku', c <= 1)
call cpy_var (ncid, 'ocean_range_numval', 'range_numval_ku')

! Retracking info: Low rate

if (mle == 'MLE4') call cpy_var (ncid, 'off_nadir_angle_wf', 'off_nadir_angle2_wf_ku')
call cpy_var (ncid, 'peakiness', 'peakiness_ku')

! Retracking info: High rate

if (.not.alt_2m) then
	call get_var (ncid, 'ocean_qual_flag_20hz', d)
	valid = (d == 0)
	call get_var (ncid, 'ocean_mqe_20hz', d)
	where (.not.valid) d = nan
	call mean_1hz (d, a, b)
	call new_var ('mqe', a)
	call get_var (ncid, 'doppler_corr delta_doppler_corr_20hz ADD', d)
	call mean_1hz (d, a, b)
	call new_var ('drange_fm', a)
	call get_var (ncid, 'sptr_jumps_corr_20hz', d)
	call new_var ('drange_sptr', d(1,:) * sec_to_m)
endif

! Sigma zero: Low rate

call cpy_var (ncid, 'ocean_sig0', 'sig0_ku')
call get_var (ncid, 'ocean_sig0_numval', c)
call cpy_var (ncid, 'ocean_sig0_rms', 'sig0_rms_ku', c <= 1)
call cpy_var (ncid, 'wind_speed_alt', 'wind_speed_alt', c == 0)

! SWH: Low rate

call cpy_var (ncid, 'swh_signed', 'swh_ku')
call get_var (ncid, 'swh_numval', c)
call cpy_var (ncid, 'swh_rms', 'swh_rms_ku', c <= 1)

! MWR: Low rate

call get_var (ncid, 'rad_state_flag_orb_prop', a)
call get_var (ncid, 'rad_state_flag_orb_init', b)
call get_var (ncid, 'rad_state_flag_l2_proc_error rad_state_validity ADD rad_state_bt_check ADD', c)
call cpy_var (ncid, 'tb_238', 'tb_238', a == 2d0 .or. b == 2d0 .or. c > 0d0)
call cpy_var (ncid, 'tb_365', 'tb_365', a == 2d0 .or. b == 2d0 .or. c > 0d0)

! Atmospheric and geophysical: Low rate

call cpy_var (ncid, 'model_dry_tropo_corr', 'dry_tropo_era')
call cpy_var (ncid, 'inv_bar_corr', 'inv_bar_static')
call cpy_var (ncid, 'hf_fluctuations_corr', 'inv_bar_mog2d')	! This is in contrast with standard_name, but supported by comment
call cpy_var (ncid, 'model_wet_tropo_corr', 'wet_tropo_era')
call cpy_var (ncid, 'rad_wet_tropo_corr', 'wet_tropo_rad')
call cpy_var (ncid, 'rad_water_vapor', 'water_vapor_rad')
call cpy_var (ncid, 'rad_liquid_water', 'liquid_water_rad')
call cpy_var (ncid, 'wind_speed_model_u', 'wind_speed_era_u')
call cpy_var (ncid, 'wind_speed_model_v', 'wind_speed_era_v')
call cpy_var (ncid, 'iono_corr_model', 'iono_nic09')
if (start_time >= 430880400d0) then	! After 1998-08-28 01:00:00 get GIM iono
	call cpy_var (ncid, 'iono_corr_gps', 'iono_gim')
endif

! MSS on REAPER is MSS UCL04, which is now obsolete
call get_var (ncid, 'geoid', a)
call new_var ('geoid_egm2008', a+dh)

! Ocean tide on REAPER is FES2004, which is now obsolete
call cpy_var (ncid, 'solid_earth_tide', 'tide_solid')
call cpy_var (ncid, 'pole_tide', 'tide_pole')

call cpy_var (ncid, 'sea_state_bias', 'ssb_hyb')

! Atmospheric correction is a factor 100 too small (wrong scale_factor)
call get_var (ncid, 'atmos_corr_sig0', a)
call new_var ('dsig0_atmos_ku', a*1d2)

! Status flags

call get_var (ncid, 'alt_state_flag', a)
call flag_set (a <= 1, flags, 11) ! Invalid range etc.
call flag_set (a <= 1, flags, 12)
call flag_set (a <= 1, flags, 13)
call flag_set (a == 3, flags, 14) ! Ice mode

call get_var (ncid, 'surface_type', a)
call flag_set (a == 2, flags, 2)	! Bit 2: Continental ice
call flag_set (a >= 2, flags, 4)	! Bit 4: Water/land
call flag_set (a >= 1, flags, 5)	! Bit 5: Ocean/other

call get_var (ncid, 'ice_flag', a)
call flag_set (a == 1, flags, 8)
call get_var (ncid, 'rad_surf_type', a)
call flag_set (a == 1, flags, 6)
call get_var (ncid, 'rad_state_bt_check', a)
call flag_set (a >= 2, flags, 9)
call flag_set (modulo(a,2d0) == 1, flags, 10)

call new_var ('flags', flags*1d0)

deallocate (tmp)

! There may be measurements with invalid times.
! If so, weed them out.

k = 0
kerr = 0
valid(1,:) = .false.
do i = 1,nrec
	if (i > 1 .and. i < nrec) then
		t = var(1)%d(nrec_buf+i-1:nrec_buf+i+1)
	else
		t(2) = var(1)%d(nrec_buf+i)
	endif
	if (t(2) < start_time .or. t(2) > end_time) then
		kerr(1) = kerr(1) + 1
	else if (i > 1 .and. i < nrec .and. t(2) > max(t(1),t(3))+1d0) then
		kerr(2) = kerr(2) + 1
	else if (t(2) > last_time+0.5d0) then
		last_time = t(2)
		valid(1,i) = .true.
		k = k + 1
	else if (k == 0) then
		kerr(3) = kerr(3) + 1
	else
		kerr(4) = kerr(4) + 1
	endif
enddo
do i = 1,4
	if (kerr(i) > 0) write (rads_log_unit,553) strerr(i),kerr(i)
enddo
if (k == 0) then
	nrec = 0
else if (k < nrec) then
	do i = 1,nvar
		var(i)%d(nrec_buf+1:nrec_buf+k) = pack(var(i)%d(nrec_buf+1:nrec_buf+nrec),valid(1,:))
	enddo
	nrec = k
endif

nrec_buf = nrec_buf + nrec

end subroutine get_reaper_data

!-----------------------------------------------------------------------
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads
integer :: i
character(len=rads_cmdl) :: original

if (nrec_out < min_rec) return	! Skip empty data sets
if (orf(ipass)%cycle < cycles(1) .or. orf(ipass)%cycle > cycles(2)) return	! Skip chunks that are not of the selected cycle
if (orf(ipass)%eqtime < times(1) .or. orf(ipass)%eqtime > times(2)) return	! Skip equator times that are not of selected range

! Set mission phase based on equator_time

call rads_set_phase (S, orf(ipass)%eqtime)

! Store relevant info
call rads_init_pass_struct (S, P)
P%cycle = orf(ipass)%cycle
P%pass = orf(ipass)%pass
P%start_time = var(1)%d(1)
P%end_time = var(1)%d(nrec_out)
P%equator_time = orf(ipass)%eqtime
P%equator_lon = orf(ipass)%eqlon

! Check which input files pertain
if (P%start_time >= start_time) then
	original = filenm
else if (P%end_time < start_time) then
	original = old_filenm
else
	original = trim(old_filenm)//' '//filenm
endif
P%original = trim(l2_version)//' data of '//l2_proc_time(:11)//': '//trim(original)

! Check which variables are empty or all zero
do i = 1,nvar
	var(i)%empty = all(isnan_(var(i)%d(1:nrec_out)))
	var(i)%zero = all(var(i)%d(1:nrec_out) == 0d0)
enddo

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

! Open output file
call rads_create_pass (S, P, nrec_out)

! Define all variables
do i = 1,nvar
	call rads_def_var (S, P, var(i)%v)
enddo

! Fill all the data fields
do i = 1,nvar
	call rads_put_var (S, P, var(i)%v, var(i)%d(1:nrec_out))
enddo

! Close the data file
call log_records (nrec_out, P)
call rads_close_pass (S, P)

end subroutine put_rads

!-----------------------------------------------------------------------
! Copy variable to RADS
!-----------------------------------------------------------------------

subroutine cpy_var (ncid, varin, varout, invalid)
integer(fourbyteint), intent(in) :: ncid
character(len=*), intent(in) :: varin, varout
logical, optional :: invalid(:)
call get_var (ncid, varin, tmp)
if (present(invalid)) where (invalid) tmp = nan
call new_var (varout, tmp)
end subroutine cpy_var

!-----------------------------------------------------------------------
! Create new RADS variable
!-----------------------------------------------------------------------

subroutine new_var (varnm, data, ndims)
! Store variables to be written later by put_var
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(in) :: data(:)
integer, optional, intent(in) :: ndims
nvar = nvar + 1
if (nvar > mvar) stop 'Too many variables'
var(nvar)%v => rads_varptr (S, varnm)
var(nvar)%d(nrec_buf+1:nrec_buf+nrec_in) = data
if (present(ndims)) var(nvar)%v%info%ndims = ndims
end subroutine new_var

!-----------------------------------------------------------------------
! Set a bit in an array of flags
!-----------------------------------------------------------------------

subroutine flag_set (a, flags, bit)
logical, intent(in) :: a(:)
integer(twobyteint), intent(inout) :: flags(:)
integer(fourbyteint), intent(in) :: bit
integer(fourbyteint) :: i
integer(twobyteint) :: j
if (size(a) /= size(flags)) stop 'Error in flag_set'
j = int(bit,twobyteint)
do i = 1,size(a)
	if (a(i)) flags(i) = ibset(flags(i),j)
enddo
end subroutine flag_set

end program rads_gen_reaper
