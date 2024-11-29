!-----------------------------------------------------------------------
! Copyright (c) 2011-2024  Remko Scharroo
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

!*rads_add_f4a_ers -- Converts ERS FDR4ALT altimeter data to RADS
!+
program rads_add_f4a_ers

! This program reads ERS FDR4ALT altimeters files and adds them to existing
! RADS data files
!
! syntax: rads_add_f4a [options] < list_of_FDR4ALT_file_names
!
! This program handles FDR4ALT Level 2 products in NetCDF format, either
! OC (Ocean and Coastal) or WA (Waves) products.
! Their format and content is described in:
!
! [1] FDR4ALT Product User Guide
!     CLS-ENV-MU-23-0237, issue 2.2, 30 Oct 2023
!
!-----------------------------------------------------------------------
!
! Variables array fields to be added from the OC products are:
! range_ku - Ocean Ku-band range (retracked)
! range_rms_ku - Std dev of Ku-band range
! range_numval_ku - Nr of averaged range measurements
! dry_tropo_era5 - ERA5 dry tropospheric correction
! wet_tropo_rad - Radiometer wet tropo correction
! iono_nic09 - NIC09 ionospheric correction
! iono_gim - GIM ionosphetic correction
! ssb_tran2019_3d - SSB Tran 2019 3D
! inv_bar_mog2d_era - MOG2D dynamic atmospheric correction (ERA driven)
! tide_internal - Internal tide
! tide_solid - Solid earth tide
! tide_pole - Pole tide
! mdt_cnescls18 - Mean dynamic topography 2018
! surface_type - Surface type flag
! dist_coast - Distance to coast
! ssha - Sea surface height anomaly
!
! Variables array fields to be added from the WA products are:
! swh_ku - Ku-band significant wave height
! swh_rms_* - Std dev of Ku-band SWH
!-----------------------------------------------------------------------
use rads
use rads_devel
use rads_devel_misc
use rads_gen
use rads_misc
use rads_netcdf
use rads_time
use rads_geo
use netcdf

! Command line arguments

integer(fourbyteint) :: ios, cycle_number, pass_number
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: ncid, ncid_m1, ncid_x1, ncid_m20, ncid_x20, varid, &
	nrec_m1, nrec_x1, nrec_m5, nrec_m20, nrec_x20
real(eightbytereal) :: first_measurement_time, last_measurement_time
logical :: oc_product

! Data variables

integer(fourbyteint), parameter :: mrec=20000, mvar=50
integer(fourbyteint) :: nvar=0, ndata=0, nout=0
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
integer(fourbyteint), parameter :: mpass = 120000 ! Enough for 12 years
type(orfinfo) :: orf(mpass)

! Other local variables

integer(fourbyteint) :: i
real(eightbytereal), parameter :: sec1990=157766400d0	! UTC seconds from 1 Jan 1985 to 1 Jan 1990
real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000

! Initialise

call synopsis ('--head')
call rads_init (S)

! Load the ORF file and change time to sec1985

call read_orf ('ER' // S%sat(2:2), orf)
orf(:)%starttime = orf(:)%starttime + sec2000
orf(:)%eqtime = orf(:)%eqtime + sec2000
ipass = 1

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

do
	read (*,'(a)',iostat=ios) infile
	if (ios /= 0) exit
	ios = nf90_close(ncid)

! Open input file

	call log_string (basename(infile), .false.)
	if (nft(nf90_open(infile,nf90_nowrite,ncid))) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Check if input is an FDR4ALT altimeter Level 2 data set

	if (nft(nf90_get_att(ncid,nf90_global,'title',arg)) .or. &
		index (arg, 'FDR4ALT Thematic Data Product') == 0) then
		call log_string ('Error: this is not an FDR4ALT altimeter Level 2 data set', .true.)
		cycle
	endif
	oc_product = index(arg, 'Ocean and Coastal') > 0

! Get the mission name and compare with opened RADS session

	call nfs(nf90_get_att(ncid,nf90_global,'mission_name',arg))
	if ((arg == 'ERS-1' .and. S%sat == 'e1') .or. &
		(arg == 'ERS-2' .and. S%sat == 'e2')) then
		! All good options
	else
		call log_string ('Error: input file does not match selected satellite', .true.)
		cycle
	endif

! Which product is on input?

	if (oc_product) then
!-----------------------------------------------------------------------
! OCEAN AND COASTAL PRODUCT (TDP_OC)
!-----------------------------------------------------------------------
! Get NetCDF ID for 1-Hz and 20-Hz data (main and expert) in case of ocean/coastal product

		call nfs(nf90_inq_grp_full_ncid(ncid, '/main/data_01', ncid_m1))
		call nfs(nf90_inq_grp_full_ncid(ncid, '/expert/data_01', ncid_x1))
		call nfs(nf90_inq_grp_full_ncid(ncid, '/main/data_20', ncid_m20))
		call nfs(nf90_inq_grp_full_ncid(ncid, '/expert/data_20', ncid_x20))

! Note that the time variables in expert groups of the ERS-1 and ERS-2 products are incorrecly named

		call nfs(nf90_inq_dimid(ncid_m1, 'time', varid))
		call nfs(nf90_inquire_dimension(ncid_m1, varid, len=nrec_m1))
		call nfs(nf90_inq_dimid(ncid_m20, 'time', varid))
		call nfs(nf90_inquire_dimension(ncid_m20, varid, len=nrec_m20))
		call nfs(nf90_inq_dimid(ncid_x1, 'time_01', varid))
		call nfs(nf90_inquire_dimension(ncid_x1, varid, len=nrec_x1))
		call nfs(nf90_inq_dimid(ncid_x20, 'time_20', varid))
		call nfs(nf90_inquire_dimension(ncid_x20, varid, len=nrec_x20))
		if (nrec_m1 == 0) then
			call log_string ('Error: file skipped: no measurements', .true.)
			cycle
		else if (nrec_m1 /= nrec_x1) then
			call log_string ('Error: file skipped: nr of 1-Hz main and expert records not the same', .true.)
			cycle
		else if (nrec_m20 /= nrec_x20) then
			call log_string ('Error: file skipped: nr of 20-Hz main and expert records not the same', .true.)
			cycle
		endif
		! Read cycle and pass number
		call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cycle_number))
		call nfs(nf90_get_att(ncid,nf90_global,'pass_number',pass_number))

	else
!-----------------------------------------------------------------------
! WAVES PRODUCT (TDP_WA)
!-----------------------------------------------------------------------
! Less to check for wave product

		call nfs(nf90_inq_dimid(ncid, 'time', varid))
		call nfs(nf90_inquire_dimension(ncid, varid, len=nrec_m5))
		if (nrec_m5 == 0) then
			call log_string ('Error: file skipped: no measurements', .true.)
			cycle
		endif
		! In WA product, cycle_number and pass_number are strings
		call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',arg))
		read (arg, *, iostat=ios) cycle_number
		call nfs(nf90_get_att(ncid,nf90_global,'pass_number',arg))
		read (arg, *, iostat=ios) pass_number
	endif

! Read start and stop time

	call nfs(nf90_get_att(ncid,nf90_global,'first_meas_time',arg))
	first_measurement_time = strp1985f (arg)
	call nfs(nf90_get_att(ncid,nf90_global,'last_meas_time',arg))
	last_measurement_time = strp1985f (arg)

! Skip passes of which the cycle number or equator crossing time is outside the specified interval

	if (first_measurement_time < S%time%info%limits(1) .or. last_measurement_time > S%time%info%limits(2) .or. &
		cycle_number < S%cycles(1) .or. cycle_number > S%cycles(2)) then
		call nfs(nf90_close(ncid))
		call log_string ('Skipped', .true.)
		cycle
	endif

! Determine L2 processing baseline, characterisation and configuration versions

	call nfs(nf90_get_att(ncid,nf90_global,'FDR_input',arg))

! Initialize

	nvar = 0
	allocate (tmp(nrec_m1))

	if (ndata+nrec_m1 > mrec) then
		call log_string ('Error: too many input measurements', .true.)
		stop
	endif

! Open the existing pass and patch it

	if (oc_product) then
!		if (nrec_m20 /= nrec_m1*20) then
!			write (rads_log_unit,*) 'Skipped: nrec_m20 does not match nrec_m1: ', nrec_m20, nrec_m1*20
!		else
			call process_pass_oc (nrec_m1, nrec_m20)
!		endif
	else
!		if (nrec_m5 /= nrec_m1*5) then
!			write (rads_log_unit,*) 'Skipped: nrec_m5 does not match nrec_m1: ', nrec_m5, nrec_m1*5
!		else
			call process_pass_wa (nrec_m1, nrec_m5)
!		endif
	endif

	deallocate (tmp)
	call nfs(nf90_close(ncid))
	ndata = ndata + nrec_m1

	! To which passes do the data belong?
	call which_pass (var(1)%d(1))

	! Look where to split this chunk of data
	do i = 2,ndata
		if (var(1)%d(i) >= orf(ipass+1)%starttime) exit
	enddo
	! It is OK to exit this loop with i = ndata + 1. This means we dump all of the memory.

	! Write out the data
	nout = i - 1 ! Number of measurements to be written out

! Write out the pending data

	call put_rads (nout)

! Number of measurements remaining
	ndata = ndata - nout

! Move the data to be beginning
	forall (i = 1:nvar)
		var(i)%d(1:ndata) = var(i)%d(nout+1:nout+ndata)
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
! OCEAN AND COASTAL PRODUCT (TDP_OC)
!-----------------------------------------------------------------------

subroutine process_pass_oc (n1, n20)
integer(fourbyteint), intent(in) :: n1, n20
integer(fourbyteint) :: i, j, n
real(eightbytereal) :: a20(n20), t20(n20), a(n1), b(n1), t(n1), par_a, par_b, par_r
integer(fourbyteint) :: nr(n1)

! Store the FDR times

call get_var (ncid_m1, 'time', t)
t = t * 86400 + sec1990	! convert from days since 1990 to seconds since 1985
call new_var ('time', t)

! Location

call get_var (ncid_m1, 'latitude', a)
do i = 1,n1
	b(i) = dhellips(1,a(i))
enddo
call get_var (ncid_x1, 'altitude', a)
call new_var ('alt_reaper_deos', a+b)

! Range

call cpy_var (ncid_x1, 'range', 'range_ku')

call get_var (ncid_m20, 'time', t20)
t20 = t20 * 86400 + sec1990	! convert from days since 1990 to seconds since 1985
call get_var (ncid_x20, 'range altitude SUB', a20)

! Find matching 1-Hz and 20-Hz times and do regression
j = 1
n = 0
nr = 0
b = nan
do i = 1,n1
	do
	if (j > n20 .or. t20(j) > t(i)+S%dt1hz/2) then
		call regression (t20(j-n:j-1)-t(i), a20(j-n:j-1), par_a, par_b, par_r, b(i), nr(i))
		n = 0
		exit
	else if (t20(j) < t(i)-S%dt1hz/2) then
		j = j + 1
	else
		n = n + 1
		j = j + 1
	endif
	enddo
enddo

call new_var ('range_rms_ku', sqrt(b))
call new_var ('range_numval_ku', dble(nr))

! Path delay

call cpy_var (ncid_x1, 'dry_tropospheric_correction', 'dry_tropo_era5')
call cpy_var (ncid_x1, 'wet_tropospheric_correction', 'wet_tropo_rad')
if (S%sat == 'e1' .and. cycle_number < 106) then
	call cpy_var (ncid_x1, 'ionospheric_correction', 'iono_nic09')
else
	call cpy_var (ncid_x1, 'ionospheric_correction', 'iono_gim')
endif

! SSB

! Compute the combined sea state bias plus high-frequency correction
call cpy_var (ncid_x1, 'sea_state_bias', 'ssb_tran2019_3d')

! IB

call cpy_var (ncid_x1, 'dynamic_atmospheric_correction', 'inv_bar_mog2d_era')

! Tides

call cpy_var (ncid_x1, 'internal_tide', 'tide_internal')
call cpy_var (ncid_x1, 'solid_earth_tide', 'tide_solid')
call cpy_var (ncid_x1, 'pole_tide', 'tide_pole')
! Not including ocean and load tides because load tide is not available

! MSS and MDT

call cpy_var (ncid_x1, 'mean_dynamic_topography', 'mdt_cnescls18')
! Skipping MSS because it is outdated

! Surface type and coastal proximity

call get_var (ncid_x1, 'surface_type', tmp)
where (tmp == 1) tmp = 3
where (tmp == 4) tmp = 2
call new_var ('surface_type', tmp)
call cpy_var (ncid_m1, 'distance_to_coast 1e-3 MUL', 'dist_coast') ! Convert m to km

! SSHA

call cpy_var (ncid_m1, 'sea_level_anomaly', 'ssha')

end subroutine process_pass_oc

!-----------------------------------------------------------------------
! WAVES PRODUCT (TDP_WA)
!-----------------------------------------------------------------------

subroutine process_pass_wa (n1, n5)
integer(fourbyteint), intent(in) :: n1, n5
real(eightbytereal) :: a5(n5)

! Match the FDR times with the already existing times

call get_var (ncid, 'time', a5)
a5 = a5 * 86400 + sec1990	! convert from days since 1990 to seconds since 1985
call rads_get_var (S, P, 'time', tmp)
if (any(abs(a5(3:nrec_m5:5) - tmp) > 5d-2)) then
	call log_string ('Warning: times do not match')
endif

! SWH

call get_var (ncid, 'swh_adjusted_filtered', a5)
call new_var ('swh_ku', a5(3:nrec_m5:5))
call get_var (ncid, 'swh_uncertainty', a5)
call new_var ('swh_rms_ku', a5(3:nrec_m5:5))

end subroutine process_pass_wa

!-----------------------------------------------------------------------
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads (nout)
integer(fourbyteint), intent(in) :: nout
real(eightbytereal), allocatable :: t(:)
integer(fourbyteint), allocatable :: idx(:)
character(len=80) :: message
integer :: i, j

if (nout == 0) return	! Skip empty data sets

! Open output file
call rads_open_pass (S, P, orf(ipass)%cycle, orf(ipass)%pass, .true.)
if (P%ndata == 0) then
	call log_string ('Error: RADS product does not exist', .true.)
	call log_records (0)
	return
else if (P%ndata < nout) then
	call log_string ('Error: RADS product is smaller than FDR4ALT', .true.)
	call log_records (0)
	return
endif

! Work out time consistency since some measurements are missing in the FDR4ALT products
! There seems to be no apparent reason why they are missing in FDR4ALT products, but they appear to match
! REAPER data with just NaN measurement

allocate (t(0:P%ndata),idx(P%ndata))
call rads_get_var (S, P, 'time', t(1:))
idx = 0
j = 0
do i = 1,P%ndata
	do
		if (var(1)%d(j) < t(i) - 2d-6) then
			j = j + 1
		else if (var(1)%d(j) > t(i) + 2d-6) then
			exit
		else
			idx(i) = j
			exit
		endif
	enddo
enddo

if (any(idx == 0)) then
	write (message, "('Warning:',i4,' FDRALT times do not match RADS')") count(idx == 0)
	call log_string (trim(message), .false.)
endif

! Check which variables are empty or all zero
do i = 1,nvar
	var(i)%empty = all(isnan_(var(i)%d(1:nout)))
	var(i)%zero = all(var(i)%d(1:nout) == 0d0)
enddo

! Write out the empty variables to be kept
if (any(var(1:nvar)%empty)) then
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

! Write to history
call rads_put_history (S, P)

! Define all variables (start from 2 as we are not copying time)
do i = 2,nvar
	call rads_def_var (S, P, var(i)%v)
enddo

! Fill all the data fields
t(0) = nan
do i = 2,nvar
	t(1:nout) = var(i)%d(1:nout)
	call rads_put_var (S, P, var(i)%v, t(idx))
enddo
deallocate (t,idx)

! Close the data file
call log_records (nout, P)
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
var(nvar)%d(ndata+1:ndata+nrec_m1) = data
if (present(ndims)) var(nvar)%v%info%ndims = ndims
end subroutine new_var

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional, intent(in) :: flag
if (rads_version ('Write FDR4ALT data to RADS', flag=flag)) return
call synopsis_devel (' [processing_options] < list_of_FDR4ALT_file_names')
stop
end subroutine synopsis

end program rads_add_f4a_ers
