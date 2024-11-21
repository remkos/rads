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
use rads_devel_netcdf
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

! Other local variables

real(eightbytereal), parameter :: sec1990=157766400d0	! UTC seconds from 1 Jan 1985 to 1 Jan 1990
integer(fourbyteint), allocatable :: idx(:)

! Initialise

call synopsis ('--head')
call rads_init (S)

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

do
	read (*,'(a)',iostat=ios) infile
	if (ios /= 0) exit

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

! Open the existing pass and patch it

	call rads_open_pass (S, P, cycle_number, pass_number, .true.)
	nrec = P%ndata

	write (*,*) "nrec, nrec_m1, nrec_m20 = ", nrec, nrec_m1, nrec_m20

	if (nrec <= 0) then
		! Skip
	else if (oc_product) then
		if (nrec_m20 /= nrec*20) then
			write (rads_log_unit,*) 'Skipped: nrec_m20 does not match nrec: ', nrec_m20,nrec*20
		else
			call process_pass_oc (nrec, nrec_m1, nrec_m20)
		endif
	else
		if (nrec_m5 /= nrec*5) then
			write (rads_log_unit,*) 'Skipped: nrec_m5 does not match nrec: ',nrec_m5,nrec*5
		else
			call process_pass_wa (nrec, nrec_m5)
		endif
	endif
	call rads_close_pass (S, P)
	call nfs(nf90_close(ncid))

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! OCEAN AND COASTAL PRODUCT (TDP_OC)
!-----------------------------------------------------------------------

subroutine process_pass_oc (n, n1, n20)
integer(fourbyteint), intent(in) :: n, n1, n20
integer(fourbyteint) :: i, j
real(eightbytereal) :: a20(n20), t20(n20), b(n), t(n)
integer(fourbyteint) :: nr(n)

! Some TDP_OC products have duplicated 1-Hz times. In that case we need to reduce the number of 1-Hz records.
! Duplicate records have the altitude set to NaN

allocate (a(n1),idx(n))
if (n1 == n) then
	do j = 1,n
		idx(j) = j
	enddo
else
	call log_string ('Warning: compacting 1-Hz arrays from TDP_OC product')
	call get_var (ncid_x1, 'altitude', a)
	j = 0
	do i = 1,n1
		if (isan_(a(i))) then
			j = j + 1
			idx(j) = i
		endif
	enddo
endif

nvar = 0

! Match the FDR times with the already existing times

call get_var (ncid_m1, 'time', a)
a = a * 86400 + sec1990	! convert from days since 1990 to seconds since 1985
call rads_get_var (S, P, 'time', t)
if (any(abs(a(idx) - t) > 1.5d-6)) then
	call log_string ('Warning: times do not match', .true.)
	deallocate (a,idx)
	return
endif

! Location

call cpy_var_i (ncid_m1, 'latitude', 'lat')
do i = 1,n
	b(i) = dhellips(1,a(idx(i)))
enddo
call cpy_var_i (ncid_m1, 'longitude', 'lon')
call get_var (ncid_x1, 'altitude', a)
call new_var ('alt_reaper_deos', a(idx)+b)

! Range

call cpy_var_i (ncid_x1, 'range', 'range_ku')

call get_var (ncid_m20, 'time', t20)
call get_var (ncid_x20, 'range altitude SUB', a20)
call trend_1hz (reshape(t20, (/20,n/)), t, reshape(a20, (/20,n/)), a(:n), b, nr)
call new_var ('range_rms_ku', b)
call new_var ('range_numval_ku', dble(nr))

! MAYBE NOT SUCH A GOOD IDEA TO ASSIGN THIS TO qual_range!
! call rads_get_var (S, P, 'flags', a)
! flags = nint(a, twobyteint)
! call nc2f (ncid_m1, 'validation_flag', 11)				! bit 11: Quality range
! call new_var ('flags', dble(flags))
! call cpy_var_i (ncid_m1, 'validation_flag', 'qual_range')

! Path delay

call cpy_var_i (ncid_x1, 'dry_tropospheric_correction', 'dry_tropo_era5')
call cpy_var_i (ncid_x1, 'wet_tropospheric_correction', 'wet_tropo_rad')
if (S%sat == 'e1' .and. cycle_number < 106) then
	call cpy_var_i (ncid_x1, 'ionospheric_correction', 'iono_nic09')
else
	call cpy_var_i (ncid_x1, 'ionospheric_correction', 'iono_gim')
endif

! SSB

! Compute the combined sea state bias plus high-frequency correction
call cpy_var_i (ncid_x1, 'sea_state_bias', 'ssb_tran2019_3d')

! IB

call cpy_var_i (ncid_x1, 'dynamic_atmospheric_correction', 'inv_bar_mog2d_era')

! Tides

call cpy_var_i (ncid_x1, 'internal_tide', 'tide_internal')
call cpy_var_i (ncid_x1, 'solid_earth_tide', 'tide_solid')
call cpy_var_i (ncid_x1, 'pole_tide', 'tide_pole')
! Not including ocean and load tides because load tide is not available

! MSS and MDT

call cpy_var_i (ncid_x1, 'mean_dynamic_topography', 'mdt_cnescls18')
! Skipping MSS because it is outdated

! Surface type and coastal proximity

call get_var (ncid_x1, 'surface_type', a)
where (a == 1) a = 3
where (a == 4) a = 2
call new_var ('surface_type', a(idx))
call cpy_var_i (ncid_m1, 'distance_to_coast 1e-3 MUL', 'dist_coast') ! Convert m to km

! SSHA

call cpy_var_i (ncid_m1, 'sea_level_anomaly', 'ssha')

! Close pass

call put_rads (skip_create = .true.)

deallocate (a,idx)

end subroutine process_pass_oc

! Version of cpy_var with index

subroutine cpy_var_i (ncid, varin, varout)
integer(fourbyteint), intent(in) :: ncid
character(len=*), intent(in) :: varin, varout
call get_var (ncid, varin, a)
call new_var (varout, a(idx))
end subroutine cpy_var_i

!-----------------------------------------------------------------------
! WAVES PRODUCT (TDP_WA)
!-----------------------------------------------------------------------

subroutine process_pass_wa (n, n5)
integer(fourbyteint), intent(in) :: n, n5
real(eightbytereal) :: a5(n5)

! Allocate data array

allocate (a(n))
nvar = 0

! Match the FDR times with the already existing times

call get_var (ncid, 'time', a5)
a5 = a5 * 86400 + sec1990	! convert from days since 1990 to seconds since 1985
call rads_get_var (S, P, 'time', a)
if (any(abs(a5(3:nrec_m5:5) - a) > 5d-2)) then
	call log_string ('Warning: times do not match')
endif

! SWH

call get_var (ncid, 'swh_adjusted_filtered', a5)
call new_var ('swh_ku', a5(3:nrec_m5:5))
call get_var (ncid, 'swh_uncertainty', a5)
call new_var ('swh_rms_ku', a5(3:nrec_m5:5))

! Close pass

call put_rads (skip_create = .true.)

deallocate (a)

end subroutine process_pass_wa

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional, intent(in) :: flag
if (rads_version ('Write FDR4ALT data to RADS', flag=flag)) return
call synopsis_devel (' [processing_options] < list_of_FDR4ALT_file_names')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:'/ &
'  -x, --no-ext              Do not add _adaptive extension for Envisat')
stop
end subroutine synopsis

end program rads_add_f4a_ers
