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

!*rads_add_f4a -- Converts FDR4ALT altimeter data to RADS
!+
program rads_add_f4a

! This program reads FDR4ALT altimeters files and adds them to existing
! RADS data files
!
! syntax: rads_add_f4a [options] < list_of_FDR4ALT_file_names
!
! This program handles FDR4ALT Level 2 products in NetCDF format.
! The format and content is described in:
!
! [1] FDR4ALT Product User Guide
!     CLS-ENV-MU-23-0237, issue 2.2, 30 Oct 2023
!
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! range_* - Ocean range (retracked)
! range_rms_* - Std dev of range
! range_numval_* - Nr of averaged range measurements
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! wet_tropo_rad - Radiometer wet tropo correction
! iono_alt - Dual-frequency ionospheric correction
! iono_alt_smooth - Filtered dual-frequency ionospheric correction
! iono_gim - GIM ionosphetic correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! ssb_cls_* - SSB
! swh_* - Significant wave height
! swh_rms_* - Std dev of SWH
! sig0_* - Sigma0
! sig0_rms_* - Std dev of sigma0
! dsig0_atmos_* - Atmospheric attenuation of sigma0
! wind_speed_alt - Altimeter wind speed
! wind_speed_rad - Radiometer wind speed
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! tide_ocean/load_got410 - GOT4.10c ocean and load tide
! tide_ocean/load_fes14 - FES2014 ocean and load tide
! tide_pole - Pole tide
! tide_solid - Solid earth tide
! topo_ace2 - ACE2 topography
! geoid_egm2008 - EGM2008 geoid
! cnescls15 - CNES/CLS15 mean sea surface
! mss_dtu18 - DTU18 mean sea surface
! tb_238 - Brightness temperature (23.8 GHz)
! tb_365 - Brightness temperature (36.5 GHz)
! flags - Engineering flags
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! rad_liquid_water - Liquid water content
! rad_water_vapor - Water vapor content
! ssha - Sea surface height anomaly
!
! Extensions _* are:
! _ku:      Ku-band
! _s:       S-band
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
	nrec_m1, nrec_x1, nrec_m20, nrec_x20
real(eightbytereal) :: first_measurement_time, last_measurement_time
logical :: oc_product
character(len=16) :: ext

! Data variables

real(eightbytereal) :: dt
integer(fourbyteint) :: i

! Other local variables

real(eightbytereal), parameter :: sec1990=157766400d0	! UTC seconds from 1 Jan 1985 to 1 Jan 1990

! Initialise

call synopsis ('--head')
!call rads_set_options (' uso all')
call rads_init (S)
dt = S%dt1hz/2

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
		index (arg, 'FDR4ALT Thematic Data Product:') == 0) then
		call log_string ('Error: this is not an FDR4ALT altimeter Level 2 data set', .true.)
		cycle
	endif
	oc_product = index(arg, 'Ocean and Coastal') > 0

! Get the mission name and compare with opened RADS session

	call nfs(nf90_get_att(ncid,nf90_global,'mission_name',arg))
	if ((arg == 'ERS-1' .and. S%sat == 'e1') .or. &
		(arg == 'ERS-2' .and. S%sat == 'e2') .or. &
		(arg == 'ENVISAT' .and. S%sat== 'n1')) then
		! All good options
	else
		call log_string ('Error: input file does not match selected satellite', .true.)
		cycle
	endif
	if (S%sat == 'n1') then
		ext = '_adaptive'
	else
		ext = ''
	endif

! Get NetCDF ID for 1-Hz and 20-Hz data (main and expert) in case of ocean/coastal product

	if (oc_product) then
		call nfs(nf90_inq_grp_full_ncid(ncid, '/main/data_01', ncid_m1))
		call nfs(nf90_inq_grp_full_ncid(ncid, '/expert/data_01', ncid_x1))
		call nfs(nf90_inq_grp_full_ncid(ncid, '/main/data_20', ncid_m20))
		call nfs(nf90_inq_grp_full_ncid(ncid, '/expert/data_20', ncid_x20))

! Read global attributes

		call nfs(nf90_inq_dimid(ncid_m1, 'time', varid))
		call nfs(nf90_inquire_dimension(ncid_m1, varid, len=nrec_m1))
		call nfs(nf90_inq_dimid(ncid_x1, 'time', varid))
		call nfs(nf90_inquire_dimension(ncid_x1, varid, len=nrec_x1))
		call nfs(nf90_inq_dimid(ncid_m20, 'time', varid))
		call nfs(nf90_inquire_dimension(ncid_m20, varid, len=nrec_m20))
		call nfs(nf90_inq_dimid(ncid_x20, 'time', varid))
		call nfs(nf90_inquire_dimension(ncid_x20, varid, len=nrec_x20))
		if (nrec_m1 /= nrec_x1) then
			call log_string ('Error: file skipped: nr of 1-Hz main and expert records not the same', .true.)
			cycle
		else if (nrec_m20 /= nrec_x20) then
			call log_string ('Error: file skipped: nr of 20-Hz main and expert records not the same', .true.)
			cycle
		endif
	else

! Less to check for wave product

		call nfs(nf90_inq_dimid(ncid, 'time', varid))
		call nfs(nf90_inquire_dimension(ncid, varid, len=nrec_m1))
	endif

! Check for empty products (just in case)

	if (nrec_m1 == 0) then
		call log_string ('Error: file skipped: no measurements', .true.)
		cycle
	endif

! Read cycle and pass number, start and stop time, equator time and longitude

	call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cycle_number))
	call nfs(nf90_get_att(ncid,nf90_global,'pass_number',pass_number))
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
	if (nrec_m1 /= nrec) then
		write (rads_log_unit,*) 'Warning: nrec_m1 does not match nrec: ',nrec_m1,nrec
	else if (oc_product .and. nrec_m1*20 /= nrec_m20) then
		write (rads_log_unit,*) 'Warning: nrec_m20 does not match nrec_m1: ',nrec_m20,nrec_m1*20
	else if (P%ndata <= 0) then
		! Skip
	else if (oc_product) then
		call process_pass_oc (nrec_m1, nrec_m20)
	else
		call process_pass_wa (nrec_m1)
	endif
	call rads_close_pass (S, P)
	call nfs(nf90_close(ncid))

enddo

call rads_end (S)

contains

subroutine process_pass_oc (n, n20)
integer(fourbyteint), intent(in) :: n, n20
integer(fourbyteint) :: nr(n)
real(eightbytereal) :: b(n), t(n), dh(n), a20(n20), t20(n20)

! Allocate data array

allocate (a(n),flags(n))

! Match the FDR times with the already existing times

call get_var (ncid_m1, 'time', a)
a = a * 86400 + sec1990	! convert from days since 1990 to seconds since 1985
call rads_get_var (S, P, 'time', t)
if (any(abs(a - t) > 0.5d-6)) then
	call log_string('Warning: times do not match')
	return
endif

! Compute ellipsoid correction

call get_var (ncid_m1, 'latitude', a)
do i = 1,n
	dh(i) = dhellips(1,a(i))	! Convert WGS-84 to TOPEX
enddo

! Range

call cpy_var (ncid_x1, 'range', 'range_ku' // ext)

call get_var (ncid_m20, 'time', t20)
call get_var (ncid_x20, 'range altitude SUB', a20)
call trend_1hz (reshape(t20, (/20,n/)), t, reshape(a20, (/20,n/)), a, b, nr)
call new_var ('range_rms_ku' // ext, b)
call new_var ('range_numval_ku' // ext, dble(nr))


! MAYBE NOT SUCH A GOOD IDEA TO ASSIGN THIS TO qual_range!
call rads_get_var (S, P, 'flags', a)
flags = nint(a, twobyteint)
call nc2f (ncid_m1, 'validation_flag', 11)				! bit 11: Quality range
call new_var ('flags' // ext, dble(flags))
call cpy_var (ncid_m1, 'validation_flag', 'qual_range' // ext)	! Make special entry in rads.xml

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
if (S%sat == 'n1') call cpy_var (ncid_x1, 'range_ssb_hfa range SUB', 'ssb_tran2019_hfa')
call cpy_var (ncid_x1, 'sea_state_bias', 'ssb_tran2019_3d')

! IB

call cpy_var (ncid_x1, 'dynamic_atmospheric_correction', 'inv_bar_mog2d_era')

! Tides

! Load tide is not available, so we cannot split
!call cpy_var (ncid_x1, 'ocean_tide_height load_tide SUB', 'tide_ocean_fes14')
!call cpy_var (ncid_x1, 'load_tide', 'tide_load_fes14')
call cpy_var (ncid_x1, 'internal_tide', 'tide_internal')
call cpy_var (ncid_x1, 'solid_earth_tide', 'tide_solid')
call cpy_var (ncid_x1, 'pole_tide', 'tide_pole')

! MSS and MDT

call get_var (ncid_x1, 'mean_sea_surface', a)
call new_var ('mss_comb15', a + dh)
call cpy_var (ncid_x1, 'mean_dynamic_topography', 'mdt_cnescls18')

! Surface type and coastal proximity

call get_var (ncid_x1, 'surface_type', a)
where (a == 1) a = 3
where (a == 4) a = 2
call new_var ('surface_type', a)
call cpy_var (ncid_m1, 'distance_to_coast 1e-3 MUL', 'dist_coast') ! Convert m to km

! SSHA

call cpy_var (ncid_m1, 'sea_level_anomaly', 'ssha' // ext)

! Close pass

call put_rads (skip_create = .true.)

deallocate (a, flags)

end subroutine process_pass_oc

subroutine process_pass_wa (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: t(n)

! Allocate data array

allocate (a(n))

! Match the FDR times with the already existing times

call get_var (ncid_m1, 'time', a)
a = a * 86400 + sec1990	! convert from days since 1990 to seconds since 1985
call rads_get_var (S, P, 'time', t)
if (any(abs(a - t) > 0.5d-6)) then
	call log_string('Warning: times do not match')
	return
endif

! SWH

call cpy_var (ncid_m1, 'swh_ocean', 'swh_ku' // ext)
call cpy_var (ncid_m1, 'swh_ocean_rms', 'swh_rms_ku' // ext)

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
call synopsis_devel (' < list_of_FDR4ALT_file_names')
write (*,1310)
1310 format (/ &
'This program adds FDR4ALT altimeter Level 2 products to existing RADS data files.')
stop
end subroutine synopsis

end program rads_add_f4a
