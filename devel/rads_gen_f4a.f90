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

!*rads_gen_f4a -- Converts FDR4ALT altimeter data to RADS
!+
program rads_gen_f4a

! This program reads FDR4ALT altimeters files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/SS/F/SSpPPPPcCCC.nc.
!    SS = satellite (e1, e2, or n1)
!     F = phase
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_f4a [options] < list_of_FDR4ALT_file_names
!
! where [options] include:
!  --min-rec <min_rec> : Specify minimum number of records per pass to process.
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
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdrf - Orbital altitude
! alt_rate - Orbital altitude rate
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
use rads_gen
use rads_misc
use rads_netcdf
use rads_time
use rads_geo
use netcdf

! Command line arguments

integer(fourbyteint) :: ios
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: ncid, ncid_m1, ncid_x1, ncid_m20, ncid_x20, cycle, pass, varid, &
	nrec_m1, nrec_x1, nrec_m20, nrec_x20
real(eightbytereal) :: first_measurement_time, last_measurement_time

! Data variables

real(eightbytereal) :: dt
integer(fourbyteint) :: i, j, n

! Other local variables

real(eightbytereal), parameter :: sec1990=157766400d0	! UTC seconds from 1 Jan 1985 to 1 Jan 1990
integer(fourbyteint), allocatable :: nr(:), idx(:)
real(eightbytereal), allocatable :: b(:), t(:), dh(:), a20(:), t20(:), a2(:,:), t2(:,:)

! Initialise

call synopsis
call rads_gen_getopt ('', ' min-rec:')
call synopsis ('--head')
call rads_set_options (' min-rec:')
if (sat /= '') call rads_init (S, sat)

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

do
	read (*,'(a)',iostat=ios) infile
	if (ios /= 0) exit

! Open input file

	call log_string (basename(infile))
	if (nft(nf90_open(infile,nf90_nowrite,ncid))) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Check if input is an FDR4ALT altimeter Level 2 data set

	if (nft(nf90_get_att(ncid,nf90_global,'title',arg)) .or. &
		arg /= 'FDR4ALT Thematic Data Product: Ocean and Coastal') then
		call log_string ('Error: this is not an FDR4ALT altimeter Level 2 data set', .true.)
		cycle
	endif

! Get the mission name and initialise RADS (if not done before)

	call nfs(nf90_get_att(ncid,nf90_global,'mission_name',arg))
	select case (arg(:3))
	case ('ERS-1')
		arg = 'e1'
	case ('ERS-2')
		arg = 'e2'
	case ('ENVISAT')
		arg = 'n1'
	case default
		call log_string ('Error: wrong misson: must be ER1, ER2, or ENV', .true.)
		cycle
	end select
	if (sat == '') then
		call rads_init (S, arg)
	else if (S%sat /= arg(:2)) then
		call log_string ('Error: wrong misson: must be S'//S%sat, .true.)
		cycle
	endif

! Get NetCDF ID for 1-Hz and 20-Hz data (main and expert)

	call nfs(nf90_inq_ncid(ncid, 'main/data_01', ncid_m1))
	call nfs(nf90_inq_ncid(ncid, 'expert/data_01', ncid_x1))
	call nfs(nf90_inq_ncid(ncid, 'main/data_20', ncid_m20))
	call nfs(nf90_inq_ncid(ncid, 'expert/data_20', ncid_x20))

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
	else if (nrec_m1 == 0) then
		call log_string ('Error: file skipped: no measurements', .true.)
		cycle
	else if (nrec_m1 < min_rec) then
		call log_string ('Warning: file skipped: too few measurements', .true.)
		cycle
	else if (nrec_m1 > mrec) then
		call log_string ('Error: file skipped: too many measurements', .true.)
		cycle
	endif

! Read cycle and pass number, start and stop time, equator time and longitude

	call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cycle))
	call nfs(nf90_get_att(ncid,nf90_global,'pass_number',pass))
	call nfs(nf90_get_att(ncid,nf90_global,'first_meas_time',arg))
	first_measurement_time = strp1985f (arg)
	call nfs(nf90_get_att(ncid,nf90_global,'last_meas_time',arg))
	last_measurement_time = strp1985f (arg)

! Initialise pass struct

	call rads_init_pass_struct (S, P)

! NOTE: No equator time and longitude, so we need to compute these

! Set mission phase based on rough equator_time

	call rads_set_phase (S, (first_measurement_time + last_measurement_time) / 2)

! Predict the equator time and longitude

	call rads_predict_equator (S, P)

! Skip passes of which the cycle number or equator crossing time is outside the specified interval

	if (P%equator_time < times(1) .or. P%equator_time > times(2) .or. P%cycle < cycles(1) .or. P%cycle > cycles(2)) then
		call nfs(nf90_close(ncid))
		call log_string ('Skipped', .true.)
		cycle
	endif

! Determine L2 processing baseline, characterisation and configuration versions

	call nfs(nf90_get_att(ncid,nf90_global,'FDR_input',arg))

! Store input file name

	P%original = trim(basename(infile)) // ' (' // trim(arg) // ')'

! Allocate variables

	allocate (a(nrec_m1),b(nrec_m1),dh(nrec_m1),flags(nrec_m1),t(nrec_m1),a20(nrec_m20),t20(nrec_m20))
	allocate (nr(nrec_m1),idx(nrec_m1+1))
	nvar = 0

! Compile flag bits

	flags = 0
	call nc2f (ncid_x1, 'surface_type', 4, eq=1)
	call nc2f (ncid_x1, 'surface_type', 4, eq=3)				! bit  4: water/land
	call nc2f (ncid_x1, 'surface_type', 4, ge=1)				! bit  5: ocean/other
	call nc2f (ncid_m1, 'validation_flag', 11, eq=0)			! bit 11: rejected measurement

! Time and location

	call get_var (ncid_m1, 'time', t)
	t = t * 86400 + sec1990	! convert from days since 1990 to seconds since 1985
	P%start_time = t(1)  ! Using 1-Hz time stamps, not the ones from the header which are
	P%end_time = a(nrec) ! the time stamps of the 20-Hz measurements
	call new_var ('time', t)
	call get_var (ncid_m1, 'latitude', a)
	do i = 1,nrec_m1
		dh(i) = dhellips(1,a(i))	! Convert WGS-84 to TOPEX
	enddo
	call new_var ('lat', a)
	call cpy_var (ncid_m1, 'longitude', 'lon')
	call get_var (ncid_x1, 'altitude', a)
	if (S%sat == 'n1') then
		call new_var ('alt_gdrf', a + dh)
	else
		call new_var ('alt_reaper', a + dh)
	endif

! Make links between 20-Hz and 1-Hz

	call get_var (ncid_m20, 'time', t20)
	t20 = t20 * 86400 + sec1990
	dt = S%dt1hz/2
	j = 1
	do i = 1,nrec_m1
		do while (t20(j) < t(i) - dt)
			j = j + 1
		enddo
		idx(i) = j
	enddo
	do while (j <= nrec_m20 .and. t20(j) <= t(nrec_m1) + dt)
		j = j + 1
	enddo
	idx(nrec_m1+1) = j

! Range

	call get_var (ncid_x20, 'range altitude SUB', a20)
	a2 = nan ; t2 = nan
	do i = 1,nrec_m1
		n = idx(i+1) - idx(i)
		t2(:n,i) = t20(idx(i):idx(i+1)-1)
		a2(:n,i) = a20(idx(i):idx(i+1)-1)
	enddo
	call trend_1hz (t2, t, a2, a, b, nr)
	call new_var ('range_rms_ku', b)
	call new_var ('range_numval_ku', dble(nr))
	call cpy_var (ncid_x1, 'range', 'range_ku')
	call cpy_var (ncid_x1, 'validation_flag', 'qual_range')	! Make special entry in rads.xml

! SWH needs to come from wave product

!	call cpy_var (ncidk, 'swh_ocean', 'swh_ku')
!	call cpy_var (ncidk, 'swh_ocean_rms', 'swh_rms_ku')

! Backscatter is not available. Copy from GDR / REAPER

!	call cpy_var (ncidk, 'sig0_ocean', 'sig0_ku')
!	call cpy_var (ncidk, 'sig0_ocean_rms', 'sig0_rms_ku')
!	call cpy_var (ncidk, 'atm_cor_sig0', 'dsig0_atmos_ku')
!	call cpy_var (ncid1, 'climato_use_flag', 'qual_dsig0_atmos')

! Wind speed. Copy from GDR / REAPER

!	call cpy_var (ncid1, 'wind_speed_alt', 'wind_speed_alt')
!	call cpy_var (ncid1, 'wind_speed_mod_u', 'wind_speed_ecmwf_u')
!	call cpy_var (ncid1, 'wind_speed_mod_v', 'wind_speed_ecmwf_v')

! Rain or ice: not available. Copy from GDR / REAPER

!	call cpy_var (ncid1, 'rain_flag', 'qual_alt_rain_ice')

! Off-nadir angle: not available

! Path delay

	call cpy_var (ncid_x1, 'dry_tropospheric_correction', 'dry_tropo_era5')
	call cpy_var (ncid_x1, 'wet_tropospheric_correction', 'wet_tropo_rad')
	if (S%sat == 'e1' .and. cycle < 106) then
		call cpy_var (ncid_x1, 'ionospheric_correction', 'iono_nic09')
	else
		call cpy_var (ncid_x1, 'ionospheric_correction', 'iono_gim')
	endif

! SSB

	if (S%sat == 'n1') then
		! Compute the combined sea state bias plus high-frequency correction
		call cpy_var (ncid_x1, 'range_ssb_hfa range SUB', 'ssb_tran2019_hfa')
	endif
	call cpy_var (ncid_x1, 'sea_state_bias', 'ssb_tran2019_3d')

! IB

	call cpy_var (ncid_x1, 'dynamic_atmospheric_correction', 'inv_bar_mog2d_era')

! Tides

	call cpy_var (ncid_x1, 'geocentric_ocean_tide load_tide SUB', 'tide_ocean_fes14')
	call cpy_var (ncid_x1, 'load_tide', 'tide_load_fes14')
	call cpy_var (ncid_x1, 'internal_tide', 'tide_internal')
	call cpy_var (ncid_x1, 'solid_earth_tide', 'tide_solid')
	call cpy_var (ncid_x1, 'pole_tide', 'tide_pole')

! MSS

	call get_var (ncid_x1, 'mean_sea_surface', a)
	call new_var ('mss_hybrid', a + dh)

! Surface type and coastal proximity

	call get_var (ncid_x1, 'surface_type', a)
	where (a == 1) a = 3
	where (a == 4) a = 2
	call new_var ('surface_type', a)
	call cpy_var (ncid_m1, 'distance_to_coast 1e-3 MUL', 'dist_coast') ! Convert m to km

! Bit flags

	call new_var ('flags', dble(flags))

! SSHA

	call cpy_var (ncid_m1, 'sea_level_anomaly', 'ssha')

! Close input file

	call nfs(nf90_close(ncid))

! Dump the data

	call put_rads
	deallocate (a, b, t, dh, flags, a20, t20, nr, idx, a2, t2)

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional, intent(in) :: flag
if (rads_version ('Write FDR4ALT data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_FDR4ALT_file_names')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --min-rec=MIN_REC         Specify minimum number of records per pass to process' // &
'This program converts FDR4ALT altimeter Level 2 products to RADS data' / &
'files with the name $RADSDATAROOT/data/SS/F/pPPPP/SSpPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_f4a
