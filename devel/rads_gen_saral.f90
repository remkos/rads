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

!*rads_gen_saral -- Converts SARAL data to RADS
!+
program rads_gen_saral

! This program reads SARAL files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/sa/a/sapPPPPcCCC.nc.
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_saral [options] < list_of_SARAL_file_names
!
! This program handles SARAL OGDR, IGDR and GDR files in NetCDF format.
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdrd - Orbit altitude
! alt_rate - Orbit altitude rate
! range_ka - Ocean range (retracked)
! range_rms_ka - Std dev of range
! range_numval_ka - Nr of averaged range measurements
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_rad - Radiometer wet tropo correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! iono_gim - GIM ionosphetic correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! tide_solid - Solid earth tide
! tide_ocean_fes14 - FES2014 ocean tide
! tide_ocean_got48 - GOT4.8 ocean tide
! tide_load_fes14 - FES2014 load tide
! tide_load_got48 - GOT4.8 load tide
! tide_pole - Pole tide
! ssb_bm3 - SSB
! mss_cnescls11 - CLS01 MSS
! geoid_egm96 - EGM96 geoid
! topo_dtm2000 - Bathymetry
! swh_ka - Significant wave height
! swh_rms_ka - Std dev of SWH
! sig0_ka - Sigma0
! sig0_rms_ka - Std dev of sigma0
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! wind_speed_alt - Altimeter wind speed
! tb_238 - Brightness temperature (K-band)
! tb_370 - Brightness temperature (Ka-band)
! peakiness_ka - Peakiness
! flags - Engineering flags
! off_nadir_angle2_wf_ka - Mispointing from waveform squared
! liquid_water - Liquid water content
! water_vapor_content - Water vapor content
!-----------------------------------------------------------------------
use rads
use rads_devel
use rads_devel_netcdf
use rads_gen
use rads_netcdf
use rads_misc
use rads_time
use netcdf

! Command line arguments

integer(fourbyteint) :: ios, j
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: ncid, cyclenr, passnr, varid
real(eightbytereal) :: equator_time
integer :: latency = rads_nrt

! Data variables

real(eightbytereal), allocatable :: b(:), d(:,:)
logical, allocatable :: valid(:,:)

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000

! Initialise

call synopsis
call rads_gen_getopt ('sa')
call synopsis ('--head')
call rads_init (S, sat)

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
	if (arg /= 'SARAL') then
		call log_string ('Error: wrong misson-name found in header', .true.)
		cycle
	endif

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

! Phase B GDRs start with cycle 100, but have been renumbered to start with cycle 36

	if (S%phase%name(1:1) == 'b') cyclenr = cyclenr - 64

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

	allocate (a(nrec),b(nrec),d(40,nrec),valid(40,nrec),flags(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
	call nc2f (ncid, 'qual_alt_1hz_off_nadir_angle_wf',1)	! bit  1: Quality off-nadir pointing
	! Use surface classification since GDR-F
	call nc2f (ncid, 'surf_class', 2, eq=4)					! bit  2: Continental ice
	call nc2f (ncid, 'surf_class', 4, eq=1)
	call nc2f (ncid, 'surf_class', 4, ge=3)					! bit  4: Water/land
	call nc2f (ncid, 'surf_class', 5, ge=1)					! bit  5: Ocean/other
	call nc2f (ncid, 'rad_surf_type',6)						! bit  6: Radiometer land flag
	call nc2f (ncid, 'ice_flag',7)							! bit  7: Ice flag
	call nc2f (ncid, 'qual_rad_1hz_tb_k',9)					! bit  9: Quality 23.8 GHz channel
	call nc2f (ncid, 'qual_rad_1hz_tb_ka',10)				! bit 10: Quality 37.0 GHz channel
	call nc2f (ncid, 'qual_alt_1hz_range',11)				! bit 11: Quality of range
	call nc2f (ncid, 'qual_alt_1hz_swh',12)					! bit 12: Quality of SWH
	call nc2f (ncid, 'qual_alt_1hz_sig0',13)				! bit 13: Quality of sigma0
	if (latency == rads_nrt) then
		call nc2f (ncid, 'orb_state_flag_diode',15,ge=2)	! bit 15: Quality of DIODE orbit
	else
		call nc2f (ncid, 'orb_state_flag_rest',15,neq=3)	! bit 15: Quality of restituted orbit
	endif

! Convert all the necessary fields to RADS

	call get_var (ncid, 'time', a)
	call new_var ('time', a + sec2000)
	call cpy_var (ncid, 'lat')
	call cpy_var (ncid, 'lon')
	call cpy_var (ncid, 'alt', 'alt_gdrf')	! Since GDR-F
	call cpy_var (ncid, 'orb_alt_rate', 'alt_rate')
	call cpy_var (ncid, 'range', 'range_ka')
	call cpy_var (ncid, 'range_numval', 'range_numval_ka')
	call cpy_var (ncid, 'range_rms', 'range_rms_ka')
	call cpy_var (ncid, 'model_dry_tropo_corr', 'dry_tropo_ecmwf')
	call cpy_var (ncid, 'model_wet_tropo_corr', 'wet_tropo_ecmwf')
	call cpy_var (ncid, 'rad_wet_tropo_corr', 'wet_tropo_rad')
	call cpy_var (ncid, 'iono_corr_gim', 'iono_gim')
	call cpy_var (ncid, 'sea_state_bias', 'ssb_tran2019')
	if (latency .ne. rads_nrt) then
		call cpy_var (ncid, 'sea_state_bias_3d_mp2', 'ssb_tran2019_3d')	! Since GDR-F
	endif
	call cpy_var (ncid, 'swh', 'swh_ka')
	call cpy_var (ncid, 'swh_rms', 'swh_rms_ka')
	call cpy_var (ncid, 'sig0', 'sig0_ka')
	call cpy_var (ncid, 'atmos_corr_sig0', 'dsig0_atmos_ka')
	call cpy_var (ncid, 'off_nadir_angle_wf', 'off_nadir_angle2_wf_ka')
	call cpy_var (ncid, 'tb_k', 'tb_238')
	call cpy_var (ncid, 'tb_ka', 'tb_370')
	call cpy_var (ncid, 'mean_sea_surface_sol1', 'mss_cnescls15')	! Since GDR-F
	call cpy_var (ncid, 'mean_sea_surface_sol2', 'mss_dtu15')		! Since GDR-F
	call cpy_var (ncid, 'geoid', 'geoid_egm96')
	call cpy_var (ncid, 'bathymetry', 'topo_dtm2000')
	call cpy_var (ncid, 'inv_bar_corr', 'inv_bar_static')
	if (latency == rads_nrt) then
		call cpy_var (ncid, 'inv_bar_corr', 'inv_bar_mog2d')
	else
		call cpy_var (ncid, 'inv_bar_corr hf_fluctuations_corr ADD', 'inv_bar_mog2d')
	endif
	call cpy_var (ncid, 'ocean_tide_sol1 load_tide_sol1 SUB', 'tide_ocean_got410')	! Since GDR-F
	call cpy_var (ncid, 'ocean_tide_sol2 load_tide_sol2 SUB ocean_tide_non_equil ADD', 'tide_ocean_fes14')	! Since GDR-F
	call cpy_var (ncid, 'load_tide_sol1', 'tide_load_got410')	! Since GDR-F
	call cpy_var (ncid, 'load_tide_sol2', 'tide_load_fes14')	! Since GDR-F
	call cpy_var (ncid, 'ocean_tide_equil', 'tide_equil')
	call cpy_var (ncid, 'ocean_tide_non_equil', 'tide_non_equil')
	call cpy_var (ncid, 'internal_tide', 'tide_internal')
	call cpy_var (ncid, 'solid_earth_tide', 'tide_solid')
	call cpy_var (ncid, 'pole_tide', 'tide_pole')
	call cpy_var (ncid, 'wind_speed_model_u', 'wind_speed_ecmwf_u')
	call cpy_var (ncid, 'wind_speed_model_v', 'wind_speed_ecmwf_v')
	call cpy_var (ncid, 'wind_speed_alt')
	call cpy_var (ncid, 'mean_wave_period_t02', 'mean_wave_period')
	call cpy_var (ncid, 'mean_wave_direction', 'mean_wave_direction')
	call cpy_var (ncid, 'rad_water_vapor', 'water_vapor_rad')
	call cpy_var (ncid, 'rad_liquid_water', 'liquid_water_rad')
	if (latency == rads_nrt) then
		call get_var (ncid, 'orb_state_flag_diode', a)
		where (a == 9)
			a = 1
		elsewhere
			a = 0
		endwhere
	else
		call get_var (ncid, 'orb_state_flag_rest', a)
		where (a == 4)
			a = 1
		elsewhere
			a = 0
		endwhere
	endif
	call new_var ('flag_manoeuvre', a)
	call get_var (ncid, 'range_used_40hz', d)
	valid = (d == 0d0)
	call get_var (ncid, 'peakiness_40hz', d)
	where (.not.valid) d = nan
	call mean_1hz (d, a, b)
	call new_var ('peakiness_ka', a)
	a = flags
	call new_var ('flags', a)
	a = latency
	call new_var ('latency', a)

! Dump the data

	call nfs(nf90_close(ncid))
	call put_rads
	deallocate (a, b, d, valid, flags)

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Write SARAL data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_SARAL_file_names')
write (*,1310)
1310 format (/ &
'This program converts SARAL OGDR/IGDR/GDR-F files to RADS data' / &
'files with the name $RADSDATAROOT/data/sa/a/pPPPP/sapPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_saral
