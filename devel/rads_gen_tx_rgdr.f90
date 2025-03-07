!-----------------------------------------------------------------------
! Copyright (c) 2011-2025  Remko Scharroo
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

!*rads_gen_tx_rgdr -- Converts TOPEX Retracked GDR data to RADS
!+
program rads_gen_tx_rgdr

! This program reads TOPEX Retracked GDR files and converts them to the
! RADS format, written into files $RADSDATAROOT/data/tx.r50/a/txpPPPPcCCC.nc.
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_tx_rgdr [options] < list_of_RGDR_file_names
!
! This program handles TOPEX Retracked GDR files in NetCDF format,
! version 5 (released Jan 2015)
! The format is described in:
!
! [1] Retracked GDR Data Record, Release 5.0
!     File: retrk-gdr-data-rec-r50nc-r2.141230.xls
!
! [2] Retracked GDR Rel 5.0
!     File: retrk-gdr-data-rec-r50-f0.150112.pdf
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdrd - Orbit altitude
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
! Extionsions _* are:
! _ku:      Ku-band retracked with MLE4
! _ku_mle3: Ku-band retracked with MLE3
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

integer(fourbyteint) :: ios, i
character(len=rads_cmdl) :: infile, arg

! Header variables

integer(fourbyteint) :: ncid, cyclenr, passnr, varid
real(eightbytereal) :: equator_time

! Data variables

real(eightbytereal), allocatable :: b(:), d(:,:)
integer(fourbyteint) :: nvar_range_c, nvar_range_ku_retrk, nvar_range_c_retrk

! Other local variables

real(eightbytereal), parameter :: sec1992 = 220838400d0	! UTC seconds from 1 Jan 1985 to 1 Jan 1992

! Initialise

call synopsis
call rads_gen_getopt ('tx.r50')
call synopsis ('--head')
call rads_init (S, sat)

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

do
	read (*,'(a)',iostat=ios) infile
	if (ios /= 0) exit

! Open input file

	call log_string (infile)
	if (nf90_open(infile,nf90_nowrite,ncid) /= nf90_noerr) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Check if input is RGDR

	if (index(infile,'RGDR_') <= 0) then
		call log_string ('Error: this is not TOPEX RGDR', .true.)
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
	if (arg /= 'TOPEX/POSEIDON') then
		call log_string ('Error: wrong misson-name found in header', .true.)
		cycle
	endif

	call nfs(nf90_get_att(ncid,nf90_global,'title',arg))
	if (arg(:4) /= 'RGDR') then
		call log_string ('Error: wrong file type (title) found in header', .true.)
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

! Store relevant info

	call rads_init_pass_struct (S, P)
	P%cycle = cyclenr
	P%pass = passnr
	P%equator_time = equator_time
	call nfs(nf90_get_att(ncid,nf90_global,'equator_longitude',arg))
	read (arg,*) P%equator_lon
	call nfs(nf90_get_att(ncid,nf90_global,'first_meas_time',arg))
	P%start_time = strp1985f(arg)
	call nfs(nf90_get_att(ncid,nf90_global,'last_meas_time',arg))
	P%end_time = strp1985f(arg)
	call nfs(nf90_get_att(ncid,nf90_global,'version',arg))
	P%original = trim(basename(infile)) // ' (' // trim(arg) // ')'

! Allocate variables

	allocate (a(nrec), b(nrec), d(10,nrec), flags(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
	call get_var (ncid, 'alton', a)
	if (any(a /= 1d0)) call log_string ('Error: wrong altimeter found in data', .true.)
	call nc2f (ncid, 'instr_state_topex', 0, eq=8)		! bit  0: Altimeter Side A/B
!	call nc2f (ncid, 'qual_alt_1hz_off_nadir_angle_wf',1)	! bit  1: Quality off-nadir pointing
	call nc2f (ncid, 'instr_state_tmr', 2, mask=32) ! #5	! bit  2: TMR status (0=Channel 21A, 1=Channel 21B)
	call nc2f (ncid, 'alt_bad_1', 3, mask=64)	! #6		! bit  3: Quality dual-frequency iono
	call nc2f (ncid, 'geo_bad_1', 4, mask=2) ! #1			! bit  4: Water/land
	call nc2f (ncid, 'geo_bad_1', 5, mask=2) ! #1			! bit  5: Ocean/other (preliminary, actually water/land!)
	call nc2f (ncid, 'geo_bad_1', 6, mask=4) ! #2			! bit  6: Radiometer land flag
	call nc2f (ncid, 'geo_bad_1', 7, mask=8) ! #3			! bit  7: Altimeter rain or ice flag
	call nc2f (ncid, 'geo_bad_2', 8, mask=1) ! #0			! bit  8: Radiometer rain or ice flag
	call nc2f (ncid, 'tmr_bad', 9, mask=1) ! #0			! bit 9&10: 0 = TMR data is OK
	call nc2f (ncid, 'tmr_bad', 10, mask=2) ! #1			!           1 = Within 25 km from land; or
													!               where >= 2 points are missing in along-track averaging; or
													!               when interpolation of TBs is from measurements > 4 sec apart.
													!           2 = Interpolation is actually extrapolation.
													!           3 = Over land; if center point from averaging is missing;
													!               when interpolation failed.
	call nc2f (ncid, 'alt_bad_1', 11, mask=159) ! #0-4,11	! bit 11: Quality range
	call nc2f (ncid, 'alt_bad_1', 12, mask=16) ! #4		! bit 12: Quality SWH
	call nc2f (ncid, 'alt_bad_2', 13, mask=64) ! #6		! bit 13: Quality Sigma0
	call nc2f (ncid, 'alt_bad_1', 14, mask=8) ! #3
	call nc2f (ncid, 'instr_state_topex', 14, mask=2)! #1	! bit 14: Tracking mode

! Convert all the necessary fields to RADS
! Time

	call cpy_var (ncid, 'tim_moy_1 9862 SUB 86400 MUL tim_moy_2 ADD tim_moy_3 ADD', 'time')

! Lat, lon, alt

	call cpy_var (ncid, 'lat')
	call cpy_var (ncid, 'lon')
	call cpy_var (ncid, 'sat_alt_1', 'alt_cnes')
	call cpy_var (ncid, 'sat_alt_2', 'alt_gdrd')
	call get_var (ncid, 'sat_alt_hi_rate', d)
	call get_var (ncid, 'dtim_pac', a)
	call new_var ('alt_rate', (d(10,:) - d(1,:)) / a / 18d0) ! alt_rate = (orbit diff over 9*10Hz) / (time diff 20Hz) / 18
	call cpy_var (ncid, 'off_nadir_angle_wf_ku')

! Altimeter range (on-board tracker)

	call get_var (ncid, 'alt_bad_1', a)
	call get_var (ncid, 'iono_corr_alt_ku', b)
	! Set iono_alt to NaN when iono_corr_alt_ku is zero AND bit 6 of alt_bad_1 is set
	do i = 1,nrec
		if (b(i) == 0 .and. btest(nint(a(i)),6)) b(i) = nan
	enddo
	call new_var ('iono_alt', b)
	call cpy_var (ncid, 'range_ku cog_corr ADD', 'range_ku')	! Explicitly add COG correction
	call new_var ('range_c', a - b / 0.17984d0) ! Compute range_c from range_ku and iono_corr
	nvar_range_c = nvar
	call cpy_var (ncid, 'range_rms_ku')
	call cpy_var (ncid, 'cog_corr net_instr_corr_range_ku SUB', 'drange_ku')
	call cpy_var (ncid, 'cog_corr net_instr_corr_range_c SUB', 'drange_c')
	call cpy_var (ncid, 'cog_corr', 'drange_cg')
	call cpy_var (ncid, 'range_numval', 'range_numval_ku')

! Atmospheric corrections

	call cpy_var (ncid, 'model_dry_tropo_corr', 'dry_tropo_ecmwf')
	call cpy_var (ncid, 'inv_bar_corr', 'inv_bar_static')
	call cpy_var (ncid, 'iono_corr_dor_ku', 'iono_doris')

	if (equator_time < 407548800d0) then
		! Correct model wet tropo as recommended in ERS-2 Validation Report 27,
		! for all ECMWF wet tropos prior to 1 December 1997.
		! This holds also for ERS-1, TOPEX and Poseidon!
		call cpy_var (ncid, 'model_wet_tropo_corr 0.85 MUL 0.006 SUB', 'wet_tropo_ecmwf')
	else if (equator_time > 538261200d0) then
		! Data after 22-Jan-2002 are affected by a modification in the ECMWF modelling at that time.
		! The newer grids are underestimating water vapour by about 4.6%
		! The correction is:  new_wet_corr = -0.4 [mm] + 1.0442 * wet_corr
		call cpy_var (ncid, 'model_wet_tropo_corr 1.0442 MUL 0.0004 SUB', 'wet_tropo_ecmwf')
	else
		call cpy_var (ncid, 'model_wet_tropo_corr', 'wet_tropo_ecmwf')
	endif
	call new_var ('wet_tropo_ecmwf', a)

! SWH and sigma0 (on-board tracker)

	call cpy_var (ncid, 'swh_ku')	! scale_factor fixed by ncatted
	call cpy_var (ncid, 'swh_c') ! scale_factor fixed by ncatted
	call cpy_var (ncid, 'swh_rms_ku')
	call cpy_var (ncid, 'swh_rms_c')
!	call cpy_var (ncid, 'swh_numval_ku')
	call cpy_var (ncid, 'net_instr_corr_swh_ku', 'dswh_ku')
	call cpy_var (ncid, 'net_instr_corr_swh_c', 'dswh_c')
	call cpy_var (ncid, 'dr_swh_att_ku', 'drange_swh_att_ku')
	call cpy_var (ncid, 'dr_swh_att_c', 'drange_swh_att_c')
	call cpy_var (ncid, 'sea_state_bias_ku', 'ssb_bm3')
	call cpy_var (ncid, 'sea_state_bias_ku_walsh', 'ssb_walsh')
	call cpy_var (ncid, 'sig0_ku')
	call cpy_var (ncid, 'sig0_c')
!	call cpy_var (ncid, 'agc_numval_ku', 'sig0_numval_ku')

! Geophysical quantities

	call cpy_var (ncid, 'mss', 'mss_cnescls11')
	call cpy_var (ncid, 'geoid', 'geoid_egm96')
	call cpy_var (ncid, 'ocean_tide_sol2 load_tide_sol2 SUB', 'tide_ocean_got410')
	call cpy_var (ncid, 'load_tide_sol2', 'tide_load_got410')
!	call cpy_var (ncid, 'ocean_tide_non_equil', 'tide_non_equil')	! Always NaN
	call cpy_var (ncid, 'solid_earth_tide', 'tide_solid')
	call cpy_var (ncid, 'pole_tide', 'tide_pole')
	call cpy_var (ncid, 'bathymetry', 'topo_dtm2000')
!	call cpy_var (ncid, 'inv_bar_corr hf_fluctuations_corr ADD', 'inv_bar_mog2d')	! Always NaN
	call cpy_var (ncid, 'wind_speed_alt')

! Retracker

	call cpy_var (ncid, 'range_ku net_instr_corr_range_ku SUB h_retrk1_ku ADD cog_corr ADD', 'range_ku_retrk')
	nvar_range_ku_retrk = nvar
	call cpy_var (ncid, 'h_retrk1_ku_rms', 'range_rms_ku_retrk')	! scale_factor fixed by ncatted
	call cpy_var (ncid, 'att_retrk1_ku', 'off_nadir_angle2_wf_ku_retrk')
	call cpy_var (ncid, 'swh_retrk1_ku', 'swh_ku_retrk')
	call cpy_var (ncid, 'numval_retrk1_ku', 'range_numval_ku_retrk')
	call get_var (ncid, 'h_retrk1_c net_instr_corr_range_c SUB cog_corr ADD', a)
	call new_var ('range_c_retrk', var(nvar_range_c)%d(1:nrec) + a)
	nvar_range_c_retrk = nvar
	call cpy_var (ncid, 'h_retrk1_c_rms', 'range_rms_c_retrk')
	call cpy_var (ncid, 'swh_retrk1_c', 'swh_c_retrk')
	call cpy_var (ncid, 'att_retrk1_c', 'off_nadir_angle2_wf_c_retrk')
	call cpy_var (ncid, 'numval_retrk1_c', 'range_numval_c_retrk')
	call cpy_var (ncid, 'net_instr_corr_retrk_ku NEG cog_corr ADD', 'drange_ku_retrk')	! scale_factor fixed by ncatted, flipped sign
	call cpy_var (ncid, 'net_instr_corr_retrk_c NEG cog_corr ADD', 'drange_c_retrk')	! scale_factor fixed by ncatted, flipped sign
!	call cpy_var (ncid, 'iono_retrk', 'iono_alt_retrk')
	call cpy_var (ncid, 'net_instr_corr_range_c net_instr_corr_range_ku SUB h_retrk1_ku ADD h_retrk1_c SUB 0.17984 MUL ' // &
		'iono_corr_alt_ku ADD', 'iono_alt_retrk')

! TMR replacement product

	call cpy_var (ncid, 'rad_wet_tropo_corr_tmrcp', 'wet_tropo_rad')	! scale_factor fixed by ncatted
	call cpy_var (ncid, 'tb_18_corr', 'tb_180')
	call cpy_var (ncid, 'tb_21_corr', 'tb_210')
	call cpy_var (ncid, 'tb_37_corr', 'tb_370')
	call cpy_var (ncid, 'atmos_sig0_corr_ku_tmrcp', 'dsig0_atmos_ku')
	call cpy_var (ncid, 'atmos_sig0_corr_c_tmrcp', 'dsig0_atmos_c')
	call cpy_var (ncid, 'rad_water_vapor', 'water_vapor_rad')
	call cpy_var (ncid, 'rad_liquid_water', 'liquid_water_rad')

! Derived from retracking and TMR update

	call cpy_var (ncid, 'wind_speed_retrk', 'wind_speed_alt_retrk')
	call cpy_var (ncid, 'sea_state_bias_retrk_ku', 'ssb_retrk')

! Flags

	a = flags
	call new_var ('flags', a)

! Dump the data

	call nfs(nf90_close(ncid))
	call put_rads
	deallocate (a, b, d, flags)

enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Write TOPEX Retracked GDR data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_RGDR_file_names')
write (*,1310)
1310 format (/ &
'This program converts TOPEX Retracked GDR files to RADS data' / &
'files with the name $RADSDATAROOT/data/tx.r50/a/pPPPP/txpPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

end program rads_gen_tx_rgdr
