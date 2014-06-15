!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2014  Remko Scharroo (Altimetrics LLC)
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
! This program handles only the REAPER ERS_ALT_2 files in netCDF format.
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
! tide_ocean_fes04 - FES2008 ocean tide
! tide_ocean_got47 - GOT4.7 ocean tide
! tide_load_fes04 - FES2008 load tide
! tide_load_got47 - GOT4.7 load tide
! tide_pole - Pole tide
! ssb_bm3 - SSB
! mss_cls01 - CLS01 MSS
! geoid_egm2008 - EGM2008 geoid
! mss_ucl04 - UCL04 MSS
! swh_ku - Significant wave height
! sig0_ku - Sigma0
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! range_rms_ku - Std dev of range
! range_numval_ku - Nr of averaged range measurements
! topo_macess - MACESS topography
! tb_238 - Brightness temperature (23.8 GHz)
! tb_365 - Brightness temperature (36.5 GHz)
! peakiness_ku - Peakiness
! flags - Engineering flags
! swh_rms_ku - Std dev of SWH
! sig0_rms_ku - Std dev of sigma0
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! liquid_water - Liquid water content
! water_vapor_content - Water vapor content
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
use netcdf
use rads
use rads_misc
use rads_time
use rads_netcdf
use rads_devel

! Command line arguments

integer(fourbyteint) :: verbose=0, c0=0, c1=999, ios
real(eightbytereal) :: t0, t1
character(len=rads_cmdl) :: infile, filenm, old_filenm = ''
character(len=rads_varl) :: optopt, optarg

! Header variables

character(len=1) :: phasenm(2)
character(len=rads_varl) :: l2_proc_time, l2_version
character(len=4) :: mle
logical :: alt_2m
real(eightbytereal) :: tnode(2), lnode(2)
integer(fourbyteint) :: orbitnr(2), cyclenr(2), passnr(2), varid

! Data variables

integer(fourbyteint), parameter :: mrec=20000, mvar=50
integer(fourbyteint) :: nvar, ndata=0, nrec=0, nout=0, ncid, ers=0
real(eightbytereal) :: start_time, end_time, last_time = 0d0
real(eightbytereal), allocatable :: tmp(:)
type(rads_sat) :: S
type(rads_pass) :: P
type :: var_
	type(rads_var), pointer :: v ! Pointer to rads_var struct
	real(eightbytereal) :: d(mrec) ! Data array
	logical :: empty ! .true. if all NaN
endtype
type(var_) :: var(mvar)

! Other local variables

real(eightbytereal), parameter :: sec1990=157766400d0	! UTC seconds from 1 Jan 1985 to 1 Jan 1990
real(eightbytereal), parameter :: sec_to_m=0.5d0*299792458d0	! seconds of 2-way range to meters 1-way
real(eightbytereal), parameter :: uso_wap=15000000.05d0	! Nominal USO frequency used in WAP products
integer :: i
logical :: new

! Initialise

t0 = nan
t1 = nan
550 format (a)

! Scan command line for options

do
	call getopt ('vC: debug: cycle: t: mjd: sec: ymd: doy:', optopt, optarg)
	select case (optopt)
	case ('!')
		exit
	case ('v')
		verbose = 1
	case ('debug')
		read (optarg,*) verbose
	case ('C', 'cycle')
		c1 = -1
		read (optarg,*,iostat=ios) c0,c1
		if (c1 < c0) c1 = c0
	case default
		if (.not.dateopt (optopt, optarg, t0, t1)) then
			call synopsis ('--help')
			stop
		endif
	end select
enddo

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

! Start reading at least first file

read (*,550,iostat=ios) infile
if (ios /= 0) then
	call synopsis ('--help')
else
	call synopsis ('--head')
endif
call get_reaper

do
	! Read the next file as long as buffer is empty or less than one orbit in memory

	do while (ndata == 0 .or. var(1)%d(ndata) - var(1)%d(1) < 6100d0)
		read (*,550,iostat=ios) infile
		if (ios /= 0) exit
		call get_reaper
	enddo

	! Look where to split this chunk of data
	new = erspass (ers, var(1)%d(1), orbitnr(1), phasenm(1), cyclenr(1), passnr(1), tnode(1), lnode(1))
	do i = 2,ndata
		if (erspass (ers, var(1)%d(i), orbitnr(2), phasenm(2), cyclenr(2), passnr(2), tnode(2), lnode(2))) exit
	enddo
	! It is OK to exit this loop with i = ndata + 1. This means we dump all of the memory.

	! Write out the data
	nout = i - 1 ! Number of measurements to be written out
	call put_rads

	! Number of measurements remaining
	ndata = ndata - nout
	if (ios /= 0 .and. ndata == 0) exit ! We are out of data

	! Move the data to be beginning
	do i = 1,nvar
		var(i)%d(1:ndata) = var(i)%d(nout+1:nout+ndata)
	enddo
enddo

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Write REAPER data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_REAPER_file_names')
write (*,1310)
1310 format (/ &
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
character(len=2) :: mission

550 format (a)
551 format (a,' ...')
552 format (i5,' records ...')
553 format (a,i6,3f18.3)

! Reduce file name to basename only

i = index(infile,'/',.true.) + 1
old_filenm = filenm
filenm = infile(i:)

! Check input file name

write (*,551,advance='no') trim(infile)
i = index(filenm, '_ERS_ALT_2')
if (i <= 0) then
	write (*,550) 'Error: Wrong input file'
	return
endif
alt_2m = (filenm(i+10:i+10) == 'M')
mission = filenm(i-7:i-6)

i = index(infile,'.nc')
if (mission == 'E1') then
	if (ers == 0) call rads_init (S, 'e1.' // strtolower(infile(i-4:i-1)), (/'reaper'/), verbose)
	ers = 1
else if (mission == 'E2') then
	if (ers == 0) call rads_init (S, 'e2.' // strtolower(infile(i-4:i-1)), (/'reaper'/), verbose)
	ers = 2
else
	write (*,550) 'Error: Unknown file type: '//mission
	return
endif

! Open input file

i = nf90_open(infile,nf90_nowrite,ncid)
write (*,*) 'ncid =',ncid
if (i /= nf90_noerr) then
	write (*,550) nf90_strerror (i)
	stop
	return
endif

! Read header records

call nfs(nf90_inq_dimid(ncid,'time',varid))
call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
write (*,552) nrec
if (nrec == 0) then	! Skip empty input files
	call nfs(nf90_close(ncid))
	return
else if (ndata+nrec > mrec) then
	write (*,553) 'Error: Too many input measurements: ', ndata+nrec
	stop
endif
call nfs(nf90_get_att(ncid,nf90_global,'l2_proc_time',l2_proc_time))
call nfs(nf90_get_att(ncid,nf90_global,'l2_software_ver',l2_version))
call nfs(nf90_get_att(ncid,nf90_global,'ocean_retracker_version_for_ocean',mle))
!call nfs(nf90_get_att(ncid,nf90_global,'mission',mission)) ! Actually still is not correct -> bypass

call get_reaper_data (nrec)
call nfs(nf90_close(ncid))
end subroutine get_reaper

subroutine get_reaper_data (nrec)
integer(fourbyteint), intent(inout) :: nrec
real(eightbytereal) :: a(nrec), b(nrec), c(nrec), d(20,nrec), dh(nrec)
integer(twobyteint) :: flags(nrec)
logical :: valid(20,nrec)
integer(fourbyteint) :: k
real(eightbytereal) :: dhellips, t(3)

553 format (a,i6,3f18.3)

! Time and orbit: Low rate

call get_var (ncid, 'time', a)
a = a + sec1990
k = min(3,nrec)

! Because first and last time can be wrong, we use the minimum of the first 3 and
! last 3 measurements as boundaries.
t(1) = last_time
t(2) = minval(a(1:k))
t(3) = maxval(a(nrec-k+1:nrec))

! There are significant overlaps between files
! Here we remove all new files that fall entirely before the end of the previous file, and
if (t(3) < t(1) + 1) then
	write (*,553) 'Warning: Removed file because of time reversal   :', nrec, t
	return
endif
start_time = t(2)
end_time = t(3)

! Discard measurements at the end of the stack that are newer than the beginning of the
! new file
k = ndata
do while (k > 0 .and. var(1)%d(k) > start_time - 0.5d0)
	k = k - 1
enddo
if (k < ndata) then
	write (*,553) 'Warning: Removed at end of buffer, time reversal :', ndata-k, &
		var(1)%d(k+1), var(1)%d(ndata), start_time
	ndata = k
endif

! Initialize

flags = 0
nvar = 0
allocate (tmp(nrec))

! Time and orbit: Low rate (cont'd)

call new_var ('time', a)
call cpy_var ('lat', 'lat')
! Compute ellipsoid corrections
do i = 1,nrec
	dh(i) = dhellips(1,tmp(i))
enddo
call cpy_var ('lon', 'lon')
call get_var (ncid, 'alt', a)
call new_var ('alt_reaper', a+dh)
call cpy_var ('orb_alt_rate', 'alt_rate')

! Range data: Low rate

call get_var (ncid, 'ocean_range_numval', c)
call cpy_var ('ocean_range', 'range_ku')
call cpy_var ('ocean_range_rms', 'range_rms_ku', c <= 1)
call cpy_var ('ocean_range_numval', 'range_numval_ku')

! Retracking info: Low rate

if (mle == 'MLE4') call cpy_var ('off_nadir_angle_wf', 'off_nadir_angle2_wf_ku')
call cpy_var ('peakiness', 'peakiness_ku')

! Retracking info: High rate

if (.not.alt_2m) then
	call get_var (ncid, 'ocean_qual_flag_20hz', d)
	valid = (d == 0)
	call get_var (ncid, 'ocean_mqe_20hz', d)
	where (.not.valid) d = nan
	call mean_1hz (d, a, b)
	call new_var ('mqe', a)
	call get_var (ncid, 'doppler_corr+delta_doppler_corr_20hz', d)
	call mean_1hz (d, a, b)
	call new_var ('drange_fm', a)
	call get_var (ncid, 'sptr_jumps_corr_20hz', d)
	call new_var ('drange_sptr', d(1,:) * sec_to_m)
endif

! Sigma zero: Low rate

call cpy_var ('ocean_sig0', 'sig0_ku')
call get_var (ncid, 'ocean_sig0_numval', c)
call cpy_var ('ocean_sig0_rms', 'sig0_rms_ku', c <= 1)
call cpy_var ('wind_speed_alt', 'wind_speed_alt', c == 0)

! SWH: Low rate

call cpy_var ('swh_signed', 'swh_ku')
call get_var (ncid, 'swh_numval', c)
call cpy_var ('swh_rms', 'swh_rms_ku', c <= 1)

! MWR: Low rate

call get_var (ncid, 'rad_state_flag_orb_prop', a)
call get_var (ncid, 'rad_state_flag_orb_init', b)
call get_var (ncid, 'rad_state_flag_l2_proc_error+rad_state_validity+rad_state_flag+rad_state_bt_check', c)
call cpy_var ('tb_238', 'tb_238', a == 2d0 .or. b == 2d0 .or. c > 0d0)
call cpy_var ('tb_365', 'tb_365', a == 2d0 .or. b == 2d0 .or. c > 0d0)

! Atmospheric and geophysical: Low rate

call cpy_var ('model_dry_tropo_corr', 'dry_tropo_era')
call cpy_var ('inv_bar_corr', 'inv_bar_static')
call cpy_var ('hf_fluctuations_corr', 'inv_bar_mog2d')	! This is in contrast with standard_name, but supported by comment
call cpy_var ('model_wet_tropo_corr', 'wet_tropo_era')
call cpy_var ('rad_wet_tropo_corr', 'wet_tropo_rad')
call cpy_var ('rad_water_vapor', 'water_vapor_rad')
call cpy_var ('rad_liquid_water', 'liquid_water_rad')
call cpy_var ('wind_speed_model_u', 'wind_speed_ecmwf_u')
call cpy_var ('wind_speed_model_v', 'wind_speed_ecmwf_v')
call cpy_var ('iono_corr_model', 'iono_nic09')
if (start_time >= 430880400d0) then	! After 1998-08-28 01:00:00 get GIM iono
	call cpy_var ('iono_corr_gps', 'iono_gim')
endif
call get_var (ncid, 'mean_sea_surface_1', a)
call new_var ('mss_cls01', a+dh)
call get_var (ncid, 'mean_sea_surface_2', a)
call new_var ('mss_ucl04', a+dh)
call get_var (ncid, 'geoid', a)
call new_var ('geoid_egm2008', a+dh)
! Need to recombine to OT+LPT
call cpy_var ('ocean_tide_sol1+ocean_tide_equil', 'tide_ocean_got47')
call cpy_var ('load_tide_sol1', 'tide_load_got47')
call cpy_var ('ocean_tide_sol2+ocean_tide_equil', 'tide_ocean_fes04')
call cpy_var ('load_tide_sol2', 'tide_load_fes04')
call cpy_var ('ocean_tide_equil', 'tide_equil')
call cpy_var ('ocean_tide_non_equil', 'tide_non_equil')
call cpy_var ('solid_earth_tide', 'tide_solid')
call cpy_var ('pole_tide', 'tide_pole')

call cpy_var ('bathymetry', 'topo_macess')
call cpy_var ('sea_state_bias', 'ssb_hyb')

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

! There may be measurements with invalid times.
! If so, weed them out.

k = 0
valid(1,:) = .true.
last_time = var(1)%d(ndata+1)
do i = 2,nrec
	t = var(1)%d(ndata+i-1:ndata+i+1)
	if (t(2) < start_time .or. t(2) > end_time) then
		write (*,553) 'Warning: Removed measurement out of time range   :', i, t
	else if (i < nrec .and. t(2) > max(t(1),t(3))+1d0) then
		write (*,553) 'Warning: Removed measurement out of time sequence:', i, t
	else if (t(2) < last_time+0.5d0) then
		write (*,553) 'Warning: Removed measurement with time reversal  :', i, t
	else
		last_time = t(2)
		cycle
	endif
	valid(1,i) = .false.
	k = k + 1
enddo
if (k > 0) then
	do i = 1,nvar
		var(i)%d(ndata+1:ndata+nrec-k) = pack(var(i)%d(ndata+1:ndata+nrec),valid(1,:))
	enddo
	nrec = nrec - k
endif

ndata = ndata + nrec

deallocate (tmp)

end subroutine get_reaper_data

!-----------------------------------------------------------------------
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads
integer :: i
character(len=rads_cmdl) :: original

if (nout == 0) return	! Skip empty data sets
if (cyclenr(1) < c0 .or. cyclenr(1) > c1) return	! Skip chunks that are not of the selected cycle
if (tnode(1) < t0 .or. tnode(1) > t1) return	! Skip equator times that are not of selected range

! Update phase name if required
phasenm(1) = strtolower(phasenm(1))
if (S%phase%name /= phasenm(1)) S%phase => rads_get_phase (S, phasenm(1))

! Store relevant info
call rads_init_pass_struct (S, P)
P%cycle = cyclenr(1)
P%pass = passnr(1)
P%start_time = var(1)%d(1)
P%end_time = var(1)%d(nout)
P%equator_time = tnode(1)
P%equator_lon = lnode(1)

! Check which input files pertain
if (P%start_time >= start_time) then
	original = filenm
else if (P%end_time < start_time) then
	original = old_filenm
else
	original = trim(old_filenm)//rads_linefeed//filenm
endif
P%original = trim(l2_version)//' data of '//l2_proc_time(:11)//rads_linefeed//trim(original)

! Check which variables are empty
do i = 1,nvar
	var(i)%empty = all(isnan_(var(i)%d(1:nout)))
enddo
if (any(var(1:nvar)%empty)) then
	write (*,550,advance='no') '... No'
	do i = 1,nvar
		if (var(i)%empty) write (*,550,advance='no') trim(var(i)%v%name)
	enddo
endif

! Open output file
call rads_create_pass (S, P, nout)

! Define all variables
do i = 1,nvar
	call rads_def_var (S, P, var(i)%v)
enddo

! Fill all the data fields
do i = 1,nvar
	call rads_put_var (S, P, var(i)%v, var(i)%d(1:nout))
enddo

! Close the data file
write (*,552) nout,trim(P%filename(len_trim(S%dataroot)+2:))
write (*,*) 'ncout =',P%ncid
call rads_close_pass (S, P)

! Formats
550 format (a,1x)
552 format ('...',i5,' records written to ',a)

end subroutine put_rads

!-----------------------------------------------------------------------
! Copy variable to RADS
!-----------------------------------------------------------------------

subroutine cpy_var (varin, varout, invalid)
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
var(nvar)%d(ndata+1:ndata+nrec) = data
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
if (size(a) /= size(flags)) stop "Error in flag_set"
j = int(bit,twobyteint)
do i = 1,size(a)
	if (a(i)) flags(i) = ibset(flags(i),j)
enddo
end subroutine flag_set

end program rads_gen_reaper
