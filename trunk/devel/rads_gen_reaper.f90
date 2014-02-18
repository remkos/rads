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
! written into files $RADSDATAROOT/data/eE/F.r/eEpPPPPcCCC.nc.
!     E = 1 or 2
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
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_rad - Radiometer wet tropo correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
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
logical :: alt_2m
real(eightbytereal) :: tnode(2), lnode(2)
integer(fourbyteint) :: orbitnr(2), cyclenr(2), passnr(2), varid, com=999

! Data variables

integer(fourbyteint), parameter :: mrec=15000, mvar=50
integer(fourbyteint) :: nvar, ndata=0, nrec=0, nout=0, ncid, ers=0
real(eightbytereal) :: start_time, end_time
real(eightbytereal), allocatable :: a(:), b(:), c(:), d(:,:), dh(:), sum_c_applied(:), sum_d_applied(:)
integer(twobyteint), allocatable :: flags(:)
integer(fourbyteint), allocatable :: f_error(:), f_applied(:)
logical, allocatable :: valid(:,:)
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
real(eightbytereal), parameter :: picosec_to_m=0.5d-12*299792458d0	! picoseconds of 2-way range to mm 1-way
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
	! If the start time is within the last file, get another file first
	if (ios /= 0) then
		if (ndata == 0) exit ! No more data left in memory and no more new files
	else if (ndata == 0 .or. var(1)%d(1) >= start_time) then
		! Read the next file
		read (*,550,iostat=ios) infile
		if (ios == 0) call get_reaper
		if (ndata == 0) cycle ! If still no data, try again
	endif

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
'files with the name $RADSDATAROOT/data/eE/F.r/pPPPP/eEpPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Store the content of a REAPER file into memory
!-----------------------------------------------------------------------

subroutine get_reaper
real(eightbytereal) :: dhellips, t(3), last_time
integer(fourbyteint) :: i, k, flag, ivar0, ivar1
character(len=80) :: string

550 format (a)
551 format (a,' ...')
552 format (i5,' records ...')
553 format (a,i5,3f18.3)

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
i = index(filenm, '_COM')
if (i > 0) read (filenm(i+4:i+4), *) com

! Open input file

if (nf90_open(infile,nf90_nowrite,ncid) /= nf90_noerr) then
	write (*,550) 'Error opening file'
	return
endif

! Check for ERS-1 or -2
! Do not trust 'mission' attribute. It is always 'E1'.

if (filenm(:2) == 'E1') then
	if (ers == 0) call rads_init (S, 'e1/a.r', verbose)
	ers = 1
else if (filenm(:2) == 'E2') then
	if (ers == 0) call rads_init (S, 'e2/a.r', verbose)
	ers = 2
else
	write (*,550) 'Error: Unknown file type'
	return
endif

! Read header records

call nfs(nf90_inq_dimid(ncid,'Record',varid))
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

! Allocate arrays

allocate (a(nrec),b(nrec),c(nrec),d(20,nrec),valid(20,nrec),f_error(nrec),f_applied(nrec), &
	flags(nrec),dh(nrec),sum_c_applied(nrec),sum_d_applied(nrec))
flags = 0
sum_d_applied = 0d0
nvar = 0

! Time and orbit: Low rate

call get_var (ncid, 'time_day_1hz', a)
call get_var (ncid, 'time_milsec_1hz', b)
call get_var (ncid, 'time_micsec_1hz', c)
a = a * 86400d0 + b * 1d-3 + c * 1d-6 + sec1990
k = min(3,nrec)
! Because first and last time can be wrong, we use the minimum of the first 3 and
! last 3 measurements as boundaries.
start_time = minval(a(1:k))
end_time = maxval(a(nrec-k+1:nrec))
! Discard measurements at the end of the stack that are newer than the beginning of this file
do while (ndata > 0 .and. var(1)%d(ndata) > start_time - 0.5d0)
	ndata = ndata - 1
enddo
call new_var ('time', a)
call get_var (ncid, 'latitude_1hz', a)
a = a*1d-6
call new_var ('lat', a)
! Compute ellipsoid corrections
do i = 1,nrec
	dh(i) = dhellips(1,a(i))
enddo
call get_var (ncid, 'longitude_1hz', a)
call new_var ('lon', a*1d-6)
call get_var (ncid, 'altitude_1hz', a)
call new_var ('alt_reaper', a*1d-3+dh)
call get_var (ncid, 'altitude_rate_1hz', a)
call new_var ('alt_rate', a*1d-3)
call get_var (ncid, 'wf_attitude_1hz', a)
call new_var ('off_nadir_angle2_wf_ku', a*a*1d-4)

! Range data: Low rate

call get_var (ncid, 'ocean_range_1hz', a)
call get_var (ncid, 'ocean_stdev_1hz', b)
call get_var (ncid, 'ocean_valid_num_1hz', c)
call invalidate (c == 0, a)
call invalidate (c <= 1, b)
call new_var ('range_ku', a*1d-3)
call new_var ('range_rms_ku', b*1d-3)
call new_var ('range_numval_ku', c)
if (alt_2m) then ! Different name in Meteo product
	call get_var (ncid, 'ocean_valid_bitmap_1hz', a)
else
	call get_var (ncid, 'f_ocean_valid_bitmap_1hz', a)
endif

! Get USO frequency and convert to a range correction assuming WAP's
! nominal USO frequency is 15000000.05 Hz
call nfs(nf90_get_att(ncid,nf90_global,'uso_applied_1', string))
i = index(string, '@')
read (string(:i-1),*) a(1) ! Get frequency from headers
a = -1d-4 * nint (795d7 * (a(1) - uso_wap) / uso_wap) ! Convert to meters and truncate at 4 decimals
call new_var ('drange_uso', a, ndims=0)

! Set "valid" array based on bitmap
! Bits are 0 when valid, 1 when invalid
! However, when all are invalid, value is 0
do i = 1,nrec
	flag = nint(a(i))
	if (c(i) == 0) then
		valid(:,i) = .false.
	else if (c(i) == 20) then
		valid(:,i) = .true.
	else
		do k = 1,20
			valid (k,i) = .not.btest(flag,k-1)
		enddo
	endif
enddo

! For some reason Meteo product has 1-Hz peakiness
! where (S)GDR have ocean_wind_1hz
if (alt_2m) then
	call get_var (ncid, 'wf_pk_1hz', a)
	call invalidate (c == 0, a)
	call new_var ('peakiness_ku', a*1d-3)
else
	call get_var (ncid, 'ocean_wind_1hz', a)
	call invalidate (c == 0, a)
	call new_var ('wind_speed_alt', a*1d-3)
endif

! Range data: High rate

if (.not.alt_2m) then
	call get_var (ncid, 'ocean_mean_quadratic_error', d)
	call mean_1hz (d, valid, a, b)
	call new_var ('mqe', a*1d-4)
	call get_var (ncid, 'wf_pk', d)
	call mean_1hz (d, valid, a, b)
	call new_var ('peakiness_ku', a*1d-3)
	call get_var (ncid, 'dop_c+delta_dop_c', d)
	valid = .true.	! For Doppler, we use all
	call mean_1hz (d, valid, a, b)
	call new_var ('drange_fm', a*1d-3)
	call get_var (ncid, 'sptr_jumps_c', d)
	call new_var ('drange_cal', d(1,:) * picosec_to_m)
	call get_var (ncid, 'sum_c_applied', d)
	sum_c_applied = d(1,:)*1d-3
	do i = 1,nrec
		do k = 2,19
			if (d(k,i) /= d(1,i)) write (*,*) 'Error: sum_c_applied changed:',i,k,d(1,i),d(k,i)
		enddo
	enddo
endif

! Get error and and correction flags first

call get_var (ncid, 'f_corr_error_1hz', a)
f_error = nint(a)
call get_var (ncid, 'f_corr_applied_1hz', a)
f_applied = nint(a)

! Sigma zero: Low rate

call get_var (ncid, 'ocean_sig0_1hz', a)
call get_var (ncid, 'ocean_sig0_stdev_1hz', b)
call get_var (ncid, 'ocean_sig0_valid_num_1hz', c)
call invalidate (c == 0, a)
call invalidate (c <= 1, b)
call new_var ('sig0_ku', a*1d-2)
call new_var ('sig0_rms_ku', b*1d-2)

! SWH: Low rate

call get_var (ncid, 'swh_signed_1hz', a)
call get_var (ncid, 'swh_stdev_1hz', b)
call get_var (ncid, 'swh_valid_num_1hz', c)
call invalidate (c == 0, a)
call invalidate (c <= 1, b)
call new_var ('swh_ku', a*1d-3)
call new_var ('swh_rms_ku', b*1d-3)

! MWR Flags: Low rate

call get_var (ncid, 'f_sea_ice_flag_1hz', a)
call flag_set (a == 1, 8)

! MWR: Low rate

call get_var (ncid, 'tb_23_8_1hz', a)
call get_var (ncid, 'tb_36_5_1hz', b)
call get_var (ncid, 'f_mwr_srf_typ_1hz', c)
call flag_set (c == 1, 6)
call get_var (ncid, 'f_mwr_interp_qual_1hz', c)
call invalidate (c == 3, a)
call invalidate (c == 3, b)
if (alt_2m) then
	call get_var (ncid, 'f_mwr_valid_1hz', c)
else ! Wrong name in (S)GDR
	call get_var (ncid, 'f_MWR_valid_1hz', c)
endif
call invalidate (c == 1, a)
call invalidate (c == 1, b)
call new_var ('tb_238', a*1d-2)
call new_var ('tb_365', b*1d-2)

! Atmospheric and geophysical: Low rate

call get_var (ncid, 'dry_c_1hz', a)
call new_var ('dry_tropo_ecmwf', a*1d-3, 1)
ivar0 = nvar
call get_var (ncid, 'ib_c_1hz', a)
call new_var ('inv_bar_static', a*1d-3, 2)
call get_var (ncid, 'mog2d_c_1hz', a)
call new_var ('inv_bar_mog2d', a*1d-3, 3)
call get_var (ncid, 'wet_c_mod_1hz', a)
call new_var ('wet_tropo_ecmwf', a*1d-3, 5)
call get_var (ncid, 'wet_c_mwr_1hz', a)
call new_var ('wet_tropo_rad', a*1d-3, 6)
call get_var (ncid, 'water_vapor_content_1hz', a)
call new_var ('water_vapor_rad', a*1d-2, -6)
call get_var (ncid, 'liquid_water_content_1hz', a)
call new_var ('liquid_water_rad', a*1d-2, -6)
call get_var (ncid, 'u_wind_1hz', a)
call new_var ('wind_speed_ecmwf_u', a*1d-3, 7)
call get_var (ncid, 'v_wind_1hz', a)
call new_var ('wind_speed_ecmwf_v', a*1d-3, 8)
call get_var (ncid, 'iono_c_mod_1hz', a)
call new_var ('iono_nic09', a*1d-3, 9)
if (start_time >= 430880400d0) then	! After 1998-08-28 01:00:00 get GIM iono
	call get_var (ncid, 'iono_c_gps_1hz', a)
	call new_var ('iono_gim', a*1d-3, 10)
endif
call get_var (ncid, 'h_mss_cls01_1hz', a)
call new_var ('mss_cls01', a*1d-3+dh, 11)
call get_var (ncid, 'h_geo_1hz', a)
call new_var ('geoid_egm2008', a*1d-3+dh, 12)
! Prior to COM5 ocean tide is OT+OLT+LPT
! Since COM5 ocean tide is split up
! Nonetheless we ignore that here and simply copy the values
! This can be fixed with rads_fix_reaper --tide
call get_var (ncid, 'h_ot_1hz', a)
call new_var ('tide_ocean_got47', a*1d-3, 13)
call get_var (ncid, 'h_ot2_1hz', a)
call new_var ('tide_ocean_fes04', a*1d-3, 14)
call get_var (ncid, 'h_olt_1hz', a)
call new_var ('tide_load_got47', a*1d-3, 15)
call get_var (ncid, 'h_olt2_1hz', a)
call new_var ('tide_load_fes04', a*1d-3, 16)
call get_var (ncid, 'h_lpt_1hz', a)
call new_var ('tide_equil', a*1d-3, 17)
call get_var (ncid, 'h_lptne_1hz', a)
call new_var ('tide_non_equil', a*1d-3, 18)
call get_var (ncid, 'h_set_1hz', a)
call new_var ('tide_solid', a*1d-3, 19)
call get_var (ncid, 'h_pol_1hz', a)
call new_var ('tide_pole', a*1d-3, 20)

call get_var (ncid, 'h_odle_1hz', a)
call new_var ('topo_macess', a*1d-3, 22)
call get_var (ncid, 'em_bias_1hz', a)
call new_var ('ssb_bm3', a*1d-3, 23)
call get_var (ncid, 'h_mss_ucl04_1hz', a)
call new_var ('mss_ucl04', a*1d-3+dh, 27)
ivar1 = nvar

call get_var (ncid, 'f_srf_typ_1hz', a)
call flag_set (a >= 3, 2)
call flag_set (a >= 2, 4)
call flag_set (a >= 1, 5)
call new_var ('flags', flags*1d0)

if (.not.alt_2m) then ! Only on (S)GDR
	call get_var (ncid, 'sig0_attn_c_1hz', a)
	call new_var ('dsig0_atmos_ku', a*1d-2)
endif

! Prior to COM5: remove applied corrections from range

if (com < 5) then
	do i = 1,nrec
		k = ndata + i
		var(7)%d(k) = var(7)%d(k) - sum_d_applied(i)
		if (.not.alt_2m .and. abs(sum_d_applied(i) - sum_c_applied(i)) > 1d-4) &
			write (*,553) 'Warning: sum_c_applied wrong: ',i,sum_d_applied(i),sum_c_applied(i),sum_d_applied(i)-sum_c_applied(i)
	enddo
endif

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

! Another problem is that longitude at the dateline are incorrectly determined,
! hence corrections at those points are wrong too.

do i = 2,nrec-1
	t = var(3)%d(ndata+i-1:ndata+i+1)
	if (t(1) > -179d0 .or. t(3) < 179d0 .or. abs(abs(t(2))-179d0) < 1d0) cycle
	write (*,553) 'Warning: Fixed incorrect longitude at date line  :', i, t
	! Average the previous and next longitude properly
	t(2) = 0.5d0 * (t(1) + t(3))
	if (t(2) > 0d0) then
		var(3)%d(ndata+i) = t(2) - 180d0
	else
		var(3)%d(ndata+i) = t(2) + 180d0
	endif
	! Average the previous and next corrections
	var(ivar0:ivar1)%d(ndata+i) = &
		0.5d0 * (var(ivar0:ivar1)%d(ndata+i-1)+var(ivar0:ivar1)%d(ndata+i+1))
	! Copy the flag words one forward
	var(ivar1+1)%d(ndata+i) = var(ivar1+1)%d(ndata+i-1)
enddo

! Close this input file

deallocate (a,b,c,d,valid,f_error,f_applied,flags,dh,sum_c_applied,sum_d_applied)

call nfs(nf90_close(ncid))

ndata = ndata + nrec
end subroutine get_reaper

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
if (S%phase%name /= phasenm(1)) S%phase => rads_get_phase (S, phasenm(1)//'.r')

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
call rads_close_pass (S, P)

! Formats
550 format (a,1x)
552 format ('...',i5,' records written to ',a)

end subroutine put_rads

!-----------------------------------------------------------------------
! Create new RADS variable
!-----------------------------------------------------------------------

subroutine new_var (varnm, data, bit, ndims)
! Store variables to be written later by put_var
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(in) :: data(:)
integer, optional, intent(in) :: bit
integer, optional, intent(in) :: ndims
integer :: i
nvar = nvar + 1
if (nvar > mvar) stop 'Too many variables'
var(nvar)%v => rads_varptr (S, varnm)
var(nvar)%d(ndata+1:ndata+nrec) = data
if (present(ndims)) var(nvar)%v%info%ndims = ndims
if (.not.present(bit)) return
if (bit > 0) then
	do i = 1,nrec
		if (btest(f_applied(i),bit)) sum_d_applied(i) = sum_d_applied(i) + data(i)
		if (btest(f_error(i),bit)) var(nvar)%d(ndata+i) = nan
	enddo
else
	do i = 1,nrec
		if (btest(f_error(i),-bit)) var(nvar)%d(ndata+i) = nan
	enddo
endif
end subroutine new_var

!-----------------------------------------------------------------------
! Set a bit in an array of flags
!-----------------------------------------------------------------------

subroutine flag_set (a, bit)
logical, intent(in) :: a(:)
integer(fourbyteint), intent(in) :: bit
integer(fourbyteint) :: i
integer(twobyteint) :: j
if (size(a) /= size(flags)) stop "Error in flag_set"
j = int(bit,twobyteint)
do i = 1,size(a)
	if (a(i)) flags(i) = ibset(flags(i),j)
enddo
end subroutine flag_set

!-----------------------------------------------------------------------
! Set to NaN elements of an array
!-----------------------------------------------------------------------

subroutine invalidate (a, b)
logical, intent(in) :: a(:)
real(eightbytereal), intent(inout) :: b(:)
if (size(a) /= size(b)) stop "Error in invalidate"
where (a) b = nan
end subroutine invalidate

end program rads_gen_reaper
