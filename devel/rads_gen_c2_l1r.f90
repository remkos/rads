!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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

!*rads_gen_c2_l1r -- Converts CryoSat Retracked Level 1 data to RADS
!
! This program reads CryoSat-2 L1R pass files and converts it to the RADS format,
! written into files $RADSDATAROOT/data/c2/F/c2pPPPPcCCC.nc.
!     F = mission phase
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_c2_l1r [options] < list_of_L1R_file_names
!
! This program handles only the CryoSat-2 L1Rs in netCDF format.
!-----------------------------------------------------------------------
!
!  Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdrd - Orbit altitude
! alt_rate - Orbit altitude rate
! range_ku - Ocean range (retracked)
! dry_tropo_ecmwf - Dry tropospheric correction
! wet_tropo_ecmwf - Wet tropo correction
! iono_bent - Bent ionospheric correction
! iono_gim - GIM ionosphetic correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! tide_solid - Solid earth tide
! tide_ocean_got00 - GOT00.1 ocean tide
! tide_load_got00 - GOT00.1 load tide
! tide_pole - Pole tide
! swh_ku - Significant wave height (retracked)
! sig0_ku - Sigma0 (retracked)
! agc_ku - AGC
! range_rms_ku - Std dev of range
! range_numval_ku - Nr of averaged range measurements
! peakiness_ku - Peakiness
! flags - Engineering flags
! drange_ku - Retracker range correction (applied)
! drange_cal - Internal calibration correction to range (applied)
! drange_fm - Doppler correction (applied)
! sig0_rms_ku - Std dev of sigma0
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! off_nadir_angle2_wf_rms_ku - Std dev of mispointing from waveform squared
! attitude_pitch - Platform pitch angle
! attitude_roll - Platform roll angle
! attitude_yaw - Platform yaw angle
! mqe - Mean quadratic error of waveform fit
! noise_floor_ku - Noise floor
! noise_floor_rms_ku - Std dev of noise floor
! flags_star_tracker - Star tracker flags
! tide_equil - Long-period tide
!-----------------------------------------------------------------------
program rads_gen_c2_l1r

use typesizes
use netcdf
use rads
use rads_netcdf
use rads_misc
use rads_time

! Command line arguments

integer(fourbyteint) :: verbose=0, c0=0, c1=999, ios
real(eightbytereal) :: t0, t1
character(160) :: arg
character(20) :: optopt, optarg
character(80), parameter :: optlist='vC: debug: sat: cycle: t: mjd: sec: ymd: doy:'


! Header variables

character(19) :: l1b_proc_time
character(14) :: l1b_version, l1r_version
character(55) :: l1r_product
real(eightbytereal) :: tai_utc, eq_time, eq_long
integer(fourbyteint) :: passnr(2),cycnr(2),recnr(2),nrec,ncid,varid

! Data variables

integer(fourbyteint), parameter :: mrec = 6000, mvar=50
integer(fourbyteint) :: nvar=0, ndata=0
real(eightbytereal), allocatable :: a(:),b(:),c(:),d(:,:),t_1hz(:),t_20hz(:,:),alt(:),dh(:)
logical, allocatable :: valid(:,:)
integer(twobyteint), allocatable :: flags(:)
type(rads_sat) :: S
type(rads_pass) :: P
type :: var_
	type(rads_var), pointer :: v ! Pointer to rads_var struct
	real(eightbytereal) :: d(mrec) ! Data array
	logical :: empty ! .true. if all NaN
endtype
type(var_) :: var(mvar)


! Other local variables

integer(fourbyteint), parameter :: maxint4=2147483647
real(eightbytereal), parameter :: sec2000=473299200d0, rev_time = 5953.45d0, rev_long = -24.858d0, &
	murad = 1d6*rad
real(eightbytereal), parameter :: pitch_bias = 0.096d0, roll_bias = 0.086d0, yaw_bias = 0d0	! Attitude biases to be added
real(eightbytereal) :: uso_corr, dhellips
integer(fourbyteint) :: i, j, m, oldcyc=0, oldpass=0, mle=3
logical :: version_a, sar

! Initialise

t0 = nan
t1 = nan
550 format (a)
551 format (a,' ...')

! Scan command line for options

call synopsis
do
	call getopt (optlist, optopt, optarg)
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

! Initialise

call synopsis ('--head')
call rads_init (S, 'c2', verbose)

!----------------------------------------------------------------------
! Read all file names for standard input
!----------------------------------------------------------------------

files: do
	read (*,550,iostat=ios) arg
	if (ios /= 0) exit files

! Open input file

	write (*,551) trim(arg)
	if (nf90_open(arg,nf90_nowrite,ncid) /= nf90_noerr) then
		write (*,550) 'Error opening file'
		cycle files
	endif

! Read cycle and pass number and determine if we need to dump the data

	call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cycnr))
	call nfs(nf90_get_att(ncid,nf90_global,'pass_number',passnr))
	if (passnr(1) /= oldpass .or. cycnr(1) /= oldcyc) then
		call put_rads (oldcyc, oldpass, ndata)
		ndata = 0
	endif
	nvar = 0

! Read header records

	call nfs(nf90_inq_dimid(ncid,'time',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
	if (nrec > mrec) then
		write (*,'("Error: Too many measurements:",i5)') nrec
		cycle files
	endif
	call nfs(nf90_get_att(ncid,nf90_global,'product',l1r_product))
	call nfs(nf90_get_att(ncid,nf90_global,'title',arg))
	if (arg /= 'CryoSat-2 Level-1 Retracked') then
		write (*,550) 'Error: Wrong input file'
		cycle files
	endif

	call nfs(nf90_get_att(ncid,nf90_global,'l1b_proc_time',l1b_proc_time))
	call nfs(nf90_get_att(ncid,nf90_global,'l1b_version',l1b_version))
	call nfs(nf90_get_att(ncid,nf90_global,'l1r_version',l1r_version))
	call nfs(nf90_get_att(ncid,nf90_global,'equator_longitude',eq_long))
	call nfs(nf90_get_att(ncid,nf90_global,'equator_time',eq_time))
	call nfs(nf90_get_att(ncid,nf90_global,'tai_utc',tai_utc))
	call nfs(nf90_get_att(ncid,nf90_global,'record_number',recnr))
	if (nf90_get_att(ncid,nf90_global,'mle_params',mle) /= nf90_noerr) mle = 3
	eq_time = eq_time + sec2000	! Equator time is already in UTC, other times are in TAI
	if (ndata + nrec > mrec) then
		write (*,'("Error: Too many accumulated measurements:",i5)') ndata+nrec
		cycle files
	endif

! Determine if this is version A (filenames ending in A001.nc)
! Determine if this originated from SAR

	version_a = index(l1r_product, '_A00') > 0
	sar = index(l1r_product, '_SIR_SA') > 0 .or. index(l1r_product, '_SIR_FBR') > 0

! Allocate arrays

	allocate (a(nrec),b(nrec),c(nrec),d(20,nrec), &
		t_1hz(nrec),t_20hz(20,nrec),alt(nrec),dh(nrec),valid(20,nrec),flags(nrec))

! Load time and location records

	call get_var (ncid, 'time_20hz', t_20hz)
	call get_var (ncid, 'time', t_1hz)
	call new_var ('time', t_1hz + sec2000 - tai_utc)
	call cpy_var ('lat', 'lat')
	! Compute ellipsoid corrections
	do i = 1,nrec
		dh(i) = dhellips(1,a(i))
	enddo
	call cpy_var ('lon', 'lon')
	call get_var (ncid, 'alt', alt)
	call new_var ('alt_gdrd', alt+dh)

! Compile flag bits

	call cpy_var ('mqe_20hz', 'mqe')
	valid = (t_20hz /= 0d0 .and. d <= 20d0)

	call get_var (ncid, 'retrack_flag_20hz', d)
	valid = (valid .and. d == 0)

	call get_var (ncid, 'surface_type', a)
	do i = 1,nrec
		b(i) = count(valid(:,i))
	enddo
	call new_var ('range_numval_ku', b)

	if (sar) then
		flags = 1
	else
		flags = 0
	endif
	call flag_set (a == 2, 2)
	call flag_set (a >= 2, 4)
	call flag_set (a >= 1, 5)
	call flag_set (b <= 10, 11)
	call flag_set (b <= 10, 12)
	call flag_set (b <= 10, 13)
	call new_var ('flags', dble(flags))

! Determine which star tracker is active. Bits 13,12,11 refer to star trackers 1,2,3
! Tests on subcycles 13-17 showed that:
! - All 20-Hz values are the same
! - A maximum of 1 star tracker is active

	call get_var (ncid, 'instr_config_flags_20hz', d)
	do i = 1,nrec
		flags = 0
		do m = 0,2
			do j = 1,20
				if (btest(nint(d(j,i)),13-m)) then
					flags = ibset (flags, m)
					exit
				endif
			enddo
		enddo
	enddo
	call new_var ('flags_star_tracker', dble(flags))

! Load other variables

	! USO factor = (nominal USO freq) / (measured USO freq)
	! Different scale factor in Version A
	call get_var (ncid, 'uso_corr_20hz', d)
	if (version_a) then
		uso_corr = 730d3 * 1d-9 * d(1,1)
	else
		uso_corr = 730d3 * 1d-15 * d(1,1)
	endif

	call cpy_var ('instr_range_corr_20hz', 'drange_cal')
	call cpy_var ('doppler_corr_20hz', 'drange_fm')

	call get_var (ncid, 'range_20hz+drange_20hz-alt_20hz', d)
	call trend_1hz (t_20hz, t_1hz, d, valid, a, b)
	call new_var ('range_ku', a + alt + uso_corr)
	call new_var ('range_rms_ku', b)

	call get_var (ncid, 'drange_20hz', d)
	call trend_1hz (t_20hz, t_1hz, d, valid, a, b)	! Temporary
	call new_var ('drange_ku', a)

	call cpy_var ('alt_rate_20hz', 'alt_rate')
	call cpy_var ('swh_20hz', 'swh_ku', 'swh_rms_ku')
	call cpy_var ('agc_20hz', 'agc_ku')
	call cpy_var ('agc_amp_20hz+dagc_eta_20hz+dagc_alt_20hz+dagc_xi_20hz+dagc_swh_20hz', 'sig0_ku', 'sig0_rms_ku')
	call cpy_var ('peakiness_20hz', 'peakiness_ku')

	if (mle == 4) call cpy_var ('xi_sq_20hz', 'off_nadir_angle2_wf_ku', 'off_nadir_angle2_wf_rms_ku')

! Convert pitch, roll, yaw from microradian to degrees and remove bias when MLE3

	call get_var (ncid, 'attitude_pitch_20hz', d)
	call mean_1hz (d/murad, valid, a, b)
	if (mle /= 4) a = a - pitch_bias
	call new_var ('attitude_pitch', a)
	c = a*a

	call get_var (ncid, 'attitude_roll_20hz', d)
	call mean_1hz (d/murad, valid, a, b)
	if (mle /= 4) a = a - roll_bias
	call new_var ('attitude_roll', a)
	c = c + a*a

	call get_var (ncid, 'attitude_yaw_20hz', d)
	call mean_1hz (d/murad, valid, a, b)
	if (mle /= 4) a = a - yaw_bias
	call new_var ('attitude_yaw', a)

	call new_var ('off_nadir_angle2_pf', c)

	call cpy_var ('noise_20hz', 'noise_floor_ku', 'noise_floor_rms_ku')
	call cpy_var ('dry_tropo', 'dry_tropo_ecmwf')
	call cpy_var ('wet_tropo', 'wet_tropo_ecmwf')
	call cpy_var ('iono_model', 'iono_bent')
	call cpy_var ('iono_gim', 'iono_gim')
	call cpy_var ('inv_baro', 'inv_bar_static')

	! Version A: DAC already includes inv_baro, which is not the intent
	! In L1R version 1.23 and later inv_baro is already taken out in the L1R file
	if (version_a .and. l1r_version < '1.23') then
		call cpy_var ('dac', 'inv_bar_mog2d')
	else
		call cpy_var ('inv_baro+dac','inv_bar_mog2d')
	endif
	call cpy_var ('tide_solid', 'tide_solid')
	call cpy_var ('tide_ocean', 'tide_ocean_got00')
	call cpy_var ('tide_load', 'tide_load_got00')
	call cpy_var ('tide_pole', 'tide_pole')
	call cpy_var ('tide_lp', 'tide_equil')

	! If input file is split between ascending/descending, dump the first chunk,
	! move the second chunk down and update the equator crossing to the new pass

	if (passnr(2) /= passnr(1)) then
		ndata = ndata + recnr(1)
		call put_rads (cycnr(1), passnr(1), ndata)

		! Move the data to be beginning
		do i = 1,nvar
			var(i)%d(1:recnr(2)) = var(i)%d(ndata+1:ndata+recnr(2))
		enddo

		! Update equator crossing info to the next pass
		eq_long = modulo (eq_long + 0.5d0 * rev_long + 180d0, 360d0)
		eq_time = eq_time + 0.5d0 * rev_time
		ndata = recnr(2)
	else
		ndata = ndata + nrec
	endif
	oldcyc = cycnr(2)
	oldpass = passnr(2)

	deallocate (a,b,c,d,t_1hz,t_20hz,alt,dh,valid,flags)

	call nfs(nf90_close(ncid))

enddo files ! Each file

! Dump whatever remains

call put_rads (oldcyc, oldpass, ndata)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Write CryoSat-2 L1R data to RADS', flag=flag)) return
write (*,1310)
1310 format (/ &
'syntax: rads_gen_c2_l1r [options] < list_of_L1R_file_names'// &
'This program converts CryoSat-2 L1R files to RADS data'/ &
'files with the name $RADSDATAROOT/data/c2/F/pPPPP/c2pPPPPcCCC.nc.'/ &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Copy variable to RADS
!-----------------------------------------------------------------------

subroutine cpy_var (varin, varout, varrms)
character(len=*), intent(in) :: varin, varout
character(len=*), intent(in), optional :: varrms
if (index(varin,'_20hz') == 0) then ! 1-Hz variable
	call get_var (ncid, varin, a)
	call new_var (varout, a)
else ! 20-Hz variable
	call get_var (ncid, varin, d)
	call mean_1hz (d, valid, a, b)
	call new_var (varout, a)
	if (present(varrms)) call new_var (varrms, b)
endif
end subroutine cpy_var

!-----------------------------------------------------------------------
! Create new RADS variable
!-----------------------------------------------------------------------

subroutine new_var (varnm, data)
! Write variables one after the other to the output file
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(in) :: data(:)
nvar = nvar + 1
if (nvar > mvar) stop 'Too many variables'
var(nvar)%v => rads_varptr (S, varnm)
var(nvar)%d = data
end subroutine new_var

!-----------------------------------------------------------------------
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads (cycnr, passnr, ndata)
integer(fourbyteint), intent(in) :: cycnr, passnr, ndata
integer(fourbyteint) :: i

if (ndata == 0) return	! Skip empty data sets
if (cycnr < c0 .or. cycnr > c1) return	! Skip chunks that are not of the selected cycle
if (eq_time < t0 .or. eq_time > t1) return	! Skip equator times that are not of selected range

! Store relevant info
call rads_init_pass_struct (S, P)
P%cycle = cycnr
P%pass = passnr
P%start_time = var(1)%d(1)
P%end_time = var(1)%d(ndata)
P%equator_time = eq_time
P%equator_lon = eq_long
P%original = 'L1R ('//trim(l1r_version)//') from L1B ('// &
	trim(l1b_version)//') data of '//l1b_proc_time

! Open output file
call rads_create_pass (S, P, ndata)

! Define all variables
do i = 1,nvar
	call rads_def_var (S, P, var(i)%v)
enddo

! Fill all the data fields
do i = 1,nvar
	call rads_put_var (S, P, var(i)%v, var(i)%d(1:ndata))
enddo

! Close the data file
write (*,552) nrec,trim(P%filename(len_trim(S%dataroot)+2:))
call rads_close_pass (S, P)

! Formats
552 format ('...',i5,' records written to ',a)

end subroutine put_rads

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

end program rads_gen_c2_l1r
