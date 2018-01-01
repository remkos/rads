!-----------------------------------------------------------------------
! Copyright (c) 2011-2018  Remko Scharroo
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
! Variables to be written to RADS are:
! time - Time since 1 Jan 85 (see note on timing bias below)
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
!
! tbias:
! - Apply timing bias according to Marco's table (different for SAR and LRM data)
! - Account for the change in timing bias in SIR1FDM/2.4 (and later) (See Ruby's e-mail of 22 Apr 2013)
! - Adjust altitude from the product accordingly (using altitude rate)
! - Note that this does NOT change the equator time or longitude!
! - This is done here, because FDM and LRM can later no longer be distinguished from eachother
! - No longer needed for Baseline C
!
! range:
! - Supposedly we need to add 7 mm (1/64 of a range gate) to FDM/LRM L1 and
!   FDM/LRM L2 range because of error in CAL1 (fixed after L1R version 2.03)
! - Subtract 673 mm from range in Baseline B products
!-----------------------------------------------------------------------
program rads_gen_c2_l1r

use rads
use rads_devel
use rads_gen
use rads_misc
use rads_netcdf
use rads_time
use netcdf

! Command line arguments

integer(fourbyteint) :: ios
character(len=rads_cmdl) :: filename, arg
character(len=rads_strl) :: filenames = ''

! Header variables

character(len=19) :: l1b_proc_time
character(len=14) :: l1b_version, l1r_version
character(len=55) :: l1r_product
character(len=1) :: baseline
real(eightbytereal) :: tai_utc, eq_time, eq_long
integer(fourbyteint) :: passnr(2), cycnr(2), recnr(2), nrec, ncid, varid, doris_nav
logical :: sar, fdm, lrm

! Data variables

integer(fourbyteint), parameter :: mrec=6000, mvar=60
integer(fourbyteint) :: nvar=0, ndata=0
real(eightbytereal), allocatable :: a(:),b(:),c(:),d(:,:),w(:,:,:),t_1hz(:),t_20hz(:,:),alt(:),alt_20hz(:,:),dh(:)
logical, allocatable :: t_valid(:,:),valid(:,:)
integer(fourbyteint), allocatable :: nvalid(:)
integer(twobyteint), allocatable :: flags(:)
type(rads_sat) :: S
type(rads_pass) :: P
type :: var_
	type(rads_var), pointer :: v ! Pointer to rads_var struct
	real(eightbytereal), allocatable :: d1(:), d2(:,:), d3(:,:,:) ! Data arrays
	logical :: skip ! .true. if to be skipped
endtype
type(var_) :: var(mvar)

! Other local variables

integer(fourbyteint), parameter :: maxint4=2147483647
real(eightbytereal), parameter :: fai = 7.3d-3
real(eightbytereal), parameter :: sec2000=473299200d0, rev_time = 5953.45d0, rev_long = -24.858d0
real(eightbytereal), parameter :: pitch_bias = 0.096d0, roll_bias = 0.086d0, yaw_bias = 0d0	! Attitude biases to be added
real(eightbytereal) :: uso_corr, dhellips, tbias, range_bias
integer(fourbyteint) :: i, j, m, oldcyc=0, oldpass=0, mle=3

! Initialise

call synopsis
call rads_gen_getopt ('c2', 'mw with-20hz with-wvf')
call synopsis ('--head')
call rads_init (S, sat)

!----------------------------------------------------------------------
! Read all file names for standard input
!----------------------------------------------------------------------

do
	read (*,'(a)',iostat=ios) filename
	if (ios /= 0) exit

! Open input file

	call log_string (filename)
	if (nf90_open(filename,nf90_nowrite,ncid) /= nf90_noerr) then
		call log_string ('Error: failed to open input file', .true.)
		cycle
	endif

! Check input file type

	call nfs(nf90_get_att(ncid,nf90_global,'title',arg))
	if (arg /= 'CryoSat-2 Level-1 Retracked') then
		call log_string ('Error: wrong input file type', .true.)
		cycle
	endif

! Read cycle and pass number

	call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cycnr))
	call nfs(nf90_get_att(ncid,nf90_global,'pass_number',passnr))

! Read length of time dimension.
! Throw out some long files (in cycle 7 only?) that span 3 passes

	call nfs(nf90_inq_dimid(ncid,'time',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
	if ((cycnr(1) == cycnr(2) .and. passnr(1) > passnr(2)) .or. &
		(passnr(1) == passnr(2) .and. nrec > 4000)) then
		call log_string ('Error: file too long (covers 3 passes)', .true.)
		cycle
	else if (nrec > mrec) then
		call log_string ('Error: too many measurements', .true.)
		cycle
	endif

! Determine if we need to dump pending data

	if (passnr(1) /= oldpass .or. cycnr(1) /= oldcyc) then
		call put_rads (oldcyc, oldpass, ndata)
		ndata = 0
		filenames = ''
	endif
	nvar = 0

! Read remaining header records

	call nfs(nf90_get_att(ncid,nf90_global,'product',l1r_product))
	call nfs(nf90_get_att(ncid,nf90_global,'l1b_proc_time',l1b_proc_time))
	call nfs(nf90_get_att(ncid,nf90_global,'l1b_version',l1b_version))
	call nfs(nf90_get_att(ncid,nf90_global,'l1r_version',l1r_version))
	call nfs(nf90_get_att(ncid,nf90_global,'doris_nav',doris_nav))
	call nfs(nf90_get_att(ncid,nf90_global,'equator_longitude',eq_long))
	call nfs(nf90_get_att(ncid,nf90_global,'equator_time',eq_time))
	call nfs(nf90_get_att(ncid,nf90_global,'tai_utc',tai_utc))
	call nfs(nf90_get_att(ncid,nf90_global,'record_number',recnr))
	if (nf90_get_att(ncid,nf90_global,'mle_params',mle) /= nf90_noerr) mle = 3
	eq_time = eq_time + sec2000	! Equator time is already in UTC, other times are in TAI
	if (ndata + nrec > mrec) then
		call log_string ('Error: too many accumulated measurements', .true.)
		cycle
	endif

! Determine baseline version

	baseline = l1r_product(52:52)

! Determine if this originated from SAR, FDM or LRM

	sar = (l1b_version(:7) == 'SIR1SAR')
	fdm = (l1b_version(:7) == 'SIR1FDM')
	lrm = (l1b_version(:7) == 'SIR1LRM')
	if (.not.sar .and. .not.fdm .and. .not.lrm) then
		call log_string ('Error: unrecognised file format', .true.)
		cycle
	endif

! Allocate arrays

	allocate (a(nrec),b(nrec),c(nrec),d(20,nrec),w(256,20,nrec), &
		t_1hz(nrec),t_20hz(20,nrec),alt(nrec),alt_20hz(20,nrec),dh(nrec), &
		t_valid(20,nrec),valid(20,nrec),nvalid(nrec),flags(nrec))


! Apply timing bias here, because they are different between between FDM, LRM, SAR
! See IPF1_datation_biases_v4.xlsx by Marco Fornari for Baseline B data.
! These have been accounted for in Baseline C.

	if (baseline >= 'C') then
		! 18-Mar-2015: The TDS_Nov2010 suggests that there is no more timing bias in Baseline C; Marco supports that idea.
		tbias = 0d0
	else if (sar) then
		! Timing bias for baseline B, in SAR mode (also PLRM)
		tbias = +0.520795d-3
	else
		! Timing bias for baseline B, in LRM mode
		tbias = -4.699112d-3
		! Partial correction of timing bias (See Ruby's e-mail of 22 Apr 2013)
		if (fdm .and. l1b_version(9:11) >= '2.4') tbias = tbias + 4.4436d-3
		!  1-Aug-2013: Additional timing bias from my own research
		!  5-Nov-2014: This applies ONLY to FDM/LRM data; CP4O demonstrated that this does not apply to PLRM/SAR!
		tbias = tbias + 0.4d-3
	endif

! Time information

	call get_var (ncid, 'time', t_1hz)
	call new_var ('time', t_1hz + sec2000 - tai_utc + tbias)
	call get_var (ncid, 'time_20hz', t_20hz)
	t_valid = (t_20hz /= 0d0)
	valid = t_valid
	where (.not.t_valid) t_20hz = nan
	call new_var_2d ('time_20hz', t_20hz + sec2000 - tai_utc + tbias)

! Location information

	call cpy_var ('lat', 'lat')

	! Compute ellipsoid corrections
	do i = 1,nrec
		dh(i) = dhellips(1,a(i))
	enddo
	call cpy_var ('lat_20hz', 'lat_20hz')
	call cpy_var ('lon', 'lon')
	call cpy_var ('lon_20hz', 'lon_20hz')
	call get_var (ncid, 'alt', alt)
	call get_var (ncid, 'alt_20hz', alt_20hz)

	! If input is FDM and there is no DORIS Navigator orbit (i.e. predicted orbit)
	! we blank the orbit out entirely: it would be useless anyhow
	if (fdm .and. doris_nav == 0) dh = nan

	! Update altitude, taking into account ellipsoid correction and timing bias
	call cpy_var ('alt_rate_20hz', '', 'alt_rate')
	call new_var ('alt_cnes', alt + dh + tbias * a)
	if (nhz /= 0) then
		forall (i = 1:20) d(i,:) = alt_20hz(i,:) + dh(:) + tbias * d(i,:)
		call new_var_2d ('alt_cnes_20hz', d)
	endif

! 20-Hz corrections

	call cpy_var ('instr_range_corr_20hz', '', 'drange_cal')
	call cpy_var ('doppler_corr_20hz', '', 'drange_fm')

! Compile flag bits; needs to be done BEFORE any averaging of the measurements

	call get_var (ncid, 'mqe_20hz', d)
	valid = (t_valid .and. d <= 20d0)
	call get_var (ncid, 'retrack_flag_20hz', d)
	valid = (valid .and. d == 0)
	do i = 1,nrec
		nvalid(i) = count(valid(:,i))
	enddo

	call get_var (ncid, 'surface_type', a)
	if (l1r_version <= '1.26') a = a * 1d3 ! Error in scale_factor
	if (sar) then
		flags = 1 ! Set bit 1 for SAR
	else
		flags = 0
	endif
	call flag_set (nint(a) == 2, 2)
	call flag_set (nint(a) >= 2, 4)
	call flag_set (nint(a) >= 1, 5)
	call flag_set (nvalid <= 10, 11)
	call flag_set (nvalid <= 10, 12)
	call flag_set (nvalid <= 10, 13)

	call new_var ('flags', dble(flags))

! Determine range bias, prior to Baseline C only
! 1) According to Marco Fornari:
!    Before Baselince C, all FDM/LRM L1/L2 data have a CAL1 which is off by 1 FAI (1/64 gate).
!    Marco suggested to ADD 7.3 mm to range, which is done in retracking at the moment, but I feel it needs to be
!    subtracted from range afterall. Hence subtract 2 FAI here before the change was made in cs2_l1b_to_l1r
! 2) Baseline C is corrected for 673 mm bias. We do that here to Baseline B as well.

	if (baseline >= 'C') then
		range_bias = 0d0
	else
		range_bias = -673d-3
		if (.not.sar .and. l1r_version < "2.03") range_bias = range_bias - 2*fai
	endif

! Range measurements

	! USO factor = (nominal USO freq) / (measured USO freq)
	call get_var (ncid, 'uso_corr_20hz', d)
	uso_corr = 730d3 * d(1,1)

	call get_var (ncid, 'range_20hz drange_20hz ADD', d)
	where (.not.valid) d = nan
	call trend_1hz (t_20hz, t_1hz, d - alt_20hz, a, b)
	call new_var ('range_ku', a + alt + uso_corr + range_bias)
	if (nhz /= 0) then
		call new_var_2d ('range_20hz_ku', d + uso_corr + range_bias)
		d = 1
		where (valid) d = 0
		call new_var_2d ('range_used_20hz_ku', d)
	endif
	call new_var ('range_rms_ku', b)
	call new_var ('range_numval_ku', dble(nvalid))

! Waves and backscatter

	call cpy_var ('swh_20hz', 'swh_20hz_ku', 'swh_ku', 'swh_rms_ku')
	call cpy_var ('agc_20hz', '', 'agc_ku')
	call cpy_var ('agc_amp_20hz dagc_eta_20hz ADD dagc_alt_20hz ADD dagc_xi_20hz ADD dagc_swh_20hz ADD', &
		'sig0_20hz_ku', 'sig0_ku', 'sig0_rms_ku')

	if (mle == 4) call cpy_var ('xi_sq_20hz', '', 'off_nadir_angle2_wf_ku', 'off_nadir_angle2_wf_rms_ku')

! Convert pitch, roll, yaw from microradian to degrees and remove bias when MLE3

	call get_var (ncid, 'attitude_pitch_20hz', d)
	where (.not.t_valid) d = nan
	d = d / rad
	if (mle /= 4) d = d - pitch_bias
	call mean_1hz (d, a, b)
	call new_var ('attitude_pitch', a)
	call new_var_2d ('attitude_pitch_20hz', d)
	c = a*a

	call get_var (ncid, 'attitude_roll_20hz', d)
	where (.not.t_valid) d = nan
	d = d / rad
	if (mle /= 4) d = d - roll_bias
	call mean_1hz (d, a, b)
	call new_var ('attitude_roll', a)
	call new_var_2d ('attitude_roll_20hz', d)
	c = c + a*a

	call get_var (ncid, 'attitude_yaw_20hz', d)
	where (.not.t_valid) d = nan
	d = d / rad
	if (mle /= 4) d = d - yaw_bias
	call mean_1hz (d, a, b)
	call new_var ('attitude_yaw', a)
	call new_var_2d ('attitude_yaw_20hz', d)

	call new_var ('off_nadir_angle2_pf', c)

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

! Waveform-related info

	call cpy_var ('peakiness_20hz', 'peakiness_20hz_ku', 'peakiness_ku')
	call cpy_var ('mqe_20hz', 'mqe_20hz_ku', 'mqe')
	call cpy_var ('noise_20hz', 'noise_floor_20hz_ku', 'noise_floor_ku', 'noise_floor_rms_ku')

! We go back to using very limited editing now to dump waveforms

	valid = t_valid
	if (nwvf > 0) then
		call get_var (ncid, 'range_20hz', d)
		call new_var_2d ('range_tracker_20hz_ku', d + uso_corr + range_bias)
		call cpy_var ('agc_20hz', 'agc_20hz_ku')
		call cpy_var ('echo_scale_20hz', 'waveform_scale_20hz')
		call cpy_var ('waveform_20hz', 'waveform_20hz')
	endif

! Geophysical corrections

	call cpy_var ('dry_tropo', 'dry_tropo_ecmwf')
	call cpy_var ('wet_tropo', 'wet_tropo_ecmwf')
	call cpy_var ('iono_model', 'iono_bent')
	call cpy_var ('iono_gim', 'iono_gim')
	call cpy_var ('inv_baro', 'inv_bar_static')

	if (.not.fdm) call cpy_var ('inv_baro dac ADD','inv_bar_mog2d')
	call cpy_var ('tide_solid', 'tide_solid')
	call cpy_var ('tide_ocean', 'tide_ocean_got00')
	call cpy_var ('tide_load', 'tide_load_got00')
	call cpy_var ('tide_pole', 'tide_pole')
	call cpy_var ('tide_lp', 'tide_equil')

	! Add current filename to list of input files

	j = index(filename, '/', .true.) + 1
	filenames = trim(filenames) // rads_linefeed // filename(j:)

	! In some cases the last 1-Hz measurement is invalid, so we need to skip it

	if (t_1hz(nrec) > 0) then
		! This is OK
	else if (recnr(2) == 0) then
		recnr(1) = recnr(1) - 1
	else
		recnr(2) = recnr(2) - 1
	endif

	! If input file is split between ascending/descending, dump the first chunk,
	! move the second chunk down and update the equator crossing to the new pass

	if (passnr(2) /= passnr(1)) then
		ndata = ndata + recnr(1)
		call put_rads (cycnr(1), passnr(1), ndata)

		! Move the data to be beginning
		do i = 1,nvar
			select case (var(i)%v%info%ndims)
			case (1)
				var(i)%d1(1:recnr(2)) = var(i)%d1(ndata+1:ndata+recnr(2))
			case (2)
				var(i)%d2(:,1:recnr(2)) = var(i)%d2(:,ndata+1:ndata+recnr(2))
			case default
				var(i)%d3(:,:,1:recnr(2)) = var(i)%d3(:,:,ndata+1:ndata+recnr(2))
			end select
		enddo
		j = index(filename, '/', .true.) + 1
		filenames = rads_linefeed // filename(j:)

		! Update equator crossing info to the next pass
		eq_long = modulo (eq_long + 0.5d0 * rev_long + 180d0, 360d0)
		eq_time = eq_time + 0.5d0 * rev_time
		ndata = recnr(2)
	else
		ndata = ndata + recnr(1)
	endif
	oldcyc = cycnr(2)
	oldpass = passnr(2)

	deallocate (a,b,c,d,w,t_1hz,t_20hz,alt,alt_20hz,dh,t_valid,valid,nvalid,flags)

	call nfs(nf90_close(ncid))

enddo ! Each file

! Dump whatever remains

call put_rads (oldcyc, oldpass, ndata)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Write CryoSat-2 L1R data to RADS', flag=flag)) return
call synopsis_devel (' < list_of_L1R_file_names')
write (*,1310)
1310 format (/ &
'Program specific [program_options] are:' / &
'  -m, --with-20hz           Include 20-Hz variables in addition to 1-Hz variables'/ &
'  -w, --with-wvf            Include waveforms (implies --with-20hz)'// &
'This program converts CryoSat-2 L1R files to RADS data' / &
'files with the name $RADSDATAROOT/data/c2/F/pPPPP/c2pPPPPcCCC.nc.' / &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Copy variable to RADS
!-----------------------------------------------------------------------

subroutine cpy_var (varin, varout, varmean, varrms)
! Copy 1-Hz data to memory buffer, or
! copy 20-Hz data or waveforms to memory buffer, optionally
! making mean and standard deviation.
! At return:
! a = 1-Hz (mean)
! b = 1-Hz stddev
! d = 20-Hz
! w = waveform
character(len=*), intent(in) :: varin, varout
character(len=*), intent(in), optional :: varmean, varrms
if (index(varin,'_20hz') == 0) then
	! 1-Hz variable
	call get_var (ncid, varin, a)
	call new_var (varout, a) ! Copy 1-Hz data
else if (varin == 'waveform_20hz') then
	! Waveform data copied varbatim
	call get_var (ncid, varin, w)
	call new_var_3d (varout, w)
else if (present(varmean)) then
	! 20-Hz variable to be averaged
	call get_var (ncid, varin, d)
	where (.not.valid) d = nan
	call new_var_2d (varout, d) ! Copy 20-Hz data
	call mean_1hz (d, a, b)
	call new_var (varmean, a)
	if (present(varrms)) call new_var (varrms, b)
else if (varout /= '' .and. nhz /= 0) then
	! 20-Hz variable copied verbatim
	call get_var (ncid, varin, d)
	where (.not.valid) d = nan
	call new_var_2d (varout, d) ! Copy 20-Hz data
endif
end subroutine cpy_var

!-----------------------------------------------------------------------
! Create new RADS variable
!-----------------------------------------------------------------------

subroutine new_var (varnm, data)
! Store 1-Hz variables to be written later by put_var
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(in) :: data(:)
if (varnm == '') return
nvar = nvar + 1
if (nvar > mvar) stop 'Too many variables'
if (.not.allocated(var(nvar)%d1)) allocate(var(nvar)%d1(mrec))
var(nvar)%v => rads_varptr (S, varnm)
var(nvar)%d1(ndata+1:ndata+nrec) = data(1:nrec)
end subroutine new_var

subroutine new_var_2d (varnm, data)
! Store 20-Hz variables to be written later by put_var
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(in) :: data(:,:)
if (varnm == '' .or. nhz == 0) return
nvar = nvar + 1
if (nvar > mvar) stop 'Too many variables'
if (.not.allocated(var(nvar)%d2)) allocate(var(nvar)%d2(20,mrec))
var(nvar)%v => rads_varptr (S, varnm)
var(nvar)%d2(:,ndata+1:ndata+nrec) = data(:,1:nrec)
end subroutine new_var_2d

subroutine new_var_3d (varnm, data)
! Store 20-Hz variables to be written later by put_var
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(in) :: data(:,:,:)
if (varnm == '' .or. nwvf==0) return
nvar = nvar + 1
if (nvar > mvar) stop 'Too many variables'
if (.not.allocated(var(nvar)%d3)) allocate(var(nvar)%d3(256,20,mrec))
var(nvar)%v => rads_varptr (S, varnm)
var(nvar)%d3(:,:,ndata+1:ndata+nrec) = data(:,:,1:nrec)
end subroutine new_var_3d

!-----------------------------------------------------------------------
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads (cycnr, passnr, ndata)
integer(fourbyteint), intent(in) :: cycnr, passnr, ndata
integer(fourbyteint) :: i, j

if (ndata == 0) return	! Skip empty data sets
if (cycnr < cycles(1) .or. cycnr > cycles(2)) return	! Skip chunks that are not of the selected cycle
if (eq_time < times(1) .or. eq_time > times(2)) return	! Skip equator times that are not of selected range

! Store relevant info
call rads_init_pass_struct (S, P)
P%cycle = cycnr
P%pass = passnr
P%start_time = var(1)%d1(1)
P%end_time = var(1)%d1(ndata)
P%equator_time = eq_time
P%equator_lon = eq_long
P%original = 'L1R ('//trim(l1r_version)//') from L1B ('// &
	trim(l1b_version)//') data of '//trim(l1b_proc_time)//filenames

! Open output file
call rads_create_pass (S, P, ndata, nhz, nwvf)

! Check for variables we want to skip because they are empty
do i = 1,nvar
	if (var(i)%v%name == 'inv_bar_mog2d') then
		var(i)%skip = all(var(i)%d1(1:ndata) == 0d0)
	else
		var(i)%skip = .false.
	endif
enddo

! Define all variables
do i = 1,nvar
	if (var(i)%skip) cycle
	call rads_def_var (S, P, var(i)%v)
	if (i == 1 .and. nhz > 0) call rads_def_var (S, P, 'meas_ind')
	if (i == 1 .and. nwvf > 0) call rads_def_var (S, P, 'wvf_ind')
enddo

! Fill all the data fields
do i = 1,nvar
	if (var(i)%skip) cycle
	select case (var(i)%v%info%ndims)
	case (1)
		call rads_put_var (S, P, var(i)%v, var(i)%d1(1:ndata))
	case (2)
		call rads_put_var (S, P, var(i)%v, var(i)%d2(:,1:ndata))
	case default
		call rads_put_var (S, P, var(i)%v, var(i)%d3(:,:,1:ndata))
	end select
	if (i == 1 .and. nhz > 0) call rads_put_var (S, P, 'meas_ind', (/(j*1d0,j=0,19)/))
	if (i == 1 .and. nwvf > 0) call rads_put_var (S, P, 'wvf_ind', (/(j*1d0,j=0,255)/))
enddo

! Close the data file
call log_records (ndata, P)
call rads_close_pass (S, P)

end subroutine put_rads

!-----------------------------------------------------------------------
! Set a bit in an array of flags
!-----------------------------------------------------------------------

subroutine flag_set (a, bit)
logical, intent(in) :: a(:)
integer(fourbyteint), intent(in) :: bit
integer(fourbyteint) :: i
integer(twobyteint) :: j
j = int(bit,twobyteint)
do i = 1,nrec
	if (a(i)) flags(i) = ibset(flags(i),j)
enddo
end subroutine flag_set

end program rads_gen_c2_l1r
