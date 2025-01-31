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

!*rads_add_dual -- Smooth dual-frequency iono data in RADS
!+
! This program adjusts the contents of RADS altimeter data files
! by smoothing the dual-frequency ionosphere correction.
!
! There are two types of smoothing implemented:
!
! * Box car filtering
!   This is a simple box car filter with a length of 35 seconds, in contrast
!   to the 21 seconds suggested in Imel [1994].
!   Editing is based on flag bits and a range relative to the mean iono
!
! * Median + Lanczos filtering
!   This is similar to the implementation of S3 and S6, except for the
!   spline interpolation at the end. Pre-editing is done based on surface
!   classification, nr values used in range compression and range rms.
!
! References:
! D. A. Imel, Evaluation of the TOPEX/Poseidon dual-frequency ionosphere
! correction, J. Geophys. Res., 99, C12, 24,895-24,906, 1994.
!
! Sentinel-6/Jason-CS ALT Level 2 Product Generation Specification (L2 ALT PGS),
! EUM/LEO-JASCS/SPE/17/901321, v5A, July 2023.
!
! Usage: rads_add_dual [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_dual

use rads
use rads_devel
use rads_misc
use meteo_subs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

integer(fourbyteint) :: ntot, cyc, pass, i, ios
integer(fourbyteint), parameter :: nmax = 3000000
character(len=40) :: ext = ''
real(eightbytereal) :: bias=0d0, f
real(eightbytereal) :: time(nmax), iono_in(nmax), iono_edit(nmax), iono_out(nmax), tmp1(nmax), tmp2(nmax), tmp3(nmax)
logical :: recompute = .false., new = .false., do_box = .true.
type(rads_var), pointer :: iono_alt

! Parameter values for box filtering

type :: box_par
	integer(twobyteint) :: mask
	real(eightbytereal) :: iwin, twin
endtype
type(box_par) :: box = box_par(2027, 35d0, 8d0)

! Parameter values for Lanczos filtering

type :: lanczos_par
	integer(fourbyteint) :: half_width_median, min_nr_median, half_width_lanczos, min_nr_lanczos, cutoff_lanczos
endtype
type(lanczos_par) :: lanczos = lanczos_par(20, 5, 50, 5, 50)

! Scan command line for options

call synopsis ('--head')
call rads_set_options ('nmx:b::l::r:: new mle: box:: lanczos:: ext: recompute::')
call rads_init (S)

! Determine conversion factor from TEC units to ionospheric delay in metres
! and mean altitude.

f = 1d0/(1d0-(S%frequency(1)/S%frequency(2))**2)

! Check options

do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('n', 'new')
		new = .true.
	case ('m', 'mle')	! For backward compatibility only
		if (rads_opt(i)%arg == '3') ext = '_mle3'
	case ('x', 'ext')
		ext = '_' // rads_opt(i)%arg(:39)
	case ('b', 'box')
		do_box = .true.
		read (rads_opt(i)%arg, *, iostat=ios) box
	case ('l', 'lanczos')
		do_box = .false.
		read (rads_opt(i)%arg, *, iostat=ios) lanczos
	case ('r', 'recompute')
		recompute = .true.
		read (rads_opt(i)%arg, *, iostat=ios) bias
	end select
enddo

! Find iono_alt[_ext] variable

iono_alt => rads_varptr (S, 'iono_alt' // ext)

! For box car filter:
! Use half the time window (in seconds) and half the range window (converted to meters)
! in the rest of the program

if (do_box) then
	box%twin = box%twin/2d0
	box%iwin = box%iwin/2d2 + 1d-6
endif

! Run process for all files within a cycle, first reading all, then writing all

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	ntot = 0
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call read_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo

	if (do_box) then
		call process_cycle_box (ntot)
	else
		call process_cycle_lanczos (ntot)
	endif

	ntot = 0
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call write_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Smooth the dual-frequency ionospheric delay', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310) box, lanczos
1310 format (/ &
'Additional [processing_options] are:' / &
'  -n, --new                 Only add smoothed iono when not yet existing' / &
'  -b, --box[=PARAMS]        Do the old box car filter with given parameters (MASK,IWIN,TWIN):' / &
'                            * MASK: exclude data based on bitmap (default:',i5,')' / &
'                            * IWIN: set editing range for iono data [cm] (default:',f4.0,')' / &
'                            * TWIN: set box car filter length [sec] (default:' f4.0,')' / &
'  -l, --lanczos[=PARAMS]    Do the median+Lanczos filter with given parameters:' / &
'                            * Half-width of the median filter [points] (default:',i3,')' / &
'                            * Minimum nr of valid points for median filter (default:',i3,')' / &
'                            * Half-width of the Lanczos filter [points] (default:',i3,')' / &
'                            * Minimum nr of valid points for Lanczos filter (default:',i3,')' / &
'                            * Cutoff period for the Lanczos filter [points] (default:',i3,')' / &
'  -r, --recompute[=BIAS]    (Re)compute ionospheric correction from delta range' / &
'                            with optional BIAS (m) added' / &
'  -x, --ext EXT             Use parameters with extension _EXT (e.g. mle3 or plrm)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Read a single pass
!-----------------------------------------------------------------------

subroutine read_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: range_ku(n), range_c(n), ssb_ku(n), ssb_c(n), iono(n), tmp(n)
integer(fourbyteint) :: i
integer(twobyteint) :: flag

! Read all the required variables

call rads_get_var (S, P, 'time', time(ntot+1:ntot+n), .true.)
if (recompute) then
	call rads_get_var (S, P, 'range_ku' // ext, range_ku, .true.)
	call rads_get_var (S, P, 'range_c', range_c, .true.)
	call rads_get_var (S, P, 'ssb' // ext, ssb_ku, .true.)
	call rads_get_var (S, P, 'ssb_c', ssb_c, .true.)
	iono = f * ((range_c + ssb_c) - (range_ku + ssb_ku)) + bias
else
	call rads_get_var (S, P, iono_alt, iono, .true.)
endif

! Do editing

iono_in(ntot+1:ntot+n) = iono
where (iono < iono_alt%info%limits(1) .or. iono > iono_alt%info%limits(2)) iono = nan ! Apply default edit level
if (do_box .and. box%mask /= 0) then
	! For box car, edit based on flags
	call rads_get_var (S, P, 'flags' // ext, tmp, .true.)
	do i = 1,n
		flag = nint(tmp(i),twobyteint)
		if (iand(flag,box%mask) /= 0) iono(i) = nan
	enddo
else if (.not.do_box) then
	! For Lanczos, edit based on surface class, number of values and range RMS
	call rads_get_var (S, P, 'surface_class', tmp, .true.)
	where (tmp /= 0 .and. tmp /= 2) iono = nan
	call rads_get_var (S, P, 'range_numval_ku' // ext, tmp) ! Apply default edit level
	where (isnan_(tmp)) iono = nan
	call rads_get_var (S, P, 'range_rms_ku' // ext, tmp) ! Apply default edit level
	where (isnan_(tmp)) iono = nan
endif
iono_edit(ntot+1:ntot+n) = iono

! Up the counter

ntot = ntot + n
end subroutine read_pass

!-----------------------------------------------------------------------
! Process an entire cycle with box car filter
!-----------------------------------------------------------------------

subroutine process_cycle_box (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, j, k, j0, j1, nsum
real(eightbytereal) :: t0, t1, wsum

! Do box car filter for each measurement (i).

j0 = 1
j1 = 1
do i = 1,n

! Find first measurement in window (starts at t0)

	t0 = time(i) - box%twin
	do while (time(j0) < t0)
		j0 = j0 + 1
	enddo

! Find last measurement in window (ends at t1)

	t1 = time(i) + box%twin
	j1 = max(i,j1)
	do while (j1 < n .and. time(j1+1) < t1)
		j1 = j1 + 1
	enddo

! Now do the box car filter to compute average
! First time with only rough outlier editing
! Second time with editing based on +/-iwin
! (This is a slight change from radsp_dual, because here in the second loop
! a iono_in > iwin could be accepted. In radsp_dual those were permanently rejected.)

	t0 = -1d0
	t1 = box%iwin
	do k = 1,2
		nsum = 0
		wsum = 0d0
		do j = j0,j1
			if (iono_edit(j) >= t0 .and. iono_edit(j) <= t1) then
				nsum = nsum + 1
				wsum = wsum + iono_in(j)
			endif
		enddo
		if (nsum == 0) exit
		wsum = wsum / nsum
		t0 = wsum - box%iwin
		t1 = wsum + box%iwin
	enddo

! Assign NaN in case nsum < 4

	if (nsum < 4) then
		iono_out(i) = nan
	else
		iono_out(i) = wsum
	endif

enddo

end subroutine process_cycle_box

!-----------------------------------------------------------------------
! Process an entire cycle with median and Lanczos filter
!-----------------------------------------------------------------------

subroutine process_cycle_lanczos (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: idx(n), ns

! Spread the values equidistantly (in time).
! This is done so that the time intervals simply become index ranges.
! One step in index is a "1-Hz" time interval, which is close to 1 second.

idx(:n) = nint((time(:n)-time(1)) / S%dt1hz + 1)
ns = idx(n)	! Number of "seconds"
if (ns > nmax) call rads_exit('Too many measurements on input to filter process')
tmp1(:ns) = nan ! Initialise with NaN
tmp1(idx(:n)) = iono_edit(:n) ! Fill the array with input values

! First do a median filter, then Lanczos filter

call median_filter (tmp1(:ns), lanczos%min_nr_median, lanczos%half_width_median, tmp2(:ns))
call lanczos_filter (tmp2(:ns), lanczos%min_nr_lanczos, lanczos%half_width_lanczos, lanczos%cutoff_lanczos, &
	tmp3(:ns), 0.5d0)

! For Lanczos filter: set output values to NaN when input is NaN
! where (isnan_(tmp2(:ns))) tmp3(:ns) = nan

! Now reduce again to the original set of times

iono_out(:n) = tmp3(idx(:n))

end subroutine process_cycle_lanczos

!*  median_filter - median filter on array of values with window counted by index
!
subroutine median_filter (values, min_values, half_width, filtered_values)
use rads
use rads_misc
real(eightbytereal), intent(in) :: values(:)
integer(fourbyteint), intent(in) :: min_values, half_width
real(eightbytereal), intent(out) :: filtered_values(:)
!
! values          : array of values to be filtered
! min_values      : minimum number of valid values per filtering window
! half_width      : half width of the filtering window (unit = index numbers)
! filtered_values : array of filtered values
!-
integer(fourbyteint) :: nr_values, nw, nv, i, k
type(quicksort_pair), allocatable :: pairs(:)

! Initialise

nr_values = size(values)
filtered_values = nan

! Escape if period is too small

if (nr_values < min_values) return

! Determine (half) width of the window

nw = min(half_width, (nr_values - 1) / 2)

allocate (pairs(2*nw+1))

! Shift the window from left to right to compute median

do i = 1, nr_values
	nv = 0
	! Count and copy only the valid values
	do k = max(1,i-nw), min(i+nw,nr_values)
		if (isnan_(values(k))) cycle ! Skip NaNs
		nv = nv + 1
		pairs(nv)%value = values(k)
	enddo
	if (nv < min_values) cycle ! Skip when no median to create
	! Quicksort the values
	call quicksort (pairs(1:nv))
	! Determine the median
	if (modulo(nv,2) == 0) then
		filtered_values(i) = sum(pairs(nv/2:nv/2+1)%value)/2
	else
		filtered_values(i) = pairs((nv+1)/2)%value
	endif
enddo

deallocate (pairs)

end subroutine median_filter

!* lanczos_filter - Low-pass Lanczos filter on array of values with window counted by index
!
subroutine lanczos_filter (values, min_values, half_width, cutoff_period, filtered_values, min_weight)
use rads
use rads_misc
real(eightbytereal), intent(in) :: values(:)
integer(fourbyteint), intent(in) :: min_values, half_width, cutoff_period
real(eightbytereal), intent(out) :: filtered_values(:)
real(eightbytereal), optional, intent(in) :: min_weight
!
! values          : array of values to be filtered
! min_values      : minimum number of valid values per filtering window
! half_width      : half width of the filtering window (unit = index numbers)
! cutoff_period   : filter cutoff period (unit = index numbers)
! filtered_values : array of filtered values
! min_weight      : (optional) minimum sum of Lanczos weights
!-
real(eightbytereal) :: sum_coeffs, sum_lanczos, freq, x, min_coeffs
integer(fourbyteint) :: nr_values, nw, nv, dx, i, k
real(eightbytereal), allocatable :: coeffs(:)

! Initialise

nr_values = size(values)
filtered_values = nan
min_coeffs = 1d30
if (present(min_weight)) min_coeffs = min_weight

! Escape if period is too small

if (nr_values < min_values) return

! Determine (half) width of the window and cuttoff frequency

nw = min(half_width, (nr_values - 1) / 2)
freq = 1d0 / cutoff_period

! Determine the weighting coefficients

allocate (coeffs(0:nw))

do i = 0, nw
	if (i == 0) then
		coeffs(0) = 2 * freq
	else
		x = pi * i
		coeffs(i) = nw * sin (2 * x * freq) * sin (x / nw) / (x * x)
	endif
enddo

! Compute the Lanczos filtered values, looping over all indices of the (filtered) values arrays

do i = 1, nr_values
	nv = 0
	sum_coeffs = 0d0
	sum_lanczos = 0d0
	do k = max(1,i-nw), min(i+nw,nr_values)
		dx = abs(k - i)
		if (isnan_(values(k))) cycle ! Skip NaN
		nv = nv + 1
		sum_lanczos = sum_lanczos + values(k) * coeffs(dx)
		sum_coeffs = sum_coeffs + coeffs(dx)
	enddo
	if (sum_coeffs >= min_coeffs .and. nv >= min_values) filtered_values(i) = sum_lanczos / sum_coeffs
!	write (*,*) i, values(i), filtered_values(i), sum_lanczos, sum_coeffs
enddo

deallocate (coeffs)

end subroutine lanczos_filter

!-----------------------------------------------------------------------
! Write a single pass
!-----------------------------------------------------------------------

subroutine write_pass (n)
use netcdf
use rads_netcdf
integer(fourbyteint), intent(in) :: n
logical :: skip
integer :: ncid

call log_pass (P)

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! If "new" option is used, skip writing iono_alt_smooth when it already exists

ncid = P%fileinfo(1)%ncid
skip = new .and. nff(nf90_inq_varid(ncid,'iono_alt_smooth' // ext,i))

! Define to output variable and write out the data

if (recompute) call rads_def_var (S, P, iono_alt)
if (.not.skip) call rads_def_var (S, P, 'iono_alt_smooth' // ext)
if (recompute) call rads_put_var (S, P, iono_alt, iono_in(ntot+1:ntot+n))
if (.not.skip) call rads_put_var (S, P, 'iono_alt_smooth' // ext, iono_out(ntot+1:ntot+n))
if (skip.and..not.recompute) then
	call log_records(0)
else
	ntot = ntot+n
	call log_records (n)
endif
end subroutine write_pass

end program rads_add_dual
