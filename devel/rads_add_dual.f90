!-----------------------------------------------------------------------
! Copyright (c) 2011-2019  Remko Scharroo
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
! The default smoothing length is 35 seconds, in contrast to the 21 seconds
! suggested in Imel [1994].
!
! References:
! D. A. Imel, Evaluation of the TOPEX/Poseidon dual-frequency ionosphere
! correction, J. Geophys. Res., 99, C12, 24,895-24,906, 1994.
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

integer(fourbyteint) :: ntot, cyc, pass, i
integer(fourbyteint), parameter :: nmax = 3000000
integer(twobyteint) :: mask = 2072 ! Bits 3, 4, 11
character(len=5) :: ext = ''
real(eightbytereal) :: twin = 35d0, iwin = 8d0, f
real(eightbytereal) :: time(nmax), lat(nmax), flags(nmax), iono1(nmax), iono2(nmax)
logical :: recompute = .false., new = .false.

! Scan command line for options

call synopsis ('--head')
call rads_set_options ('nb:i:t:m:x:r new mask: iwin: twin: mle: ext: recompute')
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
		ext = '_' // rads_opt(i)%arg(:4)
	case ('i', 'iwin')
		read (rads_opt(i)%arg,*) iwin
	case ('b', 'mask')
		read (rads_opt(i)%arg,*) mask
	case ('t', 'twin')
		read (rads_opt(i)%arg,*) twin
	case ('r', 'recompute')
		recompute = .true.
	end select
enddo

! Use half the time window (in seconds) and half the range window (converted to meters)
! in the rest of the program

twin = twin/2d0
iwin = iwin/2d2 + 1d-6

! Run process for all files within a cycle, first reading all, then writing all

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	ntot = 0
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call read_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
	call process_cycle (ntot)
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
write (*,1310) mask, iwin, twin
1310 format (/ &
'Additional [processing_options] are:' / &
'  -n, --new                 Only add smoothed iono when not yet existing' / &
'  -b, --mask MASK           Exclude data based on bitmap MASK (default:',i5,')' / &
'  -i, --iwin IWIN           Set editing range for iono data [cm] (default:',f4.0,')' / &
'  -t, --twin TWIN           Set box car filter length [sec] (default:' f4.0,')' / &
'  -r, --recompute           (Re)compute ionospheric correction from delta range' / &
'  -x, --ext EXT             Use parameters with extension _EXT (e.g. mle3 or plrm)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Read a single pass
!-----------------------------------------------------------------------

subroutine read_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: range_ku(n), range_c(n), ssb_ku(n), ssb_c(n)

! Read all the required variables

call rads_get_var (S, P, 'time', time(ntot+1:ntot+n), .true.)
call rads_get_var (S, P, 'lat', lat(ntot+1:ntot+n), .true.)
call rads_get_var (S, P, 'flags' // ext, flags(ntot+1:ntot+n), .true.)
call rads_get_var (S, P, 'iono_alt' // ext, iono1(ntot+1:ntot+n), .true.)
if (recompute) then
	call rads_get_var (S, P, 'range_ku' // ext, range_ku, .true.)
	call rads_get_var (S, P, 'range_c', range_c, .true.)
	call rads_get_var (S, P, 'ssb_cls' // ext, ssb_ku, .true.)
	call rads_get_var (S, P, 'ssb_cls_c', ssb_c, .true.)
	iono1(ntot+1:ntot+n) = f * ((range_c + ssb_c) - (range_ku + ssb_ku))
endif
ntot = ntot + n
end subroutine read_pass

!-----------------------------------------------------------------------
! Process an entire cycle
!-----------------------------------------------------------------------

subroutine process_cycle (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, j, k, j0, j1, nsum
integer(twobyteint) :: flag
real(eightbytereal) :: t0, t1, wsum

! If any mask is specified: exclude data based on bit mask

if (mask /= 0) then
	do i = 1,n
		flag = nint(flags(i),twobyteint)
		if (iand(flag,mask) /= 0) iono1(i) = nan
	enddo
endif

! Do box car filter for each measurement (i).

j0 = 1
j1 = 1
do i = 1,n

! Find first measurement in window (starts at t0)

	t0 = time(i) - twin
	do while (time(j0) < t0)
		j0 = j0 + 1
	enddo

! Find last measurement in window (ends at t1)

	t1 = time(i) + twin
	j1 = max(i,j1)
	do while (j1 < n .and. time(j1+1) < t1)
		j1 = j1 + 1
	enddo

! Now do the box car filter to compute average
! First time with only rough outlier editing
! Second time with editing based on +/-iwin
! (This is a slight change from radsp_dual, because here in the second loop
! a iono1 > iwin could be accepted. In radsp_dual those were permanently rejected.)

	t0 = -1d0
	t1 = iwin
	do k = 1,2
		nsum = 0
		wsum = 0d0
		do j = j0,j1
			if (iono1(j) >= t0 .and. iono1(j) <= t1) then
				nsum = nsum + 1
				wsum = wsum + iono1(j)
			endif
		enddo
		if (nsum == 0) exit
		wsum = wsum / nsum
		t0 = wsum - iwin
		t1 = wsum + iwin ! Maybe this should be changed to min(wsum+iwin,iwin) conform radsp_dual
	enddo

! Assign NaN in case nsum < 4

	if (nsum < 4) then
		iono2(i) = nan
	else
		iono2(i) = wsum
	endif
	! write (*,'(f14.3,i6,3f9.4,i3)') time(i),nint(flags(i)),lat(i),iono1(i),iono2(i),nsum

enddo

end subroutine process_cycle

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

if (recompute) call rads_def_var (S, P, 'iono_alt' // ext)
if (.not.skip) call rads_def_var (S, P, 'iono_alt_smooth' // ext)
if (recompute) call rads_put_var (S, P, 'iono_alt' // ext, iono1(ntot+1:ntot+n))
if (.not.skip) call rads_put_var (S, P, 'iono_alt_smooth' // ext, iono2(ntot+1:ntot+n))
if (skip.and..not.recompute) then
	call log_records(0)
else
	ntot = ntot+n
	call log_records (n)
endif
end subroutine write_pass

end program rads_add_dual
