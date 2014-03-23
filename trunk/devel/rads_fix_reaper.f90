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

!*rads_fix_reaper -- Patch RADS altimeter files of ERS/REAPER for various anomalies
!
! This program makes numerous patches to the ERS/REAPER RADS data processed
! by rads_gen_reaper. These patches include:
!
! ptr:
! - Correct range for erroneous PTR
!
! usage: rads_fix_reaper [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_reaper

use rads
use rads_misc
use rads_devel
use rads_netcdf

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

character(len=rads_cmdl) :: path
integer(fourbyteint) :: i, cyc, pass
real(eightbytereal) :: time_ptr, drange_ptr
type(grid) :: issb_hyb
logical :: lptr = .false., lssb = .false., ltbias = .false., ltide = .false., luso = .false.
real(eightbytereal), parameter :: tbias = 0.68d-3

! Scan command line for options

call synopsis ('--head')
call rads_set_options ('bpstu ptr ssb tbias tide uso all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('p', 'ptr')
		lptr = .true.
	case ('s', 'ssb')
		lssb = .true.
	case ('b', 'tbias')
		ltbias = .true.
	case ('t', 'tide')
		ltide = .true.
	case ('u', 'uso')
		luso = .true.
	case ('all')
		lptr = .true.
		lssb = .true.
		ltbias = .true.
		ltide = .true.
		luso = .true.
	end select
enddo

! Load PTR data

if (lptr) then
	call parseenv ('${RADSROOT}/ext/reaper/commissioning/diff_ptrolc_ers'//S%sat(2:2)//'.dat', path)
	open (10,file=path,status='old')
endif

! Load SSB model

if (lssb) then
	call parseenv ('${ALTIM}/data/models/reaper_ssb_hyb.nc?ssb_hyb', path)
	if (grid_load(path,issb_hyb) /= 0) call rads_exit ('Error loading '//trim(path))
endif

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

if (lptr) close (10)
if (lssb) call grid_free (issb_hyb)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Patch ERS/REAPER data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  -p, --ptr                 Correct range for PTR error (use for pre-COM5 only)' / &
'  -s, --ssb                 Add hybrid SSB model' / &
'  -b, --tbias               Correct time and altitude for timing bias' / &
'  -t, --tide                Fix the tide, subtracting load tide or adding long period tide' / &
'  -u, --uso                 Correct range for USO drift' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: x, y, time(n), range_ku(n), drange_uso(n), sig0(n), swh(n), ssb(n), got47(n), fes04(n), &
	a(n), alt_reaper(n)
integer :: n_changed, i, com

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename(len_trim(S%dataroot)+2:))

n_changed = 0
com = 999
i = index(P%original, '_COM')
if (i > 0) read (P%original(i+4:i+4), *) com

if (lptr .or. ltbias) call rads_get_var (S, P, 'time', time, .true.)
if (lptr .or. luso) call rads_get_var (S, P, 'range_ku', range_ku, .true.)

! Apply PTR correction

if (lptr) then
	do i = 1,n
		do while (time_ptr + 0.5d0 < time(i))
			call next_ptr
		enddo
		if (time_ptr - 0.5d0 > time(i)) cycle
		n_changed = n_changed + 1
		range_ku(i) = range_ku(i) - drange_ptr
	enddo
endif

! Compute SSB

if (lssb) then
	call rads_get_var (S, P, 'sig0_ku', sig0, .true.)
	call rads_get_var (S, P, 'swh_ku', swh, .true.)
	do i = 1,n
		x = sig0(i)
		if (x < issb_hyb%xmin) x = issb_hyb%xmin
		if (x > issb_hyb%xmax) x = issb_hyb%xmax
		y = swh(i)
		if (y < issb_hyb%ymin) y = issb_hyb%ymin
		if (y > issb_hyb%ymax) y = issb_hyb%ymax
		ssb(i) = grid_lininter (issb_hyb, x, y)
	enddo
	n_changed = n
endif

! Fix the tide

if (ltide) then
	call rads_get_var (S, P, 'tide_ocean_got47', got47)
	call rads_get_var (S, P, 'tide_ocean_fes04', fes04)
	if (com < 5) then	! Prior to COM5: remove load tide from ocean tide
		call rads_get_var (S, P, 'tide_load_got47', a)
		got47 = got47 - a
		call rads_get_var (S, P, 'tide_load_fes04', a)
		fes04 = fes04 - a
	else	! From COM5 onward: add long-period tide to ocean tide
		call rads_get_var (S, P, 'tide_equil', a)
		got47 = got47 + a
		fes04 = fes04 + a
	endif
	n_changed = n
endif

! Apply USO correction

if (luso) then
	call rads_get_var (S, P, 'drange_uso', drange_uso)
	range_ku = range_ku + drange_uso
	n_changed = n
endif

! Apply time tag bias

if (ltbias) then
	call rads_get_var (S, P, 'alt_reaper', alt_reaper)
	call rads_get_var (S, P, 'alt_rate', a)
	time = time + tbias
	alt_reaper = alt_reaper + tbias * a
	n_changed = n
endif

! If nothing changed, stop here

if (n_changed == 0) then
	write (*,552) 0
	return
endif

! Write out all the data

call rads_put_history (S, P)
if (luso .or. lptr) call rads_put_var (S, P, 'range_ku', range_ku)
if (ltide) then
	call rads_put_var (S, P, 'tide_ocean_got47', got47)
	call rads_put_var (S, P, 'tide_ocean_fes04', fes04)
endif
if (lssb) then
	call rads_def_var (S, P, 'ssb_hyb')
	call rads_put_var (S, P, 'ssb_hyb', ssb)
endif
if (ltbias) then
	call rads_put_var (S, P, 'time', time)
	call rads_put_var (S, P, 'alt_reaper', alt_reaper)
endif

write (*,552) n_changed

end subroutine process_pass

!-----------------------------------------------------------------------
! Get next PTR value
!-----------------------------------------------------------------------

subroutine next_ptr
character(len=80) :: line
integer :: ios, n
real(eightbytereal) :: corr1, corr2
do
	read (10,'(a)',iostat=ios) line
	if (ios /= 0) then
		time_ptr = 1d20
		exit
	endif
	if (line(1:1) == '#') cycle
	read (line,*) time_ptr, n, corr1, corr2, drange_ptr
	exit
enddo
end subroutine next_ptr

end program rads_fix_reaper
