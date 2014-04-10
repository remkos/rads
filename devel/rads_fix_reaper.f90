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
! ptr (for versions up to COM3 only):
! - Correct range for erroneous PTR
!
! tbias (for all versions COM*):
! - Add 0.68 ms to the time tags
!
! tide (for all versions COM*):
! - Prior to COM5: remove load tide from ocean tide
! - From COM5 onward: add long-period tide to ocean tide
!
! uso (for versions up to COM5 only):
! - Correct range for USO drift
!
! usage: rads_fix_reaper [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_reaper

use rads
use rads_misc
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

character(len=rads_cmdl) :: path
integer(fourbyteint) :: i, cyc, pass
real(eightbytereal) :: time_ptr, drange_ptr
logical :: lptr = .false., ltbias = .false., ltide = .false., luso = .false., lwet = .false.
real(eightbytereal), parameter :: tbias = 0.68d-3

! Scan command line for options

call synopsis ('--head')
call rads_set_options ('bpstuw ptr tbias tide uso wet all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('p', 'ptr')
		lptr = .true.
	case ('b', 'tbias')
		ltbias = .true.
	case ('t', 'tide')
		ltide = .true.
	case ('u', 'uso')
		luso = .true.
	case ('w', 'wet')
		lwet = .true.
	case ('all')
		lptr = .true.
		ltbias = .true.
		ltide = .true.
		luso = .true.
		lwet = .true.
	end select
enddo

! Load PTR data

if (lptr) then
	call parseenv ('${RADSROOT}/ext/reaper/commissioning/diff_ptrolc_ers'//S%sat(2:2)//'.dat', path)
	open (10,file=path,status='old')
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
'  -b, --tbias               Correct time and altitude for timing bias' / &
'  -t, --tide                Fix the tide, subtracting load tide or adding long period tide' / &
'  -u, --uso                 Correct range for USO drift' / &
'  -w, --wet                 Fix radiometer wet tropo correction, taking into account TB biases' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
use meteo_subs
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n), range_ku(n), drange_uso(n), got47(n), fes04(n), a(n), alt_reaper(n), &
	wet(n), sig0(n), dsig0(n), tb_238(n), tb_365(n), d_tb_238, d_tb_365, d_sig0
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

! Fix wet tropospheric correction

if (lwet) then
	call rads_get_var (S, P, 'tb_238', tb_238)
	call rads_get_var (S, P, 'tb_365', tb_365)
	call rads_get_var (S, P, 'sig0_ku', sig0)
	call rads_get_var (S, P, 'dsig0_atmos_ku', dsig0)
	if (S%sat == 'e1') then
		d_tb_238 = -5.50d0; d_tb_365 = -5.10d0; d_sig0 = -0.76d0
	else
		d_tb_238 = -5.16d0; d_tb_365 = -4.28d0; d_sig0 =  0.06d0
	endif
	do i = 1,n
		wet(i) = nn_l2_mwr (tb_238(i) + d_tb_238, tb_365(i) + d_tb_365, sig0(i)-dsig0(i)+d_sig0, 3)
	enddo
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
if (ltbias) then
	call rads_put_var (S, P, 'time', time)
	call rads_put_var (S, P, 'alt_reaper', alt_reaper)
endif
if (lwet) call rads_put_var (S, P, 'wet_tropo_rad', wet)

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
