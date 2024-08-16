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

!*rads_fix_n1 -- Patch RADS altimeter files of Sentinel-3 for various anomalies
!
! This program makes numerous patches to the Envisat RADS data processed
! by rads_gen_n1_gdr. These patches include:
!
!  --sideB                   Set side B flag
!  --biasS                   Apply S-band range bias
!  --all                     All the above
!
! usage: rads_fix_s3 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_n1

use rads
use rads_misc
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

logical :: l_sideB = .false., l_biasS = .false.
integer(fourbyteint) :: cyc, pass, i

! Biases. All need to be ADDED to the measurements
! range_s_bias:    S-band range bias for BOTH sides

real(eightbytereal), parameter :: range_s_bias = 0.1734d0
real(eightbytereal) :: f

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' sideB biasS all')
call rads_init (S)
if (S%sat /= "n1") stop

do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('sideB')
		l_sideB = .true.
	case ('biasS')
		l_biasS = .true.
	case ('all')
		l_sideB = .true.
		l_biasS = .true.
	end select
enddo

if (.not.l_sideB .and. .not.l_biasS) stop

! Determine conversion factor from range difference to ionospheric correction

f = 1d0/(1d0-(S%frequency(1)/S%frequency(2))**2)

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Patch Envisat data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --sideB                   Set side B flag' // &
'  --biasS                   Apply S-band range bias' // &
'  --all                     All the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
use rads_time
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: flags(n), range_s(n), iono_alt(n), iono_alt_smooth(n)
integer(fourbyteint) :: i
logical :: do_biasS, do_sideB

call log_pass (P)

! Switch S-band off for the whole Side-B period

do_sideB = .false.
do_biasS = l_biasS

i = cyc*10000 + pass
if (i >= 470794 .and. i <= 480847) then
	! Side B period
	do_biasS = .false.
	do_sideB = l_sideB
else if (i >= 650290) then
	! S-band failure in Side A on 17 Jan 2008 23:23:40
	do_biasS = .false.
endif

! Do this routine only for Baseline < 005

if (.not.do_sideB .and. .not.do_biasS) then
	call log_records(0)
	return
endif

! For Side B: set flag and Update Ku-band range

if (do_sideB) then
	call rads_get_var (S, P, 'flags', flags, .true.)
	do i = 1,n
		flags(i) = ibset(nint(flags(i)),0)
	enddo
endif
if (do_biasS) then
	call rads_get_var (S, P, 'range_s', range_s, .true.)
	call rads_get_var (S, P, 'iono_alt', iono_alt, .true.)
	call rads_get_var (S, P, 'iono_alt_smooth', iono_alt_smooth, .true.)
endif

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Write out all the data

if (do_sideB) call rads_put_var (S, P, 'flags', flags)
if (do_biasS) then
	call rads_put_var (S, P, 'range_s', range_s + range_s_bias)
	call rads_put_var (S, P, 'iono_alt', iono_alt + f * range_s_bias)
	call rads_put_var (S, P, 'iono_alt_smooth', iono_alt_smooth + f * range_s_bias)
endif

call log_records (n)
end subroutine process_pass

end program rads_fix_n1
