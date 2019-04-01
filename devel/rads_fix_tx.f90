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

!*rads_fix_tx -- Patch RADS altimeter files of TOPEX for various anomalies
!
! This program makes a patch to the RADS data for TOPEX
! processed by RADS3 software (tpgdrmraw and radsp_tp) from AVISO
! GDR-M products.
!
! This patch consists of the UNdoing of the Cal1 altimeter range
! correction, storing it for later use in a separate variable (drange_ku)
!
! References:
!-----------------------------------------------------------------------
program rads_fix_tx

use rads
use rads_misc
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

character(len=rads_cmdl) :: path, line
real(eightbytereal) :: x, bias(500)
integer(fourbyteint) :: cyc, pass

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' all')
call rads_init (S)

! Read bias table

call parseenv ('${ALTIM}/data/bias/RngStbUp.txt', path)
open (10, file=path, status='old')
do
	read (10,550) line
	if (line(1:10) == 'Cyc  Count') exit
enddo
do
	read (10,550) line
	if (line(1:6) == '</pre>') exit
	read (line,600) cyc, x
	x = x * 1d-3
	bias(cyc) = x
	bias(cyc+1) = x	! To also generate a value for missing (i.e. Poseidon) cycles
enddo
550 format (a)
600 format (i3,5x,f11.3)

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
if (rads_version ('Patch TOPEX data for range bias', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --all                     All (default)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: range_ku(n), drange_ku(n)
real(eightbytereal), parameter :: scale_factor = 1d-4, add_offset = 1300d3

call log_pass (P)

! Add back the range bias

call rads_get_var (S, P, 'range_ku', range_ku, .true.)
drange_ku = bias(cyc)
range_ku = range_ku + bias(cyc)

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Write out all the data

call rads_def_var (S, P, 'drange_ku', ndims=0)
call rads_def_var (S, P, 'range_ku', scale_factor=scale_factor, add_offset=add_offset)

call rads_put_var (S, P, 'range_ku', range_ku)
call rads_put_var (S, P, 'drange_ku', drange_ku)

call log_records (n)
end subroutine process_pass

end program rads_fix_tx
