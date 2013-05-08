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

!*rads_fix_sa -- Patch RADS altimeter files of SARAL for various anomalies
!
! This program makes numerous patches to the SARAL RADS data processed
! by rads_gen_saral. These patches include:
!
! ssb:
! - Change SSB to 3.5% of SWH
!
! usage: rads_fix_sa [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_sa

use rads
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

integer(fourbyteint) :: i,cyc,pass
logical :: lssb

! Scan command line for options

call synopsis
call rads_set_options (' ssb all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('ssb')
		lssb = .true.
	case ('all')
		lssb = .true.
	end select
enddo

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
if (rads_version ('$Revision$', 'Patch SARAL data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --ssb                     Change SSB to 3.5% of SWH' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: swh(n),ssb(n)

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename)

call rads_get_var (S, P, 'swh_ku', swh, .true.)

! Process data records

ssb = 3.5d-2 * swh

! If nothing changed, stop here

if (.not.lssb) then
	write (*,552) 0
	return
endif

! Write out all the data

call rads_put_history (S, P)
call rads_put_var (S, P, 'ssb_ku', ssb)

write (*,552) n
end subroutine process_pass

end program rads_fix_sa
