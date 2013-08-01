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

!*rads_add_sla -- Add SST temperature and ice concentration to RADS data
!+
! This program adds the sea level anomaly, as computed based on the
! standard rules to construct the sea level anomaly, including edit
! criteria.
!
! This field might thus NOT be the same as one would get by selecting
! --var=sla, since that depends on criteria set by the user.
!
! usage: rads_add_sla [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_sla

use rads
use rads_devel

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cyc, pass

! Initialise

call synopsis ('--head')
call rads_set_options (' all')
call rads_init (S)

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
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
if (rads_version ('$Revision$', 'Add precomputed sea level anomaly field to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  --all                     (Has no effect)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: sla(n)

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename(len_trim(S%dataroot)+2:))

! Get sea level anomaly

call rads_get_var (S, P, 'sla', sla)

! Store all data fields

call rads_put_history (S, P)

call rads_def_var (S, P, 'ssha')
call rads_put_var (S, P, 'ssha', sla)

write (*,552) n
end subroutine process_pass

end program rads_add_sla
