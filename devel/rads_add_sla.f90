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
use rads_misc

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cyc, pass, j
logical :: update = .false.

! Initialise

call synopsis ('--head')
call rads_set_options ('mu multi-hz update all')
call rads_init (S)

! Check all options
do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('u', 'update')
		update = .true.
	case ('m', 'multi-hz')
		S%n_hz_output = .true.
	end select
enddo

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
'  --all                     (Has no effect)'/ &
'  -m, --multi-hz            Do multi-Hertz SLA (only)'/ &
'  -u, --update              Update files only when there are changes')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i
real(eightbytereal) :: sla(n), tmp(n)
real(eightbytereal), parameter :: dz = 1d-4
character(len=5) :: hz

call log_pass (P)

! Get sea level anomaly

if (S%n_hz_output) then
	write (hz, '("_",i2.2,"hz")') P%n_hz
else
	hz = ''
endif
call rads_get_var (S, P, 'sla'//hz, sla)

! If requested, check for changes first

if (update) then
	i = rads_verbose; rads_verbose = -1 ! Temporarily suspend warning
	call rads_get_var (S, P, 'ssha'//hz, tmp, .true.)
	rads_verbose = i
	do i = 1,n
		if (isnan_(tmp(i)) .and. isnan_(sla(i))) cycle
		if (isnan_(tmp(i))) exit
		if (nint(tmp(i)/dz) /= nint(sla(i)/dz)) exit
	enddo
	if (i > n) then	! No changes
		call log_records (0)
		return
	endif
endif

! Store all data fields

call rads_put_history (S, P)

call rads_def_var (S, P, 'ssha'//hz)
if (S%n_hz_output) then
	i = n/P%n_hz
	call rads_put_var (S, P, 'ssha'//hz, reshape(sla, (/P%n_hz,i/)))
	call log_records (i)
else
	call rads_put_var (S, P, 'ssha', sla)
	call log_records (n)
endif

end subroutine process_pass

end program rads_add_sla
