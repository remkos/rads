!-----------------------------------------------------------------------
! Copyright (c) 2011-2026  Remko Scharroo
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

!*rads_add_sla -- Add precomputed sea level anomaly
!+
! This program adds the sea surface height anomaly, as computed based on the
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
integer(fourbyteint) :: cyc, pass, i
logical :: update = .false.

! Initialise

call synopsis ('--head')
call rads_set_options ('mux:: multi-hz update ext:: all')
call rads_init (S)

! Check all options
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('u', 'update')
		update = .true.
	case ('m', 'multi-hz')
		S%n_hz_output = .true.
	case ('x', 'ext')
		if (rads_opt(i)%arg == '') then
			call rads_parse_varlist (S, 'ssha')
		else
			call rads_parse_varlist (S, 'ssha_' // rads_opt(i)%arg)
		endif
	end select
enddo
! Default to adding 'ssha' only
if (S%nsel == 0) call rads_parse_varlist (S, 'ssha')

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata, S%nsel)
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Add precomputed sea level anomaly field to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  --all                     (Has no effect)'/ &
'  -V, --var=VAR1[,VAR2,...] Add specified variables (e.g. ssha,ssha_mle3)'/ &
'  -x, --ext EXT             Produce field ssha_EXT (e.g. "-x -x mle3" for "-Vssha,ssha_mle3")' / &
'  -m, --multi-hz            Do multi-Hertz SLA (only)'/ &
'  -u, --update              Update files only when there are changes')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n, m)
integer(fourbyteint), intent(in) :: n, m
integer(fourbyteint) :: i, k
real(eightbytereal) :: sla(n,m), tmp(n)
real(eightbytereal), parameter :: dz = 1d-4
logical :: found_changes

call log_pass (P)

found_changes = .not.update ! Look for changes only when update is set; otherwise we simply assume there are changes

! Loop through all variable names
! Note: -V will specify the output ssha name (e.g. ssha_mle3), while the input will be e.g. sla_mle3

do k = 1,m

! Get sea level anomaly

	 call rads_get_var (S, P, 'sla'//S%sel(k)%name(5:), sla(:,k))

	! If no changes were found yet, check for changes first

	if (.not.found_changes) then
		i = rads_verbose; rads_verbose = -1 ! Temporarily suspend warning
		call rads_get_var (S, P, S%sel(k), tmp, .true.)
		rads_verbose = i
		do i = 1,n
			if (isnan_(tmp(i)) .and. isnan_(sla(i,k))) cycle
			if (isnan_(tmp(i))) exit
			if (nint(tmp(i)/dz) /= nint(sla(i,k)/dz)) exit
		enddo
		if (i <= n) found_changes = .true.	! There are changes
	endif
enddo

! If no changes, we bail

if (.not.found_changes) then
	call log_records (0)
	return
endif

! Store all data fields

call rads_put_history (S, P)

do k = 1,m
	call rads_def_var (S, P, S%sel(k))
enddo
if (S%n_hz_output) then
	i = n/P%n_hz
	do k = 1,m
		call rads_put_var (S, P, S%sel(k), reshape(sla(:,k), (/P%n_hz,i/)))
	enddo
	call log_records (i)
else
	do k = 1,m
		call rads_put_var (S, P, S%sel(k), sla(:,k))
	enddo
	call log_records (n)
endif

end subroutine process_pass

end program rads_add_sla
