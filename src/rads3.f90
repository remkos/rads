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

module rads3
!
! This module provides backward compatibility of RADS3 code with the RADS4 library.
!
! The following RADS3 routines are implemented herein:
!
! subroutine getraw_init (mission, verbose)
! - Fully supported
!
! subroutine getraw (cycle, pass, ncols, maxdata, select, time, dlat, dlon, &
!   data, ndata, eqtime, eqlon, metafile)
! - Fully supported
!
! subroutine getraw_options (select, option)
! - Does not support option=0
! - Option < 0 does not remove a field from the computation of SLA
! - Option > 0 does not add a field to the computation of SLA
! - It does, however, allow modification of an alias in the computation of SLA
!
! subroutine getraw_limits (select, lo, hi)
! - Adjust only the limits for the selected field, i.e. select=901 adjust the
!   limits on dual-frequency iono, select=9 adjust the limits on the current
!   default iono correction.
!
! subroutine getraw_factors (select, fact)
! - Not supported. Had no effect on the computation of SLA
!
! subroutine getraw_stat (unit)
! - RADS4 library provides different output from RADS3
!
! function radsargs1 (usage, sat, verbose)
! - Runs rads_init based on command line arguments and prevends getraw_init from
!   running a second time
! - Provides all common synopsis info, not just for sat=, debug=, -v
!
! function radsargs2 (usage, cyc0, cyc1, pass0, pass1, dpass, nsel, sel)
! - Provides all common synopsis info, including for sat=, debug=, -v
!
! Using this module in an older RADS3 program should require only the addition
! of the line 'use rads3' before the first executable line.
!-----------------------------------------------------------------------
use rads, only: rads_sat, rads_pass
type(rads_sat), save, private :: S
type(rads_pass), save, private :: P
logical, save, private :: rads_init_done = .false.

contains

subroutine getraw_init (mission, verbose)
use rads
character(len=*), intent(in) :: mission
integer(fourbyteint), intent(in) :: verbose
rads_verbose = verbose
if (.not.rads_init_done) call rads_init (S, mission)
rads_init_done = .true.
end subroutine getraw_init

subroutine getraw (cycle, pass, ncols, maxdata, select, time, dlat, dlon, &
	data, ndata, eqtime, eqlon, metafile)
use rads
integer(fourbyteint), intent(in) :: cycle, pass, ncols, maxdata, select(ncols)
integer(fourbyteint), intent(out) :: ndata
real(eightbytereal), intent(out) :: time(maxdata), dlat(maxdata), dlon(maxdata), data(maxdata,ncols), eqtime, eqlon
character(len=*), intent(out) :: metafile
integer :: i
if (.not.rads_init_done) then
	call rads_error (S, rads_err_noinit, 'Need to call getraw_init before getraw')
	return
endif
call rads_open_pass (S, P, cycle, pass)
ndata = P%ndata
if (ndata <= 0) then
else if (ndata > maxdata) then
	call rads_error (S, rads_err_memory, 'Too little memory allocated to read data from file', P)
else
	time(1:ndata) = P%tll(:,1)
	dlat(1:ndata) = P%tll(:,2)
	dlon(1:ndata) = P%tll(:,3)
	do i = 1,ncols
		call rads_get_var (S, P, select(i), data(1:ndata,i))
	enddo
	eqlon = P%equator_lon
	eqtime = P%equator_time
	metafile = P%fileinfo(1)%name
endif
call rads_close_pass (S, P)
end subroutine getraw

subroutine getraw_options (select, option)
use rads
integer(fourbyteint), intent(in) :: select, option
character(len=4) :: name
type(rads_var), pointer :: src, sla
write (name, '(i4.4)') select*100 + abs(option)
if (option == 0) then
	call rads_error (S, rads_err_incompat, 'getraw_options('//name(1:2)//',00) has no effect')
	return
endif
src => rads_varptr (S, name(1:2))
sla => rads_varptr (S, 'sla')
if (option > 0) then
	if (.not.associated(src) .or. index(sla%info%dataname,trim(src%name)) <= 0) &
		call rads_error (S, rads_err_incompat, 'getraw_options('//name(1:2)//','//name(3:4)//') will not affect SLA computation')
else
	if (associated(src) .and. index(sla%info%dataname,trim(src%name)) > 0) &
		call rads_error (S, rads_err_incompat, 'getraw_options('//name(1:2)//',-'//name(3:4)//') will not affect SLA computation')
endif
if (associated(src)) call rads_set_alias (S, src%name, name)
end subroutine getraw_options

subroutine getraw_limits (select, lo, hi)
use rads
integer(fourbyteint), intent(in) :: select
real(eightbytereal), intent(in) :: lo, hi
character(len=4) :: name
type(rads_var), pointer :: src
write (name, '(i4.4)') select
src => rads_varptr (S, name)
if (.not.associated(src)) return
src%info%limits = (/lo, hi/)
end subroutine getraw_limits

subroutine getraw_factors (select, fact)
use rads
integer(fourbyteint), intent(in) :: select
real(eightbytereal), intent(in) :: fact
real(eightbytereal) :: x
x = fact + select ! To avoid warning about unused dummy arguments
call rads_error (S, rads_err_incompat, &
	'getraw_factor() has no effect. Change the <var name="sla"> entry in the appropriate XML file instead')
end subroutine getraw_factors

subroutine getraw_stat (unit)
use rads
integer(fourbyteint), intent(in) :: unit
rads_log_unit = unit
call rads_stat (S)
end subroutine getraw_stat

function radsargs1 (usage, sat, verbose)
use rads
integer(fourbyteint), intent(in) :: usage
character(len=*), intent(out) :: sat
integer(fourbyteint), intent(out) :: verbose
logical :: radsargs1
radsargs1 = .false.
if (usage == 3) then
	call rads_synopsis
	return
endif
if (.not.rads_init_done) call rads_init (S)
rads_init_done = .true.
radsargs1 = (S%error /= rads_noerr)
if ((radsargs1 .and. usage == 1) .or. usage == 2) call rads_synopsis
sat = S%sat//'/'//trim(S%phase%name)
verbose = rads_verbose
end function radsargs1

function radsargs2 (usage, cyc0, cyc1, pass0, pass1, dpass, nsel, sel)
use rads
integer(fourbyteint), intent(in) :: usage
integer(fourbyteint), intent(inout) :: nsel
integer(fourbyteint), intent(out) :: cyc0, cyc1, pass0, pass1, dpass, sel(nsel)
logical :: radsargs2
radsargs2 = .false.
if (usage == 3) then
	call rads_synopsis
	return
endif
if (S%nsel > nsel) then
	call rads_error (S, rads_err_memory, 'Too many data selectors on sel= argument')
	radsargs2 = .true.
	if (usage == 1) call rads_synopsis
endif
if (usage == 2) call rads_synopsis
cyc0 = S%cycles(1)
cyc1 = S%cycles(1)
pass0 = S%passes(1)
pass1 = S%passes(2)
dpass = S%passes(3)
nsel = S%nsel
sel(1:nsel) = S%sel(1:nsel)%field(1)
end function radsargs2

end module rads3
