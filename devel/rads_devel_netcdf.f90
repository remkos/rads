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

module rads_devel_netcdf
use typesizes
use rads, only: rads_sat, rads_pass, rads_var

integer(fourbyteint), parameter :: mrec=3700, mvar=140
integer(fourbyteint) :: nvar, nrec=0
real(eightbytereal), allocatable :: a(:)
integer(twobyteint), allocatable :: flags(:)
type :: var_
	type(rads_var), pointer :: v ! Pointer to rads_var struct
	real(eightbytereal) :: d(mrec) ! Data array
	logical :: empty, zero, skip_empty ! .true. if all NaN or all zero
endtype
type(var_) :: var(mvar)

type(rads_sat) :: S
type(rads_pass) :: P

!-----------------------------------------------------------------------
! Copy variable to RADS
!-----------------------------------------------------------------------
!
! subroutine cpy_var (ncid, varin, varout, do)
! use rads_devel
! use rads_netcdf
! integer(fourbyteint), intent(in) :: ncid
! character(len=*), intent(in) :: varin
! character(len=*), intent(in), optional :: varout
! logical, optional :: do
! Copy variable 'varin' from input (GDR) file to 'varout' in RADS output file.
! When 'varout' is omitted, varout=varin.
! When 'do' is present and false, the action is skipped.
private :: cpy_var_1, cpy_var_2
interface cpy_var
	module procedure cpy_var_1
	module procedure cpy_var_2
end interface cpy_var

contains

subroutine cpy_var_1 (ncid, varin, do, skip_empty)
use rads_netcdf
integer(fourbyteint), intent(in) :: ncid
character(len=*), intent(in) :: varin
logical, optional, intent(in) :: do, skip_empty
if (present(do) .and. .not.do) return
call get_var (ncid, varin, a)
call new_var (varin, a, skip_empty)
end subroutine cpy_var_1

subroutine cpy_var_2 (ncid, varin, varout, do, skip_empty)
use rads_netcdf
integer(fourbyteint), intent(in) :: ncid
character(len=*), intent(in) :: varin, varout
logical, optional, intent(in) :: do, skip_empty
if (present(do) .and. .not.do) return
call get_var (ncid, varin, a)
call new_var (varout, a, skip_empty)
end subroutine cpy_var_2

!-----------------------------------------------------------------------
! Create new RADS variable
!-----------------------------------------------------------------------

subroutine new_var (varnm, data, skip_empty)
use rads
! Write variables one after the other to the output file
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(in) :: data(:)
logical, optional, intent(in) :: skip_empty
nvar = nvar + 1
if (nvar > mvar) stop 'Too many variables allocated by new_var'
var(nvar)%v => rads_varptr (S, varnm)
var(nvar)%d(1:size(data)) = data
var(nvar)%skip_empty = (present(skip_empty) .and. skip_empty)
end subroutine new_var

!-----------------------------------------------------------------------
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads (skip_create)
use rads
use rads_misc
use rads_devel
logical, optional :: skip_create
integer(fourbyteint) :: i

! Check which variables are empty or all zero
do i = 1,nvar
	if (var(i)%skip_empty) then
		var(i)%skip_empty = all(isnan_(var(i)%d(1:nrec)))
		var(i)%empty = .false.
	else
		var(i)%empty = all(isnan_(var(i)%d(1:nrec)))
	endif
	var(i)%zero = all(var(i)%d(1:nrec) == 0d0)
enddo

! Write out the empty variables to be skipped
if (any(var(1:nvar)%skip_empty)) then
	write (rads_log_unit,551,advance='no') 'Skip empty:'
	do i = 1,nvar
		if (var(i)%skip_empty) write (rads_log_unit,551,advance='no') trim(var(i)%v%name)
	enddo
	write (rads_log_unit,551,advance='no') '...'
endif

! Write out the empty variables to be kept
if (any(var(1:nvar)%empty)) then
	write (rads_log_unit,551,advance='no') 'Empty:'
	do i = 1,nvar
		if (var(i)%empty) write (rads_log_unit,551,advance='no') trim(var(i)%v%name)
	enddo
	write (rads_log_unit,551,advance='no') '...'
endif
551 format (a,1x)

! Do the same for records that are all zero
if (any(var(1:nvar)%zero)) then
	write (rads_log_unit,551,advance='no') 'All zero:'
	do i = 1,nvar
		if (var(i)%zero) write (rads_log_unit,551,advance='no') trim(var(i)%v%name)
	enddo
	write (rads_log_unit,551,advance='no') '...'
endif

! Open output file
if (.not.(present(skip_create) .and. skip_create)) call rads_create_pass (S, P, nrec)

! Define all variables
do i = 1,nvar
	if (.not.var(i)%skip_empty) call rads_def_var (S, P, var(i)%v)
enddo

! Fill all the data fields
do i = 1,nvar
	if (.not.var(i)%skip_empty) call rads_put_var (S, P, var(i)%v, var(i)%d(1:nrec))
enddo

! Close the data file
call log_records (nrec, P)
call rads_close_pass (S, P)

end subroutine put_rads

!-----------------------------------------------------------------------
! nc2f  : Load flag field, then set corresponding bit in RADS
! ncid  : source NetCDF ID
! varnm : source variable name
! bit   : RADS bit to be set when value == 1
! eq    : set bit when value == val (optional, default = 1)
! neq   : set bit when value /= val (optional)
! ge    : set bit when value >= val (optional)
! le    : set bit when value <= val (optional)
! mask  : set bit when iand(value,mask) /= 0 (optional)

subroutine nc2f (ncid, varnm, bit, eq, neq, ge, le, mask)
use rads_netcdf
use netcdf
integer(fourbyteint), intent(in) :: ncid
character(len=*), intent(in) :: varnm
integer(fourbyteint), intent(in) :: bit
integer(fourbyteint), optional, intent(in) :: eq, neq, ge, le, mask
integer(twobyteint) :: flag(mrec), flag2d(1:1,1:nrec), mask2
integer(fourbyteint) :: i, ival, ndims, varid

if (nf90_inq_varid_warn(ncid,varnm,varid) /= nf90_noerr) return
call nfs(nf90_inquire_variable(ncid,varid,ndims=ndims))
if (ndims == 2) then	! Reduce 2D flags to 1D array
	call nfs(nf90_get_var(ncid,varid,flag2d(1:1,1:nrec)))
	flag(1:nrec) = flag2d(1,1:nrec)
else
	call nfs(nf90_get_var(ncid,varid,flag(1:nrec)))
endif
if (present(ge) .and. present(le)) then ! Set flag when ge <= value <= le
	do i = 1,nrec
		if (flag(i) >= ge .and. flag(i) <= le) flags(i) = ibset(flags(i),bit)
	enddo
else if (present(ge)) then	! Set flag when value >= ge
	do i = 1,nrec
		if (flag(i) >= ge) flags(i) = ibset(flags(i),bit)
	enddo
else if (present(le)) then	! Set flag when value <= le
	do i = 1,nrec
		if (flag(i) <= le) flags(i) = ibset(flags(i),bit)
	enddo
else if (present(neq)) then	! Set flag when value /= neq
	do i = 1,nrec
		if (flag(i) /= neq) flags(i) = ibset(flags(i),bit)
	enddo
else if (present(mask)) then	! Set flag when value and mask have common set bits
	mask2 = int(mask, twobyteint)
	do i = 1,nrec
		if (iand(flag(i), mask2) /= 0) flags(i) = ibset(flags(i),bit)
	enddo
else	! Set flag when value == 1, or when value == eq (when given)
	ival = 1
	if (present(eq)) ival = eq
	do i = 1,nrec
		if (flag(i) == ival) flags(i) = ibset(flags(i),bit)
	enddo
endif
end subroutine nc2f

end module
