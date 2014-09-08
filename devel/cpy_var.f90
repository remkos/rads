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

!-----------------------------------------------------------------------
! Copy variable to RADS
!-----------------------------------------------------------------------

subroutine cpy_var (varin, varout)
character(len=*), intent(in) :: varin
character(len=*), intent(in), optional :: varout
call get_var (ncid, varin, a)
if (present(varout)) then
	call new_var (varout, a)
else
	call new_var (varin, a)
endif
end subroutine cpy_var

!-----------------------------------------------------------------------
! Create new RADS variable
!-----------------------------------------------------------------------

subroutine new_var (varnm, data)
! Write variables one after the other to the output file
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(in) :: data(:)
nvar = nvar + 1
if (nvar > mvar) stop 'Too many variables'
var(nvar)%v => rads_varptr (S, varnm)
var(nvar)%d = data
end subroutine new_var

!-----------------------------------------------------------------------
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads
integer(fourbyteint) :: i

! Check which variables are empty
do i = 1,nvar
	var(i)%empty = all(isnan_(var(i)%d(1:nrec)))
enddo
if (any(var(1:nvar)%empty)) then
	write (*,550,advance='no') '... No'
	do i = 1,nvar
		if (var(i)%empty) write (*,550,advance='no') trim(var(i)%v%name)
	enddo
endif

! Open output file
call rads_create_pass (S, P, nrec)

! Define all variables
do i = 1,nvar
	call rads_def_var (S, P, var(i)%v)
enddo

! Fill all the data fields
do i = 1,nvar
	call rads_put_var (S, P, var(i)%v, var(i)%d(1:nrec))
enddo

! Close the data file
write (*,552) nrec,trim(P%filename(len_trim(S%dataroot)+2:))
call rads_close_pass (S, P)

! Formats
550 format (a,1x)
552 format ('...',i5,' records written to ',a)

end subroutine put_rads

!-----------------------------------------------------------------------
! nc2f: Load flag field, then set corresponding bit in RADS
! varnm : source variable name
! bit   : RADS bit to be set when value is val
! lim   : set bit when value >= lim (optional)
! val   : set bit when value == val (optional, default = 1)
! neq   : set bit when value /= val (optional)

subroutine nc2f (varnm, bit, lim, val, neq)
character(*), intent(in) :: varnm
integer(fourbyteint), intent(in) :: bit
integer(fourbyteint), optional, intent(in) :: lim,val,neq
integer(twobyteint) :: flag(mrec),flag2d(1:1,1:nrec)
integer(fourbyteint) :: i,ival,ndims

if (nf90_inq_varid(ncid,varnm,varid) /= nf90_noerr) then
	write (*,'("No such variable: ",a)') trim(varnm)
	return
endif
call nfs(nf90_inquire_variable(ncid,varid,ndims=ndims))
if (ndims == 2) then
	call nfs(nf90_get_var(ncid,varid,flag2d(1:1,1:nrec)))
	flag(1:nrec) = flag2d(1,1:nrec)
else
	call nfs(nf90_get_var(ncid,varid,flag(1:nrec)))
endif
if (present(lim)) then
	do i = 1,nrec
		if (flag(i) >= lim) flags(i) = ibset(flags(i),bit)
	enddo
else if (present(neq)) then
	do i = 1,nrec
		if (flag(i) /= neq) flags(i) = ibset(flags(i),bit)
	enddo
else
	ival = 1
	if (present(val)) ival = val
	do i = 1,nrec
		if (flag(i) == ival) flags(i) = ibset(flags(i),bit)
	enddo
endif
end subroutine nc2f
