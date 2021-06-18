!-----------------------------------------------------------------------
! Copyright (c) 2011-2021  Remko Scharroo
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

! This file is included in rads_netcdf.f90 and shared between get_var_1d
! and get_var_2d.
! Both routines are used to read a single variable from a NetCDF file
! (either 1- or 2-dimensional), or load a combination of variables
! and/or constants using RPN notation. Note that only the following RPN
! commands are available:
! --------------------------------------------
! Command | Example   | Result
! --------------------------------------------
! ADD     |  a b ADD  | a + b
! SUB     |  a b SUB  | a - b
! NEG     |  a NEG    | -a
! MUL     |  a b MUL  | a * b
! DIV     |  a b DIV  | a / b
! AND     |  a b AND  | if (isnan(a)) b else a
! HYPOT   |  a b HYPOT| sqrt(a*a+b*b)
! R2      |  a b R2   | a*a + b*b
! ISNAN   |  a ISNAN  | if (isnan(a)) 1 else 0
! IOR     |  a b IOR  | bitwise OR of a and b
! --------------------------------------------
! Only 2 buffers can be used, so the computation has to remain relatively
! simple.

real(eightbytereal) :: scale_factor, add_offset, fillvalue
integer(fourbyteint) :: i0, i1, ia, ib, ios, l, varid
logical :: with_fillvalue
i1 = 0
l = len_trim(varnm)
do
	if (i1 > l) exit
	i0 = i1
	i1 = index(varnm(i0+1:), ' ') + i0
	if (i1 == i0) i1 = l + 1
	ia = i0+1
	ib = i1-1
	select case (varnm(ia:ib))
	case ('ADD')
		array = array + temp
	case ('SUB')
		array = array - temp
	case ('NEG')
		array = -array
	case ('MUL')
		array = array * temp
	case ('DIV')
		array = array / temp
	case ('AND')
		where (array /= array) array = temp
	case ('HYPOT')
		array = sqrt(array*array + temp*temp)
	case ('R2')
		array = array*array + temp*temp
	case ('ISNAN')
		where (array == array)
			array = 0d0
		elsewhere	! NaN
			array = 1d0
		endwhere
	case ('IOR')
		where (temp /= temp)
			array = nan
		elsewhere
			array = ior(nint(array),nint(temp))
		endwhere
	case default
		if (index('.+-0123456789', varnm(ia:ia)) > 0) then
			scale_factor = 0d0
			read (varnm(ia:ib), *, iostat=ios) scale_factor
			temp = scale_factor
			cycle
		endif
		if (nf90_inq_varid_warn(ncid,varnm(ia:ib),varid) /= nf90_noerr) return
		if (nf90_get_att(ncid,varid,'scale_factor',scale_factor) /= nf90_noerr) scale_factor = 1d0
		if (nf90_get_att(ncid,varid,'add_offset',add_offset) /= nf90_noerr) add_offset = 0d0
		with_fillvalue = (nf90_get_att(ncid,varid,'_FillValue',fillvalue) == nf90_noerr)
		if (i0 == 0) then
			call nfs(nf90_get_var(ncid,varid,array))
			if (with_fillvalue) where (array == fillvalue) array = nan
			array = array * scale_factor + add_offset
		else
			call nfs(nf90_get_var(ncid,varid,temp))
			if (with_fillvalue) where (temp == fillvalue) temp = nan
			temp = temp * scale_factor + add_offset
		endif
	end select
enddo
