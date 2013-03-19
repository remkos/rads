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

!*rads_netcdf -- Module with useful interfaces to netCDF
!+
! This module provides a selections of routines that make it easier
! to work with netCDF files. It accomplishes commonly used tasks during
! the creating or reading of netCDF grids and other netCDF files
!
! To access the interface, add 'use rads_netcdf' to your Fortran 90 program
!-----------------------------------------------------------------------
module rads_netcdf
use typesizes
integer, parameter, private :: stderr = 0
real(eightbytereal), parameter, private :: nan = transfer ((/not(0_fourbyteint),not(0_fourbyteint)/),0d0)
private get_var_1d, get_var_2d

!-----------------------------------------------------------------------
!*get_var -- Load variable from netCDF file into memory
!+
!subroutine get_var (ncid, varnm, array)
!integer(fourbyteint), intent(in) :: ncid
!character(len=*), intent(in) :: varnm
!real(eightbytereal), intent(out) :: array(:) <or> array(:,:)
!
! This routine looks for a variable named <varnm> in the file associated with
! <ncid> and loads it into the 1- or 2-dimensional array <array>.
! This routine takes into account the attributes scale_factor, add_offset, and
! _FillValue.
!
! One can also add or subtract a number of fields directly. For example:
! 'alt-range'
!
! Arguments:
!  ncid  : NetCDF ID
!  varnm : Variable name
!  array : Array of data values
!-
interface get_var
	procedure get_var_1d
	procedure get_var_2d
end interface

contains

!-----------------------------------------------------------------------
!*nf90_def_axis -- Define dimension as a coordinate axis
!+
subroutine nf90_def_axis (ncid, varnm, longname, units, nx, x0, x1, dimid, varid, xtype)
use netcdf
character(len=*), intent(in) :: varnm,longname,units
integer, intent(in) :: ncid,nx
real(eightbytereal), intent(in) :: x0,x1
integer, intent(out) :: dimid,varid
integer, intent(in), optional :: xtype
!
! This routine sets up a dimension and its associated variable in a netCDF
! file. Supply variable name, long name, units, number of elements, and range.
! The call needs to happen during the define stage of creating a netCDF file.
!
! To fill the coordinate array, use nf90_put_axis after closing the define
! stage.
!
! The string <longname> can either contain the 'long_name' attribute or
! both the 'long_name' and 'units' attribute. In the latter case one can use,
! for example, longname = 'height [m]' and units = '[]', which will make this
! routine extract the 'm' as unit.
!
! When creating an 'unlimited' dimension, specify nx as nf90_unlimited (or 0).
!
! ncid     : NetCDF file ID
! varnm    : (short) variable name
! longname : Longer description of coordinate variable
! units    : Units of the coordinate variable
! nx       : Number of elements in coordinate array (can be nf90_unlimited)
! x0, x1   : Start and end of array
! dimid    : NetCDF dimension ID
! varid    : NetCDF variable ID
! xtype    : Type of variable (optional, default = nf90_double)
!-----------------------------------------------------------------------
integer(fourbyteint) :: i,j
real(eightbytereal) :: xrange(2)
xrange(1) = x0; xrange(2) = x1
call nfs(nf90_def_dim(ncid,varnm,abs(nx),dimid))
if (present(xtype)) then
	call nfs(nf90_def_var(ncid,varnm,xtype,dimid,varid))
else
	call nfs(nf90_def_var(ncid,varnm,nf90_double,dimid,varid))
endif
if (units == '[]' .or. units == '()') then
	i=index(longname,units(1:1))
	j=index(longname,units(2:2))
	if (i > 0) then
		call nfs(nf90_put_att(ncid,varid,'long_name',longname(:i-2)))
	else if (longname /= ' ') then
		call nfs(nf90_put_att(ncid,varid,'long_name',trim(longname)))
	endif
	if (j > i+1) call nfs(nf90_put_att(ncid,varid,'units',longname(i+1:j-1)))
else
	if (longname /= ' ') call nfs(nf90_put_att(ncid,varid,'long_name',trim(longname)))
	if (units /= ' ') call nfs(nf90_put_att(ncid,varid,'units',trim(units)))
endif
call nfs(nf90_put_att(ncid,varid,'actual_range',xrange))
end subroutine nf90_def_axis

!-----------------------------------------------------------------------
!*nf90_put_axis -- Fill an coordinate array
!+
subroutine nf90_put_axis (ncid, varid, len)
use netcdf
integer, intent(in) :: ncid, varid
integer, intent(in), optional :: len
!
! This routine fills a coordinate array previously set up by nf90_def_axis.
! Call this routine only after nf90_enddef().
! When the dimension was specified with length 0 (or nf90_unlimited) in
! nf90_def_axis, then the optional argument len is needed in this call.
!
! ncid    : NetCDF file ID
! varid   : NetCDF variable ID
! len     : Length of the dimension (optional, only needed for unlimited
!           dimension)
!-----------------------------------------------------------------------
integer :: dimid(1),nx,node_offset=0,i
real(eightbytereal), allocatable :: x(:)
real(eightbytereal) :: xrange(2)
call nfs(nf90_inquire_variable(ncid,varid,dimids=dimid))
call nfs(nf90_inquire_dimension(ncid,dimid(1),len=nx))
if (nx == nf90_unlimited) nx = len
call nfs(nf90_get_att(ncid,varid,'actual_range',xrange))
if (nf90_get_att(ncid,nf90_global,'node_offset',node_offset) /= nf90_noerr) node_offset=0
allocate(x(nx))
!
! This kind of prolonged way of filling x() is to make sure that:
! (1) We do not get a divide by zero when nx=1
! (2) Start and end values are exact (for node_offset=0)
! (3) Intermediate values are as close as possible to intended numbers
! The alternative, to add multiples of dx is prone to increasing errors
!
if (node_offset == 0) then
	x(1)=xrange(1)
	do i=2,nx-1
		x(i)=(xrange(1)*(nx-i)+xrange(2)*(i-1))/(nx-1)
	enddo
	x(nx)=xrange(2)
else
	do i=1,nx
		x(i)=(xrange(1)*(nx-i+0.5d0)+xrange(2)*(i-0.5d0))/(nx)
	enddo
endif
call nfs(nf90_put_var(ncid,varid,x))
deallocate(x)
end subroutine nf90_put_axis

!-----------------------------------------------------------------------
!*nff -- Return .false. upon netCDF error
!+
logical function nff(ios)
use netcdf
integer, intent(in) :: ios
!
! This is a wrapper for netCDF functions. It catches the status return code of the
! netCDF routine and checks if an error occurred. The function returns:
!
! nff  : Error (.false.) or no error (.true.)
!
! Example:
! if (nff(nf90_open ('file.nc', nf90_write, ncid))) write (*,*) 'Opening successful'
!-----------------------------------------------------------------------
nff = (ios == nf90_noerr)
end function nff

!-----------------------------------------------------------------------
!*nft -- Return .true. upon netCDF error
!+
logical function nft(ios)
use netcdf
integer, intent(in) :: ios
!
! This is a wrapper for netCDF functions. It catches the status return code of the
! netCDF routine and checks if an error occurred. The function returns:
!
! nft  : Error (.true.) or no error (.false.)
!
! Example:
! if (nft(nf90_open ('file.nc', nf90_write, ncid))) write (*,*) 'Error opening file'
!-----------------------------------------------------------------------
nft = (ios /= nf90_noerr)
end function nft

!-----------------------------------------------------------------------
!*nfs -- Stop execution upon netCDF error
!+
subroutine nfs(ios)
use netcdf
integer, intent(in) :: ios
!
! This is a wrapper for netCDF functions. It catches the status return code of the
! netCDF routine and checks if an error occurred. Upon error, an error message
! is printed to standard error and execution stops.
!
! Example:
! call nfs (nf90_open ('file.nc', nf90_write, ncid))
!-----------------------------------------------------------------------
character(len=80) :: prognm
if (ios==nf90_noerr) return
call getarg (0, prognm)
write (stderr, '(a,": ",a)') trim(prognm),trim(nf90_strerror(ios))
call exit (ios)
end subroutine nfs

subroutine get_var_1d (ncid, varnm, array)
use netcdf
integer(fourbyteint), intent(in) :: ncid
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(out) :: array(:)
real(eightbytereal) :: temp(size(array)), scale_factor, add_offset, fillvalue
integer(fourbyteint) :: i0, i1, l, varid, constant
logical :: with_fillvalue
i1 = 0
l = len_trim(varnm)
do
	if (i1 > l) exit
	i0 = i1
	i1 = scan(varnm(i0+1:), '+-') + i0
	if (i1 == i0) i1 = l + 1
	if (nf90_inq_varid(ncid,varnm(i0+1:i1-1),varid) /= nf90_noerr) then
		write (*,'("No such variable: ",a)') varnm(i0+1:i1-1)
		return
	endif
	if (nf90_get_att(ncid,varid,'scale_factor',scale_factor) /= nf90_noerr) scale_factor = 1d0
	if (nf90_get_att(ncid,varid,'add_offset',add_offset) /= nf90_noerr) add_offset = 0d0
	with_fillvalue = (nf90_get_att(ncid,varid,'_FillValue',fillvalue) /= nf90_noerr)
	constant = 0
	if (i0 == 0) then
		call nfs(nf90_get_var(ncid,varid,array))
		if (with_fillvalue) where (array == fillvalue) array = nan
		array = array * scale_factor + add_offset
	else
		call nfs(nf90_get_var(ncid,varid,temp))
		if (with_fillvalue) where (temp == fillvalue) temp = nan
		if (varnm(i0:i0) == '-') constant = -1
		if (varnm(i0:i0) == '+') constant = 1
		array = array + constant * (temp * scale_factor + add_offset)
	endif
enddo
end subroutine get_var_1d

subroutine get_var_2d (ncid, varnm, array)
use netcdf
integer(fourbyteint), intent(in) :: ncid
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(out) :: array(:,:)
real(eightbytereal) :: temp(size(array,1),size(array,2)), scale_factor, add_offset, fillvalue
integer(fourbyteint) :: i0, i1, l, varid, constant
logical :: with_fillvalue
i1 = 0
l = len_trim(varnm)
do
	if (i1 > l) exit
	i0 = i1
	i1 = scan(varnm(i0+1:), '+-') + i0
	if (i1 == i0) i1 = l + 1
	if (nf90_inq_varid(ncid,varnm(i0+1:i1-1),varid) /= nf90_noerr) then
		write (*,'("No such variable: ",a)') varnm(i0+1:i1-1)
		return
	endif
	if (nf90_get_att(ncid,varid,'scale_factor',scale_factor) /= nf90_noerr) scale_factor = 1d0
	if (nf90_get_att(ncid,varid,'add_offset',add_offset) /= nf90_noerr) add_offset = 0d0
	with_fillvalue = (nf90_get_att(ncid,varid,'_FillValue',fillvalue) /= nf90_noerr)
	constant = 0
	if (i0 == 0) then
		call nfs(nf90_get_var(ncid,varid,array))
		if (with_fillvalue) where (array == fillvalue) array = nan
		array = array * scale_factor + add_offset
	else
		call nfs(nf90_get_var(ncid,varid,temp))
		if (with_fillvalue) where (temp == fillvalue) temp = nan
		if (varnm(i0:i0) == '-') constant = -1
		if (varnm(i0:i0) == '+') constant = 1
		array = array + constant * (temp * scale_factor + add_offset)
	endif
enddo
end subroutine get_var_2d

!-----------------------------------------------------------------------
end module rads_netcdf
