!*NF_SUBS -- Module with useful interfaces to netCDF
!+
! This module provides a selections of routines that make it easier
! to work with netCDF files. It accomplishes commonly used tasks during
! the creating or reading of netCDF grids and other netCDF files
!
! To excess the interface, add "use nf_subs" to your Fortran 90 program
!-
! $Log: nf_subs.f90,v $
! Revision 1.14  2018/08/21 09:26:49  rads
! - Add optional argument 'deflate'
!
! Revision 1.13  2018/08/21 08:35:06  rads
! nf90_put_axis: updated in accordance with nf90_def_axis
!
! Revision 1.12  2018/08/21 08:28:01  rads
! - nf90_def_axis: Add attributes valid_min, valid_max, grid_step
!
! Revision 1.11  2017/09/28 12:25:47  rads
! Removed dependency on node_offset in nf90_def_axis
! Hardened nf90_put_axis in case node_offset is missing
!
! Revision 1.10  2017/09/27 18:47:23  rads
! Add optional arguments axis and standard_name in nf90_def_axis
! Remove optional argument xtype
!
! Revision 1.9  2010/09/15 14:26:24  rads
! - Avoid tiny numbers in coordinate variables
!
! Revision 1.8  2010/07/09 15:19:00  rads
! - Added optional second argument to nfs
!
! Revision 1.7  2009/09/07 14:46:02  rads
! - Added function nft
!
! Revision 1.6  2009/04/29 14:30:42  rads
! - Allow unlimited dimension in nf90_def_axis
! - Added function nint1
!
! Revision 1.5  2009/02/02 02:12:17  rads
! - Added nint2 and nint4
!
! Revision 1.4  2008/08/29 18:11:38  rads
! - Allow units in [] or () as part of longname
!
! Revision 1.3  2008/06/05 18:06:56  rads
! - Trying yet another way to make filling the axis array smarter
!
! Revision 1.2  2008/05/22 14:50:06  rads
! - Use a more precise way to fill the coordinate array
!
! Revision 1.1  2008/05/22 01:43:41  rads
! - Added nf_subs.f90 to rssubs
!-----------------------------------------------------------------------
module nf_subs

contains

!-----------------------------------------------------------------------
!&NF90_DEF_AXIS -- Define dimension as a coordinate axis
!+
subroutine nf90_def_axis (ncid, varnm, long_name, units, nx, x0, x1, dimid, varid, axis, standard_name, deflate)
use typesizes
use netcdf
character(len=*), intent(in) :: varnm, long_name, units
integer, intent(in) :: ncid, nx
real(eightbytereal), intent(in) :: x0,x1
integer, intent(out) :: dimid, varid
character(len=*), intent(in), optional :: axis, standard_name
integer, intent(in), optional :: deflate
!
! This routine sets up a dimension and its associated variable in a netCDF
! file. Supply variable name, long name, units, number of elements, and range.
! The call needs to happen during the define stage of creating a netCDF file.
!
! To fill the coordinate array, use nf90_put_axis after closing the define
! stage.
!
! The string "long_name" can either contain the "long_name" attribute or
! both the "long_name" and "units" attribute. In the latter case one can use,
! for example, long_name = "height [m]" and units = "[]", which will make this
! routine extract the "m" as unit.
!
! The optional strings <axis> and <standard_name> can be added to add these
! additional attributes to the coordinate variable.
!
! The optional value <deflate> can specify the deflate level, if required.
!
! When creating an "unlimited" dimension, specify nx as nf90_unlimited (or 0).
!
! ncid     : NetCDF file ID
! varnm    : (short) variable name
! long_name: Longer description of coordinate variable
! units    : Units of the coordinate variable
! nx       : Number of elements in coordinate array (can be nf90_unlimited)
! x0, x1   : Start and end of array
! dimid    : NetCDF dimension ID
! varid    : NetCDF variable ID
! axis     : (Optional) string for axis attribute
! standard_name : (Optional) string for standard_name attribute
! deflate  : (Optional) deflate level
!****-------------------------------------------------------------------
integer(fourbyteint) :: i, j, node_offset
real(eightbytereal) :: dx = 0d0

! Create dimension and variable for this axis
call nfs(nf90_def_dim(ncid,varnm,abs(nx),dimid))
call nfs(nf90_def_var(ncid,varnm,nf90_double,dimid,varid))

! Add attributes
if (units == '[]' .or. units == '()') then
	i = index(long_name,units(1:1))
	j = index(long_name,units(2:2))
	if (i > 0) then
		call nfs(nf90_put_att(ncid,varid,'long_name',long_name(:i-2)))
	else if (long_name /= ' ') then
		call nfs(nf90_put_att(ncid,varid,'long_name',trim(long_name)))
	endif
	if (j > i+1) call nfs(nf90_put_att(ncid,varid,'units',long_name(i+1:j-1)))
else
	if (long_name /= '') call nfs(nf90_put_att(ncid,varid,'long_name',trim(long_name)))
	if (units /= '') call nfs(nf90_put_att(ncid,varid,'units',trim(units)))
endif
if (present(standard_name) .and. standard_name /= '') call nfs(nf90_put_att(ncid,varid, &
	'standard_name', trim(standard_name)))
if (present(axis) .and. axis /= '') call nfs(nf90_put_att(ncid,varid,'axis',trim(axis)))

! Get global attribute 'node_offset' to see if we are pixel oriented
if (nf90_get_att(ncid,nf90_global,'node_offset',node_offset) /= nf90_noerr) node_offset = 0

! Add 'valid_min', 'valid_max' and 'grid_step' attributes for our CLS colleagues
dx = (x1-x0) / (nx-1)
call nfs(nf90_put_att(ncid,varid,'valid_min',x0))
call nfs(nf90_put_att(ncid,varid,'valid_max',x1))
call nfs(nf90_put_att(ncid,varid,'grid_step',dx))

! Set deflate level
if (present(deflate)) call nfs(nf90_def_var_deflate(ncid,varid,1,deflate,deflate))

end subroutine nf90_def_axis

!-----------------------------------------------------------------------
!&NF90_PUT_AXIS -- Fill an coordinate array
!+
subroutine nf90_put_axis (ncid,varid,len)
use typesizes
use netcdf
integer, intent(in) :: ncid,varid
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
integer :: dimid(1), nx, i
real(eightbytereal), allocatable :: x(:)
real(eightbytereal) :: x0, x1
call nfs(nf90_inquire_variable(ncid,varid,dimids=dimid))
call nfs(nf90_inquire_dimension(ncid,dimid(1),len=nx))
if (nx == nf90_unlimited) nx = len
call nfs(nf90_get_att(ncid,varid,'valid_min',x0))
call nfs(nf90_get_att(ncid,varid,'valid_max',x1))
allocate(x(nx))
!
! This kind of prolonged way of filling x() is to make sure that:
! (1) We do not get a divide by zero when nx=1
! (2) Start and end values are exact (for node_offset=0)
! (3) Intermediate values are as close as possible to intended numbers
! The alternative, to add multiples of dx is prone to increasing errors
!
x(1) = x0
do i = 2,nx-1
	x(i) = (x0*(nx-i)+x1*(i-1))/(nx-1)
enddo
x(nx) = x1
call nfs(nf90_put_var(ncid,varid,x))
deallocate(x)
end subroutine nf90_put_axis

!-----------------------------------------------------------------------
!&NFT -- Test for netCDF error
!+
logical function nft(ios)
use netcdf
integer, intent(in) :: ios
!
! This is a wrapper for netCDF functions. It catches the status return code of the
! netCDF routine and checks if an error occurred. If so, the function returns
! a .true.
!
! Example:
! if (nft(nf90_open ('file.nc', nf90_write, ncid))) write (*,*) 'Error opening file'
!-----------------------------------------------------------------------
nft = (ios /= nf90_noerr)
end function nft

!-----------------------------------------------------------------------
!&NFS -- Stop on netCDF error
!+
subroutine nfs(ios, string)
use netcdf
integer, intent(in) :: ios
character(*), optional, intent(in) :: string
!
! This is a wrapper for netCDF functions. It catches the status return code of the
! netCDF routine and checks if an error occurred. If so, it prints an error message
! and stops execution.
!
! Example:
! call nfs (nf90_open ('file.nc', nf90_write, ncid))
!-----------------------------------------------------------------------
character(80) :: prognm
if (ios == nf90_noerr) return
call getarg(0,prognm)
if (present(string)) then
	write (0,'(a,": ",a,": ",a)') trim(prognm),trim(nf90_strerror(ios)),trim(string)
else
	write (0,'(a,": ",a)') trim(prognm),trim(nf90_strerror(ios))
endif
stop
end subroutine nfs

!-----------------------------------------------------------------------
!&NINT1 -- Round 8-byte real to 1-byte integer
!+
elemental function nint1(x)
use typesizes
integer(onebyteint) :: nint1
real(eightbytereal), intent(in) :: x
!
! This elemental function rounds an 8-byte real to a 1-byte interger.
! If the real is out of range, or NaN, the returned value is -128.
! Since this function is elemental, it can be applied to arrays as well.
!-----------------------------------------------------------------------
integer(onebyteint), parameter :: minint1 = -127-1
if (abs(x) < 127.5d0) then
	nint1 = nint(x,onebyteint)
else	! Out of range or NaN
	nint1 = minint1
endif
end function nint1

!-----------------------------------------------------------------------
!&NINT2 -- Round 8-byte real to 2-byte integer
!+
elemental function nint2(x)
use typesizes
integer(twobyteint) :: nint2
real(eightbytereal), intent(in) :: x
!
! This elemental function rounds an 8-byte real to a 2-byte interger.
! If the real is out of range, or NaN, the returned value is -32768.
! Since this function is elemental, it can be applied to arrays as well.
!-----------------------------------------------------------------------
integer(twobyteint), parameter :: minint2 = -32767-1
if (abs(x) < 32767.5d0) then
	nint2 = nint(x,twobyteint)
else	! Out of range or NaN
	nint2 = minint2
endif
end function nint2

!-----------------------------------------------------------------------
!&NINT4 -- Round 8-byte real to 4-byte integer
!+
elemental function nint4(x)
use typesizes
integer(fourbyteint) :: nint4
real(eightbytereal), intent(in) :: x
!
! This elemental function rounds an 8-byte real to a 4-byte interger.
! If the real is out of range, or NaN, the returned value is -2147483648.
! Since this function is elemental, it can be applied to arrays as well.
!-----------------------------------------------------------------------
integer(fourbyteint), parameter :: minint4 = -2147483647-1
if (abs(x) < 2147483647.5d0) then
	nint4 = nint(x,fourbyteint)
else	! Out of range or NaN
	nint4 = minint4
endif
end function nint4

!-----------------------------------------------------------------------
end module nf_subs
