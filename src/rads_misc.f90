!-----------------------------------------------------------------------
! $Id$
!
! Copyright (C) 2011  Remko Scharroo (Altimetrics LLC)
! See LICENSE.TXT file for copying and redistribution conditions.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!-----------------------------------------------------------------------

module rads_misc
use typesizes

!***********************************************************************
!*d_int -- Convert integer to double while accounting for NaNs
!+
! elemental function d_int (i)
! integer(onebyteint) <or> integer(twobyteint) <or> integer(fourbyteint) :: i
! real(eightbytereal) :: d_int
!
! Convert an integer (or integer array) to a double, while taking into
! account that the maximum value for that type of integer (127, 3276, etc)
! indicates NaN.
!
! Arguments:
!  i     : Integer to convert to double
!  d_int : Result as a double real
!-----------------------------------------------------------------------
private :: d_int1, d_int2, d_int4
interface d_int
	module procedure d_int1
	module procedure d_int2
	module procedure d_int4
end interface d_int

contains

!***********************************************************************
!*strconcat -- Concatenate strings
!+
subroutine strconcat (val, string)
character(len=*), intent(in) :: val(:)
character(len=*), intent(out) :: string
!-
! Arguments
!  val     : array of strings
!  string  : concatenation of all strings val(:), separated by spaces
!-----------------------------------------------------------------------
integer :: i
string = val(1)
do i = 2,size(val)
	string = trim(string) // ' ' // val(i)
enddo
end subroutine strconcat

!***********************************************************************
!*chartrans -- Translate characters in a string
!+
subroutine chartrans (string, from, to)
character(len=*), intent(inout) :: string
character(len=*), intent(in) :: from, to

! This program translates each matching character of <from> in the string
! <string> with the corresponding character in <to>.
!
! <from> and <to> should be of the same length. If not, only the lesser
! number of characters are considered.
!
! Example:
!     string = 'original'
!     call chartrans (string, 'oi', 'O!')
! results in:
!     string = 'Or!g!nal'
!
! Arguments:
!  string : Character string to be translated
!  from   : Set of characters to be changed into the set <to>
!  to     : Set of characters to change <from> into
!-----------------------------------------------------------------------
integer :: i, j, n
n = min(len(from),len(to))
do i = 1,len(string)
	do j = 1,n
		if (string(i:i) == from(j:j)) then
			string(i:i) = to(j:j)
			exit
		endif
	enddo
enddo
end subroutine chartrans

!***********************************************************************
!*getlun -- Get free logical unit number
!+
function getlun()
integer :: getlun
!
! This function returns a logical unit number that is not currently
! connected to any file. A value of 0 is returned if no free unit is found.
!
! Example:
! integer :: unit, getlun
! unit = getlun()
! open (unit=unit, file='filename', status='old')
!-----------------------------------------------------------------------
integer :: unit
logical :: opened

do unit = 99,7,-1
	inquire (unit=unit,opened=opened)
	if (.not.opened) then
		getlun = unit
		return
	endif
enddo
getlun = 0
end function getlun

!***********************************************************************
!*checkenv -- Get environment variable, if it exists
!+
subroutine checkenv (env, string)
character(len=*), intent(in) :: env
character(len=*), intent(inout) :: string
!
! This routine returns the contents of environment variable <env> in the
! variable <string>.
! If the environment variable <env> is not set, <string> is unchanged.
!
! Arguments:
!   env     (input) : Name of the environment variable.
!   string  (input) : Default value for string.
!          (output) : Contents of the environment variable or default.
!-
character(len=160) :: temp
call getenv (env,temp)
if (temp /= ' ') string = temp
end subroutine checkenv

!***********************************************************************
!*outofrange - Verify if value is within limits
!+
function outofrange (limits, value)
real(eightbytereal), intent(in) :: limits(2)
real(eightbytereal), intent(inout) :: value
logical :: outofrange
!
! This function checks if value is within limits(1) and limits(2), where
! limits(1) < llimits(2). If value is not within those limits, outofrange
! is set to .true.
!
! If either of the limits is NaN, no check is performed at that limit.
!
! Arguments:
!  limits     : minimum and maximum allowed value
!  value      : value to be checked
!  outofrange : .true. if value is outside limits, .false. otherwise.
!-
outofrange = (value < limits(1) .or. value > limits(2))
end function outofrange

!***********************************************************************
!*checklon - Convert and verify longitude
!+
function checklon (lon_bounds, lon)
real(eightbytereal), intent(in) :: lon_bounds(2)
real(eightbytereal), intent(inout) :: lon
logical :: checklon
!
! This function converts the logitude lon (in degrees) to within limits
! lon_bounds, where lon_bounds(1) < lon_bounds(2). If lon is not within those
! boundaries, checklon is set to .true.
! The longitude will be wrapped (advanced or reduced by a multiple of 360
! degrees) if needed.
!
! If either of the limits is NaN, no wrapping and no check is performed.
!
! Arguments:
!  lon_bounds : longitude limits (degrees)
!  lon        : longitude (degrees)
!  checklon   : .true. if lon is outside limits, .false. otherwise.
!-
if (any(isnan(lon_bounds))) then
	checklon = .false.
	return
endif
lon = lon_bounds(1) + modulo (lon-lon_bounds(1),360d0)
checklon = (lon > lon_bounds(2))
end function checklon

!***********************************************************************
!&nint1 -- Round 8-byte real to 1-byte integer
!+
elemental function nint1(x)
integer(onebyteint) :: nint1
real(eightbytereal), intent(in) :: x
!
! This elemental function rounds an 8-byte real to a 1-byte interger.
! If the real is out of range, or NaN, the returned value is 127.
! Since this function is elemental, it can be applied to arrays as well.
!-----------------------------------------------------------------------
integer(onebyteint), parameter :: imax = huge(0_onebyteint)
real(eightbytereal), parameter :: xmin = -imax-1.5d0, xmax = imax+0.5d0
if (x > xmin .and. x < xmax) then
	nint1 = nint(x,onebyteint)
else ! Out of range or NaN
	nint1 = imax
endif
end function nint1

!***********************************************************************
!&nint2 -- Round 8-byte real to 2-byte integer
!+
elemental function nint2(x)
integer(twobyteint) :: nint2
real(eightbytereal), intent(in) :: x
!
! This elemental function rounds an 8-byte real to a 2-byte interger.
! If the real is out of range, or NaN, the returned value is 32767.
! Since this function is elemental, it can be applied to arrays as well.
!-----------------------------------------------------------------------
integer(twobyteint), parameter :: imax = huge(0_twobyteint)
real(eightbytereal), parameter :: xmin = -imax-1.5d0, xmax = imax+0.5d0
if (x > xmin .and. x < xmax) then
	nint2 = nint(x,twobyteint)
else ! Out of range or NaN
	nint2 = imax
endif
end function nint2

!***********************************************************************
!&nint4 -- Round 8-byte real to 4-byte integer
!+
elemental function nint4(x)
integer(fourbyteint) :: nint4
real(eightbytereal), intent(in) :: x
!
! This elemental function rounds an 8-byte real to a 4-byte interger.
! If the real is out of range, or NaN, the returned value is 2147483647.
! Since this function is elemental, it can be applied to arrays as well.
!-----------------------------------------------------------------------
integer(fourbyteint), parameter :: imax = huge(0_fourbyteint)
real(eightbytereal), parameter :: xmin = -imax-1.5d0, xmax = imax+0.5d0
if (x > xmin .and. x < xmax) then
	nint4 = nint(x,fourbyteint)
else ! Out of range or NaN
	nint4 = imax
endif
end function nint4

!***********************************************************************
!*to16bits - Convert to individual 16 bits
!+
pure function to16bits (x)
real(eightbytereal), intent(in) :: x
integer :: to16bits(0:15)
!
! Convert a double float first to integer and then to individual 16 bits.
! Takes into account NaN and out of bounds numbers.
!-
integer(twobyteint), parameter :: imax = huge(0_twobyteint)
real(eightbytereal), parameter :: xmin = -imax-1.5d0, xmax = imax+0.5d0
integer(twobyteint) :: i, j
if (x > xmin .and. x < xmax) then
	i = nint(x, twobyteint)
	to16bits = 0
	do j = 0,15
		if (btest(i,j)) to16bits(j) = 1
	enddo
else
	to16bits = -1
endif
end function to16bits

!***********************************************************************
!*make_nan -- Create a NaN value
!+
pure function make_nan ()
real(eightbytereal) :: make_nan
make_nan = 0d0
make_nan = make_nan / make_nan
end function make_nan

!***********************************************************************
!*cross_product -- Compute cross product of two 3-D vectors
!+
function cross_product (a, b) result(c)
real(eightbytereal), intent(in) :: a(3), b(3)
real(eightbytereal) :: c(3)
!
! This function returns the vector obtained by computing the cross
! product of two 3-dimensional vectors.
!-----------------------------------------------------------------------
c(1) = a(2) * b(3) - a(3) * b(2)
c(2) = a(3) * b(1) - a(1) * b(3)
c(3) = a(1) * b(2) - a(2) * b(1)
end function cross_product

elemental function d_int1 (i)
integer(onebyteint), intent(in) :: i
real(eightbytereal) :: d_int1
if (i == huge(0_onebyteint)) then
	d_int1 = make_nan()
else
	d_int1 = i
endif
end function d_int1

elemental function d_int2 (i)
integer(twobyteint), intent(in) :: i
real(eightbytereal) :: d_int2
if (i == huge(0_twobyteint)) then
	d_int2 = make_nan()
else
	d_int2 = i
endif
end function d_int2

elemental function d_int4 (i)
integer(fourbyteint), intent(in) :: i
real(eightbytereal) :: d_int4
if (i == huge(0_fourbyteint)) then
	d_int4 = make_nan()
else
	d_int4 = i
endif
end function d_int4

end module rads_misc
