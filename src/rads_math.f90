!****-------------------------------------------------------------------
! Copyright (c) 2011-2018  Remko Scharroo
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
!****-------------------------------------------------------------------

module rads_math
use typesizes

! Define a linked list for arrays on the math stack
type :: math_ll
	real(eightbytereal), allocatable :: data(:)
	type(math_ll), pointer :: prev
end type
integer, parameter, private :: stderr = 0
real(eightbytereal), parameter, private :: nan = transfer ((/not(0_fourbyteint),not(0_fourbyteint)/),0d0)

contains

!****f* rads_math/math_push
! SUMMARY
! Put new buffer on top of stack
!
! SYNOPSIS
subroutine math_push (n, top)
integer(fourbyteint) :: n
type(math_ll), pointer :: top
!
! PURPOSE
! This routine adds a new item to the top of the linked list of math
! buffers.
!
! ARGUMENTS
! n     : Number of elements for new data buffer
! top   : Pointer to the top of the stack
!****-------------------------------------------------------------------
type(math_ll), pointer :: temp
allocate (temp)
allocate (temp%data(n))
temp%prev => top
top => temp
end subroutine math_push

!****f* rads_math/math_pop
! SUMMARY
! Remove buffer from top of stack
!
! SYNOPSIS
subroutine math_pop (top)
type(math_ll), pointer :: top
!
! PURPOSE
! This routine removes the item on the top of the linked list of math
! buffers.
!
! ARGUMENT
! top   : Pointer to the top of the stack
!****-------------------------------------------------------------------
type(math_ll), pointer :: temp
if (.not.associated(top)) return
temp => top%prev
deallocate (top%data)
deallocate (top)
top => temp
end subroutine math_pop

!****f* rads_math/math_eval
! SUMMARY
! Evaluate math command or value
!
! SYNOPSIS
function math_eval (string, n, top)
use rads_misc
character(len=*), intent(in) :: string
integer(fourbyteint), intent(in) :: n
type(math_ll), pointer, intent(inout) :: top
integer(fourbyteint) :: math_eval
!
! ARGUMENTS
! string    : math command or value
! n         : number of elements in data buffer
! top       : Pointer to top of the stack
!
! RETURN VALUE
! math_eval : 0 = no error, -1 = no matching command
!****-------------------------------------------------------------------
type(math_ll), pointer :: temp
real(eightbytereal) :: x, y, z, w
real(eightbytereal), parameter :: pi = 4d0*atan(1d0), d2r = pi/180d0, r2d = 180d0/pi
integer(fourbyteint) :: ios, i, j, k

math_eval = 0

! Skip empty strings
if (string == '') return

! Check if string could be a value
if ((string(:1) >= '0' .and. string(:1) <= '9') .or. string(:1) == '-' .or. &
	string(:1) == '+' .or. string(:1) == '.') then
	read (string,*,iostat=ios) x
	if (ios /= 0) call math_exit ('Error parsing value; iostat =',ios)
	call math_push (n, top)
	top%data = x
	return
endif

! Try any of the defined math commands
select case (string)

!! Note: a = NaN when x or y is NaN unless indicated otherwise

! Most commonly used first
!! x y SUB a : a = x - y
case ('SUB')
	call math_check (2)
	top%prev%data = top%prev%data - top%data
	call math_pop (top)
!! x y ADD a : a = x + y
case ('ADD')
	call math_check (2)
	top%prev%data = top%prev%data + top%data
	call math_pop (top)
!! x y MUL a : a = x * y
case ('MUL')
	call math_check (2)
	top%prev%data = top%prev%data * top%data
	call math_pop (top)

! - MATH x (0 arguments to 1)
!! PI a : a = pi
case ('PI')
	call math_push (n, top)
	top%data = pi
!! E a : a = exp(1)
case ('E')
	call math_push (n, top)
	top%data = exp(1d0)

! x MATH - (1 argument to 0)
!! x POP : remove last item from stack
case ('POP')
	call math_check (1)
	call math_pop (top)

! x MATH x (1 argument to 1)
!! x NEG a : a = $-$x
case ('NEG')
	call math_check (1)
	top%data = -top%data
!! x ABS a : a = $|$x$|$
case ('ABS')
	call math_check (1)
	top%data = abs(top%data)
!! x INV a : a = 1/x
case ('INV')
	call math_check (1)
	top%data = 1d0 / top%data
!! x SQRT a : a = sqrt(x)
case ('SQRT')
	call math_check (1)
	top%data = sqrt(top%data)
!! x SQR a : a = x*x
case ('SQR')
	call math_check (1)
	top%data = top%data * top%data
!! x EXP a : a = exp(x)
case ('EXP')
	call math_check (1)
	top%data = exp(top%data)
!! x LOG a : a = ln(x)
case ('LOG')
	call math_check (1)
	top%data = log(top%data)
!! x LOG10 a : a = log10(x)
case ('LOG10')
	call math_check (1)
	top%data = log10(top%data)
!! x SIN a : a = sin(x)
case ('SIN')
	call math_check (1)
	top%data = sin(top%data)
!! x COS a : a = cos(x)
case ('COS')
	call math_check (1)
	top%data = cos(top%data)
!! x TAN a : a = tan(x)
case ('TAN')
	call math_check (1)
	top%data = tan(top%data)
!! x SIND a : a = sin(x) [x in degrees]
case ('SIND')
	call math_check (1)
	top%data = sin(top%data * d2r)
!! x COSD a : a = cos(x) [x in degrees]
case ('COSD')
	call math_check (1)
	top%data = cos(top%data * d2r)
!! x TAND a : a = tan(x) [x in degrees]
case ('TAND')
	call math_check (1)
	top%data = tan(top%data * d2r)
!! x SINH a : a = sinh(x)
case ('SINH')
	call math_check (1)
	top%data = sin(top%data)
!! x COSH a : a = cosh(x)
case ('COSH')
	call math_check (1)
	top%data = cosh(top%data)
!! x TANH a : a = tanh(x)
case ('TANH')
	call math_check (1)
	top%data = tanh(top%data)
!! x ASIN a : a = arcsin(x)
case ('ASIN')
	call math_check (1)
	top%data = asin(top%data)
!! x ACOS a : a = arccos(x)
case ('ACOS')
	call math_check (1)
	top%data = acos(top%data)
!! x ATAN a : a = arctan(x)
case ('ATAN')
	call math_check (1)
	top%data = atan(top%data)
!! x ASIND a : a = arcsin(x) [a in degrees]
case ('ASIND')
	call math_check (1)
	top%data = asin(top%data) * r2d
!! x ACOSD a : a = arccos(x) [a in degrees]
case ('ACOSD')
	call math_check (1)
	top%data = acos(top%data) * r2d
!! x ATAND a : a = arctan(x) [a in degrees]
case ('ATAND')
	call math_check (1)
	top%data = atan(top%data) * r2d
!! x ASINH a : a = arcsinh(x)
case ('ASINH')
	call math_check (1)
	top%data = asinh_(top%data)
!! x ACOSH a : a = arccosh(x)
case ('ACOSH')
	call math_check (1)
	top%data = acosh_(top%data)
!! x ATANH a : a = arctanh(x)
case ('ATANH')
	call math_check (1)
	top%data = atanh_(top%data)
!! x ISNAN a : a = 1 if x is NaN; a = 0 otherwise
case ('ISNAN')
	call math_check (1)
	where (isnan_(top%data))
		top%data = 1d0
	elsewhere
		top%data = 0d0
	endwhere
!! x ISAN a : a = 0 if x is NaN; a = 1 otherwise
case ('ISAN')
	call math_check (1)
	where (isnan_(top%data))
		top%data = 0d0
	elsewhere
		top%data = 1d0
	endwhere
!! x RINT a : a is nearest integer to x
!! x NINT a : a is nearest integer to x
case ('RINT', 'NINT')
	call math_check (1)
	where (isan_(top%data)) top%data = nint(top%data)
!! x CEIL a : a is nearest integer greater or equal to x
!! x CEILING a : a is nearest integer greater or equal to x
case ('CEIL', 'CEILING')
	call math_check (1)
	where (isan_(top%data)) top%data = ceiling(top%data)
!! x FLOOR a : a is nearest integer less or equal to x
case ('FLOOR')
	call math_check (1)
	where (isan_(top%data)) top%data = floor(top%data)
!! x D2R a : convert x from degrees to radian
case ('D2R')
	call math_check (1)
	top%data = top%data * d2r
!! x R2D a : convert x from radian to degrees
case ('R2D')
	call math_check (1)
	top%data = top%data * r2d
!! x YMDHMS a : convert seconds of 1985 to format YYYYMMDDHHMMSS
case ('YMDHMS')
	call math_check (1)
	where (isan_(top%data)) top%data = ymdhms_(top%data)
!! x SUM a : a(i) = x(1) + ... + x(i) while skipping all NaN
case ('SUM')
	call math_check (1)
	x = sum(top%data, isan_(top%data))
!! x DIF a : a(i) = x(i)-x(i-1); a(1) = NaN
case ('DIF')
	call math_check (1)
	do i = n,2,-1
		top%data(i) = top%data(i) - top%data(i-1)
	enddo
	top%data(1) = nan

! x MATH x y (1 argument to 2)
!! x DUP a b : duplicate the last item on the stack
case ('DUP')
	call math_check (1)
	call math_push (n, top)
	top%data = top%prev%data

! x y MATH x (2 arguments to 1)
!! x y DIV a : a = x / y
case ('DIV')
	call math_check (2)
	top%prev%data = top%prev%data / top%data
	call math_pop (top)
!! x y POW a : a = $x^y$
case ('POW')
	call math_check (2)
	top%prev%data = top%prev%data ** top%data
	call math_pop (top)
!! x y FMOD a : a = x modulo y
case ('FMOD')
	call math_check (2)
	top%prev%data = modulo(top%prev%data, top%data)
	call math_pop (top)
!! x y MIN a : a = the lesser of x and y
case ('MIN')
	call math_check (2)
	top%prev%data = min(top%prev%data, top%data)
	call math_pop (top)
!! x y MAX a : a = the greater of x and y
case ('MAX')
	call math_check (2)
	top%prev%data = max(top%prev%data, top%data)
	call math_pop (top)
!! x y ATAN2 a : a = arctan(x) taking into account the quadrant as determined by y
case ('ATAN2')
	call math_check (2)
	top%prev%data = atan2(top%prev%data, top%data)
	call math_pop (top)
!! x y HYPOT a : a = sqrt(x*x+y*y)
case ('HYPOT')
	call math_check (2)
	top%prev%data = hypot_(top%prev%data, top%data)
	call math_pop (top)
!! x y R2 a : a = x*x + y*y
case ('R2')
	call math_check (2)
	top%prev%data = top%prev%data*top%prev%data + top%data*top%data
	call math_pop (top)
!! x y EQ a : a = 1 if x == y; a = 0 otherwise
case ('EQ')
	call math_check (2)
	where (isnan_(top%prev%data) .or. isnan_(top%data))
		top%prev%data = nan
	elsewhere (top%prev%data == top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
!! x y NE a : a = 0 if x == y; a = 1 otherwise
case ('NEQ', 'NE')
	call math_check (2)
	where (isnan_(top%prev%data) .or. isnan_(top%data))
		top%prev%data = nan
	elsewhere (top%prev%data /= top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
!! x y LT a : a = 1 if x $<$ y; a = 0 otherwise
case ('LT')
	call math_check (2)
	where (isnan_(top%prev%data) .or. isnan_(top%data))
		top%prev%data = nan
	elsewhere (top%prev%data < top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
!! x y LE a : a = 1 if x $\le$ y; a = 0 otherwise
case ('LE')
	call math_check (2)
	where (isnan_(top%prev%data) .or. isnan_(top%data))
		top%prev%data = nan
	elsewhere (top%prev%data <= top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
!! x y GT a : a = 1 if x $>$ y; a = 0 otherwise
case ('GT')
	call math_check (2)
	where (isnan_(top%prev%data) .or. isnan_(top%data))
		top%prev%data = nan
	elsewhere (top%prev%data > top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
!! x y GE a : a = 1 if x $\ge$ y; a = 0 otherwise
case ('GE')
	call math_check (2)
	where (isnan_(top%prev%data) .or. isnan_(top%data))
		top%prev%data = nan
	elsewhere (top%prev%data >= top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
!! x y NAN a : a = NaN if x == y; a = x otherwise
case ('NAN')
	call math_check (2)
	where (top%prev%data == top%data) top%prev%data = nan
	call math_pop (top)
!! x y AND a : a = y if x is NaN; a = x otherwise
case ('AND')
	call math_check (2)
	where (isnan_(top%prev%data)) top%prev%data = top%data
	call math_pop (top)
!! x y OR a : a = NaN if y is NaN; a = x otherwise
case ('OR')
	call math_check (2)
	where (isnan_(top%data)) top%prev%data = nan
	call math_pop (top)
!! x y IAND a : a = bitwise AND of x and y
case ('IAND')
	call math_check (2)
	where (isnan_(top%prev%data) .or. isnan_(top%data))
		top%prev%data = nan
	elsewhere
		top%data = iand(nint(top%prev%data),nint(top%data))
	endwhere
	call math_pop (top)
!! x y IOR a : a = bitwise OR of x and y
case ('IOR')
	call math_check (2)
	where (isnan_(top%prev%data) .or. isnan_(top%data))
		top%prev%data = nan
	elsewhere
		top%data = ior(nint(top%prev%data),nint(top%data))
	endwhere
	call math_pop (top)
!! x y BTEST a : a = 1 if bit y of x is set; a = 0 otherwise
case ('BTEST')
	call math_check (2)
	where (isnan_(top%prev%data) .or. isnan_(top%data))
		top%prev%data = nan
	elsewhere (btest(nint(top%prev%data),nint(top%data)))
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
!! x y AVG a : a = 0.5*(x+y) [when x or y is NaN a returns the other value]
case ('AVG')
	call math_check (2)
	where (isnan_(top%prev%data))
		top%prev%data = top%data
	elsewhere (isan_(top%data))
		top%prev%data = 0.5d0 * (top%prev%data + top%data)
	endwhere
	call math_pop (top)
!! x y DXDY a : a(i) = (x(i+1)-x(i-1))/(y(i+1)-y(i-1)); a(1) = a(n) = NaN
case ('DXDY')
	call math_check (2)
	x = top%prev%data(1)
	top%prev%data(1) = nan
	do i = 2,n-1
		y = (top%prev%data(i+1) - x) / (top%data(i+1) - top%data(i-1))
		x = top%prev%data(i)
		top%prev%data(i) = y
	enddo
	top%prev%data(n) = nan
	call math_pop (top)

! x y MATH x y (2 arguments to 2)
!! x y EXCH a b : exchange the last two items on the stack (NaNs have no influence)
case ('EXCH')
	call math_check (2)
	temp => top%prev
	top%prev => temp%prev
	temp%prev => top
	top => temp

! x y z MATH x (3 arguments to 1)
!! x y z INRANGE a : a = 1 if x is between y and z; a = 0 otherwise (also in case of any NaN)
case ('INRANGE')
	call math_check (3)
	where (top%prev%prev%data >= top%prev%data .and. top%prev%prev%data <= top%data)
		top%prev%prev%data = 1d0
	elsewhere
		top%prev%prev%data = 0d0
	endwhere
	call math_pop (top)
	call math_pop (top)
!! x y z BOXCAR a : a = filter x along monotonic dimension y with boxcar of length z (NaNs are skipped)
case ('BOXCAR')
	call math_check (3)
	z = 0.5d0 * top%data(1)
	do i = 1,n
		x = 0d0
		k = 0
		do j = i,1,-1
			if (isnan_(top%prev%prev%data(j))) cycle
			if (abs(top%prev%data(j)-top%prev%data(i)) > z) exit
			x = x + top%prev%prev%data(j)
			k = k + 1
		enddo
		do j = i+1,n
			if (isnan_(top%prev%prev%data(j))) cycle
			if (abs(top%prev%data(j)-top%prev%data(i)) > z) exit
			x = x + top%prev%prev%data(j)
			k = k + 1
		enddo
		top%data(i) = x/k
	enddo
	top%prev%prev%data = top%data
	call math_pop (top)
	call math_pop (top)
!! x y z GAUSS a : a = filter x along monotonic dimension y with Gauss function with sigma z (NaNs are skipped)
case ('GAUSS')
	call math_check (3)
	z = -0.5d0 / top%data(1) / top%data(1)
	do i = 1,n
		x = 0d0
		y = 0d0
		do j = i,1,-1
			if (isnan_(top%prev%prev%data(j))) cycle
			w = top%prev%data(j) - top%prev%data(i)
			w = exp (w * w * z)
			if (w < 1d-8) exit
			x = x + w * top%prev%prev%data(j)
			y = y + w
		enddo
		do j = i+1,n
			if (isnan_(top%prev%prev%data(j))) cycle
			w = top%prev%data(j) - top%prev%data(i)
			w = exp (w * w * z)
			if (w < 1d-8) exit
			x = x + w * top%prev%prev%data(j)
			y = y + w
		enddo
		top%data(i) = x/y
	enddo
	top%prev%prev%data = top%data
	call math_pop (top)
	call math_pop (top)

! Default (error)
case default
	math_eval = -1

end select

contains

subroutine math_exit (text, n)
character(len=*), intent(in) :: text
integer(fourbyteint), intent(in) :: n
character(len=80) :: prog
call getarg (0,prog)
write (stderr, '(a,": math_eval (",a,") ",a,1x,i0)') trim(prog), trim(string), trim(text), n
call exit (20)
end subroutine math_exit

subroutine math_check (n)
integer(fourbyteint), intent(in) :: n
if (n >= 1 .and. .not.associated(top)) call math_exit ('No items on stack. Requested items =',n)
if (n >= 2 .and. .not.associated(top%prev)) call math_exit ('Only 1 item on stack. Requested items =',n)
if (n >= 3 .and. .not.associated(top%prev%prev)) call math_exit ('Only 2 items on stack. Requested items =',n)
end subroutine math_check

! These functions are based on other functions

elemental function ymdhms_ (x)
use rads_time
real(eightbytereal), intent(in) :: x
real(eightbytereal) :: ymdhms_
integer(fourbyteint) :: yy, mm, dd, hh, mn
real(eightbytereal) :: ss
call sec85ymdhms (x, yy, mm, dd, hh, mn, ss)
ymdhms_ = modulo (yy,100) * 1d10 + mm * 1d8 + dd * 1d6 + hh * 1d4 + mn * 1d2 + ss
end function ymdhms_

! These functions added because they are only specified in Fortran 2008 and later

elemental function asinh_ (x)
real(eightbytereal), intent(in) :: x
real(eightbytereal) :: asinh_
asinh_ = log(x + sqrt(x*x+1))
end function asinh_

elemental function acosh_ (x)
real(eightbytereal), intent(in) :: x
real(eightbytereal) :: acosh_
acosh_ = log(x + sqrt(x*x-1))
end function acosh_

elemental function atanh_ (x)
real(eightbytereal), intent(in) :: x
real(eightbytereal) :: atanh_
atanh_ = 0.5d0 * log((x+1)/(x-1))
end function atanh_

elemental function hypot_ (x, y)
real(eightbytereal), intent(in) :: x, y
real(eightbytereal) :: hypot_
hypot_ = sqrt (x*x + y*y)
end function hypot_

end function math_eval

end module rads_math
