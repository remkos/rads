!-----------------------------------------------------------------------
! $Id$
!
! Copyright (C) 2012  Remko Scharroo (Altimetrics LLC)
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

module rads_math
use typesizes

! Define a linked list for arrays on the math stack
type :: math_ll
	real(eightbytereal), allocatable :: data(:)
	type(math_ll), pointer :: prev
end type
integer, parameter, private :: stderr = 0

contains

!***********************************************************************
!*math_push -- Put new buffer on top of stack
!+
subroutine math_push (n, top)
integer(fourbyteint) :: n
type(math_ll), pointer :: top
!
! This routine adds a new item to the top of the linked list of math
! buffers.
!
! Arguments:
!   n     : Number of elements for new data buffer
!   top   : Pointer to the top of the stack
!-----------------------------------------------------------------------
type(math_ll), pointer :: temp
allocate (temp)
allocate (temp%data(n))
temp%prev => top
top => temp
end subroutine math_push

!***********************************************************************
!*math_pop -- Remove buffer from top of stack
!+
subroutine math_pop (top)
type(math_ll), pointer :: top
!
! This routine removes the item on the top of the linked list of math
! buffers.
!
! Argument:
!   top   : Pointer to the top of the stack
!-----------------------------------------------------------------------
type(math_ll), pointer :: temp
if (.not.associated(top)) return
temp => top%prev
deallocate (top%data)
deallocate (top)
top => temp
end subroutine math_pop

!***********************************************************************
!*math_eval -- Evaluate math command or value
!+
function math_eval (string, n, top)
character(len=*), intent(in) :: string
integer(fourbyteint), intent(in) :: n
type(math_ll), pointer, intent(inout) :: top
integer(fourbyteint) :: math_eval
!
! Arguments:
!  string : math command or value
!  n      : number of elements in data buffer
!  top    : Pointer to top of the stack
!
! Returned value:
!  math_eval : 0 = no error, -1 = no matching command
!-----------------------------------------------------------------------
type(math_ll), pointer :: temp
real(eightbytereal) :: value, nan
real(eightbytereal), parameter :: pi = 4d0*atan(1d0), d2r = pi/180d0, r2d = 180d0/pi
integer(fourbyteint) :: ios

math_eval = 0

! Skip empty strings
if (string == '') return

! Init nan
nan = 0d0
nan = nan / nan

! Check if string could be a value
if ((string(:1) >= '0' .and. string(:1) <= '9') .or. string(:1) == '-' .or. &
	string(:1) == '+' .or. string(:1) == '.') then
	read (string,*,iostat=ios) value
	if (ios /= 0) call math_exit ('Error parsing value; iostat =',ios)
	call math_push (n, top)
	top%data = value
	return
endif

! Try any of the defined math commands
select case (string)

! Most commonly used first
case ('SUB')
	call math_check (2)
	top%prev%data = top%prev%data - top%data
	call math_pop (top)
case ('OR')
	call math_check (2)
	where (isnan(top%data)) top%prev%data = nan
	call math_pop (top)

! - MATH x (0 arguments to 1)
case ('PI')
	call math_push (n, top)
	top%data = pi
case ('E')
	call math_push (n, top)
	top%data = exp(1d0)

! x MATH - (1 argument to 0)
case ('POP')
	call math_check (1)
	call math_pop (top)

! x MATH x (1 argument to 1)
case ('NEG')
	call math_check (1)
	top%data = -top%data
case ('ABS')
	call math_check (1)
	top%data = abs(top%data)
case ('INV')
	call math_check (1)
	top%data = 1d0 / top%data
case ('SQRT')
	call math_check (1)
	top%data = sqrt(top%data)
case ('SQR')
	call math_check (1)
	top%data = top%data * top%data
case ('EXP')
	call math_check (1)
	top%data = exp(top%data)
case ('LOG')
	call math_check (1)
	top%data = log(top%data)
case ('LOG10')
	call math_check (1)
	top%data = log10(top%data)
case ('SIN')
	call math_check (1)
	top%data = sin(top%data)
case ('COS')
	call math_check (1)
	top%data = cos(top%data)
case ('TAN')
	call math_check (1)
	top%data = tan(top%data)
case ('SIND')
	call math_check (1)
	top%data = sin(top%data * d2r)
case ('COSD')
	call math_check (1)
	top%data = cos(top%data * d2r)
case ('TAND')
	call math_check (1)
	top%data = tan(top%data * d2r)
case ('SINH')
	call math_check (1)
	top%data = sin(top%data)
case ('COSH')
	call math_check (1)
	top%data = cosh(top%data)
case ('TANH')
	call math_check (1)
	top%data = tanh(top%data)
case ('ASIN')
	call math_check (1)
	top%data = asin(top%data)
case ('ACOS')
	call math_check (1)
	top%data = acos(top%data)
case ('ATAN')
	call math_check (1)
	top%data = atan(top%data)
case ('ASIND')
	call math_check (1)
	top%data = asin(top%data) * r2d
case ('ACOSD')
	call math_check (1)
	top%data = acos(top%data) * r2d
case ('ATAND')
	call math_check (1)
	top%data = atan(top%data) * r2d
case ('ASINH')
	call math_check (1)
	top%data = asinh_(top%data)
case ('ACOSH')
	call math_check (1)
	top%data = acosh_(top%data)
case ('ATANH')
	call math_check (1)
	top%data = atanh_(top%data)
case ('ISNAN')
	call math_check (1)
	where (isnan(top%data))
		top%data = 1d0
	elsewhere
		top%data = 0d0
	endwhere
case ('RINT', 'NINT')
	call math_check (1)
	top%data = nint(top%data)
case ('CEIL', 'CEILING')
	call math_check (1)
	top%data = ceiling(top%data)
case ('FLOOR')
	call math_check (1)
	top%data = floor(top%data)
case ('D2R')
	call math_check (1)
	top%data = top%data * d2r
case ('R2D')
	call math_check (1)
	top%data = top%data * r2d
case ('YMDHMS')
	call math_check (1)
	top%data = ymdhms_(top%data)
	
! x MATH x y (1 argument to 2)
case ('DUP')
	call math_check (1)
	call math_push (n, top)
	top%data = top%prev%data

! x y MATH x (2 arguments to 1)
case ('ADD')
	call math_check (2)
	top%prev%data = top%prev%data + top%data
	call math_pop (top)
case ('MUL')
	call math_check (2)
	top%prev%data = top%prev%data * top%data
	call math_pop (top)
case ('DIV')
	call math_check (2)
	top%prev%data = top%prev%data / top%data
	call math_pop (top)
case ('POW')
	call math_check (2)
	top%prev%data = top%prev%data ** top%data
	call math_pop (top)
case ('FMOD')
	call math_check (2)
	top%prev%data = modulo(top%prev%data, top%data)
	call math_pop (top)
case ('MIN')
	call math_check (2)
	top%prev%data = min(top%prev%data, top%data)
	call math_pop (top)
case ('MAX')
	call math_check (2)
	top%prev%data = max(top%prev%data, top%data)
	call math_pop (top)
case ('ATAN2')
	call math_check (2)
	top%prev%data = atan2(top%prev%data, top%data)
	call math_pop (top)
case ('HYPOT')
	call math_check (2)
	top%prev%data = hypot_(top%prev%data, top%data)
	call math_pop (top)
case ('R2')
	call math_check (2)
	top%prev%data = top%prev%data*top%prev%data + top%data*top%data
	call math_pop (top)
case ('EQ')
	call math_check (2)
	where (top%prev%data == top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
case ('NEQ', 'NE')
	call math_check (2)
	where (top%prev%data /= top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
case ('LT')
	call math_check (2)
	where (top%prev%data < top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
case ('LE')
	call math_check (2)
	where (top%prev%data <= top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
case ('GT')
	call math_check (2)
	where (top%prev%data > top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
case ('GE')
	call math_check (2)
	where (top%prev%data >= top%data)
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
case ('NAN')
	call math_check (2)
	where (top%prev%data == top%data) top%prev%data = nan
	call math_pop (top)
case ('AND')
	call math_check (2)
	where (isnan(top%prev%data)) top%prev%data = top%data
	call math_pop (top)
case ('BTEST')
	call math_check (2)
	where (isnan(top%prev%data) .or. isnan(top%data))
		top%prev%data = nan
	elsewhere (btest(nint(top%prev%data),nint(top%data)))
		top%prev%data = 1d0
	elsewhere
		top%prev%data = 0d0
	endwhere
	call math_pop (top)
case ('AVG')
	call math_check (2)
	where (isnan(top%prev%data))
		top%prev%data = top%data
	elsewhere (.not.isnan(top%data))
		top%prev%data = 0.5d0 * (top%prev%data + top%data)
	endwhere
	call math_pop (top)

! x y MATH x y (2 arguments to 2)
case ('EXCH')
	call math_check (2)
	temp => top%prev
	top%prev => temp%prev
	temp%prev => top
	top => temp

! x y z MATH x (3 arguments to 1)
case ('INRANGE')
	call math_check (3)
	where (top%prev%prev%data >= top%prev%data .and. top%prev%prev%data <= top%data)
		top%prev%prev%data = 1d0
	elsewhere
		top%prev%prev%data = 0d0
	endwhere
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
stop
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
