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

module rads_misc
use typesizes

! There are used by getopt
character(len=160), private, save :: getopt_arg
logical, private, save :: getopt_new = .true.
integer, private, save :: getopt_ind = 1, getopt_chr = 2

!***********************************************************************
!*d_int -- Convert integer to double while accounting for NaNs
!+
! elemental function d_int (i)
! integer(onebyteint) <or> integer(twobyteint) <or> integer(fourbyteint) :: i
! real(eightbytereal) :: d_int
!
! Convert an integer (or integer array) to a double, while taking into
! account that the maximum value for that type of integer (127, 32767, etc)
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
!*getopt -- Get option from command line argument
!
integer function getopt (opts, optopt, optarg, unit, opterr)
character(len=*), intent(in) :: opts
character(len=*), intent(out) :: optopt, optarg
integer, intent(in), optional :: unit
logical, intent(in), optional :: opterr
!
! This function mimics both the GNU C <getopt> and <getopt_long> functions,
! and provides several enhancements over those functions (see below).
!
! At every call, <getopt> reads a new option and (if required) its
! option argument. By default, <getopt> reads from the command line,
! but if <unit> is specified, it will read from a file containing options
! that is associated with the i/o unit <unit>.
!
! The argument <opts> may contain both short options (i.e. '-f') or long
! options (i.e. '--file'), in the following way.
! 1) Start with a string of single character short options, optionally
!    insert a ':' when the option requires an argument.
! 2) Append, each separated by a space the long options, ending in ':'
!    when the option requires an argument.
!
! Example:
! opts = 'f:qo:v quiet verbose file: output: outcome help'
!
! Note 1: There needs to be no association between short and long options
!    This is different from <getopt_long> and actually more practical in
!    Fortran, because we can easily parse long and short options later with
!    case ('v', 'verbose')
! Note 2: Long options are not required, simply leave them out if not needed.
! Note 3: Short options are not required, leave them out if not needed, but
!    then still start <opts> with a space.
!
! The argument <optopt> returns the short or the (full) long option. This
! routine will match any shorter unique version of the long option but still
! return the full one. At the same time, <optarg> will provide the argument
! to the option.
!
! Examples of input, optopt and optarg using opts above:
! -q           optopt = 'q', optarg = ''
! --quiet      optopt = 'quiet', optarg = ''
! --quiet=0    optopt = 'quiet', optarg = '' (Addition =0 is ignored)
! -fFILE       optopt = 'f', optarg = 'FILE'
! -f FILE      optopt = 'f', optarg = 'FILE' (Two arguments joined)
! --file=FILE  optopt = 'file', optarg = 'FILE'
! --out FILE   optopt = 'output', optarg = 'FILE' (Returning first match)
! file.txt     optopt = '', optarg = 'file.txt' (No match)
!
! The function <getopt> returns one of the following error codes:
! -3) An option with a required argument was found, but there was no
!     required argument, because of reaching the end of the input or
!     because the option was followed by another option.
!     <optopt> returns the option, <optarg> returns ''.
! -2) An argument starts with '--' and no match with an option is found.
!     <optopt> returns '?' and <optarg> the full argument.
! -1) An argument starts with '-' and no match with an option is found.
!     <optopt> returns '?' and <optarg> the full argument.
!  0) No error or warning. But if no match was found, <optopt> returns ''
!     and <optarg> the full argument.
!  1) And of input is reached. <optopt> and <optarg> remain unchanged.
!
! Arguments:
!   opts     : list of short and long options
!   optopt   : option part of the option, without -- or -
!   optarg   : argument of the option
!   unit     : (Optional) input unit number, otherwise command line
!   opterr   : (Optional) print error messages with .true.
!-----------------------------------------------------------------------
integer :: ios = 0, input, i, j, k, n
logical :: print_err

! Set default return code
getopt = 0

! Select input unit
if (.not.present(unit)) then
	input = 0
else
	input = unit
endif

! Determine if error printing is necessary
if (.not.present(opterr)) then
	print_err = .false.
else
	print_err = opterr
endif

! Read next argument when needed
if (getopt_new) call nextarg

! Return on error
if (ios /= 0) then
	getopt = 1
	return
endif

! Arguments '-' and '--' by itself should be regarded as normal arguments
if (getopt_arg(1:2) == '- ' .or. getopt_arg(1:3) == '-- ') then
	optopt = ' '
	optarg = getopt_arg
	return
endif

! Default is no option argument
optarg = ''

! Handle double-dash options
if (getopt_arg(1:2) == '--') then
	getopt_new = .true.
	i = len_trim(getopt_arg)
	j = scan (getopt_arg, ' :=')
	if (j > 3) i = j - 1
	n = index(opts, ' ' // getopt_arg(3:i))
	if (n == 0) then
		getopt = -2
		optopt = '?'
		optarg = getopt_arg
		return
	endif
	n = n + 1
	k = scan(opts(n:), ' :') + n - 2
	optopt = opts(n:k)
	if (opts(k+1:k+1) /= ':') return
	if (j > 3) then
		optarg = adjustl(getopt_arg(j+1:))
		return
	endif
	call nextarg
	if (ios /= 0) then
		if (print_err) print '(a,a,a)', &
			'Warning: option "',trim(getopt_arg),'" should contain argument, but end reached'
		getopt = -3
	else if (getopt_arg(1:1) == '-' .and. getopt_arg(2:2) /= ' ') then
		if (print_err) print '(a,a,a)', &
			'Warning: option "',trim(getopt_arg),'" should contain argument, but new option reached'
		getopt = -3
		getopt_new = .false.
	else
		optarg = getopt_arg
	endif
	return
endif

! Handle single-dash options
if (getopt_arg(1:1) == '-') then
	getopt_new = .false.
	optopt = getopt_arg(getopt_chr:getopt_chr)
	getopt_chr = getopt_chr + 1
	if (getopt_arg(getopt_chr:) == '') getopt_new = .true.
	n = 0
	do
		n = n + 1
		if (opts(n:n) == ' ') exit
		if (opts(n:n) /= optopt(:1)) cycle
		if (opts(n+1:n+1) /= ':') return
		if (.not.getopt_new) then
			optarg = adjustl(getopt_arg(getopt_chr:))
			getopt_new = .true.
			return
		endif
		call nextarg
		if (ios /= 0) then
			if (print_err) print '(a,a,a)', &
				'Warning: option "-',trim(optopt),'" should contain argument, but end reached'
			getopt = -3
		else if (getopt_arg(1:1) == '-' .and. getopt_arg(2:2) /= ' ') then
			if (print_err) print '(a,a,a)', &
				'Warning: option "-',trim(optopt),'" should contain argument, but new option reached'
			getopt = -3
			getopt_new = .false.
		else
			optarg = getopt_arg
			getopt_new = .true.
		endif
		return
	enddo
	getopt = -1
	optopt = '?'
	optarg = getopt_arg
	return
endif

! For compatibility, handle options like lat=0,180
j = index(getopt_arg, '=')
if (j > 1) then
	getopt_new = .true.
	j = scan (getopt_arg, ':=') - 1
	n = index (opts, ' ' // getopt_arg(:j))
	if (n > 0) then
		n = n + 1
		k = scan (opts(n:), ' :') + n - 2
		optopt = opts(n:k)
		optarg = getopt_arg(j+2:)
		return
	endif
endif

! Return simple argument
optopt = ' '
optarg = getopt_arg

contains

! Get the next argument from command line or file
subroutine nextarg ()
if (input > 0) then
	read (input, '(a)', iostat=ios) getopt_arg
else if (getopt_ind > iargc()) then
	ios = 1
else
	call getarg (getopt_ind, getopt_arg)
	getopt_ind = getopt_ind + 1
	ios = 0
endif
getopt_chr = 2
end subroutine nextarg

end function getopt

!***********************************************************************
!*getopt_reset -- Reset to beginning of command line options
!+
subroutine getopt_reset (unit)
integer, intent(in), optional :: unit
!
! By calling this routine the pointer for <getopt> is reset to the
! beginning of the argument list, or the file on input <unit> is rewound
! to the beginning.
!
! Argument:
!   unit  : (Optional) unit of the input file, command line otherwise
!-----------------------------------------------------------------------
if (present(unit)) then
	rewind (unit)
else
	getopt_ind = 1
endif
getopt_chr = 2
getopt_new = .true.
end subroutine getopt_reset

!***********************************************************************
!*strtolower -- Convert string to lower case
!+
function strtolower (string) result (lower)
character(len=*), intent(in) :: string
character(len=len(string)) :: lower
!
! Convert character string to all lower case letters, leaving all non-
! alphabetical characters as is.
!
! Arguments:
!  string : String to be converted
!  lower  : String returned in lower case
!-----------------------------------------------------------------------
integer :: i
do i = 1,len(string)
	if (string(i:i) >= 'A' .and. string(i:i) <= 'Z') then
		lower(i:i) = achar(iachar(string(i:i)) + 32)
	else
		lower(i:i) = string(i:i)
	endif
enddo
end function strtolower

!***********************************************************************
!*strtoupper -- Convert string to upper case
!+
function strtoupper (string) result (upper)
character(len=*), intent(in) :: string
character(len=len(string)) :: upper
!
! Convert character string to all upper case letters, leaving all non-
! alphabetical characters as is.
!
! Arguments:
!  string : String to be converted
!  upper  : String returned in upper case
!-----------------------------------------------------------------------
integer :: i
do i = 1,len(string)
	if (string(i:i) >= 'a' .and. string(i:i) <= 'z') then
		upper(i:i) = achar(iachar(string(i:i)) - 32)
	else
		upper(i:i) = string(i:i)
	endif
enddo
end function strtoupper

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
	j = index(from,string(i:i))
	if (j >= 1 .and. j <= n) string(i:i) = to(j:j)
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
!   env    : Name of the environment variable.
!   string : Upon input: default value for string.
!          : Upon output: contents of the environment variable or default.
!-----------------------------------------------------------------------
character(len=160) :: temp
call getenv (env,temp)
if (temp /= ' ') string = temp
end subroutine checkenv

!***********************************************************************
!*parseenv -- Parse a string for embedded environment variable
!+
recursive subroutine parseenv (input, output)
character(len=*), intent(in) :: input
character(len=*), intent(inout) :: output
!
! This routine parse the string <input> for embedded environment
! variables and replaces them with its value. The result is returned
! in <output>.
! Environment variables are of the form ${parameter}.
! Alternatively use ${parameter:-word} to substitude the expansion of
! word when parameter is unset or empty; or use ${parameter:+word}
! to substitute the expansion of word when parameter is not empty
! (same meaning as in the bash shell).
! Nesting of variable names is not allowed.
!
! Arguments:
!   input  : String to be parsed.
!   output : String with environment variables replaced.
!-----------------------------------------------------------------------
character(len=160) :: env
integer :: j, k = 0, l, m, n
output = input
do
	l = len_trim(output)
	j = index (output(1:), '${')
	if (j == 0) exit
	j = j + 2
	n = 1
	do k = j,l
		if (output(k:k) == '}') then
			n = n - 1
			if (n == 0) exit
		else if (output(k:k) == '{') then
			n = n + 1
		endif
	enddo
	if (n > 0) exit
	k = k - 1
	m = index(output(j:k), ':')
	if (m == 0) then
		call getenv (output(j:k),env)
	else
		m = j + m - 1
		call getenv (output(j:m-1),env)
		if (index ('-=?', output(m+1:m+1)) > 0) then
			if (env == '') env = output(m+2:k)
		else if (output(m+1:m+1) == '+') then
			if (env /= '') env = output(m+2:k)
		endif
	endif
	output = output(:j-3) // trim(env) // output(k+2:)
enddo
end subroutine parseenv

!***********************************************************************
!*outofrange - Verify if value is within limits
!+
function outofrange (limits, value)
real(eightbytereal), intent(in) :: limits(2)
real(eightbytereal), intent(inout) :: value
logical :: outofrange
!
! This function checks if value is within <limits(1:2)>, where
! <limits(1)> is less than <limits(2)>. If <value> is not within those
! limits, <outofrange> is set to .true.
!
! If either of the limits is NaN, no check is performed at that limit.
!
! Arguments:
!  limits     : minimum and maximum allowed value
!  value      : value to be checked
!  outofrange : .true. if value is outside limits, .false. otherwise.
!-----------------------------------------------------------------------
outofrange = (value < limits(1) .or. value > limits(2))
end function outofrange

!***********************************************************************
!*nint1 -- Round 8-byte real to 1-byte integer
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
!*nint2 -- Round 8-byte real to 2-byte integer
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
!*nint4 -- Round 8-byte real to 4-byte integer
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
!-----------------------------------------------------------------------
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
!
! This function returns a double float NaN scalar
!-----------------------------------------------------------------------
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

!***********************************************************************
!*splitarg -- Split command line argument in option and option-argument
!+
subroutine splitarg (arg, option, optarg, pos)
character(len=*), intent(in) :: arg
character(len=*), intent(out) :: option, optarg
integer(fourbyteint), optional :: pos
!
! This routine splits a command line argument <arg> in a form like
! arg:    '-A'  '-x1'  '--lon=0,10'  'sel=sla'  lim:sla=-1,1  'filename'
! into and option (-? or --*) and an argument. For the above examples:
! option: '-A'  '-x'   '--lon'       '--sel'    '--lim'       ''
! optarg: ''    '1'    '0,10'        'sla'      'sla=-1,1'    ''
!
! Thus the <option> is either -? (dash followed by a single character), or
! --* (double-dash followed by one or more characters).
! Note that for the old-style argument 'sel=', the prefix '--' is added.
! When this is not an recognisable option, a blank string is returned.
!
! The <optarg>, or optional argument, is what follows -? or the ':' or the
! '=' sign (whatever comes first).
! When there is no optional argument, or arg is not a recognisable option,
! a blank string is returned.
!
! Finally, <pos> returns the position of the equal sign in optarg.
!
! Arguments:
!  arg    : Command line argument
!  option : Identifier of the option
!  optarg : Argument of the option
!  pos    : Optional: position of the equal sign in optarg
!-----------------------------------------------------------------------
integer :: i, j
i = index(arg, ':')
j = index(arg, '=')
if (arg(:2) == '--') then
	if (i > 0) then
		option = arg(:i-1)
		optarg = arg(i+1:)
		j = j - i
	else if (j > 0) then
		option = arg(:j-1)
		optarg = arg(j+1:)
		j = 0
	else
		option = arg
		optarg = ''
	endif
else if (arg(:1) == '-') then
	option = arg(:2)
	optarg = arg(3:)
	j = j - 2
else if (j == 0) then
	option = ''
	optarg = ''
else if (i > 0) then
	option = '--' // arg(:i-1)
	optarg = arg(i+1:)
	j = j - i
else
	option = '--' // arg(:j-1)
	optarg = arg(j+1:)
	j = 0
endif
if (present(pos)) pos = j
end subroutine splitarg

!***********************************************************************
!*iqsort -- Make a quick sort of an INTEGER*4 array
!+
subroutine iqsort (idx, val, nr)
use typesizes
integer(fourbyteint), intent(inout) :: idx(nr)
integer(fourbyteint), intent(in) :: val(:), nr
!
! This routine quick-sorts an integer array <idx> according to the
! corresponding value in array <val> of type integer.
! The array <idx> contains <nr> pointers to the values in array <val>,
! which does not have to be equally large.
!
! Arguments:
!  idx (input) : One-dimensional array of index numbers
!     (output) : Row of index numbers sorted on corresponding value
!  val (input) : One-dimensional array of values
!  nr  (input) : Number of indices/values
!
! Example 1:
! Assume you have 6 values 100, 10, 11, 21, 17, and 90, stored in array
! <val>. These values have to be sorted in ascending order. Your input
! will be like this:
!
!   idx    1    2    3    4    5    6
!   val  100   10   11   21   17   90
!    nr    6
!
! After running iqsort, the output will look like this:
!
!   idx    2    3    5    4    6    1
!   val  100   10   11   21   17   90
!    nr    6
!
! Note that the array <val> has not been changed, but that <idx> is sorted by
! the order in which <val> should be sorted. You may print the values in
! the ascending order as follows:
!
!     do i=1,nr
!        write (*,*) val(idx(i))
!     enddo
!
! Example 2:
! It is also possible that the indices are not in logical order. As long as
! they point to the respective locations of the values in array <val>, the
! sorting will be done accordingly. Assume you have the following input:
! (values ** are irrelevant)
!
!   idx    1    2    4    5    6    8
!   val  100   10   **   21   17   90   **   1
!    nr    6
!
! Then after running iqsort, the result will be:
!
!   idx    8    2    5    4    6    1
!   val  100   10   **   21   17   90   **   1
!    nr    6
!
! Printing the values in ascending order after this goes the same as in
! Example 1.
!-----------------------------------------------------------------------
integer(fourbyteint), parameter :: m = 7, nstack = 2000
integer(fourbyteint) :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
integer(fourbyteint) :: a
jstack = 0
l = 1
ir = nr
do
	if (ir-l < m) then
		jloop: do j = l+1,ir
			indxt = idx(j)
			a = val(indxt)
			do i = j-1,1,-1
				if (val(idx(i)) < a) then
					idx(i+1) = indxt
					cycle jloop
				endif
				idx(i+1) = idx(i)
			enddo
			i = 0
			idx(i+1)=indxt
		enddo jloop
		if (jstack == 0) return
		ir = istack(jstack)
		l = istack(jstack-1)
		jstack = jstack-2
	else
		k = (l+ir) / 2
		itemp = idx(k)
		idx(k) = idx(l+1)
		idx(l+1) = itemp
		if (val(idx(l+1)) > val(idx(ir))) then
			itemp = idx(l+1)
			idx(l+1) = idx(ir)
			idx(ir) = itemp
		endif
		if (val(idx(l)) > val(idx(ir))) then
			itemp = idx(l)
			idx(l) = idx(ir)
			idx(ir) = itemp
		endif
		if (val(idx(l+1)) > val(idx(l))) then
			itemp = idx(l+1)
			idx(l+1) = idx(l)
			idx(l) = itemp
		endif
		i = l+1
		j = ir
		indxt = idx(l)
		a = val(indxt)
		do
			do
				i = i+1
				if (val(idx(i)) >= a) exit
			enddo
			do
				j = j-1
				if (val(idx(j)) <= a) exit
			enddo
			if (j < i) exit
			itemp=idx(i)
			idx(i)=idx(j)
			idx(j)=itemp
		enddo
		idx(l) = idx(j)
		idx(j) = indxt
		jstack = jstack+2
		if (jstack > nstack) stop 'NSTACK too small in IQSORT'
		if (ir-i+1 >= j-1) then
			istack(jstack) = ir
			istack(jstack-1) = i
			ir = j-1
		else
			istack(jstack) = j-1
			istack(jstack-1) = l
			l = i
		endif
	endif
enddo
end subroutine iqsort

!***********************************************************************
!*mean_variance -- Compute mean and standard deviation of a series
!
subroutine mean_variance (x, mean, variance)
real(eightbytereal), intent(in) :: x(:)
real(eightbytereal), intent(out) :: mean, variance
!
! This routine computes the mean and variance of a series
! of values <x> using the method proposed by West (1979), since it is
! more stable in the way it computes the variance.
!
! This routine does not check for NaN values.
!
! Reference:
!  West, D. H. D
!  Updating mean and variance estimates: an improved method
!  Communications of the ACM, 22(9), 532-535, 1979
!
! Arguments:
!  x        : Series of values
!  mean     : Mean of series <x>
!  variance : Variance of series <x>
!-----------------------------------------------------------------------
real(eightbytereal) :: q, r, sum2
integer(fourbyteint) :: i, n
n = size(x)
if (n == 0) then
	mean = make_nan()
	variance = mean
	return
else if (n == 1) then
	mean = x(1)
	variance = make_nan()
	return
endif
mean = 0d0
sum2 = 0d0
do i = 1,n
	q = x(i) - mean
	r = q / i
	mean = mean + r
	sum2 = sum2 + r * q * (i-1)
enddo
variance = sum2/(n-1)
end subroutine mean_variance

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
