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

!*radsreconfig -- Convert RADS3 RMF and NML files to RADS4 XML file
!+
! This program reads the getraw.rmf and getraw.nml files used in RADS3
! and converts them to the rads.xml configuration file used in RADS4.
!
! When called without arguments, radsreconfig reads all files of the form:
! getraw.{rmf,nml} and getraw_$sat.{rmf,nml} in the current directory
! and writes out a single file rads.xml.
!
! When called with one of more filenames as arguments, the program combines
! those and provides output to stdout. 
!-----------------------------------------------------------------------
program radsreconfig
use rads
use typesizes
use rads_misc
use xmlparse

type(rads_sat) :: S
type(xml_parse) :: X
character(len=80) :: arg
integer, parameter :: nsat = 10
character(len=2) :: satnm(nsat+1) = (/ 'c2', 'e1', 'e2', 'g1', 'gs', 'j1', 'j2', 'n1', 'pn', 'tx', '??' /)

integer :: passes(2),cycles(2)=0,i,j
character(len=4) :: radsfmt
real(eightbytereal) :: incl(1)=0,period(1)=0,freq(1)=0,dt1hz(1)=0
integer(fourbyteint), parameter :: maxsat=20
integer(fourbyteint) :: options(0:99)=0,maxalias
character(len=2) :: satpref(maxsat)
character(len=8) :: satname(maxsat)
character(len=80) :: texts(0:99),satnams(maxsat),attr(2,3),string(1)
character(len=80) :: flags(0:15),flagon(0:15),flagoff(0:15)
character(len=16) :: formts(0:99)
real(eightbytereal) :: limits(2,0:99),factors(0:99)=0,scales(0:99)=0,offsets(0:99)=0
character(len=80) :: pass_fmt,cycle_fmt
character(len=1) :: phase_def=''
character(len=6) :: grid_type(3) = (/'grid  ','grid_c','grid_q'/)

namelist /getraw_nml/ texts,limits,factors,options,scales,offsets, &
	satname,satpref,satnams,incl,period,freq,dt1hz,pass_fmt,cycle_fmt, &
	phase_def,flags,flagon,flagoff,formts,passes,cycles, maxalias,radsfmt

! Usage message
1300 format ('radsreconfig -- Convert RADS3 NML and RMF files to RADS4 rads.xml' // &
'usage: radsconfig [ FILENAME ... ]' // &
'When used without arguments:'/ &
'  Search directory for NML and RMF config files and create rads.xml' / &
'When used with arguments:'/ &
'  Convert files on command line and send resulting XML code to standard output')

! Check command line arguments
call getarg (1, arg)
if (arg == '--help' .or. arg == '-?') then
	write (*,1300)
	stop
else if (arg == '') then
	call xml_open (X, 'rads.xml', .false.)
	call xml_options (X, ignore_whitespace=.true.)
	call rads_setup ('')
	call convert_rmf ('getraw.rmf')
	call convert_nml ('getraw.nml')
	call rads_end (S)
	do j = 1,nsat
		call rads_setup ('getraw_'//satnm(j))
		call convert_rmf ('getraw_'//satnm(j)//'.rmf')
		call convert_nml ('getraw_'//satnm(j)//'.nml')
		call rads_end (S)
	enddo
else
	call xml_open (X, '-', .false.)
	call xml_options (X, ignore_whitespace=.true.)
	do i = 1,iargc()
		call getarg (i, arg)
		call rads_setup (arg)
		if (index(arg,'.rmf') > 0) call convert_rmf (arg)
		if (index(arg,'.nml') > 0) call convert_nml (arg)
		call rads_end (S)
	enddo
endif

call xml_close (X)

contains

subroutine rads_setup (filename)
character(len=*) :: filename
integer :: j

! Check if this file is satellite-specific
do j = 1,nsat
	if (index(filename, 'getraw_'//satnm(j)) > 0) exit
enddo

! Set up the S structure
call rads_init_sat_struct (S)
allocate (S%var(rads_var_chunk))
S%var = rads_var ('', null(), .false., rads_nofield)
call rads_read_xml (S, trim(radsdataroot) // '/conf/rads.xml')
if (S%error == rads_err_xml_file) call rads_exit ('Required XML file '//trim(radsdataroot)//'/conf/rads.xml does not exist')
call rads_read_xml (S, trim(radsuserroot) // '/.rads/rads.xml')
end subroutine rads_setup

subroutine convert_rmf (filename)
character(len=*) :: filename
integer :: i, j, ix, iy, type
type(rads_var), pointer :: var
character(len=5) :: field
character(len=640) :: line, math(1), long_name, gridnm(1)
logical :: dummy 

! Open RMF file
open (10, file=filename, status='old', iostat=i)
if (i /= 0) return

! Output XML
attr(1,1) = 'sat'
attr(2,1) = S%sat
if (S%sat /= '??') call xml_put (X, 'if', attr, 1, string, 0, 'open')

do
	read (10, '(a)', iostat=i) line
	if (i /= 0) exit
	if (line(:4) == 'MATH') then
		read (line(5:), *) i, math, long_name

		! Create new variable
		write (field, '(i4.4)') i
		var => rads_varptr (S, field, null())
		attr(1,1) = 'name'
		attr(2,1) = var%name
		attr(1,2) = 'field'
		attr(2,2) = field
		call xml_put (X, 'var', attr, 2, string, 0, 'open')
		call xml_long_name (long_name)		


		! Add math string
		call replace_string (math(1))
		call del_string (math(1), ' =', dummy)
		attr(1,1) = 'source'
		attr(2,1) = 'math'
		call xml_put (X, 'data', attr, 1, math, 1, 'elem')
		call xml_put (X, 'var', attr, 0, string, 0, 'close')

	else if (line(:4) == 'GRID') then
		if (line(5:5) == 'L') then
			type = 1
			read (line(6:), *) i, ix, iy, gridnm, long_name
		else if (line(5:5) == 'S') then
			type = 2
			read (line(6:), *) i, ix, iy, gridnm, long_name
		else
			read (line(5:), *) i, type, ix, iy, gridnm, long_name
		endif

		! Create new variable
		write (field, '(i4.4)') i
		var => rads_varptr (S, field, null())
		attr(1,1) = 'name'
		attr(2,1) = var%name
		attr(1,2) = 'field'
		attr(2,2) = field
		call xml_put (X, 'var', attr, 2, string, 0, 'open')
		call xml_long_name (long_name)
		write (field, '(i4.4)') ix
		var => rads_varptr (S, field)
		attr(1,1) = 'x'
		attr(2,1) = var%name
		write (field, '(i4.4)') iy
		var => rads_varptr (S, field)
		attr(1,2) = 'y'
		attr(2,2) = var%name
		attr(1,3) = 'source'
		attr(2,3) = trim(grid_type(type))
		call xml_put (X, 'data', attr, 3, gridnm, 1, 'elem')
		call xml_put (X, 'var', attr, 0, string, 0, 'close')

	else if (line(:5) == 'ALIAS') then
		read (line(6:), *) ix, iy, long_name

		! Create new variable
		j = S%nvar
		write (field, '(i4.4)') ix
		var => rads_varptr (S, field, null())
		attr(1,1) = 'name'
		attr(1,2) = 'field'
		attr(2,2) = field
		if (S%nvar == j) then	! Variable already existed
			attr(2,1) = var%name
			call xml_put (X, 'var', attr, 2, string, 0, 'open')
			write (field, '(i4.4)') iy
			var => rads_varptr (S, field, null())
			call xml_text (var%name, 'backup')
			call xml_put (X, 'var', attr, 0, string, 0, 'close')
		else
			attr(2,1) = 'f' // field
			write (field, '(i4.4)') iy
			var => rads_varptr (S, field, null())
			string(1) = var%name
			call xml_put (X, 'alias', attr, 2, string, 1, 'elem')
		endif
	endif
enddo

if (S%sat /= '??') call xml_put (X, 'if', attr, 0, string, 0, 'close')
close (10)
end subroutine convert_rmf

subroutine convert_nml (filename)
character(len=*) :: filename
integer :: i
logical :: newmath, newqual
type(rads_var), pointer :: var
character(len=5) :: field
character(len=640) :: qual, math(1)

! Init some variables
formts=''
limits = make_nan()
factors = 0d0
factors(4) = 1d0
factors(6:16) = -1d0
factors(38) = -1d0
options = 100

! Read NML file to be converted
open (10, file=filename, status='old', iostat=i)
if (i /= 0) return
read (10, nml=getraw_nml, iostat=i)
close (10)
if (i /= 0) return

! Output XML
attr(1,1) = 'sat'
attr(2,1) = S%sat
if (S%sat /= '??') call xml_put (X, 'if', attr, 1, string, 0, 'open')

! Parsing the variables
attr(1,1) = 'name'
do i = 0,99
	if (all(isnan(limits(:,i))) .and. formts(i) == '') cycle
	if (limits(1,i) < -1d19) limits(1,i) = make_nan()
	if (limits(2,i) >  1d19) limits(2,i) = make_nan()

	write (field, '(i2.2)') i
	var => rads_varptr (S, field)
	attr(2,1) = var%info%name

	call xml_put (X, 'var', attr, 1, string, 0, 'open')
	if (.not.all(isnan(limits(:,i)))) then
		where (isnan(limits(:,i))) limits(:,i) = var%info%limits(:)
		call xml_dble (limits(:,i), 'limits', 'f0.3')
	endif
	if (formts(i) /= '') call xml_text (formts(i), 'format')
	call xml_put (X, 'var', attr, 0, string, 0, 'close')
enddo

! Aliases
do i = 0,99
	if (options(i) == 0 .or. options(i) == 100) cycle

	write (field, '(i2.2)') i
	var => rads_varptr (S, field)
	attr(2,1) = var%name
	string(1) = var%info%name
	call xml_put (X, 'alias', attr, 1, string, 1, 'elem')
enddo

! Load math and quality strings
var => rads_varptr (S, 'sla')
qual = var%info%quality_flag
math(1) = ''
if (var%info%datasrc == rads_src_math) math(1) = var%info%dataname

newmath = .false.
newqual = .false.
attr(1,1) = 'name'
attr(2,1) = 'sla'
attr(1,2) = 'action'

do i = 6,99
	if (options(i) == 100) cycle
	write (field, '(i2.2)') i
	var => rads_varptr (S, field)
	if (options(i) <= 0) then
		call del_string (math(1), trim(var%name)//' SUB ', newmath)
		string(1) = var%name
		attr(2,2) = 'delete'
		if (.not.newqual) call xml_put (X, 'var', attr, 1, string, 0, 'open')
		call xml_put (X, 'quality_flag', attr(:,2:2), 1, string, 1, 'elem')
		newqual = .true.
	else if (factors(i) == 0d0) then
		string(1) = var%name
		attr(2,2) = 'merge'
		if (.not.newqual) call xml_put (X, 'var', attr, 1, string, 0, 'open')
		call xml_put (X, 'quality_flag', attr(:,2:2), 1, string, 1, 'elem')
		newqual = .true.
	else if (factors(i) < 0d0) then
		call add_string (math(1), trim(var%name)//' SUB ', newmath)
	else if (factors(i) > 0d0) then
		call add_string (math(1), trim(var%name)//' ADD ', newmath)
	endif
enddo

if (newmath) then
	if (.not.newqual) call xml_put (X, 'var', attr, 1, string, 0, 'open')
	attr(1,1) = 'source'
	attr(2,1) = 'math'
	call xml_put (X, 'data', attr, 1, math, 1, 'elem')
endif
if (newqual.or.newmath) call xml_put (X, 'var', attr, 0, string, 0, 'close')

if (S%sat /= '??') call xml_put (X, 'if', attr, 0, string, 0, 'close')
end subroutine convert_nml


subroutine xml_long_name (long_name)
character(len=*) :: long_name
integer :: k1, k2
k1 = index(long_name,'[')
k2 = index(long_name,']')
if (k1 > 0 .and. k2 > k1) then
	call xml_text (long_name(:k1-2), 'long_name')
	call xml_text (long_name(k1+1:k2-1), 'units')
else
	call xml_text (long_name, 'long_name')
	call xml_text ('-', 'units')
endif
end subroutine xml_long_name

subroutine xml_dble (var, name, format)
real(eightbytereal) :: var(:)
character(len=*) :: name, format
integer :: i,n
character(len=640) :: attr(2,1), string(size(var))
n = size(var)
do i = 1,n
	if (abs(var(i)-nint(var(i))) < 1d-10) then
		write (string(i), '(i20)') nint(var(i))
	else
		write (string(i), '('//format//')') var(i)
	endif
enddo
call xml_put (X, name, attr, 0, string, n, 'elem')
end subroutine xml_dble

subroutine xml_text (var, name)
character(len=*) :: var, name
character(len=640) :: attr(2,1), string(1)
if (var == '') return
string(1) = var
call xml_put (X, name, attr, 0, string, 1, 'elem')
end subroutine xml_text

subroutine add_string (string, sub, change)
character(len=*) :: string, sub
logical :: change
integer :: i, l
l = len_trim(sub)
i = index(string, sub(:l))
if (i > 0) return	! String already contains sub
change = .true.
string(l+1:) = ' '//sub
end subroutine add_string

subroutine del_string (string, sub, change)
character(len=*) :: string, sub
logical :: change
integer :: i, l
l = len(sub)
i = index(string, sub)
if (i == 0) return	! String does not contain sub
change = .true.
string(i:) = string(i+l:)
end subroutine del_string

subroutine replace_string (string)
character(len=*) :: string
integer :: i, j
type (rads_var), pointer :: var
character(len=5) :: field
character(len=1) :: insert
do
	i = index(string, '%')
	insert = ''
	if (i == 0) then
		i = index(string, '$')
		insert = rads_noedit
	endif
	if (i == 0) return
	read (string(i+1:), *) j
	write (field, '(i4.4)') j
	var => rads_varptr (S, field)
	j = index(string(i:), ' ')
	if (associated(var)) then
		string(i:) = trim(var%name) // trim(insert) // string(i+j-1:)
	else
		string(i:) = 'f' // trim(field) // trim(insert) // string(i+j-1:)
	endif
enddo
end subroutine replace_string

end program radsreconfig
