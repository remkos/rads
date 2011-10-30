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

!*rads_nml2xml -- Convert RADS3 NML file to RADS4 XML file
!-
program rads_nml2xml
use rads
use typesizes
use rads_misc
use xmlparse

type(rads_sat) :: S
type(xml_parse) :: X
character(len=80) :: filenm
logical :: init
character(len=2) :: satnm(10) = (/ 'c2', 'e1', 'e2', 'g1', 'gs', 'j1', 'j2', 'n1', 'pn', 'tx' /)

integer :: passes(2),cycles(2)=0,i,j
character(len=4) :: radsfmt
real(eightbytereal) :: incl(1)=0,period(1)=0,freq(1)=0,dt1hz(1)=0
integer(fourbyteint), parameter :: maxsat=20
integer(fourbyteint) :: options(0:99)=0,maxalias
character(len=2) :: satpref(maxsat)
character(len=8) :: satname(maxsat)
character(len=80) :: texts(0:99),satnams(maxsat),attr(2,1),string(1)
character(len=80) :: flags(0:15),flagon(0:15),flagoff(0:15)
character(len=16) :: formts(0:99)
real(eightbytereal) :: limits(2,0:99),factors(0:99)=0,scales(0:99)=0,offsets(0:99)=0,nan
character(len=80) :: pass_fmt,cycle_fmt
character(len=1) :: phase_def=''

namelist /getraw_nml/ texts,limits,factors,options,scales,offsets, &
	satname,satpref,satnams,incl,period,freq,dt1hz,pass_fmt,cycle_fmt, &
	phase_def,flags,flagon,flagoff,formts,passes,cycles, maxalias,radsfmt

! Initialize
nan=0d0; nan = nan/nan

! Open XML on stdout
call xml_open (X, '-', .false.)
call xml_options (X, ignore_whitespace=.true.)

! Do things for all files on the command line
do i = 1,iargc()
	call getarg (i, filenm)
	init = .false.
	do j = 1,10
		if (index(filenm, 'getraw_'//satnm(j)) > 0) then
			call rads_init (S, satnm(j))
			attr(1,1) = 'sat'
			attr(2,1) = satnm(j)
			init = .true.
			exit
		endif
	enddo
	if (.not.init) call rads_init (S, 'j2')
	if (init) call xml_put (X, 'if', attr, 1, string, 0, 'open')


	! Init some variables
	formts=''
	limits = nan
	factors = nan
	factors(4) = 1d0
	factors(6:16) = -1d0
	factors(38) = -1d0
	options = 100

	! Read NML file to be converted
	open (10, file=filenm)
	read (10, nml=getraw_nml)
	close (10)
	call xml_output ()
	if (init) call xml_put (X, 'if', attr, 0, string, 0, 'close')
	call rads_end (S)
enddo

call xml_close (X)

contains

subroutine xml_output ()
integer :: i
logical :: newmath, newqual
type(rads_var), pointer :: var
character(len=5) :: field
character(len=640) :: qual, math

! Initialize
attr(1,1) = 'name'

! Parsing the variables
do i = 0,99
	if (all(isnan(limits(:,i))) .and. formts(i) == '') cycle
	if (limits(1,i) < -1d19) limits(1,i) = nan
	if (limits(2,i) >  1d19) limits(2,i) = nan

	write (field, '(i2.2)') i
	var => rads_varptr (S, field)
	attr(2,1) = var%info%name

	call xml_put (X, 'var', attr, 1, string, 0, 'open')
	if (.not.all(isnan(limits(:,i)))) then
		where (isnan(limits(:,i))) limits(:,i) = var%info%limits(:)
		call xml_dble (limits(:,i), 'limits', '2f0.3')
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
math = var%info%math

newmath = .false.
newqual = .false.
do i = 6,99
	if (options(i) == 100) cycle
	write (field, '(i2.2)') i
	var => rads_varptr (S, field)
	if (options(i) <= 0) then
		call del_string (math, trim(var%name)//' SUB ', newmath)
		call del_string (qual, trim(var%info%name)//' ', newqual)
	else if (factors(i) == 0d0) then
		write (field, '(i4.4)') i*100+options(i)
		var => rads_varptr (S, field)
		call add_string (qual, trim(var%name)//' ', newqual)
	else if (factors(i) < 0d0) then
		call add_string (math, trim(var%name)//' SUB ', newmath)
	else if (factors(i) > 0d0) then
		call add_string (math, trim(var%name)//' ADD ', newmath)
	endif
enddo

if (newmath .or. newqual) then
	attr(2,1) = 'sla'
	call xml_put (X, 'var', attr, 1, string, 0, 'open')
	if (newmath) call xml_text (math, 'math')
	if (newqual) call xml_text (qual, 'quality_flag')
	call xml_put (X, 'var', attr, 0, string, 0, 'close')
endif

end subroutine xml_output

subroutine xml_dble (var, name, format)
real(eightbytereal) :: var(:)
character(len=*) :: name, format
integer :: i,n
character(len=640) :: attr(2,1), string(size(var))
n = size(var)
if (sum(abs(var)) /= 0d0) then
	do i = 1,n
		write (string(i), '('//format//')') var(i)
	enddo
	call xml_put (X, name, attr, 0, string, n, 'elem')
endif
end subroutine xml_dble

subroutine xml_text (var, name)
character(len=*) :: var, name
character(len=640) :: attr(2,1), string(1)
if (var /= '') then
	string(1) = var
	call xml_put (X, name, attr, 0, string, 1, 'elem')
endif
end subroutine xml_text

subroutine xml_int (var, name, format)
integer(fourbyteint) :: var(:)
character(len=*) :: name,format
character(len=640) :: attr(2,1), string(1)
if (sum(abs(var)) /= 0) then
	write (string(1), '('//format//')') var
	call xml_put (X, name, attr, 0, string, 1, 'elem')
endif
end subroutine xml_int

subroutine add_string (string, sub, change)
character(len=*) :: string, sub
logical :: change
integer :: i, l
i = index(string, sub)
if (i > 0) return	! String already contains sub
l = len_trim(string)
change = .true.
string(l+1:) = ' '//sub
end subroutine add_string

subroutine del_string (string, sub, change)
character(len=*) :: string, sub
logical :: change
integer :: i, l
i = index(string, sub)
if (i == 0) return	! String does not contain sub
l = len(sub)
change = .true.
string(i:) = string(i+l:)
end subroutine del_string

end program rads_nml2xml
