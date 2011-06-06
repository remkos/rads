!*rads_nml2xml -- Convert RADS3 NML file to RADS4 XML file
!+
program rads_nml2xml
use typesizes
use rads_misc
use xmlparse

include 'config.inc'

integer :: passes(2),cycles(2)=0,i,j,k,k1,k2
character(len=4) :: radsfmt
real(eightbytereal) :: incl(1)=0,period(1)=0,freq(1)=0,dt1hz(1)=0
integer(fourbyteint), parameter :: maxsat=20
integer(fourbyteint) :: options(0:99)=0,maxalias
character(len=2) :: satpref(maxsat)
character(len=8) :: satname(maxsat)
character(len=80) :: texts(0:99),satnams(maxsat),fldnm(0:9999),attr(2,1)='name',string(1)
character(len=80) :: flags(0:15),flagon(0:15),flagoff(0:15)
character(len=16) :: formts(0:99)
real(eightbytereal) :: limits(2,0:99),factors(0:99)=0,scales(0:99)=0,offsets(0:99)=0,nan
character(len=80) :: pass_fmt,cycle_fmt
character(len=1) :: phase_def=''
type(XML_PARSE) :: X

namelist /getraw_nml/ texts,limits,factors,options,scales,offsets, &
	satname,satpref,satnams,incl,period,freq,dt1hz,pass_fmt,cycle_fmt, &
	phase_def,flags,flagon,flagoff,formts,passes,cycles, maxalias,radsfmt

namelist /raw2nc_nml/ fldnm

! Init some variables
texts=''; fldnm=''; formts=''
nan=0d0; nan = nan/nan
limits = nan

! Read NML file to be converted from standard input
read (*,nml=getraw_nml)

! Now parse all the global settings in this file
call xml_open (X, '-', .false.)
call xml_options (X, ignore_whitespace=.true.)
call xml_dble (incl, 'inclination', 'f7.3')
call xml_dble (period, 'period', 'f6.1')
call xml_dble (freq, 'frequency', 'f6.3')
call xml_dble (dt1hz, 'dt1hz', 'f5.3')
call xml_text (phase_def, 'phase')
call xml_int (cycles, 'cycles', '2i4')
call xml_int (passes, 'passes', '2i5')

! Read raw2nc.nml file
call checkenv ('RADSDATAROOT', radsdataroot)
open (10, file=trim(radsdataroot)//'/nml/raw2nc.nml')
read (10, nml=raw2nc_nml)
close (10)

! Continue parsing the variables
do i = 0,99
	if (all(isnan(limits(:,i))) .and. texts(i) == '' .and. formts(i) == '') cycle
	if (limits(1,i) < -1d19) limits(1,i) = nan
	if (limits(2,i) >  1d19) limits(2,i) = nan

	! Variables
	do j = 1,99
		if (fldnm(i*100+j) == '') cycle
		attr(2,1) = fldnm(i*100+j)
		if (j == 0 .and. options(i) /= 0) then
			string(1) = fldnm(i*100+abs(options(i)))
			call xml_put (X, 'alias', attr, 1, string, 1, 'elem')
			cycle
		endif
		call xml_put (X, 'var', attr, 1, string, 0, 'open')
		if (texts(i) /= '') then
			k1 = index(texts(i),'[')
			k2 = index(texts(i),']')
			if (k1 > 0 .and. k2 > k1) then
				call xml_text (texts(i)(:k1-2), 'long_name')
				call xml_text (texts(i)(k1+1:k2-1), 'units')
			else
				call xml_text (texts(i), 'long_name')
				call xml_text ('-', 'units')
			endif
		endif
		if (.not.all(isnan(limits(:,i)))) call xml_dble (limits(:,i), 'limits', '2f10.1')
		call xml_text (formts(i), 'format')
		call xml_text (fldnm(i*100+j), 'netcdf')
		call xml_put (X, 'var', attr, 0, string, 0, 'close')
	enddo
enddo

! Aliases
do i = 0,99
	if (fldnm(i*100) /= '' .and. options(i) /= 0) then
		attr(2,1) = fldnm(i*100)
		string(1) = fldnm(i*100+abs(options(i)))
		call xml_put (X, 'alias', attr, 1, string, 1, 'elem')
	endif
enddo

! Write out all the aliases for the numbered items
do i = 0,99
	k = abs(options(i))
	if (k /= 0) then
		write (attr(2,1), '(i2)') i
		string(1) = fldnm(i*100+k)
		call xml_put (X, 'alias', attr, 1, string, 1, 'elem')
	endif
	cycle
	do j = 1,99
		if (fldnm(i*100+j) == '') cycle
		write (attr(2,1), '(i2,i2.2)') i,j
		string(1) = fldnm(i*100+j)
		call xml_put (X, 'alias', attr, 1, string, 1, 'elem')
	enddo
enddo

call xml_close (X)

contains

subroutine xml_dble (var, name, format)
real(eightbytereal) :: var(:)
character(len=*) :: name, format
integer :: i,n
character(len=80) :: attr(2,1), string(size(var))
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
character(len=80) :: attr(2,1), string(1)
if (var /= '') then
	string(1) = var
	call xml_put (X, name, attr, 0, string, 1, 'elem')
endif
end subroutine xml_text

subroutine xml_int (var, name, format)
integer(fourbyteint) :: var(:)
character(len=*) :: name,format
character(len=80) :: attr(2,1), string(1)
if (sum(abs(var)) /= 0) then
	write (string(1), '('//format//')') var
	call xml_put (X, name, attr, 0, string, 1, 'elem')
endif
end subroutine xml_int

end program rads_nml2xml
