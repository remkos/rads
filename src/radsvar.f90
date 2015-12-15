!-----------------------------------------------------------------------
! Copyright (c) 2011-2015  Remko Scharroo
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

!*radsvar -- List RADS variables
!+
program radsvar
!
! This program lists all the variables in RADS for a given satellite
! and mission. To make sure that only the available
!
! usage: radsvar [RADS_options] [options]
!-----------------------------------------------------------------------
use typesizes
use rads
use rads_misc
use rads_time
use rads_netcdf
use netcdf

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cycle, pass, i, j, ncid, varid
type(rads_varinfo), pointer :: info
logical :: show_x = .false.

! Initialize RADS or issue help
call synopsis
call rads_set_options ('x')
call rads_init (S)
if (S%error /= rads_noerr) call rads_exit ('Fatal error')

! Scan program specific options
do i = 1,rads_nopt
	if (rads_opt(i)%opt == 'x') show_x = .true.
enddo

! Look for the first existing file
loop: do cycle = S%cycles(1),S%cycles(2)
	do pass = S%passes(1),S%passes(2)
		call rads_open_pass (S, P, cycle, pass)
		if (P%ndata > 0) exit loop
		call rads_close_pass (S, P)
	enddo
enddo loop

write (*,600) timestamp(), trim(S%command)
if (show_x) write (*,610)
600 format ('# RADS variable list'/'# Created: ',a,' UTC: ',a/ &
'#'/'# Possible output records:'/ &
'# A  field  alias     var_name             long_name'/ &
'# D  field  var_name  default_value        long_name'/ &
'# G  field  var_name  grid_file_name       long_name'/ &
'# M  field  var_name  math_rpn_statment    long_name'/ &
'# N  field  var_name  netcdf_var_name      long_name')
610 format ( &
'# X  field  var_name  non_existent_netcdf  long_name')

! If no file, let it be known
if (P%ndata == 0) write (*,'(a)') '# Could not find any corresponding file, continuing without ...'
write (*,'(a)') '#'

! List the variables
ncid = P%finfo(1)%ncid
do i = 1,S%nvar
	info => S%var(i)%info
	if (S%var(i)%name /= info%name) then
		call list ('A',S%var(i))
	else if (info%datasrc == rads_src_none) then
		call list ('U',S%var(i))
	else if (info%datasrc == rads_src_math) then
		call list ('M',S%var(i))
	else if (info%datasrc == rads_src_constant) then
		call list ('C',S%var(i))
	else if (info%datasrc == rads_src_nc_att) then
		j = index(info%dataname,':')
		if (j == 1) then
			varid = nf90_global
		else if (nft(nf90_inq_varid(ncid,info%dataname(:j-1),varid))) then
			call list ('X',S%var(i))
			cycle
		endif
		if (nft(nf90_inquire_attribute(ncid,varid,info%dataname(j+1:),xtype=j))) then
			call list ('X',S%var(i))
		else
			call list ('N',S%var(i))
		endif
	else if (info%datasrc == rads_src_nc_var) then
		if (.not.nft(nf90_inq_varid(ncid,info%dataname,varid))) then
			call list ('N',S%var(i))
		else if (info%default /= huge(0d0)) then
			call list ('D',S%var(i))
		else
			call list ('X',S%var(i))
		endif
	else
		call list ('G',S%var(i))
	endif
	if (associated(S%var(i)%inf1)) call list ('1',S%var(i))
	if (associated(S%var(i)%inf2)) call list ('2',S%var(i))
enddo

! Close the session
call rads_close_pass (S, P)
call rads_end (S)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('Make list of RADS variables')) return
call rads_synopsis
write (stderr,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -x                        Also list defined but unavailable variables')
stop
end subroutine synopsis

!***********************************************************************

subroutine list (type, var)
character(len=1), intent(in) :: type
type(rads_var), intent(in) :: var
character(len=40) :: string
select case (type)
case ('A')
	write (*, 600) type,fields(var),var%name,var%info%name,trim(var%long_name)
case ('1')
	write (*, 600)  'A',fields(var),var%name,var%inf1%name,trim(var%long_name)
case ('2')
	write (*, 600)  'A',fields(var),var%name,var%inf2%name,trim(var%long_name)
case ('D')
	write (string,*) var%info%default
	write (*, 600) type,fields(var),var%name,adjustl(string),trim(var%long_name)
case ('N', 'M', 'G')
	write (*, 600) type,fields(var),var%name,var%info%dataname,trim(var%long_name)
case default
	if (show_x) write (*, 600) type,fields(var),var%name,var%info%dataname,trim(var%long_name)
end select
600 format (a1,1x,a7,1x,a25,1x,a40,1x,a)
end subroutine list

function fields (var)
type(rads_var), intent(in) :: var
character(len=7) :: fields
if (var%field(1) == rads_nofield) then
	fields = ''
else if (var%field(2) /= rads_nofield) then
	write (fields, '(i2,i5)') var%field(2:1:-1)
else if (var%field(1) < 100) then
	write (fields, '(i2,5x)') var%field(1)
else
	write (fields, '(2x,i5)') var%field(1)
endif
end function fields

end program radsvar
