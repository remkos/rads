!-----------------------------------------------------------------------
! Copyright (c) 2011-2016  Remko Scharroo
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
integer(fourbyteint) :: cycle, pass, i, ncid, ios, intervals = 0
type(rads_var), pointer :: var
logical :: show_x = .false., round = .false.
real(eightbytereal) :: factor = 1d0
character(len=rads_varl) :: prefix = ''

! Initialize RADS or issue help
call synopsis
call rads_set_options ('d:i:j:p:x')
call rads_init (S)
if (S%error /= rads_noerr) call rads_exit ('Fatal error')

! Scan program specific options
nullify (var)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('V', 'var', 'sel')
		var => rads_varptr(S, rads_opt(i)%arg)
	case ('d')
		read (rads_opt(i)%arg, *, iostat=ios) factor
	case ('i')
		read (rads_opt(i)%arg, *, iostat=ios) intervals
	case ('j')
		read (rads_opt(i)%arg, *, iostat=ios) intervals
		round = .true.
	case ('p')
		prefix = rads_opt(i)%arg(:rads_varl)
	case ('x')
		show_x = .true.
	end select
enddo

! Look for the first existing file
loop: do cycle = S%cycles(1),S%cycles(2)
	do pass = S%passes(1),S%passes(2)
		call rads_open_pass (S, P, cycle, pass)
		if (P%ndata > 0) exit loop
		call rads_close_pass (S, P)
	enddo
enddo loop

if (associated(var)) then
	! If specific variable if given, print variables for insertion in bash script
	call list_variable
else
	! If no variable is given, list all variables
	call list_variables
endif

! Close the session
call rads_close_pass (S, P)
call rads_end (S)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('Make list of RADS variables')) return
call rads_synopsis
write (*,1300)
1300 format (/ &
'With the option "-V VAR", program specific [program_options] are:'/ &
'  -d FACTOR                 Scale the plotting range for differences by FACTOR (def: 1)'/ &
'  -i INTERVALS              Select nr of intervals to determine a plotting bin size'/ &
'  -j INTERVALS              Same as -i INTERVAL and round the bin size'/ &
'  -p PREFIX                 Set prefix for script variables'// &
'Without the option "-V VAR", program specific [program_options] are:'/ &
'  -x                        Also list defined but unavailable variables')
stop
end subroutine synopsis

!***********************************************************************
! List attributes of single variable, to be used in bash script
subroutine list_variable
type(rads_varinfo), pointer :: info
real(eightbytereal) :: plot_bin, diff_range(2), diff_bin
info => var%info

! If no plot_range then copy limits
if (all(isnan_(info%plot_range))) info%plot_range = info%limits

! Centre and scale plot_range to get diff_range for plotting differences
diff_range = factor * (info%plot_range - 0.5d0 * sum(info%plot_range))

! Compute bin size and adjust range if necessary
call adjust_range (info, info%plot_range, plot_bin)
call adjust_range (info, diff_range, diff_bin)

! Print the general satellite information
write (*,600) timestamp(), trim(S%command)
write (*,610) trim(prefix), 'userroot', trim(S%userroot)
write (*,610) trim(prefix), 'dataroot', trim(S%dataroot)
write (*,610) trim(prefix), 'sat', trim(S%sat)
write (*,610) trim(prefix), 'branch', trim(S%branch(1))
write (*,610) trim(prefix), 'phase', trim(S%phase%name)
write (*,610) trim(prefix), 'satellite', trim(S%satellite)
write (*,623) trim(prefix), 'cycles', S%cycles
write (*,623) trim(prefix), 'passes', S%passes

! Print the variable information
write (*,610) trim(prefix), 'alias', trim(var%name)
write (*,610) trim(prefix), 'var', trim(info%name)
write (*,610) trim(prefix), 'long_name', trim(var%long_name)
write (*,610) trim(prefix), 'units', trim(info%units)
write (*,610) trim(prefix), 'data', trim(info%dataname)
if (info%source /= '') write (*,610) trim(prefix), 'source', trim(info%source)
if (info%comment /= '') write (*,610) trim(prefix), 'comment', trim(info%comment)
if (info%flag_meanings /= '') write (*,610) trim(prefix), 'flag_meanings', trim(info%flag_meanings)
if (all(isan_(info%limits))) write (*,632) trim(prefix), 'limits', info%limits
if (any(isnan_(info%plot_range))) then
	! Nothing
else if (intervals == 0) then
	write (*,632) trim(prefix), 'plot_range', info%plot_range
	write (*,632) trim(prefix), 'diff_range', diff_range
else
	write (*,633) trim(prefix), 'plot_range', info%plot_range, plot_bin
	write (*,633) trim(prefix), 'diff_range', diff_range, diff_bin
endif

! Formats
600 format ('# RADS variable info'/'# Created: ',a,' UTC: ',a/'#')
610 format (a,a,'="',a,'"')
623 format (a,a,'=(',i0,2(1x,i0),')')
632 format (a,a,'=(',f0.4,1x,f0.4,')')
633 format (a,a,'=(',f0.4,2(1x,f0.4),')')
end subroutine list_variable

!***********************************************************************
! List all the variables
subroutine list_variables
integer :: i, j, varid
type(rads_varinfo), pointer :: info

write (*,600) timestamp(), trim(S%command)
if (show_x) write (*,610)
600 format ('# RADS variable list'/'# Created: ',a,' UTC: ',a/'#'/ &
'# Possible output records:'/ &
'# A  field  alias     var_name             long_name'/ &
'# D  field  var_name  default_value        long_name'/ &
'# G  field  var_name  grid_file_name       long_name'/ &
'# M  field  var_name  math_rpn_statment    long_name'/ &
'# N  field  var_name  netcdf_var_name      long_name')
610 format ( &
'# X  field  var_name  non_existent_netcdf  long_name')

! If no file, let it be known
if (P%ndata == 0) write (*,'(a)') '# Could not find any corresponding file, continuing without it.'
write (*,'(a)') '#'

! Now start the list
ncid = P%fileinfo(1)%ncid
do i = 1,S%nvar
	info => S%var(i)%info
	if (S%var(i)%name /= info%name) then
		call list ('A', S%var(i))
	else if (info%datasrc == rads_src_none) then
		call list ('U', S%var(i))
	else if (info%datasrc == rads_src_math) then
		call list ('M', S%var(i))
	else if (info%datasrc == rads_src_constant) then
		call list ('C', S%var(i))
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
			call list ('N', S%var(i))
		endif
	else if (info%datasrc == rads_src_nc_var) then
		if (.not.nft(nf90_inq_varid(ncid,info%dataname,varid))) then
			call list ('N',S%var(i))
		else if (info%default /= huge(0d0)) then
			call list ('D', S%var(i))
		else
			call list ('X', S%var(i))
		endif
	else
		call list ('G', S%var(i))
	endif
	if (associated(S%var(i)%inf1)) call list ('1', S%var(i))
	if (associated(S%var(i)%inf2)) call list ('2', S%var(i))
enddo
end subroutine list_variables

subroutine adjust_range (info, rng, bin)
type(rads_varinfo), pointer, intent(in) :: info
real(eightbytereal), intent(inout) :: rng(2)
real(eightbytereal), intent(out) :: bin
! Compute bin size; will be NaN when not requested
bin = (rng(2)-rng(1)) / intervals
if (round) call round_up (bin)	! Round if requested
! Adjust bin size if it is below the data resolution
if (info%nctype /= nf90_real .and. info%nctype /= nf90_double .and. bin <= info%scale_factor) then
	bin = info%scale_factor
endif
	rng(1) = (nint(rng(1)/bin)-0.5d0)*bin
	rng(2) = (nint(rng(2)/bin)+0.5d0)*bin
end subroutine adjust_range

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
