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

module rads
use typesizes
use rads_grid

implicit none

include 'config.inc'

! Dimensions
integer(fourbyteint), parameter :: rads_var_chunk = 100, rads_varl = 40
! RADS4 data types
integer(fourbyteint), parameter :: rads_type_other = 0, rads_type_sla = 1, rads_type_flagword = 2, &
	rads_type_time = 11, rads_type_lat = 12, rads_type_lon = 13
! RADS4 data sources
integer(fourbyteint), parameter :: rads_src_none = 0, rads_src_nc_var = 10, rads_src_nc_att = 11, rads_src_math = 20, &
	rads_src_bilinear = 30, rads_src_cspline = 31, rads_src_query = 32
! RADS4 warnings
integer(fourbyteint), parameter :: rads_warn_xml_file = -1, rads_warn_alias_circular = -2, rads_warn_nc_file = -3
! RADS4 errors
integer(fourbyteint), parameter :: rads_noerr = 0, &
	rads_err_nc_file = 1, rads_err_nc_parse = 2, rads_err_nc_close = 3, rads_err_memory = 4, &
	rads_err_var = 5, rads_err_source = 6, rads_err_nc_var = 7, rads_err_nc_get = 8, &
	rads_err_xml_parse = 9, rads_err_xml_file = 10, rads_err_alias = 11, rads_err_math = 12, &
	rads_err_cycle = 13, rads_err_nc_create = 14, rads_err_nc_put = 15
! RADS3 errors or incompatibilities
integer(fourbyteint) :: rads_err_incompat = 101, rads_err_noinit = 102
integer(twobyteint), parameter :: rads_nofield = -1
real(eightbytereal), parameter :: pi = 4d0*atan(1d0), rad = pi/180d0
character(len=1), parameter, private :: char_linefeed = char(10), char_noedit = '%'
integer, parameter, private :: stderr = 0, stdout = 6

private :: traxxing, rads_get_phase

type :: rads_varinfo
	character(len=rads_varl) :: name                 ! Short name used by RADS
	character(len=80) :: long_name                   ! Long name (description) of variable
	character(len=80) :: standard_name               ! Optional pre-described CF-compliant 'standard' name ('' if not available)
	character(len=80) :: source                      ! Optional data source ('' if none)
	character(len=rads_varl) :: units                ! Optional units of variable ('' if none)
	character(len=320) :: flag_meanings              ! Optional meaning of flag values ('' if none)
	character(len=320) :: comment                    ! Optional comment ('' if none)
	character(len=rads_varl) :: format               ! Fortran format for output
	character(len=rads_varl) :: netcdf               ! NetCDF variable name
	character(len=rads_varl) :: backup               ! Optional backup RADS variable name (or value) ('' if none)
	character(len=320) :: quality_flag               ! Quality flag associated with this variable ('' if none)
	character(len=640) :: math                       ! Math string (alternative to netcdf)
	character(len=rads_varl) :: gridx, gridy         ! RADS variable names of the grid x and y coordinates
	type(grid), pointer :: grid                      ! Pointer to grid for interpolation (alternative to netcdf)
	real(eightbytereal) :: limits(2)                 ! Lower and upper limit for editing
	real(eightbytereal) :: add_offset, scale_factor  ! Offset and scale factor in case of netCDF
	integer(fourbyteint) :: ndims                    ! Number of dimensions of variable
	integer(fourbyteint) :: nctype                   ! netCDF data type (nf90_int, etc.)
	integer(fourbyteint) :: varid                    ! netCDF variable ID
	integer(fourbyteint) :: datatype                 ! Type of data (rads_type_other|flagword|time|lat|lon)
	integer(fourbyteint) :: datasrc                  ! Retrieval source (rads_src_nc_var|nc_att|math|bilinear|cspline)
	integer(fourbyteint) :: cycle, pass              ! Last processed cycle and pass
	integer(fourbyteint) :: selected, rejected       ! Number of selected or rejected measurements
	real(eightbytereal) :: xmin, xmax, mean, sum2    ! Minimum, maximum, mean, sum squared deviation
endtype

type :: rads_var
	character(len=rads_varl) :: name                 ! Variable name (or alias thereof)
	type(rads_varinfo), pointer :: info              ! Link to struct of type(rads_varinfo)
	logical(twobyteint) :: noedit                    ! .true. if editing is suspended
	integer(twobyteint) :: field                     ! RADS3 field number (rads_nofield = none)
endtype

type :: rads_phase
	character(len=rads_varl) :: name, mission        ! Name (1-letter), and mission description
	character(len=160) :: dataroot                   ! Root directory of satellite and phase
	integer(fourbyteint) :: cycles(2),passes(2)      ! Cycle and pass limits
	real(eightbytereal) :: start_time, end_time      ! Start time and end time of this phase
	real(eightbytereal) :: ref_time, ref_lon         ! Time and longitude of equator crossing of "reference pass"
	integer(fourbyteint) :: ref_cycle, ref_pass      ! Cycle and pass number of "reference pass"
	real(eightbytereal) :: repeat_days               ! Length of repeat period in days
	real(eightbytereal) :: repeat_shift              ! Eastward shift of track pattern for near repeats
	integer(fourbyteint) :: repeat_passes            ! Number of passes per repeat period
endtype

type :: rads_sat
	character(len=160) :: dataroot                   ! Root directory of RADS data directory
	character(len=160) :: command                    ! Command line
	character(len=8) :: satellite                    ! Satellite name
	real(eightbytereal) :: nan                       ! NaN value
	real(eightbytereal) :: dt1hz                     ! "1 Hz" sampling interval
	real(eightbytereal) :: frequency(2)              ! Frequency (GHz) of primary and secondary channel
	real(eightbytereal) :: eqlonlim(0:1,2)           ! Equator longitude limits for ascending and descending passes
	real(eightbytereal) :: inclination               ! Satellite inclination and orbital period
	integer(fourbyteint) :: cycles(3),passes(3)      ! Cycle and pass limits and steps
	integer(fourbyteint) :: error                    ! Error code (positive is fatal, negative is warning)
	integer(fourbyteint) :: nvar, nsel               ! Number of available and selected variables and aliases
	integer(fourbyteint) :: debug                    ! Silent (0), verbose (1), or debug level
	integer(fourbyteint) :: pass_stat(7)             ! Statistics of rejection at start of rads_open_pass
	integer(fourbyteint) :: total_read, total_inside ! Total number of measurements read and inside region
	character(len=2) :: sat                          ! 2-Letter satellite abbreviation
	integer(twobyteint) :: satid                     ! Numerical satellite identifier
	type(rads_var), pointer :: var(:)                ! List of all (possibly) available variables and aliases
	type(rads_var), pointer :: sel(:)                ! List of selected variables and aliases
	type(rads_var), pointer :: time, lat, lon        ! Pointers to time, lat, lon variables
	type(rads_phase), pointer :: phases(:), phase    ! Definitions of all mission phases and pointer to current phase
endtype

type :: rads_pass
	character(len=160) :: filename                   ! Name of the netCDF pass file
	character(len=160) :: original                   ! Name of the original (GDR) pass file
	character(len=2048), allocatable :: history(:)   ! File creation history
	real(eightbytereal) :: equator_time, equator_lon ! Equator time and longitude
	real(eightbytereal) :: start_time, end_time      ! Start and end time of pass
	real(eightbytereal), pointer :: tll(:,:)         ! Time, lat, lon matrix
	integer(fourbyteint) :: cycle, pass              ! Cycle and pass number
	integer(fourbyteint) :: ncid, ndims              ! NetCDF ID of pass file and number of dimensions
	integer(fourbyteint) :: ndata                    ! Number of data points
	integer(fourbyteint) :: first_meas, last_meas    ! Measurement index of first and last point in region
	integer(fourbyteint) :: trkid                    ! Numerical track identifiers
	character(len=2) :: sat                          ! 2-Letter satellite abbreviation
	integer(twobyteint) :: satid                     ! Numerical satellite identifiers
	type (rads_pass), pointer :: next                ! Pointer to next pass in linked list
endtype

!***********************************************************************
!*rads_init -- Initialize RADS4
!+
! subroutine rads_init (S, sat)
! type(rads_sat), intent(inout) :: S or S(:)
! character(len=*), intent(in), optional :: sat or sat(:)
!
! This routine initializes the <S> struct with the information pertaining
! to given satellite/mission phase <sat>, which is to be formed as 'e1',
! or 'e1g', or 'e1/g'. If no phase is specified, all mission phases will be
! queried.
!
! The routine will read the satellite/mission specific setup xml-files and store
! all the information in the stuct <S>. The xml-files polled are:
!   $RADSDATAROOT/xml/rads.xml
!   ~/.rads/rads.xml
!   rads.xml
!
! The <S> and <sat> arguments can either a single element or an array. In the
! latter case, one <S> struct will be initialized for each <sat>.
! To parse command line options after this, use rads_parse_cmd.
!
! Only if the <sat> argument is omitted, then the routine will parse
! the command line for arguments in the form:
! sat=<sat> cycle=<lo>,<hi>,<step> pass=<lo>,<hi>,<step> lim:<var>=<lo>,<hi>
! lat=<lo>,<hi> lon=<lo>,<hi> alias:<var>=<var> opt:<value>=<value>
! opt=<value>,... fmt:<var>=<value>
!
! If more than one sat= argument is given, then all further options following
! this argument until the next sat= argument, plus all options prior to the
! first sat= argument will pertain to this mission.
!
! Execution will be halted when the dimension of <S> is insufficient to
! store information of multiple missions, or when required xml-files are
! missing.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  sat      : Optional: Satellite/mission abbraviation
!-----------------------------------------------------------------------
private :: rads_init_sat_0d, rads_init_sat_1d, rads_init_cmd_0d, rads_init_cmd_1d, rads_load_options, rads_parse_options
interface rads_init
	module procedure rads_init_sat_0d
	module procedure rads_init_sat_1d
	module procedure rads_init_cmd_0d
	module procedure rads_init_cmd_1d
end interface rads_init

!***********************************************************************
!*rads_get_var -- Read variable (data) from RADS4 file
!+
! recursive subroutine rads_get_var (S, P, var, data, noedit)
! type(rads_sat), intent(inout) :: S
! type(rads_pass), intent(inout) :: P
! character(len=*) <or> integer(fourbyteint) <or> type(rads_var), intent(in) :: var
! real(eightbytereal), intent(out) :: data(:)
! logical, intent(in), optional :: noedit
!
! This routine loads the data from a single variable <var> into the
! buffer <data>. This command must be preceeded by <rads_open_pass>.
! The variable <var> can be addressed as a variable name, a RADS3-type
! field number or a varlist item.
!
! The array <data> must be at the correct size to contain the entire
! pass of data, i.e., it must have the dimension P%ndata.
! If no data are available (and no backup variable) then NaN is
! returned in the array <data>.
!
! Arguments:
!  S      : Satellite/mission dependent structure
!  P      : Pass dependent structure
!  var    : (string) Name of the variable to be read.
!                    If <var> ends with % editing is skipped.
!           (integer) Field number.
!           (type(rads_var)) Variable struct (e.g. S%sel(i))
!  data   : Data returned by this routine
!  noedit : Optional: set to .true. to skip editing;
!           set to .false. to allow editing (default)
!
! Error code:
!  S%error  : rads_noerr, rads_err_var, rads_err_memory, rads_err_source
!-----------------------------------------------------------------------
private rads_get_var_by_name, rads_get_var_by_var, rads_get_var_by_number, rads_get_var_common
interface rads_get_var
	module procedure rads_get_var_by_name
	module procedure rads_get_var_by_var
	module procedure rads_get_var_by_number
end interface rads_get_var

!***********************************************************************
!*rads_def_var -- Define variable(s) to be written to RADS data file
!+
! subroutine rads_def_var (S, P, var, nctype, scale_factor, add_offset)
!
! type(rads_sat), intent(inout) :: S
! type(rads_pass), intent(inout) :: P
! type(rads_var), intent(in) :: var <or> var(:)
! integer(fourbyteint), intent(in), optional :: nctype
! real(eightbytereal), intent(in), optional :: scale_factor, add_offset
!
! This routine defines the variables to be written to a RADS data file by
! the reference to the variable struct of type(rads_var), or a
! number of variables referenced by an array of structures.
! This will create the netCDF variable and its attributes in the file
! previously opened with rads_create_pass or rads_open_pass.
!
! The optional arguments <nctype>, <scale_factor>, <add_offset>, <ndims>  can
! be used to overrule those value in the <var%info> struct.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  P        : Pass structure
!  var      : Structure(s) of variable(s) of type(rads_var)
!  nctype   : (Optional) Data type in netCDF file
!  scale_factor : (Optional) Value of the scale_factor attribute
!  add_offset : (Optional) Value of the add_offset attribute
!  ndims    : (Optional) Number of dimensions of the variable
!
! Error codes:
!  S%error  : rads_noerr, rads_err_nc_var
!-----------------------------------------------------------------------
private :: rads_def_var_0d, rads_def_var_1d
interface rads_def_var
	module procedure rads_def_var_0d
	module procedure rads_def_var_1d
end interface rads_def_var

!***********************************************************************
!*rads_stat -- Print the RADS statistics for a given satellite
!+
! subroutine rads_stat (S, unit)
! type(rads_sat), intent(in) :: S <or> S(:)
! integer(fourbyteint), intent(in), optional :: unit
!
! This routine prints out the statistics of all variables that were
! processed per mission (indicated by scalar or array <S>), to the output
! on unit <unit>. If argument <unit> is not given, then standard output
! is used.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  unit     : Fortran output unit (6 = stdout (default), 0 = stderr)
!-----------------------------------------------------------------------
private :: rads_stat_0d, rads_stat_1d
interface rads_stat
	module procedure rads_stat_0d
	module procedure rads_stat_1d
end interface rads_stat

contains

!***********************************************************************
!*rads_init_sat_0d -- Initialize RADS4 by satellite
!+
subroutine rads_init_sat_0d (S, sat, debug)
use rads_misc
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: sat
integer(fourbyteint), intent(in), optional :: debug
!
! This routine initializes the <S> struct with the information pertaining to satellite
! and mission phase <sat>, which is to be formed as 'e1', or 'e1g', or 'e1/g'.
! If no phase is specified, a default is assumed (usually 'a').
! The routine will read the satellite and mission specific setup xml-files and store
! all the information in the stuct <S>.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  sat      : Satellite/mission abbreviation
!  debug    : Optional: debug level (default = 0)
!
! Error code:
!  S%error  : rads_noerr, rads_err_xml_parse, rads_err_var, rads_warn_xml_file
!-----------------------------------------------------------------------
integer(fourbyteint) :: i
character(len=160) :: userroot

! Decipher the satellite and phase
if (len_trim(sat) < 2) call rads_exit ('Satellite/phase has fewer than 2 characters')
S%sat = sat(:2)
S%satellite = S%sat
S%satid = 0

! Set some global variables
S%dataroot = radsdataroot
call checkenv ('RADSDATAROOT',S%dataroot)
call checkenv ('HOME',userroot)
call get_command (S%command, status=i)
if (i < 0) S%command (len(S%command)-2:) = '...'

! Set all values in <S> struct to default
S%debug = 0
if (present(debug)) S%debug = debug
S%nan = make_nan()
S%dt1hz = 1d0
S%inclination = 90d0
S%frequency = (/13.8d0, S%nan/)
S%error = rads_noerr
nullify (S%sel, S%phases)
S%nvar = 0
S%nsel = 0
allocate (S%var(rads_var_chunk))
S%var = rads_var ('', null(), .false., rads_nofield)
S%total_read = 0
S%total_inside = 0
S%pass_stat = 0
S%cycles(3) = 1
S%passes(3) = 1
S%error = rads_noerr

! Read the global rads.xml setup file, the one in ~/.rads and the one in the current directory
call rads_read_xml (S, trim(S%dataroot) // '/xml/rads.xml', .true.)
call rads_read_xml (S, trim(userroot) // '/.rads/rads.xml', .false.)
call rads_read_xml (S, 'rads.xml', .false.)

! When a phase/mission is specifically given, load the appropriate settings
i = 3
if (sat(i:i) == '/' .or. sat(i:i) == ':') i = i + 1
if (sat(i:) /= '') then
	S%phase => rads_get_phase(S, sat(i:i))
	if (.not.associated(S%phase)) call rads_exit ('No such mission phase '//sat(i:i)//' of satellite '//sat(:2))
	S%phase%dataroot = trim(S%dataroot) // '/' // S%sat // '/' // trim(sat(i:))
	S%cycles(1:2) = S%phase%cycles
	S%passes(1:2) = S%phase%passes
else
	! By default, use the largest possible cycle and pass range and set the first (default) phase
	S%phase => S%phases(1)
	S%cycles(1) = minval(S%phases%cycles(1))
	S%cycles(2) = maxval(S%phases%cycles(2))
	S%passes(1) = minval(S%phases%passes(1))
	S%passes(2) = maxval(S%phases%passes(2))
endif

! Link the time, lat, lon variables for later use
S%time => rads_varptr (S, 'time')
S%lat => rads_varptr (S, 'lat')
S%lon => rads_varptr (S, 'lon')

! Now limit cycles based on time limits (if any)
if (.not.any(isnan(S%time%info%limits))) then
	S%cycles(1) = max(S%cycles(1), rads_time_to_cycle (S, S%time%info%limits(1)))
	S%cycles(2) = min(S%cycles(2), rads_time_to_cycle (S, S%time%info%limits(2)))
endif

! Compute equator longitude limits
call traxxing (S)

! List the variables
if (S%debug >= 3) then
	do i=1,S%nvar
		write (*,*) i,S%var(i)%name,S%var(i)%info%name,S%var(i)%field
	enddo
endif

end subroutine rads_init_sat_0d

subroutine rads_init_sat_1d (S, sat, debug)
type(rads_sat), intent(inout) :: S(:)
character(len=*), intent(in) :: sat(:)
integer(fourbyteint), intent(in), optional :: debug
integer :: i
if (size(S) /= size(sat)) call rads_exit ('Size of "S" and "sat" in rads_init should be the same')
if (size(S) < 1) call rads_exit ('Size of "S" in rads_init should be at least 1')
do i = 1,size(S)
	call rads_init_sat_0d (S(i), sat(i), debug)
enddo
end subroutine rads_init_sat_1d

subroutine rads_init_cmd_0d (S)
type(rads_sat), intent(inout) :: S
integer :: isatopt(1), n, debug
character(len=640), pointer :: options(:), selopt
nullify (selopt)
call rads_load_options (options, isatopt, n, debug)
if (n < 1) call rads_exit ('Failed to find "sat=" on command line')
call rads_init_sat_0d (S, options(isatopt(1))(5:), debug)
call rads_parse_options (S, options(:isatopt(1)-1), selopt)
call rads_parse_options (S, options(isatopt(1)+1:), selopt)
if (associated(selopt)) call rads_parse_varlist (S, selopt(5:))
deallocate (options)
end subroutine rads_init_cmd_0d

subroutine rads_init_cmd_1d (S)
type(rads_sat), intent(inout) :: S(:)
integer :: isatopt(size(S)), i, n, debug
character(len=640), pointer :: options(:), selopt
S%sat = ''
S%error = rads_noerr
call rads_load_options (options, isatopt, n, debug)
if (n < 1) call rads_exit ('Failed to find "sat=" on command line')
do i = 1,n
	nullify (selopt)
	call rads_init_sat_0d (S(i), options(isatopt(i))(5:), debug)
	call rads_parse_options (S(i), options(:isatopt(1)-1), selopt)
	if (i < n) then
		call rads_parse_options (S(i), options(isatopt(i)+1:isatopt(i+1)-1), selopt)
	else
		call rads_parse_options (S(i), options(isatopt(i)+1:), selopt)
	endif
	if (associated(selopt)) call rads_parse_varlist (S(i), selopt(5:))
enddo
deallocate (options)
end subroutine rads_init_cmd_1d

!***********************************************************************
!*rads_load_options -- Extract options from command line
!+
subroutine rads_load_options (options, isatopt, n, debug)
use rads_misc
character(len=640), pointer, intent(out) :: options(:)
integer, intent(out) :: isatopt(:), n, debug
!
! Scan the command line for common options (sat=, debug=, -v) and
! store these, and all others in the array <options> after allocating it
! to the right size.
!
! The <options> array will thus have the same size as the number of options
! on the command line, unless the first option is 'args=<filename>'.
! In the latter case options are first read from the file and then
! from the command line.
!
! The array <isatopt> will contain the indices of all the options that
! start with 'sat='. Finally, <debug> gives the result of the 'debug='
! or '-v' options.
!
! Arguments:
!  options : List of command line options
!  isatopt : List of indices of options starting with 'sat='
!  n       : Number of 'sat=' options
!  debug   : Verbose level
!-----------------------------------------------------------------------
integer :: m, iarg, nopts, ios, iunit
character(len=80) :: arg
m = size(isatopt)
debug = 0

! If first option is 'args=' then load the options from file
nopts = 0
call getarg(1,arg)
if (arg(:5) == 'args=') then
	iunit = getlun()
	open (iunit, file=arg(6:), status='old', iostat=ios)
	nopts = -1
	do while (ios /= 0)
		nopts = nopts + 1
		read (iunit, *, iostat=ios)
	enddo
	allocate (options(nopts+iargc()-1))
	rewind (iunit)
	do iarg = 1,nopts
		read (iunit, '(a)') options(iarg)
	enddo
	close (iunit)
	do iarg = 2,iargc()
		nopts = nopts + 1
		call getarg(iarg, options(nopts))
	enddo
else ! Load the options from the command line
	nopts = iargc()
	allocate (options(nopts))
	do iarg = 1,nopts
		call getarg(iarg, options(iarg))
	enddo
endif

! Now hunt for 'sat=' and 'debug=' in command line options
n = 0
do iarg = 1,nopts
	if (options(iarg)(:4) == 'sat=') then
		n = n + 1
		if (n > m) call rads_exit ('Too many "sat=" options on command line')
		isatopt(n) = iarg
	else if (options(iarg)(:6) == 'debug=') then
		read (options(iarg)(7:),*,iostat=ios) debug
	else if (options(iarg)(:2) == '-v') then
		debug = debug + 1
	endif
enddo
end subroutine rads_load_options

!***********************************************************************
!*rads_parse_options -- Parse RADS options
!+
subroutine rads_parse_options (S, options, selopt)
use rads_time
use rads_misc
type(rads_sat), intent(inout) :: S
character(len=*) :: options(:)
character(len=*), pointer :: selopt
!
! This routine fills the <S> struct with the information pertaining
! to a given satellite and mission as read from the options in the
! array <options>.
! The sat= option is ignored (it should be dealt with separately).
! The sel= option is not parsed, but a pointer is returned to the element
! of <options> starting with sel=
!
! Arguments:
!  S       : Satellite/mission dependent structure
!  options : Array of options (usually command line arguments)
!  selopt  : Argument containing sel= (if found)
!-----------------------------------------------------------------------
integer :: i
do i = 1,size(options)
	call rads_parse_option (options(i))
enddo
if (S%error /= rads_noerr) call rads_exit ('Error parsing command line arguments')

contains

subroutine rads_parse_option (arg)
character(len=*), target, intent(inout) :: arg
integer :: ios, i, j, k, k0, k1
real(eightbytereal) :: val(2)
! Scan a single command line argument (option)
i = index(arg,':')
j = index(arg,'=')
val = S%nan
if (arg(:6) == 'cycle=') then
	call chartrans(arg(6:),'/-x',',,,')
	S%cycles(2) = -1
	S%cycles(3) = 1
	read (arg(7:),*,iostat=ios) S%cycles
	if (S%cycles(2) < 0) S%cycles(2) = S%cycles(1)
else if (arg(:5) == 'pass=') then
	call chartrans(arg(5:),'/-x',',,,')
	S%passes(2) = -1
	S%passes(3) = 1
	read (arg(6:),*,iostat=ios) S%passes
	if (S%passes(2) < 0) S%passes(2) = S%passes(1)
else if (arg(:4) == 'lat=' .or. arg(:4) == 'lon=' .or. arg(:5) == 'time=' .or. arg(:4) == 'sla=' &
	.or. arg(:4) == 'lim:') then
	call rads_set_limits (S, arg(i+1:j-1), string=arg(j+1:))
else if (arg(:6) == 'alias:') then
	call rads_set_alias (S, arg(i+1:j-1), arg(j+1:))
else if (arg(:4) == 'fmt:') then
	call rads_set_format (S, arg(i+1:j-1), arg(j+1:))
else if (datearg(arg, val(1), val(2))) then
	call rads_set_limits (S, 'time', val(1), val(2))
else if (arg(:6) == 'debug=') then
	read (arg(7:),*,iostat=ios) S%debug
else if (arg(:4) == 'sel=') then
	selopt => arg

! The next are only for compatibility with RADS 3
else if (arg(:2) == 'h=') then
	call rads_set_limits (S, 'sla', string=arg(j+1:))
else if (arg(:4) == 'opt:') then
	call rads_set_alias (S, arg(i+1:j-1), arg(i+1:j-1)//arg(j+1:))
else if (arg(:4) == 'opt=') then
	k0 = j + 1
	k1 = len_trim(arg)
	do k = j+1,k1
		if (k == k1 .or. arg(k+1:k+1) == ',' .or. arg(k+1:k+1) == '/') then
			call rads_set_alias (S, arg(k0:k-2), arg(k0:k))
			k0 = k + 2
		endif
	enddo
endif
end subroutine rads_parse_option

end subroutine rads_parse_options

!***********************************************************************
!*rads_parse_varlist -- Parse a string of variables
!+
subroutine rads_parse_varlist (S, string)
type (rads_sat), intent(inout) :: S
character(len=*) :: string
!
! Parse the string <string> of comma- or slash-separated names of variables.
! This will allocate the array S%sel and update counter S%nsel.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  string   : String of variables
!
! Error codes:
!  S%error  : rads_err_var
!-----------------------------------------------------------------------
integer(fourbyteint) :: i0, i, l, n, noedit
type(rads_var), pointer :: temp(:), var

! Avoid parsing empty string
l = len_trim(string)
if (l == 0) then
	call rads_error (S, rads_err_var, 'No variables selected')
	return
endif

! Count the number of commas or slashes first to determine number of variables
n = 1
do i = 1,l
	if (string(i:i) == ',' .or. string(i:i) == '/') n = n + 1
enddo

! (Re)allocate memory
if (associated(S%sel)) then
	allocate (temp(S%nsel+n))
	temp(1:S%nsel) = S%sel(1:S%nsel)
	deallocate (S%sel)
	S%sel => temp
else
	allocate (S%sel(n))
endif

! Parse the variable list
i0 = 1
do i = 1,l
	if (i == l .or. string(i+1:i+1) == ',' .or. string(i+1:i+1) == '/') then
		noedit = 0
		if (string(i:i) == char_noedit) noedit = 1
		var => rads_varptr(S, string(i0:i-noedit))
		if (associated(var)) then
			S%nsel = S%nsel + 1
			S%sel(S%nsel) = var
			if (noedit == 1) S%sel(S%nsel)%noedit = .true.
		else
			call rads_error (S, rads_err_var, 'Unknown variable "'//string(i0:i)//'" removed from list of variables')
		endif
		i0 = i + 2
	endif
enddo
end subroutine rads_parse_varlist

!***********************************************************************
!*rads_parse_cmd -- Parse command line options
!+
subroutine rads_parse_cmd (S)
type(rads_sat), intent(inout) :: S
!
! After initialising RADS for a single satellite using rads_init_sat(S,sat)
! this routine can be used to update the <S> struct with information from
! the command line options. In contrast to rads_init_sat(S) this does not
! require that sat= is one of the command line options.
!
! Argument:
!  S : Satellite/mission dependent structure
!-----------------------------------------------------------------------
integer :: isatopt(1), n, debug
character(len=640), pointer :: options(:), selopt
nullify (selopt)
call rads_load_options (options, isatopt, n, debug)
call rads_parse_options (S, options, selopt)
if (associated(selopt)) call rads_parse_varlist (S, selopt(5:))
deallocate (options)
end subroutine rads_parse_cmd

!***********************************************************************
!*rads_end -- End RADS4
!+
elemental subroutine rads_end (S)
type(rads_sat), intent(inout) :: S
!
! This routine ends RADS by freeing up all <S> space
!
! Argument:
!  S   : Satellite/mission dependent struct or array of structs
!-----------------------------------------------------------------------
integer(fourbyteint) :: i,ios
if (S%sat == '') return
do i = S%nvar,1,-1
	if (associated(S%var(i)%info)) then
		if (S%var(i)%name == S%var(i)%info%name) then
			if (associated(S%var(i)%info%grid)) call grid_free(S%var(i)%info%grid)
			deallocate (S%var(i)%info, stat=ios)
		endif
		nullify (S%var(i)%info)
	endif
enddo
deallocate (S%var, S%sel, S%phases, stat=ios)
end subroutine rads_end

!***********************************************************************
!*rads_open_pass -- Open RADS pass file
!+
subroutine rads_open_pass (S, P, cycle, pass)
use netcdf
use rads_netcdf
use rads_time
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
integer(fourbyteint), intent(in) :: cycle, pass
!
! This routine opens a netCDF file for access to the RADS machinery.
! However, prior to opening the file, three tests are performed to speed
! up data selection:
! (1) All passes outside the preset cycle and pass limits are rejected.
! (2) Based on the time of the reference pass, the length of the repeat
!     cycle and the number of passes per cycle, a rough estimate is
!     made of the temporal extent of the pass. If this is outside the
!     selected time window, then the pass is rejected.
! (3) Based on the equator longitude and the pass number of the reference
!     pass, the length of the repeat cycle and the number of passes in
!     the repeat cycle, an estimate is made of the equator longitude of
!     the current pass. If this is outside the limits set in S%eqlonlim
!     then the pass is rejected.
!
! If the pass is rejected based on the above critetia or when no netCDF
! file exists, S%error returns the warning value rads_warn_nc_file.
! If the file cannot be read properly, rads_err_nc_parse is returned.
! Also, in both cases, P%ndata will be set to zero.
!
! The file opened with this routine should be closed by using rads_close_pass.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  P        : Pass structure
!  cycle    : Cycle number
!  pass     : Pass number
!
! Error codes:
!  S%error  : rads_noerr, rads_warn_nc_file, rads_err_nc_parse
!-----------------------------------------------------------------------
character(len=40) :: date
character(len=640) :: string
integer(fourbyteint) :: i,k,ascdes
real(eightbytereal) :: d
real(eightbytereal), allocatable :: tll(:,:)

! Initialise
S%error = rads_warn_nc_file
P = rads_pass ('', '', null(), S%nan, S%nan, S%nan, S%nan, null(), cycle, pass, 0, 1, 0, 0, 0, 0, S%sat, S%satid, null())
ascdes = modulo(pass,2)	! 1 if ascending, 0 if descending

if (S%debug >= 2) write (*,*) 'Checking cycle/pass : ',cycle,pass

! Do checking on cycle limits
if (cycle < S%cycles(1) .or. cycle > S%cycles(2)) then
	S%pass_stat(1) = S%pass_stat(1) + 1
	return
endif

! If the cycle is out of range for the current phase, look for a new phase
call rads_set_phase

! Do checking on pass limits (which may include new phase limits)
if (pass < S%passes(1) .or. pass > S%passes(2) .or. pass > S%phase%passes(2)) then
	S%pass_stat(2) = S%pass_stat(2) + 1
	return
endif

! Make estimates of the equator time and longitude and time range, which helps screening
d = S%phase%repeat_days * 86400d0 / S%phase%repeat_passes
P%equator_time = S%phase%ref_time + ((cycle - S%phase%ref_cycle) * S%phase%repeat_passes + (pass - S%phase%ref_pass)) * d
P%start_time = P%equator_time - 0.5d0 * d
P%end_time = P%equator_time + 0.5d0 * d
d = -nint(S%phase%repeat_days) * 360d0 / S%phase%repeat_passes
P%equator_lon = modulo(S%phase%ref_lon + (pass - S%phase%ref_pass) * d + modulo(pass - S%phase%ref_pass,2) * 180d0, 360d0)
if (S%debug >= 4) write (*,*) 'Estimated start/end/equator time/longitude = ', &
	P%equator_time, P%start_time, P%end_time, P%equator_lon

! Do checking on the time criteria (only when such are given)
if (.not.all(isnan(S%time%info%limits))) then
	if (P%end_time + 300d0 < S%time%info%limits(1) .or. P%start_time - 300d0 > S%time%info%limits(2)) then ! Allow 5 minute slop
		if (S%debug >= 2) write (*,*) 'Bail out on estimated time:', S%time%info%limits, P%start_time, P%end_time
		S%pass_stat(3) = S%pass_stat(3) + 1
		return
	endif
endif

! Do checking on the longitude criteria (only when those are limiting)
d = P%equator_lon
if (S%eqlonlim(ascdes,2) - S%eqlonlim(ascdes,1) < 360d0) then
	if (checklon(S%eqlonlim(ascdes,:),d)) then
		if (S%debug >= 2) write (*,*) 'Bail out on estimated equator longitude:', S%eqlonlim(ascdes,:), P%equator_lon
		S%pass_stat(4+ascdes) = S%pass_stat(4+ascdes) + 1
		return
	endif
endif

! Open pass file
S%pass_stat(6+ascdes) = S%pass_stat(6+ascdes) + 1
write (P%filename, '(a,"/c",i3.3,"/",a2,"p",i4.4,"c",i3.3,".nc")') trim(S%phase%dataroot), cycle, S%sat, pass, cycle
if (S%debug >= 2) write (*,*) 'Opening : '//trim(P%filename)
if (nft(nf90_open(P%filename,nf90_nowrite,P%ncid))) return

! Read global attributes
S%error = rads_err_nc_parse
if (nft(nf90_inquire(P%ncid,ndimensions=P%ndims))) return
if (nft(nf90_inquire_dimension(P%ncid,1,len=P%ndata))) return
if (nft(nf90_get_att(P%ncid,nf90_global,'equator_longitude',P%equator_lon))) return
P%equator_lon = S%lon%info%limits(1) + modulo (P%equator_lon - S%lon%info%limits(1), 360d0)
if (nft(nf90_get_att(P%ncid,nf90_global,'equator_time',date))) return
P%equator_time = strp1985f(date)
if (nft(nf90_get_att(P%ncid,nf90_global,'first_meas_time',date))) return
P%start_time = strp1985f(date)
if (nft(nf90_get_att(P%ncid,nf90_global,'last_meas_time',date))) return
P%end_time = strp1985f(date)
if (nft(nf90_get_att(P%ncid,nf90_global,'cycle_number',i)) .or. i /= cycle) return
if (nft(nf90_get_att(P%ncid,nf90_global,'pass_number',i)) .or. i /= pass) return
if (S%debug >= 3) write (*,*) 'Start/end/equator time/longitude = ', P%equator_time, P%start_time, P%end_time, P%equator_lon

S%error = rads_noerr

! Return if the file is empty
if (P%ndata == 0) return

! Update number of "read" measurements
S%total_read = S%total_read + P%ndata

! Read history
if (nff(nf90_inquire_attribute(P%ncid,nf90_global,'history',attnum=k))) then
	allocate (P%history(1))
	call nfs(nf90_get_att(P%ncid,nf90_global,'history',P%history(1)))
else if (nff(nf90_inquire_attribute(P%ncid,nf90_global,'log01',attnum=k))) then ! Read logs
	allocate (P%history(1))
	i = 0
	do
		i = i + 1
		write (date, '("log",i2.2)') i
		if (nft(nf90_get_att(P%ncid,nf90_global,date,string))) exit
		if (i == 1) then
			P%history(1) = string
			! If log01 has original file name, save it
			k = index(string, 'RAW data from file')
			if (k > 0) then
				P%original = string(k+19:)
			else
				k = index(string, 'RAW data from')
				P%original = string(k+14:)
			endif
		else
			P%history(1) = trim(string) // char_linefeed // P%history(1)
		endif
	enddo
endif

! Read time, lon, lat
allocate (tll(P%ndata,3))
call rads_get_var_by_var (S, P, S%time, tll(:,1))
call rads_get_var_by_var (S, P, S%lat, tll(:,2))
call rads_get_var_by_var (S, P, S%lon, tll(:,3))

! Look for first non-NaN measurement
do i = 1,P%ndata
	if (.not.any(isnan(tll(i,:)))) exit
enddo
if (i > P%ndata) then ! Got no non-NaNs
	P%ndata = 0
	deallocate (tll)
	return
endif
P%first_meas = i

! Now look for last non-NaN measurement
do i = P%ndata, P%first_meas, -1
	if (.not.any(isnan(tll(i,:)))) exit
enddo
P%last_meas = i

! Allocate appropriately sized time, lat, lon arrays
P%ndata = P%last_meas - P%first_meas + 1
allocate (P%tll(P%ndata,3))
P%tll = tll(P%first_meas:P%last_meas,:)
deallocate (tll)

! Successful ending; update number of measurements in region
S%error = rads_noerr
S%total_inside = S%total_inside + P%ndata

contains

subroutine rads_set_phase
integer :: i
if (cycle >= S%phase%cycles(1) .and. cycle <= S%phase%cycles(2)) return
do i = 1,size(S%phases)
	if (cycle >= S%phases(i)%cycles(1) .and. cycle <= S%phases(i)%cycles(2)) then
		S%phase => S%phases(i)
		return
	endif
enddo
S%phase => S%phases(1)
end subroutine rads_set_phase

end subroutine rads_open_pass

!***********************************************************************
!*rads_close_pass -- Open RADS pass file
!+
subroutine rads_close_pass (S, P, keep)
use netcdf
use rads_netcdf
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
logical, intent(in), optional :: keep
!
! This routine closes a netCDF file previously opened by rads_open_pass.
! The routine will reset the ncid element of the <P> structure to
! indicate that the passfile is closed.
! If <keep> is set to .true., then the time, lat and lon elements of
! the <P> structure are kept. Otherwise, they are deallocated along with
! the log entries.
! A second call to rads_close_pass without the keep argment can subsequently
! deallocate the time, lat and lon elements of the <P> structure.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  P        : Pass structure
!  keep     : Keep the P%tll matrix (destroy by default)
!
! Error code:
!  S%error  : rads_noerr, rads_err_nc_close
!-----------------------------------------------------------------------
integer :: ios
S%error = rads_noerr
if (P%ncid > 0 .and. nft(nf90_close(P%ncid))) S%error = rads_err_nc_close
P%ncid = 0
deallocate (P%history, stat=ios)
if (.not.present(keep) .or. .not.keep) deallocate (P%tll, stat=ios)
end subroutine rads_close_pass

!***********************************************************************
!*rads_get_var_by_number -- Read variable (data) from RADS4 file by integer
!+
recursive subroutine rads_get_var_by_number (S, P, field, data, noedit)
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
integer(fourbyteint), intent(in) :: field
real(eightbytereal), intent(out) :: data(:)
logical, intent(in), optional :: noedit
!
! This routine is added for backward compatibility with the field numbers
! in RADS3. It functions the same as rads_get_var, except that it
! uses the old field number as indicator for the variable.
!-----------------------------------------------------------------------
character(len=4) :: name
integer :: i
logical :: skip_edit

S%error = rads_noerr
if (P%ndata <= 0) return ! Skip empty files

! Do we need to skip editing?
if (present(noedit)) then
	skip_edit = noedit
else
	skip_edit = .false.
endif

! Look for field number in variable list
do i = 1,S%nvar
	if (S%var(i)%field == field) then
		call rads_get_var_common (S, P, S%var(i)%name, S%var(i)%info, data(:P%ndata), skip_edit)
		return
	endif
enddo
write (name,'(i4)') field
call rads_error (S, rads_err_var, 'Variable with field number '//name//' not found')
data(:P%ndata) = S%nan
end subroutine rads_get_var_by_number

!***********************************************************************
!*rads_get_var_by_var -- Read variable (data) from RADS by type(rads_var)
!+
recursive subroutine rads_get_var_by_var (S, P, var, data, noedit)
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(inout) :: var
real(eightbytereal), intent(out) :: data(:)
logical, intent(in), optional :: noedit
!
! This routine is an alternative to rads_get_var_by_name in the sense
! that it addresses the variable by a varlist item, like S%sel(i), thus
! shortcutting the need to (re)do the search for varinfo.
!-----------------------------------------------------------------------
logical :: skip_edit

S%error = rads_noerr
if (P%ndata <= 0) return ! Skip empty files

! Do we need to skip editing?
if (var%noedit) then
	skip_edit = .true.
else if (present(noedit)) then
	skip_edit = noedit
else
	skip_edit = .false.
endif

call rads_get_var_common (S, P, var%name, var%info, data(:P%ndata), skip_edit)
end subroutine rads_get_var_by_var

!***********************************************************************
!*rads_get_var_by_name -- Read variable (data) from RADS by character
!+
recursive subroutine rads_get_var_by_name (S, P, varname, data, noedit)
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
real(eightbytereal), intent(out) :: data(:)
logical, intent(in), optional :: noedit
!
! This routine loads the data from a single variable <varname>, addressed
! by a character string into the buffer <data>. This version of rads_get_var
! has an additional optional argument <varinfo>, which, if specified,
! by-passes the need to search the list of variables.
!-----------------------------------------------------------------------
type(rads_var), pointer :: var
integer(fourbyteint) :: l
logical :: skip_edit

S%error = rads_noerr
if (P%ndata <= 0) return ! Skip empty files

! If varname ends with %, suspend editing, otherwise follow noedit, or default = .false.
l = len_trim(varname)
if (varname(l:l) == char_noedit) then
	l = l - 1
	skip_edit = .true.
else if (present(noedit)) then
	skip_edit = noedit
else
	skip_edit = .false.
endif

var => rads_varptr (S, varname(:l))
if (.not.associated(var)) then
	data(:P%ndata) = S%nan
	return
endif
call rads_get_var_common (S, P, var%name, var%info, data(:P%ndata), skip_edit)
end subroutine rads_get_var_by_name


!***********************************************************************
!*rads_get_var_common -- Read variable (data) from RADS database (common to all)
!+
recursive subroutine rads_get_var_common (S, P, varname, info, data, skip_edit)
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
type(rads_varinfo), pointer, intent(inout) :: info
real(eightbytereal), intent(out) :: data(:)
logical, intent(in) :: skip_edit
!-
integer(fourbyteint) :: ios
type(rads_var), pointer :: var

! Check size of array
if (size(data) < P%ndata) then
	call rads_error (S, rads_err_memory, 'Too little memory allocated to read data')
	return
endif

! Load data depending on type of data source
select case (info%datasrc)
case (rads_src_nc_var)
	call rads_get_var_nc
case (rads_src_nc_att)
	call rads_get_att_nc
case (rads_src_math)
	call rads_get_var_math
case (rads_src_bilinear, rads_src_cspline, rads_src_query)
	call rads_get_var_grid
case default
	data = S%nan
	return
end select

! Editing or error handling
if (S%error == rads_noerr) then ! No error, edit if requested
	if (.not.skip_edit) call rads_edit_data
else if (info%backup == '') then ! No backup alternative, print error
	call rads_error (S, rads_err_var, 'Error loading variable "'//trim(varname)//'"')
	data = S%nan
else if (info%backup == 'nan' .or. info%backup(:1) == '-' .or. info%backup(:1) == '.' .or. &
	(info%backup(:1) >= '0' .and. info%backup(:1) <= '9')) then ! Numerical backup
	read (info%backup,*,iostat=ios) data(1)
	if (ios == 0) then
		S%error = rads_noerr
		data = data(1)
	else
		S%error = rads_err_var
		data = S%nan
	endif
else ! Go look for alternative variable
	S%error = rads_noerr
	var => rads_varptr (S, info%backup)
	if (.not.associated(var)) then
		data = S%nan
	else
		call rads_get_var_by_var (S, P, var, data, skip_edit)
		! Copy a few variables upward
		info%long_name = var%info%long_name
		info%standard_name = var%info%standard_name
		info%comment = var%info%comment
		! When this returns, editing and statistics will have been done, so we should directly return
		return
	endif
endif

call rads_quality_check	! Check quality flags (if provided)

call rads_update_stat	! Update statistics on variable

contains

recursive subroutine rads_get_var_nc ! Get data variable from RADS netCDF file
use netcdf
use rads_netcdf
integer(fourbyteint) :: start(1), e, nctype, ndims
real(eightbytereal) :: fillvalue, scale_factor, add_offset

! If time, lat, lon are already read, return those arrays upon request
if (P%first_meas == 0) then
else if (info%datatype == rads_type_time) then
	data = P%tll(:,1)
	return
else if (info%datatype == rads_type_lat) then
	data = P%tll(:,2)
	return
else if (info%datatype == rads_type_lon) then
	data = P%tll(:,3)
	return
endif

! Look for the variable name in the netCDF file (or take the stored one)
if (P%cycle == info%cycle .and. P%pass == info%pass) then
	! Keep old varid
else if (nff(nf90_inq_varid(P%ncid, info%netcdf, info%varid))) then
	! Read variable attributes if not yet set
	if (info%nctype == 0) info%nctype = nctype
	if (info%long_name(:1) == ' ') e =nf90_get_att(P%ncid, info%varid, 'long_name', info%long_name)
	if (info%units(:1) == ' ') e =nf90_get_att(P%ncid, info%varid, 'units', info%units)
	if (info%standard_name(:1) == ' ') e =nf90_get_att(P%ncid, info%varid, 'standard_name', info%standard_name)
	if (info%comment(:1) == ' ') e =nf90_get_att(P%ncid, info%varid, 'comment', info%comment)
else
	! Failed to find variable
	S%error = rads_err_nc_var
	return
endif
e = nf90_inquire_variable (P%ncid, info%varid, xtype=nctype, ndims=ndims)

! Load the data
start = max(1,P%first_meas)
select case (ndims)
case (0) ! Constant
	if (nft(nf90_get_var(P%ncid, info%varid, data(1)))) then
		call rads_error (S, rads_err_nc_get, 'Error reading netCDF constant "'//trim(info%netcdf)//'"')
		return
	endif
	data = data(1)
case (1) ! Array
	if (nft(nf90_get_var(P%ncid, info%varid, data(1:P%ndata), start))) then
		call rads_error (S, rads_err_nc_get, 'Error reading netCDF array "'//trim(info%netcdf)//'"')
		return
	endif
case default
	call rads_error (S, rads_err_nc_get, 'Wrong dimensions of variable "'//trim(info%netcdf)//'" (not 0 or 1)')
	return
end select

! Set NaN values and apply optional scale_factor and add_offset
if (nff(nf90_get_att(P%ncid, info%varid, '_FillValue', fillvalue))) where (data == fillvalue) data = S%nan
if (nff(nf90_get_att(P%ncid, info%varid, 'scale_factor', scale_factor))) data = data * scale_factor
if (nff(nf90_get_att(P%ncid, info%varid, 'add_offset', add_offset))) data = data + add_offset
end subroutine rads_get_var_nc

recursive subroutine rads_get_att_nc ! Get data attribute from RADS netCDF file
use netcdf
use rads_netcdf
use rads_time
integer(fourbyteint) :: varid, i
character(len=26) :: date

! First locate the colon in the name
i = index(info%netcdf, ':')

! If name starts with colon, then we have a global attribute, else a variable attribute
if (i == 1) then
	varid = nf90_global
else if (nft(nf90_inq_varid(P%ncid, info%netcdf(:i-1), varid))) then
	S%error = rads_err_nc_var
	return
endif
if (nft(nf90_inquire_attribute (P%ncid, varid, info%netcdf(i+1:), xtype=info%nctype))) info%nctype = 0

! Read the attribute
if (info%nctype == nf90_char) then
	! This is likely a date string
	if (nft(nf90_get_att(P%ncid, varid, info%netcdf(i+1:), date))) then
		call rads_error (S, rads_err_nc_get, 'Error reading netCDF attribute "'//trim(info%netcdf)//'"')
		return
	endif
	data = strp1985f (date)
else
	! Load an integer or float value
	if (nft(nf90_get_att(P%ncid, varid, info%netcdf(i+1:), data(1)))) then
		call rads_error (S, rads_err_nc_get, 'Error reading netCDF attribute "'//trim(info%netcdf)//'"')
		return
	endif
	data = data(1)
endif
end subroutine rads_get_att_nc

recursive subroutine rads_get_var_math ! Get data variable from MATH statement
use rads_math
type(math_ll), pointer :: top
integer(fourbyteint) :: i, j, l, istat

! Start with a nullified 'top'
nullify(top)

! Process the math commands left to right
l = len_trim(info%math)
i = 0
do
	j = index(info%math(i+1:),' ')
	istat = math_eval(info%math(i+1:i+j),P%ndata,top)
	if (istat /= 0) then  ! No command or value, likely to be a variable
		call math_push (P%ndata,top)
		call rads_get_var_by_name (S, P, info%math(i+1:i+j), top%data)
	endif
	if (S%error /= rads_noerr) exit
	i = i+j
	if (i >= l) exit
enddo

! When no error, copy top of stack to output
if (S%error == rads_noerr) data = top%data

! See if we have leftovers on the stack
i = -1
do while (associated(top)) 
	i = i + 1
	call math_pop(top)
enddo
if (i > 0) then
	write (*,*) i,' remaining items on stack'
	call rads_error (S, rads_noerr, 'Cleaned up')
endif
end subroutine rads_get_var_math

subroutine rads_get_var_grid ! Get data by interpolating a grid
real (eightbytereal) :: x(P%ndata), y(P%ndata)
integer(fourbyteint) :: i

! Load grid if not yet done
if (info%grid%ntype /= 0) then	! Already loaded
else if (grid_load (info%grid%filenm, info%grid) /= 0) then
	info%datasrc = rads_src_none	! Will not attempt to read grid again
	S%error = rads_err_source
	return
endif

! Get x and y coordinates
call rads_get_var_by_name (S, P, info%gridx, x, .true.)
if (S%error <= rads_noerr) call rads_get_var_by_name (S, P, info%gridy, y, .true.)

if (S%error > rads_noerr) then
else if (info%datasrc == rads_src_bilinear) then
	forall (i = 1:P%ndata) data(i) = grid_lininter (info%grid, x(i), y(i))
else if (info%datasrc == rads_src_cspline) then
	forall (i = 1:P%ndata) data(i) = grid_splinter (info%grid, x(i), y(i))
else
	forall (i = 1:P%ndata) data(i) = grid_query (info%grid, x(i), y(i))
endif
end subroutine rads_get_var_grid

subroutine rads_edit_data ! Edit the data base on limits
integer(fourbyteint) :: mask(2)

if (all(isnan(info%limits))) then
	! Skip editing when both limits are NaN
	return
else if (info%datatype == rads_type_lon) then
	! When longitude, shift the data first above the lower (west) limit
	data = info%limits(1) + modulo (data - info%limits(1), 360d0)
	! Then only check the upper (east) limit
	where (data > info%limits(2)) data = S%nan
else if (info%datatype == rads_type_flagword) then
	! Do editing of flagword
	mask = nint(info%limits)
	! Reject when any of the set bits of the lower limit are set in the data
	if (mask(1) /= 0) where (iand(nint(data),mask(1)) /= 0) data = S%nan
	! Reject when any of the set bits of the upper limit are not set in the data
	if (mask(2) /= 0) where (iand(nint(data),mask(2)) /= mask(2)) data = S%nan
else
	! Do editing of any other variable
	where (data < info%limits(1) .or. data > info%limits(2)) data = S%nan
endif
end subroutine rads_edit_data

subroutine rads_quality_check
real(eightbytereal) :: qual(P%ndata)
integer(fourbyteint) :: i, j, l
if (info%quality_flag == '') return

! Process all quality checks left to right
l = len_trim(info%quality_flag)
i = 0
do
	j = index(info%quality_flag(i+1:),' ')
	if (info%quality_flag(i+1:i+j) /= 'and') then	! Skip 'and' string
		call rads_get_var_by_name (S, P, info%quality_flag(i+1:i+j), qual)
		if (S%error /= rads_noerr) exit
		where (isnan(qual)) data = qual
	endif
	i = i+j
	if (i >= l) exit
enddo
end subroutine rads_quality_check

subroutine rads_update_stat	! Update the statistics for given var
real(eightbytereal) :: q, r
integer(fourbyteint) :: i

! If stats are already done for this cycle and pass, skip it
if (P%cycle == info%cycle .and. P%pass == info%pass) return
info%cycle = P%cycle
info%pass = P%pass

! Count selected and rejected measurements and update pass statistics
! Use method of West (1979)
do i = 1,P%ndata
	if (isnan(data(i))) then
		info%rejected = info%rejected + 1
	else
		info%selected = info%selected + 1
		info%xmin = min(info%xmin, data(i))
		info%xmax = max(info%xmax, data(i))
		q = data(i) - info%mean
		r = q / info%selected
		info%mean = info%mean + r
		info%sum2 = info%sum2 + r * q * (info%selected - 1)
	endif
enddo
end subroutine rads_update_stat

end subroutine rads_get_var_common

!***********************************************************************
!*rads_read_xml -- Read RADS4 XML file
!+
subroutine rads_read_xml (S, filename, fatal)
use netcdf
use xmlparse
use rads_time
use rads_misc
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: filename
logical, intent(in) :: fatal
!
! This routine parses a RADS4 XML file and fills the <S> struct with
! information pertaining to the given satellite and all variable info
! encountered in that file.
!
! The execution terminates on any error, and also on any warning if
! fatal = .true.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  filename : XML file name
!  fatal    : If .true., then all warnings are fatal.
!
! Error code:
!  S%error  : rads_noerr, rads_err_xml_parse, rads_warn_xml_file
!-----------------------------------------------------------------------
type(xml_parse) :: X
integer, parameter :: max_lvl = 20
character(len=rads_varl) :: tag, attr(2,max_lvl), name, tags(max_lvl)
character(len=160) :: val(max_lvl)
integer :: nattr, nval, i, ios, skip, skip_lvl
integer(twobyteint) :: field
logical :: endtag
type(rads_varinfo), pointer :: info => null()
type(rads_var), pointer :: var => null()
type(rads_phase), pointer :: phase => null()

S%error = rads_noerr
skip_lvl = 0

! Open XML file
call xml_open (X, filename, .true.)
if (.not.X%error) then
else if (fatal) then
	S%error = rads_err_xml_file
	call rads_exit ('Required xml-file '//trim(filename)//' does not exist')
	return
else
	S%error = rads_warn_xml_file
	return
endif
if (S%debug >= 2) write (*,*) 'Parsing XML file '//trim(filename)
call xml_options (X, ignore_whitespace = .true.)

! Parse XML file, store information in S struct
do
	if (X%eof) exit
	call xml_get (X, tag, endtag, attr, nattr, val, nval)
	if (X%error) then
		S%error = rads_err_xml_parse
		exit
	endif

	! Process closing tags
	if (endtag) then
		if (tag /= tags(X%level+1)) &
			call xmlparse_error ('Closing tag "'//trim(tag)//'" follows opening tag "'//trim(tags(X%level+1))//'"')
		if (X%level < skip_lvl) skip_lvl = 0	! Stop skipping when descended back below the starting level
		cycle  ! Ignore all other end tags
	endif
	
	! Process opening tags
	tags(X%level) = tag
	if (skip_lvl > 0 .and. X%level >= skip_lvl) cycle	! Skip all data at level equal or larger than skip level

	! Check if we need to skip this level
	skip = 0
	skip_lvl = 0
	do i = 1,nattr
		if (attr(1,i) == 'sat') then
			if (skip == 0) skip = 1
			if (index(attr(2,i),S%sat) > 0) skip = -1
		endif
	enddo
	if (skip == 1) then
		skip_lvl = X%level
		cycle
	endif

	select case (tag)
	case ('satellite')
		S%satellite = val(1)(:8)
	case ('satid')
		read (val(1),*,iostat=ios) S%satid
	case ('phase')
		if (has_name()) then
			phase => rads_get_phase (S, name, .true.)
			phase%name = attr(2,1)
			phase%dataroot = trim(S%dataroot) // '/' // S%sat // '/' // trim(name)
		endif
	case ('mission')
		phase%mission = val(1)(:rads_varl)
	case ('passes')
		read (val(1),*,iostat=ios) phase%passes
	case ('cycles')
		read (val(1),*,iostat=ios) phase%cycles
	case ('ref_pass') ! First part is a date string, at least 19 characters long
		i = index(val(1)(20:),' ')+19
		phase%ref_time = strp1985f(val(1)(:i))
		phase%ref_pass = 1
		read (val(1)(i+1:),*,iostat=ios) phase%ref_lon, phase%ref_cycle, phase%ref_pass
	case ('start_time')
		phase%start_time = strp1985f(val(1))
	case ('end_time')
		phase%end_time = strp1985f(val(1))
	case ('repeat')
		read (val(1),*,iostat=ios) phase%repeat_days, phase%repeat_passes, phase%repeat_shift
	case ('dt1hz')
		read (val(1),*,iostat=ios) S%dt1hz
	case ('inclination')
		read (val(1),*,iostat=ios) S%inclination
	case ('frequency')
		read (val(1),*,iostat=ios) S%frequency
	case ('alias')
		if (has_name (field)) call rads_set_alias (S, name, val(1), field)
	case ('var')
		if (has_name (field)) then
			var => rads_varptr (S, name, null())
			info => var%info
			if (field > rads_nofield) var%field = field
		endif
	case ('grid')
		info%datasrc = rads_src_bilinear
		info%gridx = '' ; info%gridy = ''
		do i = 1,nattr
			select case (attr(1,i))
			case ('x')
				info%gridx = attr(2,i)
			case ('y')
				info%gridy = attr(2,i)
			case ('inter')
				if (attr(2,i)(:1) == 'c') info%datasrc = rads_src_cspline
				if (attr(2,i)(:1) == 'q') info%datasrc = rads_src_query
			end select
		enddo
		if (info%gridx == '' .or. info%gridy == '') then
			call xmlparse_error ('Grid tag requires both x and y attributes')
			cycle
		endif
		allocate (info%grid)
		info%grid%filenm = val(1)(:128)
		info%grid%ntype = 0	! This signals that the grid was not loaded yet
	case ('long_name')
		call assign_or_append (info%long_name)
	case ('standard_name')
		info%standard_name = val(1)(:80)
		select case (val(1))
		case ('time')
			info%datatype = rads_type_time
		case ('latitude')
			info%datatype = rads_type_lat
		case ('longitude')
			info%datatype = rads_type_lon
		case ('sea_surface_height_above_sea_level')
			info%datatype = rads_type_sla
		case default
			info%datatype = rads_type_other
		end select
	case ('source')
		info%source = val(1)(:80)
	case ('units')
		info%units = val(1)(:rads_varl)
	case ('flag_values', 'flag_meanings')
		call assign_or_append (info%flag_meanings)
	case ('flag_masks')
		call assign_or_append (info%flag_meanings)
		info%datatype = rads_type_flagword
	case ('comment')
		call assign_or_append (info%comment)
	case ('netcdf')
		info%datasrc = rads_src_nc_var
		info%netcdf = val(1)(:rads_varl)
		if (index(info%netcdf, ':') > 0) info%datasrc = rads_src_nc_att
	case ('backup')
		info%backup = val(1)(:rads_varl)
	case ('quality_flag')
		call assign_or_append (info%quality_flag)
	case ('format')
		info%format = val(1)(:rads_varl)
	case ('compress')
		i = index(val(1), ' ')
		select case (val(1)(:i-1))
		case ('int1', 'byte')
			info%nctype = nf90_int1
		case ('int2', 'short')
			info%nctype = nf90_int2
		case ('int4', 'int')
			info%nctype = nf90_int4
		case ('real', 'real4', 'float')
			info%nctype = nf90_real
		case default
			info%nctype = nf90_double
		end select
		read (val(1)(i+1:), *, iostat=ios) info%scale_factor, info%add_offset
	case ('limits') ! Do not use routine rads_set_limits here!
		read (val(1), *, iostat=ios) info%limits
	case ('math')
		call assign_or_append (info%math)
		info%datasrc = rads_src_math
	case ('if', '!--')
		! Dummy and comment tags
	case default
		call xmlparse_error ('Unknown tag "'//trim(tag)//'"')
	end select
enddo

! Close XML file
call xml_close (X)

if (X%level > 0) call xmlparse_error ('Did not close tag "'//trim(tags(X%level))//'"')

if (S%error > rads_noerr) call rads_exit ('Fatal errors occurred while parsing xml-file '//trim(filename))

contains

function has_name (field)
integer(twobyteint), optional :: field
logical :: has_name
integer :: i, ios
name = ''
if (present(field)) field = rads_nofield
do i = 1,nattr
	if (attr(1,i) == 'name') then
		name = attr(2,i)
	else if (present(field) .and. attr(1,i) == 'field') then
		read (attr(2,i), *, iostat=ios) field
	endif
enddo
has_name = (name /= '')
if (.not.has_name) call xmlparse_error ('Tag <'//trim(tag)//'> requires name attribute')
end function has_name

subroutine assign_or_append (string)
! Assign or append concatenation of val(1:nval) to string
character(*), intent(inout) :: string
integer :: i, j
j = 2
do i = 1,nattr
	if (attr(1,i) == 'action' .and. attr(2,i) == 'append') j = 1
enddo
if (j == 1 .and. string /= '') then
else if (nval == 0) then
	string = ''
else
	string = val(1)
	j = 2
endif
do i = j,nval
	string = trim(string) // ' ' // val(i)
enddo
end subroutine assign_or_append

subroutine xmlparse_error (string)
! Issue error message with filename and line number
character(*), intent(in) :: string
character(160) :: text
write (text, 1300) trim(filename), X%lineno, string
call rads_error (S, rads_err_xml_parse, text)
1300 format ('Error parsing ',a,' at line',i5,': ',a)
end subroutine xmlparse_error

end subroutine rads_read_xml

!***********************************************************************
!*rads_varptr -- Returns the pointer to a given variable
!+
function rads_varptr (S, varname, tgt) result (ptr)
use netcdf
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: varname
type(rads_var), intent(in), pointer, optional :: tgt
type(rads_var), pointer :: ptr
!
! This function returns a pointer to the structure of type(rads_var)
! pointing to an element of the list of variables S%var.
!
! If the variable with the name <varname> is not available, it will either
! return a NULL pointer or create a new variable and return its pointer,
! depending on whether the optional argument <tgt> is given.
!
! The presence of <tgt> determines how the routine deals with
! non-existent variables. 
! 1. When no <tgt> argument is given:
!    - When <varname> exists: Return pointer to its structure
!    - Otherwise: Return a NULL pointer
! 2. When <tgt> is a NULL pointer:
!    - When <varname> exists: Return pointer to its structure
!    - Otherwise: Initialize a new variable and return its pointer
! 3. When <tgt> is an associated pointer:
!    - When <varname> exists: Remove its current info structure (unless
!      it is an alias), redirect it to the info structure of <tgt>,
!      and return its pointer
!    - Otherwise: Initialize a new variable, copy the info pointer
!      from <tgt> and return its pointer
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  varname  : Name of the RADS variable
!  tgt      : Optional target pointer (see above)
!  rads_varptr : Pointer to the structure for <varname>
!
! Error code:
!  S%error  : rads_noerr, rads_err_var
!-----------------------------------------------------------------------
integer(fourbyteint) :: i, n
type(rads_var), pointer :: temp(:)
integer(twobyteint) :: field

S%error = rads_noerr

! If variable name is numerical, look for field number
if (varname(1:1) >= '0' .and. varname(1:1) <= '9') then
	field = -999
	read (varname, *, iostat=i) field
	do i = 1,S%nvar
		if (S%var(i)%field == field) exit
	enddo
else
! Look for the matching variable name
	do i = 1,S%nvar
		if (S%var(i)%name(1:1) /= varname(1:1)) then
			cycle
		else if (S%var(i)%name == varname) then	! Match found: assign info and return
			exit
		endif
	enddo
endif

! If match found, return pointer
if (i <= S%nvar) then
	ptr => S%var(i)
	if (present(tgt) .and. associated(tgt)) then
		if (ptr%name == ptr%info%name) then
			if (associated(ptr%info%grid)) call grid_free(ptr%info%grid)
			deallocate (ptr%info)
		endif
		ptr%info => tgt%info
	endif
	return
! No match found, and none should be created: return null pointer and error
else if (.not.present(tgt)) then
	nullify (ptr)
	call rads_error (S, rads_err_var, 'Variable "'//trim(varname)//'" not found')
	return
endif

! If we got here, we need to make a new variable. Do we also need to allocate more space?
n = size(S%var)
if (i > n) then
	allocate (temp(n + rads_var_chunk))
	temp(1:n) = S%var
	temp(n+1:n+rads_var_chunk) = rads_var ('', null(), .false., rads_nofield)
	deallocate (S%var)
	S%var => temp
	if (S%debug >= 3) write (*,*) 'Increased S%var:',n,n+rads_var_chunk
endif
S%nvar = i
ptr => S%var(i)
ptr%name = varname

! Finally assign the info struct to that of tgt or make a new one
if (associated(tgt)) then
	ptr%info => tgt%info
else
	allocate (ptr%info)
	ptr%info = rads_varinfo (varname, varname, '', '', '', '', '', 'f7.3', '', '', '', '', '', '', null(), &
	S%nan, 0d0, 1d0, 1, nf90_double, 0, 0, 0, 0, 0, 0, 0, S%nan, S%nan, 0d0, 0d0)
endif

end function rads_varptr

!***********************************************************************
!*rads_set_alias -- Set alias to an already defined variable
!+
subroutine rads_set_alias (S, alias, varname, field)
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: alias, varname
integer(twobyteint), intent(in), optional :: field
!
! This routine defines an alias to an existing variable.
! If alias is already defined as an alias or variable, it will be overruled.
! The alias will need to point to an already existing variable or alias.
! Circular aliases are ignored.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  alias    : New alias for an existing variable
!  varname  : Existing variable name
!  field    : Optional: new field number to associate with alias
!
! Error code:
!  S%error  : rads_noerr, rads_err_alias, rads_err_var,
!             rads_warn_alias_circular
!-----------------------------------------------------------------------
type(rads_var), pointer :: tgt, src
S%error = rads_noerr
if (alias == '') then
	call rads_error (S, rads_err_alias, 'Alias is empty')
	return
else if (varname == '') then
	call rads_error (S, rads_err_var, 'Variable name is empty')
	return
else if (alias == varname) then
	call rads_error (S, rads_warn_alias_circular, 'Circular alias "'//trim(varname)//'" ignored')
	return
endif
tgt => rads_varptr (S, varname)
if (.not.associated(tgt)) then
	call rads_error (S, rads_err_var, 'Alias target "'//trim(varname)//'" of "'//trim(alias)//'" not found')
	return
endif
src => rads_varptr (S, alias, tgt)
if (present(field) .and. field /= rads_nofield) src%field = field
end subroutine rads_set_alias

!***********************************************************************
!*rads_set_limits -- Set limits on given variable
!+
subroutine rads_set_limits (S, varname, lo, hi, string)
use rads_misc
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: varname
real(eightbytereal), intent(in), optional :: lo, hi
character(len=*), intent(inout), optional :: string
!
! This routine set the lower and upper limits for a given variable in
! RADS.
! The limits can either be set by giving the lower and upper limits
! as double floats <lo> and <hi> or as a character string <string> which
! contains the two number separated by whitespace, a comma or a slash.
! In case only one number is given, only <lo> or <hi> (following the
! separator) is set.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  varname  : Variable name
!  lo, hi   : Lower and upper limit
!  string   : String of up to two values, with separating whitespace
!             or comma or slash.
!
! Error code:
!  S%error  : rads_noerr, rads_err_var
!-----------------------------------------------------------------------
type(rads_var), pointer :: var
integer :: ios
var => rads_varptr (S, varname)
if (associated(var)) then
	if (present(lo)) var%info%limits(1) = lo
	if (present(hi)) var%info%limits(2) = hi
	if (present(string)) then
		call chartrans(string, '/', ',')
		read (string, *, iostat=ios) var%info%limits
	endif
	if (var%info%datatype == rads_type_lat .or. var%info%datatype == rads_type_lon) then
		! If latitude or longitude limits are changed, recompute equator longitude limits
		call traxxing (S)
	else if (var%info%datatype == rads_type_time) then
		! If time limits are changed, also limit the cycles
		S%cycles(1) = max(S%cycles(1), rads_time_to_cycle (S, lo))
		S%cycles(2) = min(S%cycles(2), rads_time_to_cycle (S, hi))
	endif
endif
end subroutine rads_set_limits

!***********************************************************************
!*rads_set_format -- Set output format for given variable
!+
subroutine rads_set_format (S, varname, format)
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: varname, format
!
! This routine set the FORTRAN format specifier of output of a given
! variable in RADS.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  varname  : Variable name
!  format   : FORTRAN format specifier (e.g. 'f10.3')
!
! Error code:
!  S%error  : rads_noerr, rads_err_var
!-----------------------------------------------------------------------
type(rads_var), pointer :: var
var => rads_varptr(S, varname)
if (associated(var)) var%info%format = format
end subroutine rads_set_format

!***********************************************************************
!*rads_stat_0d -- Print the RADS statistics for a given satellite
!+
subroutine rads_stat_0d (S, unit)
type(rads_sat), intent(in) :: S
integer(fourbyteint), intent(in), optional :: unit
!
! This is the scalar version of rads_stat.
!-----------------------------------------------------------------------
integer(fourbyteint) :: iunit, i
iunit = 6
if (present(unit)) iunit = unit
write (iunit, 700)
write (iunit, 710) trim(S%satellite), S%sat

write (iunit, 720) 'PASSES QUERIED', sum(S%pass_stat)
write (iunit, 722)
write (iunit, 730) 'Cycle number limits', S%pass_stat(1), sum(S%pass_stat(2:7)), S%cycles
write (iunit, 730) 'Pass number limits' , S%pass_stat(2), sum(S%pass_stat(3:7)), S%passes
write (iunit, 731) 'Time limits' , S%pass_stat(3), sum(S%pass_stat(4:7)), datestring (S%time%info%limits)
write (iunit, 732) 'Equator longitude limits (asc)' , S%pass_stat(5:7:2), S%eqlonlim(1,:)
write (iunit, 732) 'Equator longitude limits (des)' , S%pass_stat(4:6:2), S%eqlonlim(0,:)

write (iunit, 721) 'PASSES AND MEASUREMENTS READ', sum(S%pass_stat(6:7)), S%total_read
write (iunit, 723)
call rads_stat_line (S%time%info)
call rads_stat_line (S%lat%info)
call rads_stat_line (S%lon%info)

write (iunit, 720) 'MEASUREMENTS IN REQUESTED PERIOD AND REGION', S%total_inside
write (iunit, 723)
do i = 1,S%nvar
	if (.not. associated(S%var(i)%info)) then ! Skip undefined variables
	else if (S%var(i)%info%selected + S%var(i)%info%rejected == 0) then ! Skip "unused" variables
	else if (S%var(i)%info%datatype >= rads_type_time) then ! Skip time, lat, lon
	else if (S%var(i)%name /= S%var(i)%info%name) then ! Skip aliases
	else ! Print statistics line for whatever remains
		call rads_stat_line (S%var(i)%info)
	endif
enddo

700 format (134('#'))
710 format ('# Editing statistics for ',a,' (',a,')')
720 format ('#'/'# ',a,t53,i10/'#')
721 format ('#'/'# ',a,t43,2i10/'#')
722 format ('#',t43,'  REJECTED  SELECTED       LOWER       UPPER        STEP')
723 format ('#',t43,'  REJECTED  SELECTED       LOWER       UPPER         MIN         MAX        MEAN      STDDEV')
730 format ("# ",a,t43,2i10,3i12)
731 format ("# ",a,t43,2i10,2a)
732 format ("# ",a,t43,2i10,2f12.3)

contains

subroutine rads_stat_line (info)
type(rads_varinfo), intent(in) :: info
real(eightbytereal), parameter :: sec2000 = 473299200d0
write (iunit, '("# ",a," [",a,"]",t43,2i10)', advance='no') trim(info%long_name), trim(info%units), &
	info%rejected, info%selected
if (info%units(:13) /= 'seconds since') then
	write (iunit, '(6f12.3)') info%limits, info%xmin, info%xmax, info%mean, sqrt(info%sum2/(info%selected-1))
else if (info%units(15:18) == '1985') then
	write (iunit, '(5a,f12.0)') datestring (info%limits), datestring(info%xmin), datestring(info%xmax), &
		datestring(info%mean), sqrt(info%sum2/(info%selected-1))
else
	write (iunit, '(5a,f12.0)') datestring (info%limits+sec2000), datestring(info%xmin+sec2000), datestring(info%xmax+sec2000), &
		datestring(info%mean+sec2000), sqrt(info%sum2/(info%selected-1))
endif
end subroutine rads_stat_line

elemental function datestring (sec)
use rads_time
real(eightbytereal), intent(in) :: sec
character(len=12) :: datestring
integer(fourbyteint) :: yy, mm, dd, hh, mn
real(eightbytereal) :: ss
if (isnan(sec)) then
	datestring = '         NaN'
else
	call sec85ymdhms(sec, yy, mm, dd, hh, mn, ss)
	write (datestring, '(1x,3i2.2,"/",2i2.2)') mod(yy,100), mm, dd, hh, mn
endif
end function datestring

end subroutine rads_stat_0d

subroutine rads_stat_1d (S, unit)
type(rads_sat), intent(in) :: S(:)
integer(fourbyteint), intent(in), optional :: unit
integer :: i
do i = 1,size(S)
	if (S(i)%sat /= '') call rads_stat_0d (S(i), unit)
enddo
end subroutine rads_stat_1d

!***********************************************************************
!*rads_exit -- Exit RADS with error message
!+
subroutine rads_exit (string)
character(len=*), intent(in) :: string
!
! This routine terminates RADS after printing an error message.
!
! Argument:
!  string   : Error message
!-----------------------------------------------------------------------
character(len=80) :: progname
call getarg (0, progname)
write (stderr, '(a,": ",a)') trim(progname),trim(string)
stop
end subroutine rads_exit

!***********************************************************************
!*rads_error -- Print error message
!+
subroutine rads_error (S, ierr, string)
type(rads_sat), intent(inout) :: S
integer(fourbyteint), intent(in) :: ierr
character(len=*), intent(in) :: string
!
! This routine prints an error message.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  ierr     : Error code
!  string   : Error message
!
! Error code:
!  S%error  : Will be set to ierr when not rads_noerr
!-----------------------------------------------------------------------
character(len=80) :: progname
call getarg (0, progname)
write (stderr, '(a,": ",a)') trim(progname),trim(string)
if (ierr /= rads_noerr) S%error = ierr
end subroutine rads_error

!***********************************************************************
!*rads_rev -- Get the library or program revision number
!+
function rads_rev (string)
character(len=*), optional :: string
integer(fourbyteint) :: rads_rev
!
! This function returns the revision number given by string (which needs
! to be an SVN rev-tag) or the revision number of the RADS library (when
! string is not given).
!
! Argument:
!  string  : SVN revision tag
!
! Return value:
!  rads_rev: Revision number of SVN revision tag, or of library
!-----------------------------------------------------------------------
character(len=20) :: libversion = '$Revision$'
integer :: l, ios
save libversion
ios = 1
if (present(string)) then
	l = len_trim(string)-1
	if (string(2:10) == 'Revision:') read (string(11:l),*,iostat=ios) rads_rev
	if (ios == 0) return
endif
l = len_trim(libversion)-1
read (libversion(11:l),*,iostat=ios) rads_rev
end function rads_rev

!***********************************************************************
!*rads_version -- Print message about current program and version
!+
function rads_version (revision, description, unit)
character(len=*), intent(in) :: revision
character(len=*), intent(in), optional :: description
integer(fourbyteint), intent(in), optional :: unit
logical :: rads_version
!
! This routine prints out a message in one of the following forms:
! 1) When first command line argument is empty or --help:
!    "rads_program (r<number>): <description>"
! 2) When the first command line argument is --version:
!    "rads_program: revision <number>, library revision <number>"
!    The program then terminates here.
! 3) When no <description> is given:
!    "rads_program (r<number>)"
! 4) No output
!
! Arguments:
!  revision    : SVN revision tag
!  description : One-line description of program
!  unit        : Fortran output unit (6 = stdout (default), 0 = stderr)
!
! Return value:
!  rads_version: .false. if output is of type 1, otherwise .true.
!-----------------------------------------------------------------------
integer :: iunit
character(len=80) :: progname,arg
call getarg (0, progname)
call getarg (1, arg)
rads_version = .true.
if (present(unit)) then
	iunit = unit
else
	iunit = stdout
endif
if (arg == '--version') then
	write (iunit, 1320) trim(progname), rads_rev(revision), rads_rev()
	stop
else if (.not.present(description)) then
	write (iunit, 1300) trim(progname), max(rads_rev(),rads_rev(revision))
else if (arg == '--help' .or. arg == '') then
	write (iunit, 1310) trim(progname), max(rads_rev(),rads_rev(revision)), trim(description)
	rads_version = .false.
else
endif
1300 format (a,' (r',i2,')')
1310 format (a,' (r',i2,'): ',a)
1320 format (a,': revision ',i2,', library revision ',i2)
end function rads_version

!***********************************************************************
!*rads_synopsis -- Print general usage information for all RADS programs
!+
subroutine rads_synopsis (unit)
integer(fourbyteint), intent(in), optional :: unit
!
! This routine prints out the usage information that is common to all
! RADS programs.
!
! Arguments:
!  unit        : Fortran output unit (6 = stdout (default), 0 = stderr)
!-----------------------------------------------------------------------
integer :: iunit
character(len=80) :: progname
call getarg (0, progname)
if (present(unit)) then
	iunit = unit
else
	iunit = stdout
endif
write (iunit, 1300) trim(progname)
1300 format (/ &
'usage: ',a,' sat=sat[/phase] [rads_dataselectors] [rads_options] [program_options]' // &
'Required argument is:'/ &
'  sat=sat[/phase]   : specify satellite [and phase] (e.g. e1/g, tx)'// &
'Optional [rads_dataselectors] are:'/ &
'  sel=var1,...      : select variables to be read'/ &
'  cycle=c0[,c1[,dc]]: specify first and last cycle and modulo'/ &
'  pass=p0[,p1[,dp]] : specify first and last pass and modulo'/ &
'  lon=lon0,lon1     : specify longitude boundaries (deg)'/ &
'  lat=lat0,lat1     : specify latitude  boundaries (deg)'/ &
'  t=t0,t1           : specify time selection'/ &
'                      (optionally use ymd=, doy=, sec= for [YY]YYMMDD[HHMMSS], YYDDD, or SEC85)'/ &
'  sla=sla0,sla1     : specify range for SLA (m)'/ &
'  alias:var1=var2   : use variable var2 when var1 is requested'/ &
'  lim:var=min,max   : specify edit data range for variable var'/ &
'  fact:var=f        : set multiplication factor for option var to f'/ &
'  xml=xmlfile       : load xml-file in addition to defaults'// &
'Still working for backwards compatibility with RADS 3 are options:'/ &
'  opt=j             : use selection code j when j/100 requested (now alias:var1=var2)'/ &
'  opt:i=j           : set option for data item i to j (now alias:var1=var2)'/ &
'  h=h0,h1           : specify range for SLA (m) (now sla=h0,h1)'// &
'Common [rads_options] are:'/ &
'  -v                : be verbose'/ &
'  debug=level       : set debug level'/ &
'  args=filename     : get any of the above arguments from filename (one argument per line)'/ &
'  --help            : this syntax massage'/ &
'  --version         : version info')
end subroutine rads_synopsis

!***********************************************************************
!*rads_get_phase -- Get pointer to satellite phase info
!+
function rads_get_phase (S, name, allow_new) result (phase)
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: name
logical, intent(in), optional :: allow_new
type(rads_phase), pointer :: phase
!-----------------------------------------------------------------------
integer(fourbyteint) :: i, n
logical :: new
type(rads_phase), pointer :: temp(:)
nullify (phase)
new = .false.
if (present(allow_new)) new = allow_new
n = 0
if (associated(S%phases)) n = size(S%phases)

! Search for the correct phase name
do i = 1,n
	if (S%phases(i)%name(1:1) /= name(1:1)) then
		cycle
	else if (S%phases(i)%name == name) then
		phase => S%phases(i)
		return
	endif
enddo

! No matching name found. Only continue when new = .true.
if (.not.new) return

! Allocate S%phases for the first time, or reallocate more space
if (n == 0) then
	allocate (S%phases(1))
	n = 1
else
	allocate (temp(n+1))
	temp(1:n) = S%phases(1:n)
	deallocate (S%phases)
	S%phases => temp
	n = n + 1
endif

! Initialize the new phase information and direct the pointer
S%phases(n) = rads_phase (name, '', '', (/999,0/), (/9999,0/), S%nan, S%nan, S%nan, S%nan, 0, 0, S%nan, S%nan, 0)
phase => S%phases(n)

end function rads_get_phase

!***********************************************************************
!*rads_time_to_cycle -- Determine cycle number for given epoch
!+
function rads_time_to_cycle (S, time)
type(rads_sat), intent(inout) :: S
real(eightbytereal), intent(in) :: time
integer(fourbyteint) :: rads_time_to_cycle
!
! Given an epoch <time> in seconds since 1985, determine the cycle in
! which that epoch falls.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  time     : Time in seconds since 1985
!
! Return value:
!  rads_time_to_cycle : Cycle number in which <time> falls
!-----------------------------------------------------------------------
integer :: i
real(eightbytereal) :: d, t0
i = 1
S%error = rads_noerr
do i = 1,size(S%phases)-1
	if (time < S%phases(i+1)%start_time) exit
enddo
d = S%phases(i)%repeat_days * 86400d0 / S%phases(i)%repeat_passes ! Length in seconds of single pass
t0 = S%phases(i)%ref_time - (S%phases(i)%ref_pass - 0.5d0) * d ! Time of start of ref_cycle
rads_time_to_cycle = floor((time - t0) / (S%phases(i)%repeat_days * 86400d0)) + S%phases(i)%ref_cycle
end function rads_time_to_cycle

!***********************************************************************
!*rads_cycle_to_time -- Determine start or end time of cycle
!+
function rads_cycle_to_time (S, cycle)
type(rads_sat), intent(inout) :: S
integer(fourbyteint), intent(in) :: cycle
real(eightbytereal) :: rads_cycle_to_time
!
! Given a cycle number, estimate the start time of that cycle.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  cycle    : Cycle number
!
! Return value:
!  rads_cycle_to_time : Start time of <cycle> in seconds since 1985
!-----------------------------------------------------------------------
integer :: i
real(eightbytereal) :: d, t0
i = 1
do i = 1,size(S%phases)-1
	if (cycle < S%phases(i+1)%cycles(1)) exit
enddo
d = S%phases(i)%repeat_days * 86400d0 / S%phases(i)%repeat_passes ! Length in seconds of single pass
t0 = S%phases(i)%ref_time - (S%phases(i)%ref_pass - 0.5d0) * d ! Time of start of ref_cycle
rads_cycle_to_time = max(S%phases(i)%start_time, t0 + (cycle - S%phases(i)%ref_cycle) * S%phases(i)%repeat_days * 86400d0)
end function rads_cycle_to_time

!***********************************************************************
!*traxxing -- Determine which tracks cross an area
!+
subroutine traxxing (S)
type(rads_sat), intent(inout) :: S
!
! Given an area and a certain satellite in a low circular orbit with
! known inclination and orbital period, the question is: which tracks
! cross the area?
!
! The area is specified by longitude and latitude limits defined in
! <S%lon%info%limits> and <S%lat%info%limits>.
! The satellite orbit is described by the inclination inclination,
! repeat period and number of repeat passes: <S%inclination>,
! <S%phase%repeat_days>, and <S%phase%repeat_passes>.
!
! The answer to the question is returned in <S%eqlonlim>: 4 values
! representing the lower and upper bounds for the longitude of the
! equator passage for ascending and descending tracks, resp.
!
! Argument:
!  S        : Satellite/mission dependent structure
!-----------------------------------------------------------------------
real(eightbytereal) :: u(2),l(2),latmax
integer(fourbyteint) :: i

! Initialize
latmax = min(S%inclination,180d0-S%inclination)

do i = 1,2
	! Compute argument of latitude for lower and upper latitude boundary
	if (S%lat%info%limits(i) < -latmax) then
		u(i) = -90d0
		l(i) = -90d0
		if (S%inclination > 90d0) l(i) = -l(i)
	else if (S%lat%info%limits(i) > latmax) then
		u(i) = 90d0
		l(i) = 90d0
		if (S%inclination > 90d0) l(i) = -l(i)
	else
		u(i) = asin(sin(S%lat%info%limits(i)*rad)/sin(S%inclination*rad))/rad
		l(i) = asin(tan(S%lat%info%limits(i)*rad)/tan(S%inclination*rad))/rad
	endif

	! Compute the longitude advance corrected for time
	l(i) = l(i) - 2d0*u(i)*S%phase%repeat_days/S%phase%repeat_passes
enddo

! Compute the equator crossing longitudes for ascending tracks
! Add 2-degree margin
S%eqlonlim(1,1) = S%lon%info%limits(1) - max(l(1),l(2)) - 2d0
S%eqlonlim(1,2) = S%lon%info%limits(2) - min(l(1),l(2)) + 2d0

! Compute the equator crossing longitudes for descending tracks
! Add 2-degree margin
S%eqlonlim(0,1) = S%lon%info%limits(1) + min(l(1),l(2)) - 2d0
S%eqlonlim(0,2) = S%lon%info%limits(2) + max(l(1),l(2)) + 2d0

if (S%debug >= 3) write (*,*) "Eqlonlim = ",S%eqlonlim

end subroutine traxxing

!***********************************************************************
!*rads_progress_bar -- Print and update progress of scanning cycles/passes
!+
subroutine rads_progress_bar (S, P, nselpass, unit)
type(rads_sat), intent(in) :: S
type(rads_pass), intent(in) :: P
integer(fourbyteint), intent(in) :: nselpass
integer(fourbyteint), intent(in), optional :: unit
!
! This routine prints information on the progress of cycling through
! cycles and passes of satellite data. If used, it should preferably be
! called before every call of rads_close_pass.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  P        : Pass structure
!  nselpass : Number of measurements selected in this pass
!  unit     : Optional Fortran output unit (6 = stdout (default), 0 = stderr)
!-----------------------------------------------------------------------
integer :: pos_old = 50, lin_old = -1, cycle_old = -1, pos, lin, logunit, i

! Formats for printing progress report

700 format(79('*')// &
'Data selection for satellite ',a,' phase ',a// &
'x = pass has no data in period and area'/ &
'- = pass file does not exist'/ &
'o = pass file has no valid data'/ &
'# = pass file has valid data')
710 format(/'Cycle Pass  ', &
'....+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8....+....9....+....0')
720 format(/2i5,'  ')
750 format(a)

logunit = 6
if (present(unit)) logunit = unit
if (cycle_old < 0) write (logunit,700) trim(S%satellite),trim(S%phase%name)
pos = mod(P%pass-1,100)+1
lin = int((P%pass-1)/100)
if (P%cycle /= cycle_old) then
	write (logunit,710,advance='no')     ! Print header
	pos_old = 100
	cycle_old = P%cycle
endif
if (pos_old == 100 .or. lin_old /= lin) then
	write (logunit,720,advance='no') P%cycle,int((P%pass-1)/100)*100+1
	! Start new line, print cycle nr and pass nr
	pos_old=0
endif
do i=pos_old+1,pos-1
	write (logunit,750,advance='no') ' ' ! Fill up to current pass with ' '
enddo
if (S%error /= rads_noerr) then
	write (logunit,750,advance='no') '-' ! Print - for non existing pass
else if (P%ndata == 0) then
	write (logunit,750,advance='no') 'x' ! Print x for empty pass
else if (nselpass == 0) then
	write (logunit,750,advance='no') 'o' ! Print o when no data is selected
else
	write (logunit,750,advance='no') '#' ! Print # when data available
endif
pos_old=pos
lin_old=lin

end subroutine rads_progress_bar

!***********************************************************************
!*rads_create_pass -- Create RADS pass (data) file
!+
subroutine rads_create_pass (S, P, ndata, name)
use netcdf
use rads_netcdf
use rads_time
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
integer(fourbyteint), intent(in), optional :: ndata
character(len=*), intent(in), optional :: name
!
! This routine creates a new RADS netCDF data file. If one of the same
! file name already exists, it is removed.
! The file is initialized with the appropriate global attributes, and
! the primary dimension ('time') will be set up. This dimension can
! either be fixed (ndata > 0) or unlimited (ndata == 0).
! If an unlimited dimension is selected, then the pass-related global
! attributes, like cycle, pass, equator_time, equator_lon will not be written
! to the file.
! When the <ndata> argument is omitted, the dimension <P%ndata> is used.
!
! The optional argument <name> can have one of three forms:
! - Left empty it specifies the current directory
! - Ending in '/' it specifies a named directory in which to store the pass
!   file of the form <sat>p<pass>c<cycle>.nc
! - Otherwise it specifies the filename in full.
! If <name> is omitted, the default directory is used:
! $RADSDATAROOT/<sat>/<phase>/c<cycle>/
! Any directory that does not yet exist will be created.
!
! Upon entry, the <P> structure needs to contain the relevant information on
! the pass (cycle, pass, equator_time, equator_lon). Upon return, the
! P%ncid and P%dimid will be updated.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  P        : Pass structure
!  ndata    : Length of the primary ('time') dimension (use 0 for unlimited)
!  name     : Name of directory in which to store pass file, or file name
!
! Error codes:
!  S%error  : rads_noerr, rads_err_nc_create
!-----------------------------------------------------------------------
integer(fourbyteint) :: l, e
logical :: exist
character(len=26) :: date(3)
character(len=10) :: yymmdd
integer :: values(9)
S%error = rads_noerr

! Build the filename, make directory if needed
if (.not.present(name)) then
	write (P%filename, '(a,"/c",i3.3,"/",a2,"p",i4.4,"c",i3.3,".nc")') trim(S%phase%dataroot), S%sat, P%cycle, P%pass, P%cycle
	l = len_trim(P%filename)-15
	inquire (file = P%filename(:l), exist = exist)
	if (.not.exist) call system ('mkdir -p ' // P%filename(len_trim(P%filename)-13:))
else if (name == '') then
	write (P%filename, '(a2,"p",i4.4,"c",i3.3,".nc")') S%sat, P%pass, P%cycle
else if (name(len_trim(name):) == '/') then
	write (P%filename, '(a,a2,"p",i4.4,"c",i3.3,".nc")') trim(name), S%sat, P%pass, P%cycle
	inquire (file = name, exist = exist)
	if (.not.exist) call system ('mkdir -p ' // name)
else
	P%filename = name
endif

! Create the (new) data file
if (present(ndata)) P%ndata = ndata
if (S%debug >= 2) write (*,*) 'Creating ',trim(P%filename),P%ndata
if (nft(nf90_create(P%filename, nf90_write+nf90_nofill, P%ncid))) then
	call rads_error (S, rads_err_nc_create, 'Error creating ' // trim(P%filename))
	return
endif

! Define the principle dimension
if (nft(nf90_def_dim (P%ncid, 'time', P%ndata, values(1)))) then
	call rads_error (S, rads_err_nc_create, 'Error creating time dimension of ' // trim(P%filename))
	return
endif

! Specify the global attributes
e = nf90_put_att(P%ncid, nf90_global, 'Conventions', 'CF-1.5') + &
	nf90_put_att(P%ncid, nf90_global, 'title', 'RADS 4.0 pass file') + &
	nf90_put_att(P%ncid, nf90_global, 'institution', 'Altimetrics / NOAA / TU Delft') + &
	nf90_put_att(P%ncid, nf90_global, 'source', 'radar altimeter') + &
	nf90_put_att(P%ncid, nf90_global, 'reference_document', 'RADS Data Manual, Issue 4.0') + &
	nf90_put_att(P%ncid, nf90_global, 'mission_name', trim(S%satellite)) + &
	nf90_put_att(P%ncid, nf90_global, 'mission_phase', S%phase%name(:1))
if (ndata > 0) then
	date = strf1985f ((/P%equator_time,P%start_time,P%end_time/))
	e = e + &
	nf90_put_att(P%ncid, nf90_global, 'cycle_number', P%cycle) + &
	nf90_put_att(P%ncid, nf90_global, 'pass_number', P%pass) + &
	nf90_put_att(P%ncid, nf90_global, 'equator_longitude', P%equator_lon) + &
	nf90_put_att(P%ncid, nf90_global, 'equator_time', date(1)) + &
	nf90_put_att(P%ncid, nf90_global, 'first_meas_time', date(2)) + &
	nf90_put_att(P%ncid, nf90_global, 'last_meas_time', date(3))
endif
e = e + nf90_put_att(P%ncid, nf90_global, 'original', P%original)

call gmtime(time(),values)
write (yymmdd, '(i4.4,2("-",i2.2))') values(6)+1900,values(5)+1,values(4)
if (allocated(P%history)) then
	e = e + nf90_put_att(P%ncid, nf90_global, 'history', yymmdd//': '//trim(S%command)//char_linefeed//trim(P%history(1)))
else
	e = e + nf90_put_att(P%ncid, nf90_global, 'history', yymmdd//': '//trim(S%command))
endif
if (e > 0) call rads_error (S, rads_err_nc_create, 'Error writing global attributes to '//trim(P%filename))
end subroutine rads_create_pass

!***********************************************************************
!*rads_def_var_0d -- Define variable to be written to RADS data file
!
subroutine rads_def_var_0d (S, P, var, nctype, scale_factor, add_offset)
use netcdf
use rads_netcdf
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
integer(fourbyteint), intent(in), optional :: nctype
real(eightbytereal), intent(in), optional :: scale_factor, add_offset
!
! Zero-dimensional version of rads_def_var.
!-----------------------------------------------------------------------
type(rads_varinfo), pointer :: info
integer(fourbyteint) :: e, i, n
integer, parameter :: dimid(1:4) = (/ 1, 2, 3, 4 /)
integer(onebyteint), parameter :: flag_values(0:9) = int((/ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 /), onebyteint)
S%error = rads_noerr
info => var%info
if (present(nctype)) info%nctype = nctype
if (present(scale_factor)) info%scale_factor = scale_factor
if (present(add_offset)) info%add_offset = add_offset
	
if (nft(nf90_def_var(P%ncid, var%name, info%nctype, dimid(1:info%ndims), info%varid))) then
	call rads_error (S, rads_err_nc_var, 'Error creating variable '//trim(var%name)//' in '//trim(P%filename))
	return
endif
e = nf90_put_att (P%ncid, info%varid, 'long_name', trim(info%long_name))
if (info%standard_name /= '') e = e + nf90_put_att (P%ncid, info%varid, 'standard_name', trim(info%standard_name))
if (info%source /= '') e = e + nf90_put_att (P%ncid, info%varid, 'source', trim(info%source))
if (info%units /= '') e = e + nf90_put_att (P%ncid, info%varid, 'units', trim(info%units))
if (info%datatype >= rads_type_time) then
	! No _FillValue, flag_values or flag_masks
else if (info%flag_meanings /= '') then
	! Count the number of spaces in flag_meanings
	n = 0
	do i = 1,len_trim(info%flag_meanings)
		if (info%flag_meanings(i:i) == ' ') n = n + 1
	enddo
	if (info%datatype == rads_type_flagword) then
		if (info%nctype == nf90_int1) then
			e = e + nf90_put_att (P%ncid, info%varid, 'flag_masks', int(2**flag_values(0:n),onebyteint))
		else
			e = e + nf90_put_att (P%ncid, info%varid, 'flag_masks', int(2**flag_values(0:n),twobyteint))
		endif
	else
		if (info%nctype == nf90_int1) then
			e = e + nf90_put_att (P%ncid, info%varid, 'flag_values', int(flag_values(0:n),onebyteint))
		else
			e = e + nf90_put_att (P%ncid, info%varid, 'flag_values', int(flag_values(0:n),twobyteint))
		endif
	endif
	e = e + nf90_put_att (P%ncid, info%varid, 'flag_meanings', info%flag_meanings)
else if (info%nctype == nf90_int1) then
	e = e + nf90_put_att (P%ncid, info%varid, '_FillValue', huge(0_onebyteint))
else if (info%nctype == nf90_int2) then
	e = e + nf90_put_att (P%ncid, info%varid, '_FillValue', huge(0_twobyteint))
else if (info%nctype == nf90_int4) then
	e = e + nf90_put_att (P%ncid, info%varid, '_FillValue', huge(0_fourbyteint))
endif
if (info%scale_factor /= 1d0) e = e + nf90_put_att (P%ncid, info%varid, 'scale_factor', info%scale_factor)
if (info%add_offset /= 0d0)  e = e + nf90_put_att (P%ncid, info%varid, 'add_offset', info%add_offset)
if (info%datatype < rads_type_time) e = e + nf90_put_att (P%ncid, info%varid, 'coordinates', 'lon lat')
if (info%comment /= '') e = e + nf90_put_att (P%ncid, info%varid, 'comment', info%comment)
if (var%field /= rads_nofield) e = e + nf90_put_att (P%ncid, info%varid, 'field', var%field)
if (e > 0) call rads_error (S, rads_err_nc_var, &
	'Error writing attributes for variable '//trim(var%name)//' in '//trim(P%filename))
info%cycle = P%cycle
info%pass = P%pass
end subroutine rads_def_var_0d

subroutine rads_def_var_1d (S, P, var, nctype, scale_factor, add_offset)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var(:)
integer(fourbyteint), intent(in), optional :: nctype
real(eightbytereal), intent(in), optional :: scale_factor, add_offset
integer :: i
do i = 1,size(var)
	call rads_def_var_0d (S, P, var(i), nctype, scale_factor, add_offset)
	if (S%error /= rads_noerr) return
enddo
end subroutine rads_def_var_1d

!***********************************************************************
!*rads_put_var -- Write data for variable to RADS file
!+
subroutine rads_put_var (S, P, var, data, start)
use netcdf
use rads_netcdf
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
real(eightbytereal), intent(in) :: data(:)
integer(fourbyteint), intent(in), optional :: start
!
! This routine writes the data array <data> for the variable <var>
! (referenced by the structure of type(rads_var)) to the netCDF file
! previously opened with rads_create_pass or rads_open_pass.
!
! The data in <data> are in the original SI units (like [m] or [s])
! and will be converted to the internal units based on the values of
! <var%info%nctype>, <var%info%scale_factor>, <var%info%add_offset>.
!
! The optional argument <start> can be added to indicate an offset
! for storing the data from the first available position in the file.
! For example: start=101 first skips 100 records.
!
! Arguments:
!  S        : Satellite/mission dependent structure
!  P        : Pass structure
!  var      : Structure of variable of type(rads_var)
!  data     : Data to be written (in original (SI) units)
!  start    : (Optional) Position of the first data point in the file
!
! Error codes:
!  S%error  : rads_noerr, rads_err_nc_put
!-----------------------------------------------------------------------
type(rads_varinfo), pointer :: info
integer(fourbyteint) :: e, first(1)
S%error = rads_noerr
e = nf90_enddef(P%ncid) ! Make sure to get out of define mode
info => var%info
if (P%cycle == info%cycle .and. P%pass == info%pass) then
	! Keep old varid
else if (nft(nf90_inq_varid (P%ncid, var%name, info%varid))) then
	call rads_error (S, rads_err_nc_var, 'No variable '//trim(var%name)//' in '//trim(P%filename))
	return
endif
if (present(start)) then
	first = start
else
	first = 1
endif

select case (info%nctype)
case (nf90_int1)
	e = nf90_put_var (P%ncid, info%varid, nint1((data - info%add_offset) / info%scale_factor), first)
case (nf90_int2)
	e = nf90_put_var (P%ncid, info%varid, nint2((data - info%add_offset) / info%scale_factor), first)
case (nf90_int4)
	e = nf90_put_var (P%ncid, info%varid, nint4((data - info%add_offset) / info%scale_factor), first)
case default
	e = nf90_put_var (P%ncid, info%varid, (data - info%add_offset) / info%scale_factor, first)
end select
if (e > 0) call rads_error (S, rads_err_nc_put, 'Error writing data for variable '//trim(var%name)//' to '//trim(P%filename))
end subroutine rads_put_var

end module rads
