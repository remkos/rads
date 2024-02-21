!****h* module/rads
! SUMMARY
! RADS main module
!
! SYNOPSIS
module rads
use typesizes
use rads_grid, only: grid

! * Parameters
! Dimensions
integer(fourbyteint), parameter :: rads_var_chunk = 100, rads_varl = 40, &
	rads_naml = 160, rads_cmdl = 640, rads_strl = 1600, rads_hstl = 3200, &
	rads_cyclistl = 50, rads_optl = 50, rads_max_branches = 5
! RADS4 data types
integer(fourbyteint), parameter :: rads_type_other = 0, rads_type_sla = 1, &
	rads_type_flagmasks = 2, rads_type_flagvalues = 3, rads_type_time = 11, &
	rads_type_lat = 12, rads_type_lon = 13, rads_type_dim = 14
! RADS4 data sources
integer(fourbyteint), parameter :: rads_src_none = 0, rads_src_nc_var = 10, &
	rads_src_nc_att = 11, rads_src_math = 20, rads_src_grid_lininter = 30, &
	rads_src_grid_splinter = 31, rads_src_grid_query = 32, rads_src_grid_linphase = 33, &
	rads_src_constant = 40, rads_src_flags = 50, rads_src_tpj = 60
! RADS4 warnings
integer(fourbyteint), parameter :: rads_warn_nc_file = -3
! RADS4 errors
integer(fourbyteint), parameter :: rads_noerr = 0, &
	rads_err_nc_file = 1, rads_err_nc_parse = 2, rads_err_nc_close = 3, rads_err_memory = 4, &
	rads_err_var = 5, rads_err_source = 6, rads_err_nc_var = 7, rads_err_nc_get = 8, &
	rads_err_xml_parse = 9, rads_err_xml_file = 10, rads_err_alias = 11, rads_err_math = 12, &
	rads_err_cycle = 13, rads_err_nc_create = 14, rads_err_nc_put = 15, rads_err_nophase = 16
! Additional RADS4 helpers
character(len=1), parameter :: rads_linefeed = char(10), rads_tab = char(9), rads_noedit = '_'
! RADS3 errors or incompatibilities
integer(fourbyteint), parameter :: rads_err_incompat = 101, rads_err_noinit = 102
integer(twobyteint), parameter :: rads_nofield = -1
! Math constants
real(eightbytereal), parameter :: pi = 3.1415926535897932d0, rad = pi/180d0
! I/O parameters
integer, parameter :: stderr = 0, stdin = 5, stdout = 6

! * Variables
! I/O variables
integer(fourbyteint), save :: rads_verbose = 0       ! Verbosity level
integer(fourbyteint), save :: rads_log_unit = stdout ! Unit number for statistics logging

! * RADS4 variable structures
type :: rads_varinfo                                 ! Information on variable used by RADS
	character(len=rads_varl) :: name                 ! Short name of variable used by RADS
	character(len=rads_naml) :: long_name            ! Long name (description) of variable
	character(len=rads_naml) :: standard_name        ! Optional CF 'standard' name ('' if none)
	character(len=rads_naml) :: source               ! Optional data source ('' if none)
	character(len=rads_naml) :: parameters           ! Optional link to model parameters ('' if none)
	character(len=rads_strl) :: dataname             ! Name associated with data (e.g. NetCDF var name)
	character(len=rads_cmdl) :: flag_meanings        ! Optional meaning of flag values ('' if none)
	character(len=rads_cmdl) :: quality_flag         ! Quality flag(s) associated with var ('' if none)
	character(len=rads_cmdl) :: comment              ! Optional comment ('' if none)
	character(len=rads_varl) :: units                ! Optional units of variable ('' if none)
	character(len=rads_varl) :: format               ! Fortran format for output
	character(len=rads_varl) :: gridx, gridy         ! RADS variable names of the grid x and y coords
	type(grid), pointer :: grid                      ! Pointer to grid (if data source is grid)
	real(eightbytereal) :: default                   ! Optional default value (Inf if not set)
	real(eightbytereal) :: limits(2)                 ! Lower and upper limit for editing
	real(eightbytereal) :: plot_range(2)             ! Suggested range for plotting
	real(eightbytereal) :: add_offset, scale_factor  ! Offset and scale factor in case of NetCDF
	real(eightbytereal) :: xmin, xmax, mean, sum2    ! Minimum, maximum, mean, sum squared deviation
	logical :: boz_format                            ! Format starts with B, O or Z.
	integer(fourbyteint) :: ndims                    ! Number of dimensions of variable
	integer(fourbyteint) :: brid                     ! Branch ID (default 1)
	integer(fourbyteint) :: nctype, varid            ! NetCDF data type (nf90_int, etc.) and input var ID
	integer(fourbyteint) :: datatype                 ! Type of data (one of rads_type_*)
	integer(fourbyteint) :: datasrc                  ! Retrieval source (one of rads_src_*)
	integer(fourbyteint) :: cycle, pass              ! Last processed cycle and pass
	integer(fourbyteint) :: selected, rejected       ! Number of selected or rejected measurements
endtype

type :: rads_var                                     ! Information on variable or alias
	character(len=rads_varl), pointer :: name        ! Pointer to short name of variable (or alias)
	character(len=rads_naml), pointer :: long_name   ! Pointer to long name (description) of variable
	type(rads_varinfo), pointer :: info, inf1, inf2  ! Links to structs of type(rads_varinfo)
	logical(twobyteint) :: noedit                    ! .true. if editing is suspended
	integer(twobyteint) :: field(2)                  ! RADS3 field numbers (rads_nofield = none)
endtype

type :: rads_cyclist                                 ! List of cycles
	integer(fourbyteint) :: n, i                     ! Number of elements in list, additional value
	integer(fourbyteint) :: list(rads_cyclistl)      ! List of values
endtype

type :: rads_phase                                   ! Information about altimeter mission phase
	character(len=rads_varl) :: name, mission        ! Name (1-letter), and mission description
	integer(fourbyteint) :: cycles(2), passes        ! Cycle range and maximum number of passes
	real(eightbytereal) :: start_time, end_time      ! Start time and end time of this phase
	real(eightbytereal) :: ref_time, ref_lon         ! Time and lon of equator crossing of "ref. pass"
	integer(fourbyteint) :: ref_cycle, ref_pass      ! Cycle and pass number of "reference pass"
	integer(fourbyteint) :: ref_orbit                ! Absolute orbit nr of "reference pass" (at eq.)
	real(eightbytereal) :: pass_seconds              ! Length of pass in seconds
	real(eightbytereal) :: repeat_days               ! Length of repeat period in days
	real(eightbytereal) :: repeat_shift              ! Eastward shift of track pattern per cycle for near repeats
	integer(fourbyteint) :: repeat_nodal             ! Length of repeat period in nodal days
	integer(fourbyteint) :: repeat_passes            ! Number of passes per repeat period
	type(rads_cyclist), pointer :: subcycles         ! Subcycle definition (if requested)
endtype

type :: rads_sat                                     ! Information on altimeter mission
	character(len=rads_naml) :: userroot             ! Root directory of current user (i.e. $HOME)
	character(len=rads_naml) :: dataroot             ! Root directory of RADS data (i.e. $RADSDATAROOT)
	character(len=rads_varl) :: branch(rads_max_branches) ! Name of optional branches
	character(len=rads_varl) :: spec                 ! Temporary holding space for satellite specs
	character(len=rads_cmdl) :: command              ! Command line
	character(len=rads_naml), pointer :: glob_att(:) ! Global attributes
	character(len=8) :: satellite                    ! Satellite name
	real(eightbytereal) :: dt1hz                     ! "1 Hz" sampling interval
	real(eightbytereal) :: frequency(2)              ! Frequency (GHz) of primary and secondary channel
	real(eightbytereal) :: inclination               ! Satellite inclination (deg)
	real(eightbytereal) :: eqlonlim(0:1,2)           ! Equator lon limits for asc. and desc. passes
	real(eightbytereal) :: centroid(3)               ! Lon, lat, distance (in rad) selection criteria
	real(eightbytereal) :: xover_params(2)           ! Crossover parameters used in radsxoconv
	integer(fourbyteint) :: cycles(3),passes(3)      ! Cycle and pass limits and steps
	integer(fourbyteint) :: error                    ! Error code (positive = fatal, negative = warning)
	integer(fourbyteint) :: pass_stat(7)             ! Stats of rejection at start of rads_open_pass
	integer(fourbyteint) :: total_read, total_inside ! Total nr of measurements read and inside region
	integer(fourbyteint) :: nvar, nsel               ! Nr of available and selected vars and aliases
	logical :: n_hz_output                           ! Produce multi-Hz output
	character(len=2) :: sat                          ! 2-Letter satellite abbreviation
	integer(twobyteint) :: satid                     ! Numerical satellite identifier
	type(rads_cyclist), pointer :: excl_cycles       ! Excluded cycles (if requested)
	type(rads_var), pointer :: var(:)                ! List of available variables and aliases
	type(rads_var), pointer :: sel(:)                ! List of selected variables and aliases
	type(rads_var), pointer :: time, lat, lon        ! Pointers to time, lat, lon variables
	type(rads_phase), pointer :: phases(:)           ! Definitions of all mission phases
	type(rads_phase), pointer :: phase               ! Pointer to current phase
endtype

type :: rads_file                                    ! Information on RADS data file
	integer(fourbyteint) :: ncid                     ! NetCDF ID of pass file
	character(len=rads_cmdl) :: name                 ! Name of the NetCDF pass file
endtype

type :: rads_pass                                    ! Pass structure
	character(len=rads_strl) :: original             ! Name of the original (GDR) pass file(s)
	character(len=rads_hstl), pointer :: history     ! File creation history
	real(eightbytereal) :: equator_time, equator_lon ! Equator time and longitude
	real(eightbytereal) :: start_time, end_time      ! Start and end time of pass
	real(eightbytereal), pointer :: tll(:,:)         ! Time, lat, lon matrix
	integer(twobyteint), pointer :: flags(:)         ! Array of engineering flags
	logical :: rw                                    ! NetCDF file opened for read/write
	integer(fourbyteint) :: cycle, pass              ! Cycle and pass number
	integer(fourbyteint) :: nlogs                    ! Number of RADS3 log entries
	integer(fourbyteint) :: ndata                    ! Number of data points (1-Hz)
	integer(fourbyteint) :: n_hz, n_wvf              ! Size second/third dimension (0=none)
	integer(fourbyteint) :: first_meas, last_meas    ! Index of first and last point in region
	integer(fourbyteint) :: time_dims                ! Dimensions of time/lat/lon stored
	integer(fourbyteint) :: trkid                    ! Numerical track identifiers
	type (rads_file) :: fileinfo(rads_max_branches)  ! File information for pass files
	type (rads_sat), pointer :: S                    ! Pointer to satellite/mission structure
	type (rads_pass), pointer :: next                ! Pointer to next pass in linked list
endtype

type :: rads_option                                  ! Information on command line options
	character(len=rads_varl) :: opt                  ! Option (without the - or --)
	character(len=rads_cmdl) :: arg                  ! Option argument
	integer :: id                                    ! Identifier in form 10*nsat + i
endtype

! These command line options can be accessed by RADS programs
type(rads_option), allocatable, target, save :: &
	rads_opt(:)                                      ! List of command line options
integer(fourbyteint), save :: rads_nopt = 0          ! Number of command line options saved

!
! PURPOSE
! This module provides the main functionalities for the RADS4 software.
! To use any of the following subroutines and functions, add the following
! line in your Fortran 90 (or later) code:
!	use rads
!-----------------------------------------------------------------------
! COPYRIGHT
! Copyright (c) 2011-2024  Remko Scharroo
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

real(eightbytereal), parameter, private :: nan = transfer ((/not(0_fourbyteint),not(0_fourbyteint)/),0d0)
integer(onebyteint), parameter, private :: flag_values(0:15) = &
	int((/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/), onebyteint)
integer(twobyteint), parameter, private :: flag_masks (0:15) = &
	int((/1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,-32768/),twobyteint)
integer(fourbyteint), parameter, private :: id_sat = 1, id_xml = 2, id_alias = 3, id_misc = 4, id_var = 5

include 'config.f90'

! Some private variables to keep

character(len=*), parameter, private :: default_short_optlist = 'S:X:vqV:C:P:A:F:R:L:Q:Z:', &
	default_long_optlist = ' t: h: args: sat: xml: debug: verbose log: quiet var: sel: cycle: pass: alias:' // &
	' cmp: compress: fmt: format: lat: lon: time: sla: limits: opt: mjd: sec: ymd: doy: quality-flag: quality_flag: region:'
character(len=rads_strl), save, private :: rads_optlist = default_short_optlist // default_long_optlist

private :: rads_traxxing, rads_free_sat_struct, rads_free_pass_struct, rads_free_var_struct, rads_set_limits_by_flagmask

!****f* rads/rads_init
! SUMMARY
! Initialize RADS4
!
! SYNTAX
! subroutine rads_init (S, sat, xml)
! type(rads_sat), intent(inout) :: S <or> S(:)
! character(len=*), intent(in), optional :: sat <or> sat(:)
! character(len=*), intent(in), optional :: xml(:)
!
! PURPOSE
! This routine initializes the <S> struct with the information pertaining
! to given satellite/mission phase <sat>, which is to be formed as 'e1',
! or 'e1g', or 'e1/g'. If no phase is specified, all mission phases will be
! queried.
!
! The <S> and <sat> arguments can either a single element or an array. In the
! latter case, one <S> struct will be initialized for each <sat>.
! To parse command line options after this, use rads_parse_cmd.
!
! Only if the <sat> argument is omitted, then the routine will parse
! the command line for arguments in the form:
! --sat=<sat> --cycle=<lo>,<hi>,<step> --pass=<lo>,<hi>,<step>
! --lim:<var>=<lo>,<hi> --lat=<lo>,<hi> --lon=<lo>,<hi> --alias:<var>=<var>
! --opt:<value>=<value> --opt=<value>,... --fmt:<var>=<value>
! or their equivalents without the = or : separators after the long name,
! or their equivalents without the initial --, or the short options -S, -C,
! -P, -L, -F
!
! The routine will read the satellite/mission specific setup XML files and
! store all the information in the stuct <S>. The XML files polled are:
!   $RADSDATAROOT/conf/rads.xml
!   ~/.rads/rads.xml
!   rads.xml
!   <xml> (from the optional array of file names)
!
! If more than one -S option is given, then all further options following
! this argument until the next -S option, plus all options prior to the
! first -S option will pertain to this mission.
!
! Execution will be halted when the dimension of <S> is insufficient to
! store information of multiple missions, or when required XML files are
! missing.
!
! The verbosity level can be controlled by setting rads_verbose before
! calling this routine (default = 0). The output unit for log info can
! be controlled by setting rads_log_unit up front (default = stdout).
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! sat      : (optional) Satellite/mission abbreviation
! xml      : (optional) Array of names of additional XML files to be loaded
!****-------------------------------------------------------------------
private :: rads_init_sat_0d, rads_init_sat_1d, &
	rads_init_cmd_0d, rads_init_cmd_1d
interface rads_init
	module procedure rads_init_sat_0d
	module procedure rads_init_sat_1d
	module procedure rads_init_cmd_0d
	module procedure rads_init_cmd_1d
end interface rads_init

!****f* rads/rads_end
! SUMMARY
! End RADS4
!
! SYNTAX
! subroutine rads_end (S)
! type(rads_sat), intent(inout) :: S <or> S(:)
!
! PURPOSE
! This routine ends RADS by freeing up all <S> space and other allocated
! global arrays.
!
! ARGUMENT
! S        : Satellite/mission dependent struct or array of structs
!****-------------------------------------------------------------------
private :: rads_end_0d, rads_end_1d
interface rads_end
	module procedure rads_end_0d
	module procedure rads_end_1d
end interface rads_end

!****f* rads/rads_get_var
! SUMMARY
! Read variable (data) from RADS4 file
!
! SYNTAX
! recursive subroutine rads_get_var (S, P, var, data, noedit)
! type(rads_sat), intent(inout) :: S
! type(rads_pass), intent(inout) :: P
! character(len=*) :: var
! <or> integer(fourbyteint) :: var
! <or> type(rads_var), intent(in) :: var
! real(eightbytereal), intent(out) :: data(:)
! logical, intent(in), optional :: noedit
!
! PURPOSE
! This routine loads the data from a single variable <var> into the
! buffer <data>. This command must be preceeded by <rads_open_pass>.
! The variable <var> can be addressed as a variable name, a RADS3-type
! field number or a varlist item.
!
! The array <data> must be at the correct size to contain the entire
! pass of data, i.e., it must have the dimension P%ndata.
! If no data are available and no default value and no secondary aliases
! then NaN is returned in the array <data>.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass dependent structure
! var      : (string) Name of the variable to be read.
!                     If <var> ends with % editing is skipped.
!            (integer) Field number.
!            (type(rads_var)) Variable struct (e.g. S%sel(i))
! data     : Data returned by this routine
! noedit   : (optional) Set to .true. to skip editing on limits and/or
!            quality flags; set to .false. to allow editing (default)
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_var, rads_err_memory, rads_err_source
!****-------------------------------------------------------------------
private :: rads_get_var_by_name, rads_get_var_by_var, rads_get_var_by_number, &
	rads_get_var_by_name_2d, rads_get_var_helper, rads_get_var_common
interface rads_get_var
	module procedure rads_get_var_by_name
	module procedure rads_get_var_by_var
	module procedure rads_get_var_by_number
	module procedure rads_get_var_by_name_2d
end interface rads_get_var

!****if* rads/rads_def_var
! SUMMARY
! Define variable(s) to be written to RADS data file
!
! SYNTAX
! subroutine rads_def_var (S, P, var, nctype, scale_factor, add_offset, ndims, coordinates, varid)
! type(rads_sat), intent(inout) :: S
! type(rads_pass), intent(inout) :: P
! type(rads_var), intent(in) :: var <or> var(:)
! <or> character(len=*), intent(in) :: var
! integer(fourbyteint), intent(in), optional :: nctype, ndims
! real(eightbytereal), intent(in), optional :: scale_factor, add_offset
! integer(fourbyteint), intent(out), optional :: varid
!
! PURPOSE
! This routine defines the variables to be written to a RADS data file by
! the reference to the variable struct of type(rads_var), or a
! number of variables referenced by an array of structures.
! This will create the NetCDF variable and its attributes in the file
! previously opened with rads_create_pass or rads_open_pass.
!
! The optional arguments <nctype>, <scale_factor>, <add_offset>, <ndims>  can
! be used to overrule those value in the <var%info> struct.
! The optional argument <coordinates> can be used to switch off the addition of the
! 'coordinates' attribute. Default is 'on' for most variables.
!
! ARGUMENTS
! S            : Satellite/mission dependent structure
! P            : Pass structure
! var          : Structure(s) of variable(s) of type(rads_var) or name of variable
! nctype       : (optional) Data type in NetCDF file
! scale_factor : (optional) Value of the scale_factor attribute
! add_offset   : (optional) Value of the add_offset attribute
! ndims        : (optional) Number of dimensions of the variable
! coordinates  : (optional) No not add coordinates attribute if .false. (default is .true.)
! varid        : variable ID of created NetCDF variable
!
! ERROR CODE
! S%error      : rads_noerr, rads_err_nc_var
!****-------------------------------------------------------------------
private :: rads_def_var_by_var_0d, rads_def_var_by_var_1d, rads_def_var_by_name
interface rads_def_var
	module procedure rads_def_var_by_var_0d
	module procedure rads_def_var_by_var_1d
	module procedure rads_def_var_by_name
end interface rads_def_var

!****if* rads/rads_put_var
! SUMMARY
! Write data for variable to RADS file
!
! SYNTAX
! subroutine rads_put_var (S, P, var, data, start)
! use netcdf
! use rads_netcdf
! use rads_misc
! type(rads_sat), intent(inout) :: S
! type(rads_pass), intent(inout) :: P
! type(rads_var), intent(in) :: var <or>
! character(len=*), intent(in) :: varname
! real(eightbytereal), intent(in) :: data(:) <or> data(:,:)
! integer(fourbyteint), optional, intent(in) :: start(:)
!
! PURPOSE
! This routine writes the data array <data> for the variable <var>
! (referenced by the structure of type(rads_var)) to the NetCDF file
! previously opened with rads_create_pass or rads_open_pass.
! Alternatively, the variable can be referenced by name <varname>.
!
! The data in <data> are in the original SI units (like [m] or [s])
! and will be converted to the internal units based on the values of
! <var%info%nctype>, <var%info%scale_factor>, <var%info%add_offset>.
!
! The argument <start> can be added to indicate an offset
! for storing the data from the first available position in the file.
! For example: start=101 first skips 100 records.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass structure
! var      : Structure of variable of type(rads_var)
! data     : Data to be written (in original (SI) units)
! start    : Position of the first data point in the file
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_nc_put
!****-------------------------------------------------------------------
private :: rads_put_var_helper, &
	rads_put_var_by_var_0d, rads_put_var_by_name_0d, &
	rads_put_var_by_var_1d, rads_put_var_by_var_1d_start, rads_put_var_by_name_1d, &
	rads_put_var_by_var_2d, rads_put_var_by_var_2d_start, rads_put_var_by_name_2d, &
	rads_put_var_by_var_3d, rads_put_var_by_var_3d_start, rads_put_var_by_name_3d
interface rads_put_var
	module procedure rads_put_var_by_var_0d
	module procedure rads_put_var_by_name_0d
	module procedure rads_put_var_by_var_1d
	module procedure rads_put_var_by_var_1d_start
	module procedure rads_put_var_by_name_1d
	module procedure rads_put_var_by_var_2d
	module procedure rads_put_var_by_var_2d_start
	module procedure rads_put_var_by_name_2d
	module procedure rads_put_var_by_var_3d
	module procedure rads_put_var_by_var_3d_start
	module procedure rads_put_var_by_name_3d
end interface rads_put_var

!****f* rads/rads_stat
! SUMMARY
! Print the RADS statistics for a given satellite
!
! SYNTAX
! subroutine rads_stat (S)
! type(rads_sat), intent(in) :: S <or> S(:)
! integer(fourbyteint), intent(in), optional :: unit
!
! PURPOSE
! This routine prints out the statistics of all variables that were
! processed per mission (indicated by scalar or array <S>), to the output
! on unit <rads_log_unit>.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
!****-------------------------------------------------------------------
private :: rads_stat_0d, rads_stat_1d
interface rads_stat
	module procedure rads_stat_0d
	module procedure rads_stat_1d
end interface rads_stat

!****if* rads/rads_parse_varlist
! SUMMARY
! Parse a string of variables
!
! SYNTAX
! subroutine rads_parse_varlist (S, string)
! type (rads_sat), intent(inout) :: S
! character(len=*), intent(in) :: string <or> string(:)
!
! PURPOSE
! Parse the string <string> of comma- or slash-separated names of variables.
! This will allocate the array S%sel and update counter S%nsel.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! string   : String of variables
!
! ERROR CODE
! S%error  : rads_err_var
!****-------------------------------------------------------------------
private :: rads_parse_varlist_0d, rads_parse_varlist_1d
interface rads_parse_varlist
	module procedure rads_parse_varlist_0d
	module procedure rads_parse_varlist_1d
end interface rads_parse_varlist

!****if* rads/rads_set_phase
! SUMMARY
! Set the pointer to satellite phase info within the S struct
!
! SYNTAX
! subroutine rads_set_phase (S, name <or> cycle <or> time, error)
! type(rads_sat), intent(inout) :: S
! character(len=*), intent(in) :: name <or>
! integer(fourbyteint), intent(in) :: cycle <or>
! real(eightbytereal), intent(in) :: time
! logical, optional, intent(out) :: error
!
! PURPOSE
! Set the pointer <S%phase> to the proper phase definitions for the mission
! phase. The selection of the phase can be done by phase name <name>,
! cycle number <cycle> or time <time>.
!
! In matching the phase name with the database, only the first letter is
! used. However the path to the directory with NetCDF files will use the
! entire phase name. In case multiple phases have the same name, the first
! one will be pointed to.
!
! For efficiency, the routine first checks if the current phase is correct.
!
! If no matching phase can be found the outcome depends on the use of the
! optional argument <error>:
! - If <error> is present, then <error> will refurn .true. (otherwise .false.)
! - If <error> is not present, the routine prints an error message and the
!   calling program stops.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! name     : Name of phase
! cycle    : Cycle number
! time     : UTC time in seconds since 1985
!****-------------------------------------------------------------------
private :: rads_set_phase_by_name, rads_set_phase_by_cycle, rads_set_phase_by_time
interface rads_set_phase
	module procedure rads_set_phase_by_name
	module procedure rads_set_phase_by_cycle
	module procedure rads_set_phase_by_time
end interface rads_set_phase

contains

!****if* rads/rads_init_sat_0d
! SUMMARY
! Initialize RADS4 by satellite
!
! SYNOPSIS
subroutine rads_init_sat_0d (S, sat, xml)
use rads_misc
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: sat
character(len=*), intent(in), optional :: xml(:)
!
! PURPOSE
! This routine initializes the <S> struct with the information pertaining to satellite
! and mission phase <sat>, which is to be formed as 'e1', or 'e1g', or 'e1/g'.
! If no phase is specified, a default is assumed (usually 'a').
! The routine will read the satellite and mission specific setup XML files and store
! all the information in the stuct <S>.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! sat      : Satellite/mission abbreviation
! xml      : Array of additional XML files to be loaded
! verbose  : (optional) verbosity level (default = 0)
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_xml_file, rads_err_xml_parse, rads_err_var
!****-------------------------------------------------------------------
integer(fourbyteint) :: i, n
logical, allocatable :: mask(:)

call rads_init_sat_struct (S)

! Store satellite specification
! S%sat, S%branch, S%satellite will be populated/overwritten later by rads_read_xml
! S%spec will then be replaced by the phase name (if specified)
if (len_trim(sat) < 2) call rads_exit ('Satellite/phase has fewer than 2 characters')
S%spec = sat
S%sat = sat(:2)
S%satellite = S%sat

! Set some global variables
S%dataroot = radsdataroot
call checkenv ('RADSDATAROOT',S%dataroot)
call checkenv ('HOME',S%userroot)
call get_command (S%command, status=i)
if (i < 0) S%command (len(S%command)-2:) = '...'

! Set all values in <S> struct to default
allocate (S%var(rads_var_chunk))
S%var = rads_var (null(), null(), null(), null(), null(), .false., rads_nofield)

! Read the global rads.xml setup file, the one in ~/.rads and the one in the current directory
call rads_read_xml (S, trim(S%dataroot)//'/conf/rads.xml')
if (S%error == rads_err_xml_file) call rads_exit ('Required XML file '//trim(S%dataroot)//'/conf/rads.xml does not exist')
call rads_read_xml (S, trim(S%userroot)//'/.rads/rads.xml')
if (index(S%command, 'radsreconfig') > 0) return ! Premature bailing from here in radsreconfig
call rads_read_xml (S, 'rads.xml')

! Now read the xml files specified on the command line, or the xml argument.
! If it is not an absolute path name, also look for it in $RADSDATAROOT/conf as well as in $HOME/.rads
if (present(xml)) then
	do i = 1,size(xml)
		call rads_read_xml (S, xml(i))
		if (S%error /= rads_err_xml_file) cycle
		if (xml(i)(1:1) == '/') call rads_exit ('Requested XML file '//trim(xml(i))//' does not exist')
		call rads_read_xml (S, trim(S%dataroot)//'/conf/'//xml(i))
		if (S%error /= rads_err_xml_file) cycle
		call rads_read_xml (S, trim(S%userroot)//'/.rads/'//xml(i))
		if (S%error /= rads_err_xml_file) cycle
		call rads_message ('Requested XML file '//trim(xml(i))//' does not exist in current directory,')
		call rads_exit ('nor in '//trim(S%dataroot)//'/conf, nor in '//trim(S%userroot)//'/.rads')
	enddo
endif

! If no phases are defined, then the satellite is not known
if (.not.associated(S%phases)) call rads_exit ('Satellite "'//S%sat//'" unknown')

! When NaN, set end time of a phase to the start time of the next phase
! When a phase contains only one cycle, compute the maximum pass number
n = size(S%phases)
do i = 1,n-1
	if (isnan_(S%phases(i)%end_time)) S%phases(i)%end_time = S%phases(i+1)%start_time
	if (S%phases(i)%cycles(1) == S%phases(i)%cycles(2)) &
		S%phases(i)%passes = S%phases(i)%ref_pass + 2 * (S%phases(i+1)%ref_orbit-S%phases(i)%ref_orbit) - 1
enddo
if (isnan_(S%phases(n)%end_time)) S%phases(n)%end_time = 2051222400d0 ! 2050-01-01

! Default phase is the first mission phase
S%phase => S%phases(1)

! When no phase/mission is specified: use the largest possible cycle and pass range
! When a phase/mission is specified: load the appropriate settings
! This can handle multiple "phase" specifications with the same name
allocate (mask(size(S%phases)))
mask = (S%spec == '' .or. S%phases%name(1:1) == S%spec(1:1))
if (S%spec /= '') call rads_set_phase (S, S%spec)
S%cycles(1) = minval(S%phases%cycles(1), mask)
S%cycles(2) = maxval(S%phases%cycles(2), mask)
S%passes(2) = maxval(S%phases%passes, mask)
deallocate (mask)

! Check if variables include the required <data> identifier
do i = 1,S%nvar
	if (S%var(i)%info%dataname(:1) == ' ') then
		call rads_message ('Variable "'//trim(S%var(i)%name)// &
			'" has no data source specified; assuming NetCDF variable "'//trim(S%var(i)%name)//'"')
		S%var(i)%info%dataname = S%var(i)%name
		S%var(i)%info%datasrc = rads_src_nc_var
	endif
enddo

! Link the time, lat, lon variables for later use
S%time => rads_varptr (S, 'time')
S%lat => rads_varptr (S, 'lat')
S%lon => rads_varptr (S, 'lon')

! Limit cycles based on time range and compute equator longitude limits
call rads_set_limits (S, 'time')
call rads_set_limits (S, 'lon')

! List the variables
if (rads_verbose >= 3) then
	do i = 1,S%nvar
		write (*,'(i3,1x,a,a,2i5,2f14.4)') i,S%var(i)%name,S%var(i)%info%name,S%var(i)%field,S%var(i)%info%limits
	enddo
endif
end subroutine rads_init_sat_0d

subroutine rads_init_sat_1d (S, sat, xml)
type(rads_sat), intent(inout) :: S(:)
character(len=*), intent(in) :: sat(:)
character(len=*), intent(in), optional :: xml(:)
integer :: i
if (size(S) /= size(sat)) call rads_exit ('Size of arguments "S" and "sat" in routine rads_init should be the same')
if (size(S) < 1) call rads_exit ('Size of argument "S" in routine rads_init should be at least 1')
do i = 1,size(S)
	call rads_init_sat_0d (S(i), sat(i), xml)
enddo
end subroutine rads_init_sat_1d

subroutine rads_init_cmd_0d (S)
type(rads_sat), intent(inout) :: S
integer :: sopt(1), nsat
!
call rads_load_options (nsat)
if (nsat < 1) call rads_exit ('Failed to find any "-S" or "--sat" options on command line')
if (nsat > 1) call rads_exit ('Too many "-S" or "--sat" options on command line')
sopt = minloc (rads_opt%id, rads_opt%id==10+id_sat) ! Position of the first (only) -S option
!
! The next three pack functions isolate the -X option arguments, the command line
! arguments associcated with the -S option, and  all -V option arguments.
call rads_init_sat_0d (S, rads_opt(sopt(1))%arg, &
	pack(rads_opt%arg, mod(rads_opt%id,10)==id_xml))
call rads_parse_options (S, pack(rads_opt, mod(rads_opt%id,10)==id_alias))
call rads_parse_options (S, pack(rads_opt, mod(rads_opt%id,10)==id_misc))
call rads_parse_varlist (S, pack(rads_opt%arg, mod(rads_opt%id,10)==id_var))
end subroutine rads_init_cmd_0d

subroutine rads_init_cmd_1d (S)
type(rads_sat), intent(inout) :: S(:)
integer :: i, sopt(1), nsat
!
call rads_load_options (nsat)
if (nsat < 1) call rads_exit ('Failed to find any "-S" or "--sat" options on command line')
if (nsat > size(S)) call rads_exit ('Too many "-S" or "--sat" options on command line')
!
do i = 1,nsat
	sopt = minloc (rads_opt%id, rads_opt%id==i*10+id_sat) ! Position of i-th -S option
	! The next three pack functions isolate the -X option arguments, the command line
	! arguments associcated with the i-th -S option, and  all -V option arguments.
	call rads_init_sat_0d (S(i), rads_opt(sopt(1))%arg, &
		pack(rads_opt%arg, rads_opt%id==id_xml .or. rads_opt%id==i*10+id_xml))
	call rads_parse_options (S(i), pack(rads_opt, rads_opt%id==id_alias .or. rads_opt%id==i*10+id_alias))
	call rads_parse_options (S(i), pack(rads_opt, rads_opt%id==id_misc .or. rads_opt%id==i*10+id_misc))
	call rads_parse_varlist (S(i), pack(rads_opt%arg, mod(rads_opt%id,10)==id_var))
enddo
!
! Blank out the rest
do i = nsat+1,size(S)
	call rads_init_sat_struct (S(i))
enddo
end subroutine rads_init_cmd_1d

!*****f* rads/rads_init_sat_struct
! SUMMARY
! Initialize empty rads_sat struct
!
! SYNOPSIS
subroutine rads_init_sat_struct (S)
type(rads_sat), intent(inout) :: S
!
! PURPOSE
! This routine initializes the <S> struct with the bare minimum.
! It is later updated in rads_init_sat_0d.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
!****-------------------------------------------------------------------
! gfortran 4.3.4 to 4.4.1 segfault on the next line if this routine is made pure or elemental,
! so please leave it as a normal routine.
S = rads_sat ('', '', '', '', '', null(), '', 1d0, (/13.8d0, nan/), 90d0, nan, nan, nan, 1, 1, rads_noerr, &
	0, 0, 0, 0, 0, .false., '', 0, null(), null(), null(), null(), null(), null(), null(), null())
end subroutine rads_init_sat_struct

!****if* rads/rads_free_sat_struct
! SUMMARY
! Free all allocated memory from rads_sat struct and reinitialise
!
! SYNOPSIS
subroutine rads_free_sat_struct (S)
type(rads_sat), intent(inout) :: S
!
! PURPOSE
! This routine frees the <S> struct of type rads_sat and all allocated
! memory in it. It then reinitialises a clean struct.
!
! ARGUMENT
! S        : Satellite/mission dependent structure
!****-------------------------------------------------------------------
integer(fourbyteint) :: i, ios
if (S%sat == '') return
do i = S%nvar,1,-1
	call rads_free_var_struct (S, S%var(i), .false.)
enddo
do i = 1,size(S%phases)
	deallocate (S%phases(i)%subcycles, stat=ios)
enddo
deallocate (S%glob_att, S%var, S%sel, S%phases, S%excl_cycles, stat=ios)
call rads_init_sat_struct (S)
end subroutine rads_free_sat_struct

!*****f* rads/rads_init_pass_struct
! SUMMARY
! Initialize empty rads_pass struct
!
! SYNOPSIS
subroutine rads_init_pass_struct (S, P)
type(rads_sat), target, intent(in) :: S
type(rads_pass), intent(inout) :: P
!
! PURPOSE
! This routine initializes the <P> struct with the bare minimum.
! This is only really necessary prior to calling rads_create_pass.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass dependent structure
!****-------------------------------------------------------------------
! gfortran 4.4.1 segfaults on the next line if this routine is made pure or elemental,
! so please leave it as a normal routine.
P = rads_pass ('', null(), nan, nan, nan, nan, null(), null(), .false., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
	rads_file (0, ''), null(), null())
P%S => S
end subroutine rads_init_pass_struct

!****if* rads/rads_free_pass_struct
! SUMMARY
! Free all allocated memory from rads_pass struct and reinitialise
!
! SYNOPSIS
subroutine rads_free_pass_struct (S, P)
type(rads_sat), intent(in) :: S
type(rads_pass), intent(inout) :: P
!
! PURPOSE
! This routine frees the <P> struct of type rads_pass.
! Afterwards the routine reinitialises a clean struct.
!
! ARGUMENT
! S        : Satellite/mission dependent structure
! P        : Pass dependent structure
!****-------------------------------------------------------------------
integer(fourbyteint) :: ios
deallocate (P%history, P%flags, P%tll, stat=ios)
call rads_init_pass_struct (S, P)
end subroutine rads_free_pass_struct

!****if* rads/rads_free_var_struct
! SUMMARY
! Free all allocated memory from rads_var struct
!
! SYNOPSIS
pure subroutine rads_free_var_struct (S, var, alias)
use rads_grid
type(rads_sat), intent(inout) :: S
type(rads_var), intent(inout) :: var
logical, intent(in) :: alias
!
! PURPOSE
! This routine frees the <var> struct of type rads_var and all allocated
! memory in it. Use alias = .true. when called from rads_set_alias, so
! that it prevents removing info structs used elsewhere.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! var      : Variable struct to be freed
! alias    : .true. if called to create alias
!****-------------------------------------------------------------------
integer(fourbyteint) :: i, n, ios
type(rads_varinfo), pointer :: info
info => var%info ! This is needed to make the associated () work later on
if (.not.associated(var%info)) return
if (alias .and. associated(var%name,info%name)) then
	! This is an original variable. Do not alias it if it is used more than once
	n = 0
	do i = 1,S%nvar
		if (associated(S%var(i)%info,var%info)) n = n + 1
		if (associated(S%var(i)%info,var%inf1)) n = n + 1
		if (associated(S%var(i)%info,var%inf2)) n = n + 1
	enddo
	if (n > 1) then
		S%error = rads_err_alias
		return
	endif
endif

! Unassociate or deallocate long_name
if (associated(var%long_name,info%long_name)) then
	nullify (var%long_name)
else
	deallocate (var%long_name, stat=ios)
endif

! Unassociate or deallocate name and info struct
if (associated(var%name,info%name)) then
	if (associated(var%info%grid)) call grid_free(var%info%grid)
	nullify (var%name)
	deallocate (var%info, stat=ios)
else
	deallocate (var%name, stat=ios)
	nullify (var%info)
endif

! Clean out rest of var
nullify (var%inf1,var%inf2)
end subroutine rads_free_var_struct

!****f* rads/rads_set_options
! SUMMARY
! Specify the list of command specific options
!
! SYNPOSIS
subroutine rads_set_options (optlist)
character(len=*), intent(in), optional :: optlist
!
! PURPOSE
! Add the command specific options to the list of common RADS options.
! The argument <optlist> needs to have the same format as in the routine
! <getopt> in the <rads_misc> module. The short options will be placed
! before the common ones, the long options will be placed after them.
!
! ARGUMENT
! optlist  : (optional) list of command specific short and long options
!****-------------------------------------------------------------------
integer :: i
if (.not.present(optlist)) return

! Find first space in optlist
i = index(optlist, ' ')
! Now merge it appropriately with default list
if (i == 0) then	! Only short options specified in oplist
	rads_optlist = optlist // default_short_optlist // default_long_optlist
else if (i == 1) then	! Only long options specified in optlist
	rads_optlist = default_short_optlist // optlist // default_long_optlist
else ! Short options in front, long options in the middle
	rads_optlist = optlist(1:i-1) // default_short_optlist // optlist(i:) // default_long_optlist
endif
end subroutine rads_set_options

!****if* rads/rads_load_options
! SUMMARY
! Extract options from command line
!
! SYNOPSIS
subroutine rads_load_options (nsat)
use rads_misc
integer, intent(out) :: nsat
!
! PURPOSE
! Scan the command line for all common RADS 4 command line options
! like --cycle, --sat, -P, etc, and store only these, not any other
! command line options, in an external array <rads_option>.
! In addition any arguments from a file pointed to by the --args <file>
! option are loaded as well.
!
! The array <opt%id> will contain the identifiers of the various options.
! They are given as <nsat>*10+<type>, where <nsat> is 0 before the first
! -S option, 1 starting with the first -S option, 2 starting with the second
! -S option, etc, and <type> is 1 for -S, 2 for -X, 4 for -V, 3 for all
! other known options, 0 for all other.
!
! The external variable <rads_nopt> will be set to the number of command
! line options saved.
!
! The value <nsat> returns the total amount of -S options.
!
! ARGUMENTS
! nsat     : Number of '-S' or '--sat' options
!****-------------------------------------------------------------------
integer :: ios, iunit, i, nopt
type(rads_option), pointer :: opt(:)

! Initialize
nopt = 1
nsat = 0
iunit = 0
allocate (rads_opt(rads_optl))
opt => rads_opt
opt = rads_option ('!', '', 9999)
call getopt_reset

! If any option is 'args=' then load the options from file
! Otherwise load the options from the command line

do
	! Check if there are not too many options
	if (nopt > size(opt)) call rads_exit ('Too many command line options')

	! Get next option and its argument
	call getopt (rads_optlist, opt(nopt)%opt, opt(nopt)%arg, iunit)

! Now hunt for '-S', --sat', '-v', '--verbose', '--debug', '-X', and '--xml' among the command line options
! The array opt%id is filled to indicate where to find -S, -X, -v|q, -V and other options
! Also look for '--args' option to load arguments from file

	select case (opt(nopt)%opt)
	case ('!') ! End of argument list
		if (iunit == 0) then
			nopt = nopt - 1 ! Because we count always one ahead
			exit ! Done with command line
		else
			close (iunit) ! Done with argument file, continue with command line
			iunit = 0
		endif

	! We do not save any of the next options
	case (':', '::')
		call rads_opt_error (opt(nopt)%opt, opt(nopt)%arg)
	case ('q', 'quiet')
		rads_verbose = -1
	case ('v', 'verbose')
		rads_verbose = rads_verbose + 1
	case ('debug')
		read (opt(nopt)%arg, *, iostat=ios) rads_verbose
		if (ios /= 0) call rads_opt_error (opt(nopt)%opt, opt(nopt)%arg)
	case ('log')
		if (opt(nopt)%arg == '-' .or. opt(nopt)%arg == 'stdout') then
			rads_log_unit = stdout
		else if (opt(nopt)%arg == '+' .or. opt(nopt)%arg == 'stderr') then
			rads_log_unit = stderr
		else
			rads_log_unit = getlun()
			open (rads_log_unit, file=opt(nopt)%arg, status='replace', iostat=ios)
			if (ios /= 0) rads_log_unit = stdout
		endif
	case ('args')
		! Open file with command line arguments. If not available, ignore
		iunit = getlun()
		open (iunit, file=opt(nopt)%arg, status='old', iostat=ios)
		if (ios /= 0) iunit = 0

	! The following options we memorise. They will later be parsed in the order:
	! id_sat, id_xml, id_alias, id_misc, id_var
	case (' ') ! Non-option arguments
		opt(nopt)%id = nsat*10
		nopt = nopt + 1
	case ('S', 'sat')
		nsat = nsat + 1
		opt(nopt)%id = nsat*10+id_sat
		nopt = nopt + 1
	case ('X', 'xml')
		opt(nopt)%id = nsat*10+id_xml
		nopt = nopt + 1
	case ('A', 'alias')
		opt(nopt)%id = nsat*10+id_alias
		nopt = nopt + 1
	case ('V', 'var', 'sel')
		opt(nopt)%id = nsat*10+id_var
		nopt = nopt + 1
	case default ! Here for other known options to be parsed in rads_parse_options
		opt(nopt)%id = nsat*10+id_misc
		nopt = nopt + 1
	end select
enddo

if (rads_verbose >= 2) then
	write (*,'(i6,1x,a)') nopt, ' command line options:'
	do i = 1, nopt
		write (*,'(i6,5a,i6)') i,' optopt = ',trim(opt(i)%opt),'; optarg = ',trim(opt(i)%arg),'; optid = ',opt(i)%id
	enddo
endif
rads_nopt = nopt
end subroutine rads_load_options

!****if* rads/rads_parse_options
! SUMMARY
! Parse RADS options
!
! SYNOPSIS
subroutine rads_parse_options (S, opt)
use rads_time
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_option), intent(in) :: opt(:)
!
! PURPOSE
! This routine fills the <S> struct with the information pertaining
! to a given satellite and mission as read from the options in the
! array <opt>.
! The -S or --sat option is ignored (it should be dealt with separately).
! The -V or --var option is not parsed, but a pointer is returned to the element
! of the list of variables following -V or --var. Same for --sel, and the
! backward compatible options sel= and var=.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! opt      : Array of options (usually command line arguments)
!****-------------------------------------------------------------------
integer :: i
do i = 1,size(opt)
	call rads_parse_option (opt(i))
enddo
if (S%error /= rads_noerr) call rads_exit ('Error parsing command line arguments')

contains

subroutine rads_parse_option (opt)
type(rads_option), intent(in) :: opt
integer :: j, k0, k1, ios
real(eightbytereal) :: val(2)
! Scan a single command line argument (option)
val = nan
j = index(opt%arg, '=')
select case (opt%opt)
case ('C', 'cycle')
	S%cycles(2) = -1
	S%cycles(3) = 1
	call read_val (opt%arg, S%cycles, '/-', iostat=ios)
	if (ios > 0) call rads_opt_error (opt%opt, opt%arg)
	if (S%cycles(2) < 0) S%cycles(2) = S%cycles(1)
case ('P', 'pass')
	if (opt%arg(:1) == 'a') then ! Only ascending passes
		S%passes(1) = S%passes(1)/2*2+1
		S%passes(3) = 2
	else if (opt%arg(:1) == 'd') then ! Only descending passes
		S%passes(1) = (S%passes(1)+1)/2*2
		S%passes(3) = 2
	else
		S%passes(2) = -1
		S%passes(3) = 1
		call read_val (opt%arg, S%passes, '/-', iostat=ios)
		if (ios > 0) call rads_opt_error (opt%opt, opt%arg)
		if (S%passes(2) < 0) S%passes(2) = S%passes(1)
	endif
case ('A', 'alias')
	if (opt%arg(j+1:j+1) == ',' .or. opt%arg(j+1:j+1) == '/') then
		! -Aalias=,var is processed as -Aalias=alias,var
		call rads_set_alias (S, opt%arg(:j-1), opt%arg(:j-1) // opt%arg(j+1:))
	else
		call rads_set_alias (S, opt%arg(:j-1), opt%arg(j+1:))
	endif
case ('F', 'fmt', 'format')
	call rads_set_format (S, opt%arg(:j-1), opt%arg(j+1:))
case ('Z', 'cmp', 'compress')
	call rads_set_compress (S, opt%arg(:j-1), opt%arg(j+1:))
case ('R', 'region')
	call rads_set_region (S, opt%arg)
case ('lat', 'lon', 'sla')
	call rads_set_limits (S, opt%opt, string=opt%arg, iostat=ios)
	if (ios > 0) call rads_opt_error (opt%opt, opt%arg)
case ('L', 'limits')
	call rads_set_limits (S, opt%arg(:j-1), string=opt%arg(j+1:), iostat=ios)
	if (ios > 0) call rads_opt_error (opt%opt, opt%arg)
case ('Q', 'quality-flag', 'quality_flag')
	if (j > 0) then
		call rads_set_quality_flag (S, opt%arg(:j-1), opt%arg(j+1:))
	else
		call rads_set_quality_flag (S, 'sla', opt%arg)
	endif

! The next are only for compatibility with RADS 3
case ('h')
	call rads_set_limits (S, 'sla', string=opt%arg, iostat=ios)
	if (ios > 0) call rads_opt_error (opt%opt, opt%arg)
case ('opt')
	if (j > 0) then
		if (len_trim(opt%arg) == j+1) then
			call rads_set_alias (S, opt%arg(:j-1), opt%arg(:j-1)//'0'//opt%arg(j+1:))
			call rads_set_quality_flag (S, 'sla', opt%arg(:j-1)//'0'//opt%arg(j+1:))
		else
			call rads_set_alias (S, opt%arg(:j-1), opt%arg(:j-1)//opt%arg(j+1:))
			call rads_set_quality_flag (S, 'sla', opt%arg(:j-1)//opt%arg(j+1:))
		endif
	else
		k1 = 0
		do
			if (.not.next_word (opt%arg, k0, k1)) exit
			if (k1 == k0) cycle
			call rads_set_alias (S, opt%arg(k0:k1-3), opt%arg(k0:k1-1))
			call rads_set_quality_flag (S, 'sla', opt%arg(k0:k1-1))
		enddo
	endif
case ('time', 't:', 'sec', 'mjd', 'doy', 'ymd') ! Finally try date arguments
	call dateopt(opt%opt, opt%arg, val(1), val(2), iostat=ios)
	if (ios > 0) call rads_opt_error (opt%opt, opt%arg)
	call rads_set_limits (S, 'time', val(1), val(2))
end select
end subroutine rads_parse_option

end subroutine rads_parse_options

!****if* rads/rads_parse_varlist_0d
! SUMMARY
! Parse a string of variables
!
! SYNOPSIS
subroutine rads_parse_varlist_0d (S, string)
use rads_misc
type (rads_sat), intent(inout) :: S
character(len=*), intent(in) :: string
!
! Parse the string <string> of comma- or slash-separated names of variables.
! This will allocate the array S%sel and update counter S%nsel.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! string   : String of variables
!
! ERROR CODE
! S%error  : rads_err_var
!****-------------------------------------------------------------------
integer(fourbyteint) :: i0, i1, n, noedit
type(rads_var), pointer :: temp(:), var

! Count the number of words first
i1 = 0
n = 0
do
	if (.not.next_word (string, i0, i1)) exit
	if (i1 == i0) cycle
	n = n + 1
enddo

! If no variables, exit
if (n == 0) then
	call rads_error (S, rads_err_var, 'No variables selected')
	return
endif

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
i1 = 0
do
	if (.not.next_word (string, i0, i1)) exit
	if (i1 == i0) cycle
	noedit = 1
	if (string(i1-1:i1-1) == rads_noedit) noedit = 2
	var => rads_varptr(S, string(i0:i1-noedit))
	if (.not.associated(var)) then
		call rads_error (S, rads_err_var, 'Unknown variable "'//string(i0:i1-noedit)//'" removed from list of variables')
		cycle
	endif
	S%nsel = S%nsel + 1
	S%sel(S%nsel) = var
	if (noedit == 2) S%sel(S%nsel)%noedit = .true.
	! Switch on multi-Hz output when a variable is 2-dimensional
	if (var%info%ndims > 1) S%n_hz_output = .true.
enddo

end subroutine rads_parse_varlist_0d

subroutine rads_parse_varlist_1d (S, string)
type (rads_sat), intent(inout) :: S
character(len=*), intent(in) :: string(:)
integer :: i
do i = 1,size(string)
	call rads_parse_varlist_0d (S, string(i))
enddo
end subroutine rads_parse_varlist_1d

!****if* rads/rads_parse_cmd
!*rads_parse_cmd -- Parse command line options
!
! SYNOPSIS
subroutine rads_parse_cmd (S)
type(rads_sat), intent(inout) :: S
!
! PURPOSE
! After initialising RADS for a single satellite using rads_init_sat (S,sat)
! this routine can be used to update the <S> struct with information from
! the command line options. In contrast to rads_init_sat (S) this does not
! require that -S or --sat is one of the command line options.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! optlist  : (optional) list of command specific short and long options
!****-------------------------------------------------------------------
integer :: nsat
call rads_load_options (nsat)
call rads_parse_options (S, rads_opt(1:rads_nopt))
call rads_parse_varlist (S, pack(rads_opt%arg, mod(rads_opt%id,10)==id_var))
end subroutine rads_parse_cmd

subroutine rads_end_0d (S)
type(rads_sat), intent(inout) :: S
integer :: ios
deallocate (rads_opt, stat=ios)
call rads_free_sat_struct (S)
if (rads_log_unit > stdout) close (rads_log_unit)
end subroutine rads_end_0d

subroutine rads_end_1d (S)
type(rads_sat), intent(inout) :: S(:)
integer :: ios, i
deallocate (rads_opt, stat=ios)
do i = 1,size(S)
	call rads_free_sat_struct (S(i))
enddo
if (rads_log_unit > stdout) close (rads_log_unit)
end subroutine rads_end_1d

!****f* rads/rads_open_pass
! SUMMARY
! Open RADS pass file
!
! SYNOPSIS
subroutine rads_open_pass (S, P, cycle, pass, rw, echofilepaths)
use netcdf
use rads_netcdf
use rads_time
use rads_misc
use rads_geo
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
integer(fourbyteint), intent(in) :: cycle, pass
logical, intent(in), optional :: rw
logical, intent(in), optional :: echofilepaths
!
! PURPOSE
! This routine opens a NetCDF file for access to the RADS machinery.
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
! If the pass is rejected based on the above critetia or when no NetCDF
! file exists, S%error returns the warning value rads_warn_nc_file.
! If the file cannot be read properly, rads_err_nc_parse is returned.
! Also, in both cases, P%ndata will be set to zero.
!
! By default the file is opened for reading only. Specify perm=nf90_write to
! open for reading and writing.
! The file opened with this routine should be closed by using rads_close_pass.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass structure
! cycle    : Cycle number
! pass     : Pass number
! rw       : (optional) Set read/write permission (def: read only)
! echofilepaths : (optional) Print to standard output the file path before
!            it is being checked for existence.
!
! ERROR CODE
! S%error  : rads_noerr, rads_warn_nc_file, rads_err_nc_parse
!****-------------------------------------------------------------------
character(len=40) :: date
character(len=5) :: hz
character(len=rads_strl) :: string
integer(fourbyteint) :: i, j1, j2, k, ascdes, ncid
real(eightbytereal), pointer :: temp(:,:)
logical :: error

! Initialise
S%error = rads_warn_nc_file
call rads_init_pass_struct (S, P)
P%cycle = cycle
P%pass = pass
ascdes = modulo(pass,2)	! 1 if ascending, 0 if descending
P%ndata = 0

if (rads_verbose >= 2) write (*,'(a,a3,2i5)') 'Checking sat/cycle/pass : ', S%sat, cycle, pass

! If the cycle is out of range for the current phase, look for a new phase
call rads_set_phase (S, cycle, error)

! Do checking on cycle limits
if (error) then
	S%pass_stat(1) = S%pass_stat(1) + 1
	return
endif

! Do check on excluded cycles
if (associated(S%excl_cycles)) then
	if (any(S%excl_cycles%list(:S%excl_cycles%n) == cycle)) then
		S%pass_stat(1) = S%pass_stat(1) + 1
		return
	endif
endif

! Do checking on pass limits (which may include new phase limits)
if (pass < S%passes(1) .or. pass > S%passes(2) .or. pass > S%phase%passes) then
	S%pass_stat(2) = S%pass_stat(2) + 1
	return
endif

! Predict equator crossing info
call rads_predict_equator (S, P)

! Do checking of pass ends on the time criteria (only when such are given)
if (.not.all(isnan_(S%time%info%limits))) then
	if (P%end_time + 300d0 < S%time%info%limits(1) .or. P%start_time - 300d0 > S%time%info%limits(2)) then ! Allow 5 minute slop
		if (rads_verbose >= 2) write (*,'(a,4f15.3)') 'Bail out on estimated time:', S%time%info%limits, P%start_time, P%end_time
		S%pass_stat(3) = S%pass_stat(3) + 1
		return
	endif
endif

! Do checking of equator longitude on the longitude criteria (only when those are limiting)
if (S%eqlonlim(ascdes,2) - S%eqlonlim(ascdes,1) < 360d0) then
	if (checklon(S%eqlonlim(ascdes,:),P%equator_lon)) then
		if (rads_verbose >= 2) write (*,'(a,3f12.6)') 'Bail out on estimated equator longitude:', &
			S%eqlonlim(ascdes,:), P%equator_lon
		S%pass_stat(4+ascdes) = S%pass_stat(4+ascdes) + 1
		return
	endif
endif

! Update pass opening stats
S%pass_stat(6+ascdes) = S%pass_stat(6+ascdes) + 1

! Open pass file
600 format (a,'/',a,'/',a,'/c',i3.3,'/',a2,'p',i4.4,'c',i3.3,'.nc')
write (P%fileinfo(1)%name, 600) trim(S%dataroot), trim(S%branch(1)), trim(S%phase%name), cycle, S%sat, pass, cycle

! optional echo of the file in which we are looking for data
! ...intentionally do this prior to testing if file exists
601 format ('# checkfile: ',a)
if (present(echofilepaths)) then
    if (echofilepaths) then
        write (*,601) trim(P%fileinfo(1)%name)
    endif
endif

if (present(rw)) then
	P%rw = rw
else
	P%rw = .false.
endif
if (P%rw) then
	if (rads_verbose >= 2) write (*,'(2a)') 'Opening for read/write: ',trim(P%fileinfo(1)%name)
	if (nft(nf90_open(P%fileinfo(1)%name,nf90_write,ncid))) return
else
	if (rads_verbose >= 2) write (*,'(2a)') 'Opening for read only: ',trim(P%fileinfo(1)%name)
	if (nft(nf90_open(P%fileinfo(1)%name,nf90_nowrite,ncid))) return
endif
P%fileinfo(1)%ncid = ncid

! Get the time dimension
S%error = rads_err_nc_parse
k = index(S%time%info%dataname, ' ') - 1
if (nft(nf90_inq_dimid(ncid,S%time%info%dataname(:k),i))) return
if (nft(nf90_inquire_dimension(ncid,i,len=P%ndata))) return

! Get the optional second (n-Hz) dimension
if (nf90_inq_dimid(ncid,'meas_ind',i) + nf90_inquire_dimension(ncid,i,len=P%n_hz) /= nf90_noerr) P%n_hz = 0

! Read global attributes
if (nft(nf90_get_att(ncid,nf90_global,'equator_longitude',P%equator_lon))) return
P%equator_lon = S%lon%info%limits(1) + modulo (P%equator_lon - S%lon%info%limits(1), 360d0)
if (nft(nf90_get_att(ncid,nf90_global,'equator_time',date))) return
P%equator_time = strp1985f(date)
if (nft(nf90_get_att(ncid,nf90_global,'first_meas_time',date))) return
P%start_time = strp1985f(date)
if (nft(nf90_get_att(ncid,nf90_global,'last_meas_time',date))) return
P%end_time = strp1985f(date)
if (nft(nf90_get_att(ncid,nf90_global,'cycle_number',i)) .or. i /= cycle) return
if (nft(nf90_get_att(ncid,nf90_global,'pass_number',i)) .or. i /= pass) return
if (rads_verbose >= 3) write (*,'(a,3f15.3,f12.6)') 'Start/end/equator time/longitude = ', &
	P%equator_time, P%start_time, P%end_time, P%equator_lon

S%error = rads_noerr

! Return if the file is empty
if (P%ndata == 0) return

! Update number of "read" measurements
S%total_read = S%total_read + P%ndata

! Read history
if (nff(nf90_inquire_attribute(ncid,nf90_global,'history',attnum=k))) then
	allocate (P%history)
	call nfs(nf90_get_att(ncid,nf90_global,'history',P%history))
	if (nft(nf90_get_att(ncid,nf90_global,'original',P%original))) P%original = ''
else if (nff(nf90_inquire_attribute(ncid,nf90_global,'log01',attnum=k))) then ! Read logs (RADS3)
	allocate (P%history)
	i = 0
	do
		i = i + 1
		write (date, '("log",i2.2)') i
		if (nft(nf90_get_att(ncid,nf90_global,date,string))) exit
		j1 = index(string, '|', .true.) ! Index of last '|'
		j2 = index(string, ':') ! Index of first ':'
		if (i == 1) then
			! If log01 has original file name, save it
			k = index(string, 'RAW data from file')
			if (k > 0) then
				P%original = string(k+19:)
			else
				k = index(string, 'RAW data from')
				if (k > 0) then
					P%original = string(k+14:)
				else
					P%original = string(j2+1:)
				endif
			endif
			! Write date and command to history
			P%history = string(:11) // ':' // string(j1+1:j2-1)
		else
			P%history = string(:11) // ':' // string(j1+1:j2-1) // rads_linefeed // P%history
		endif
		P%nlogs = i
	enddo
endif

! Read time, lon, lat
! If read-only: also check their limits
allocate (P%tll(P%ndata,3))
call rads_get_var_common (S, P, S%time, P%tll(:,1), P%rw)
call rads_get_var_common (S, P, S%lat, P%tll(:,2), P%rw)
call rads_get_var_common (S, P, S%lon, P%tll(:,3), P%rw)
P%time_dims = 1

! If requested, check for distance to centroid
if (.not.P%rw .and. S%centroid(3) > 0d0) then
	do i = 1,P%ndata
		if (any(isnan_(P%tll(i,2:3)))) cycle
		if (sfdist(P%tll(i,2)*rad, P%tll(i,3)*rad, S%centroid(2), S%centroid(1)) > S%centroid(3)) P%tll(i,2:3) = nan
	enddo
endif

! Look for first non-NaN measurement
do i = 1, P%ndata
	if (.not.any(isnan_(P%tll(i,:)))) exit
enddo
if (i > P%ndata) then ! Got no non-NaNs
	P%ndata = 0
	return
endif
P%first_meas = i

! Now look for last non-NaN measurement
do i = P%ndata, P%first_meas, -1
	if (.not.any(isnan_(P%tll(i,:)))) exit
enddo
P%last_meas = i

! Open files in other (optional) branches
do i = 2, rads_max_branches
	if (S%branch(i) == '') exit
	write (P%fileinfo(i)%name, 600) &
		trim(S%dataroot), trim(S%branch(i)), trim(S%phase%name), cycle, S%sat, pass, cycle
	if (P%rw) then
		if (rads_verbose >= 2) write (*,'(2a)') 'Opening for read/write: ',trim(P%fileinfo(i)%name)
		if (nft(nf90_open(P%fileinfo(i)%name,nf90_write,ncid))) return
	else
		if (rads_verbose >= 2) write (*,'(2a)') 'Opening for read only: ',trim(P%fileinfo(i)%name)
		if (nft(nf90_open(P%fileinfo(i)%name,nf90_nowrite,ncid))) return
	endif
	if (nft(nf90_inquire_dimension(ncid,1,len=k))) return
	if (k /= P%ndata) call rads_message ('Dimension in ' // trim(P%fileinfo(i)%name) // ' does not match ' // &
		trim(P%fileinfo(1)%name))
	P%fileinfo(i)%ncid = ncid
enddo

! If subset is requested, reallocate appropriately sized time, lat, lon arrays
! If multi-Hertz data: load multi-Hertz fields
i = P%last_meas - P%first_meas + 1
if (S%n_hz_output .and. P%n_hz > 0) then
	P%ndata = i * P%n_hz
	deallocate (P%tll)
	allocate (P%tll(P%ndata,3))
	write (hz, '("_",i2.2,"hz")') P%n_hz
	call rads_get_var (S, P, 'time'//hz, P%tll(:,1))
	call rads_get_var (S, P, 'lat' //hz, P%tll(:,2))
	call rads_get_var (S, P, 'lon' //hz, P%tll(:,3))
	P%time_dims = 2
else if (i /= P%ndata) then
	P%ndata = i
	allocate (temp(P%ndata,3))
	temp = P%tll(P%first_meas:P%last_meas,:)
	deallocate (P%tll)
	P%tll => temp
endif

! Successful ending; update number of measurements in region
S%error = rads_noerr
S%total_inside = S%total_inside + P%ndata

end subroutine rads_open_pass

!****f* rads/rads_close_pass
! SUMMARY
! Close RADS pass file
!
! SYNOPSIS
subroutine rads_close_pass (S, P, keep)
use netcdf
use rads_netcdf
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
logical, intent(in), optional :: keep
!
! PURPOSE
! This routine closes a NetCDF file previously opened by rads_open_pass.
! The routine will reset the ncid element of the <P> structure to
! indicate that the passfile is closed.
! If <keep> is set to .true., then the history, flags, time, lat, and lon
! elements the <P> structure are kept. I.e., they are not deallocated
! but only their links are removed. Otherwise, they are deallocated along
! with the log entries.
! A second call to rads_close_pass without the keep argment can subsequently
! deallocate the time, lat and lon elements of the <P> structure.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass structure
! keep     : Keep the P%tll matrix (destroy by default)
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_nc_close
!****-------------------------------------------------------------------
integer :: i
integer(fourbyteint) :: ncid
S%error = rads_noerr
do i = 1,rads_max_branches
	ncid = P%fileinfo(i)%ncid
	if (ncid > 0 .and. nft(nf90_close(ncid))) S%error = rads_err_nc_close
enddo
P%fileinfo = rads_file (0, '')
if (.not.present(keep) .or. .not.keep) call rads_free_pass_struct (S, P)
end subroutine rads_close_pass

!****if* rads/rads_get_var_by_number
! SUMMARY
! Read variable (data) from RADS4 file by integer
!
! SYNOPSIS
recursive subroutine rads_get_var_by_number (S, P, field, data, noedit)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
integer(fourbyteint), intent(in) :: field
real(eightbytereal), intent(out) :: data(:)
logical, intent(in), optional :: noedit
!
! PURPOSE
! This routine is added for backward compatibility with the field numbers
! in RADS3. It functions the same as rads_get_var, except that it
! uses the old field number as indicator for the variable.
!****-------------------------------------------------------------------
character(len=4) :: name
integer :: i
logical :: skip_edit

S%error = rads_noerr
if (rads_get_var_helper (S, P, data)) return ! Check data sizes

! Do we need to skip editing?
if (present(noedit)) then
	skip_edit = noedit
else
	skip_edit = .false.
endif

! Look for field number in variable list
do i = 1,S%nvar
	if (any(S%var(i)%field == field)) then
		call rads_get_var_common (S, P, S%var(i), data(:P%ndata), skip_edit)
		return
	endif
enddo
write (name,'(i0)') field
call rads_error (S, rads_err_var, 'No variable with field number "'//name//'" was defined for "'//trim(S%branch(1))//'"')
data(:P%ndata) = nan
end subroutine rads_get_var_by_number

!****if* rads/rads_get_var_by_var
! SUMMARY
! Read variable (data) from RADS by type(rads_var)
!
! SYNOPSIS
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
!****-------------------------------------------------------------------
logical :: skip_edit

S%error = rads_noerr
if (rads_get_var_helper (S, P, data)) return ! Check data sizes

! Do we need to skip editing?
if (var%noedit) then
	skip_edit = .true.
else if (present(noedit)) then
	skip_edit = noedit
else
	skip_edit = .false.
endif

call rads_get_var_common (S, P, var, data(:P%ndata), skip_edit)
end subroutine rads_get_var_by_var

!****if* rads/rads_get_var_by_name
! SUMMARY
! Read variable (data) from RADS by character
!
! SYNOPSIS
recursive subroutine rads_get_var_by_name (S, P, varname, data, noedit)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
real(eightbytereal), intent(out) :: data(:)
logical, intent(in), optional :: noedit
!
! PURPOSE
! This routine loads the data from a single variable <varname>, addressed
! by a character string into the buffer <data>.
!****-------------------------------------------------------------------
type(rads_var), pointer :: var
integer(fourbyteint) :: l
logical :: skip_edit

S%error = rads_noerr
if (rads_get_var_helper (S, P, data)) return ! Check data sizes

! If varname ends with %, suspend editing, otherwise follow noedit, or default = .false.
l = len_trim(varname)
if (varname(l:l) == rads_noedit) then
	l = l - 1
	skip_edit = .true.
else if (present(noedit)) then
	skip_edit = noedit
else
	skip_edit = .false.
endif

var => rads_varptr (S, varname(:l))
if (.not.associated(var)) then
	data(:P%ndata) = nan
	return
endif
call rads_get_var_common (S, P, var, data(:P%ndata), skip_edit)
end subroutine rads_get_var_by_name

recursive subroutine rads_get_var_by_name_2d (S, P, varname, data, noedit)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
real(eightbytereal), intent(out) :: data(:,:)
logical, intent(in), optional :: noedit
real(eightbytereal), allocatable :: temp(:)
integer(fourbyteint) :: k
k = P%ndata
P%ndata = P%ndata*P%n_hz
allocate (temp(P%ndata))
call rads_get_var_by_name (S, P, varname, temp, noedit)
data(1:P%n_hz,1:k) = reshape(temp,(/P%n_hz,k/))
deallocate(temp)
P%ndata = k
end subroutine rads_get_var_by_name_2d

!****if* rads/rads_get_var_helper
! SUMMARY
! Helper routine for all rads_get_var_by_* routines
!
! SYNOPSIS
logical function rads_get_var_helper (S, P, data)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(in) :: P
real(eightbytereal), intent(in) :: data(:)
!
! PURPOSE
! This routine checks two things:
! - If P%ndata <= 0, then return .true.
! - If allocated memory is too small, then return .true. and print error
! - Else return .false.
!****-------------------------------------------------------------------
if (P%ndata <= 0) then
	rads_get_var_helper = .true.
else if (size(data) < P%ndata) then
	call rads_error (S, rads_err_memory, 'Too little memory allocated to read data from file', P)
	rads_get_var_helper = .true.
else
	rads_get_var_helper = .false.
endif
end function rads_get_var_helper

!****if* rads/rads_get_var_common
! SUMMARY
! Read variable (data) from RADS database (common to all)
!
! SYNOPSIS
recursive subroutine rads_get_var_common (S, P, var, data, noedit)
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(inout) :: var
real(eightbytereal), intent(out) :: data(:)
logical, intent(in) :: noedit
!
! PURPOSE
! This routine is common to all rads_get_var_* routines. It should only
! be called from those routines.
!****-------------------------------------------------------------------
type(rads_varinfo), pointer :: info
integer :: i

if (rads_verbose >= 4) write (*,'(2a)') 'rads_get_var_common: ',trim(var%name)

! Check size of array
if (size(data) < P%ndata) then
	call rads_error (S, rads_err_memory, 'Too little memory allocated to read data from file', P)
	return
endif

! Set pointer to info struct
info => var%info

do i = 1,3 ! This loop is here to allow processing of aliases
	S%error = rads_noerr

	! Load data depending on type of data source
	select case (info%datasrc)
	case (rads_src_nc_var)
		call rads_get_var_nc
	case (rads_src_nc_att)
		call rads_get_var_nc_att
	case (rads_src_math)
		call rads_get_var_math
	case (rads_src_grid_lininter, rads_src_grid_splinter, rads_src_grid_query, rads_src_grid_linphase)
		call rads_get_var_grid
	case (rads_src_constant)
		call rads_get_var_constant
	case (rads_src_flags)
		call rads_get_var_flags
	case (rads_src_tpj)
	    call rads_get_var_tpj
	case default
		data = nan
		return
	end select

	! Error handling
	if (S%error == rads_noerr) then ! No error, go to editing stage
		exit
	else if (info%default /= huge(0d0)) then ! Default value
		data = info%default
		S%error = rads_noerr
		exit

	! Look for alternative variables and run loop again
	else if (i == 1 .and. associated(var%inf1)) then
		info => var%inf1
	else if (i == 2 .and. associated(var%inf2)) then
		info => var%inf2

	! Ran out of options
	else
		call rads_error (S, rads_err_var, 'Could not find any data for variable "'//trim(var%name)//'" in file', P)
		data = nan
		exit
	endif

enddo ! End alias loop

! Edit this data if required and no error occurred
if (.not.noedit .and. S%error == rads_noerr) then
	if (any(info%limits == info%limits)) call rads_check_limits
	if (info%quality_flag(:1) /= ' ') call rads_check_quality_flag
endif
call rads_update_stat	! Update statistics on variable

contains

include "rads_tpj.f90"

recursive subroutine rads_get_var_nc ! Get data variable from RADS NetCDF file
use netcdf
use rads_netcdf
integer(fourbyteint) :: start(2), count(2), e, i, nf_get_vara_double, ncid
real(eightbytereal) :: x

! If time, lat, lon are already read, return those arrays upon request
if (P%time_dims /= info%ndims) then
	! Read from NetCDF file
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

! Look for the variable name in the NetCDF file (or take the stored one)
ncid = P%fileinfo(info%brid)%ncid
if (P%cycle == info%cycle .and. P%pass == info%pass) then
	! Keep old varid, but produce error when already tried and failed
	if (info%varid == 0) then
		S%error = rads_err_nc_var
		return
	endif
else if (nff(nf90_inq_varid(ncid, info%dataname, info%varid))) then
	! Read variable attributes if not yet set, or if we read/write
	if (P%rw .or. info%nctype == 0) e = nf90_inquire_variable (ncid, info%varid, xtype=info%nctype)
	if (P%rw .or. info%long_name(:1) == ' ') e = nf90_get_att(ncid, info%varid, 'long_name', info%long_name)
	if (P%rw .or. info%units(:1) == ' ') e = nf90_get_att(ncid, info%varid, 'units', info%units)
	if (P%rw .or. info%standard_name(:1) == ' ') e = nf90_get_att(ncid, info%varid, 'standard_name', info%standard_name)
	if (P%rw .or. info%comment(:1) == ' ') e = nf90_get_att(ncid, info%varid, 'comment', info%comment)
else
	! Failed to find variable
	S%error = rads_err_nc_var
	info%varid = 0
	return
endif
e = nf90_inquire_variable (ncid, info%varid, ndims=info%ndims)

! Load the data
start(1) = max(1,abs(P%first_meas))
if (info%ndims == 0) then
	! Constant to be converted to 1-dimensional array
	if (nft(nf90_get_var(ncid, info%varid, data(1)))) then
		call rads_error (S, rads_err_nc_get, 'Error reading NetCDF constant "'//trim(info%dataname)//'" in file', P)
		return
	endif
	data = data(1)
	info%ndims = 1
else if (info%ndims == 1 .and. S%n_hz_output .and. P%n_hz > 0 .and. P%first_meas > 0) then
	! 1-dimensional array with duplicated 1-Hz values
	if (nft(nf90_get_var(ncid, info%varid, data(1:P%ndata:P%n_hz), start))) then
		call rads_error (S, rads_err_nc_get, 'Error reading NetCDF array "'//trim(info%dataname)//'" in file', P)
		return
	endif
	forall (i = 1:P%ndata:P%n_hz) data(i+1:i+P%n_hz-1) = data(i)
else if (info%ndims == 1) then
	! 1-dimensional array of 1-Hz values
	if (nft(nf90_get_var(ncid, info%varid, data, start))) then
		call rads_error (S, rads_err_nc_get, 'Error reading NetCDF array "'//trim(info%dataname)//'" in file', P)
		return
	endif
else if (info%ndims == 2) then
	! 2-dimensional array, read as single column
	start(2) = start(1)
	start(1) = 1
	count(1) = P%n_hz
	count(2) = P%ndata / P%n_hz
	! We use the Fortran 77 routine here so that we can easily read a 2D field into a 1D array
	if (nft(nf_get_vara_double(ncid, info%varid, start, count, data))) then
		call rads_error (S, rads_err_nc_get, 'Error reading NetCDF array "'//trim(info%dataname)//'" in file', P)
		return
	endif
else
	call rads_error (S, rads_err_nc_get, 'Too many dimensions for variable "'//trim(info%dataname)//'" in file', P)
	return
endif

! Set NaN values and apply optional scale_factor and add_offset
! If we read/write, we also store the scale factor and add_offset
if (nff(nf90_get_att(ncid, info%varid, '_FillValue', x))) where (data == x) data = nan

if (nff(nf90_get_att(ncid, info%varid, 'scale_factor', x))) then
	data = data * x
else
	x = 1d0
endif
if (P%rw) info%scale_factor = x
if (nff(nf90_get_att(ncid, info%varid, 'add_offset', x))) then
	data = data + x
else
	x = 0d0
endif
if (P%rw) info%add_offset = x
end subroutine rads_get_var_nc

subroutine rads_get_var_nc_att ! Get data attribute from RADS NetCDF file
use netcdf
use rads_netcdf
use rads_time
use rads_misc
integer(fourbyteint) :: varid, i, ncid
character(len=26) :: date

ncid = P%fileinfo(info%brid)%ncid

! First locate the colon in the name
i = index(info%dataname, ':')

! If name starts with colon, then we have a global attribute, else a variable attribute
if (i == 1) then
	varid = nf90_global
else if (nft(nf90_inq_varid(ncid, info%dataname(:i-1), varid))) then
	S%error = rads_err_nc_var
	return
endif
if (nft(nf90_inquire_attribute (ncid, varid, info%dataname(i+1:), xtype=info%nctype))) info%nctype = 0

! Read the attribute
if (info%nctype == nf90_char) then
	! This is likely a date string
	if (nft(nf90_get_att(ncid, varid, info%dataname(i+1:), date))) then
		call rads_error (S, rads_err_nc_get, 'Error reading NetCDF attribute "'//trim(info%dataname)//'" in file', P)
		return
	endif
	data = strp1985f (date)
else
	! Load an integer or float value
	if (nft(nf90_get_att(ncid, varid, info%dataname(i+1:), data(1)))) then
		call rads_error (S, rads_err_nc_get, 'Error reading NetCDF attribute "'//trim(info%dataname)//'" in file', P)
		return
	endif
	data = data(1)
endif
end subroutine rads_get_var_nc_att

subroutine rads_get_var_flags ! Get value from flag word or vice versa
use netcdf
use rads_netcdf
integer(fourbyteint) :: start(1), i, j, k, bits(2), ncid

ncid = P%fileinfo(1)%ncid
if (info%dataname /= 'flags') then
	! Extract single flags from flagword
	if (.not.associated(P%flags)) then
		! Flags need to be loaded first
		if (nft(nf90_inq_varid (ncid, 'flags', info%varid))) then
			! Failed to find variable
			S%error = rads_err_nc_var
			return
		endif
		start = max(1,P%first_meas)
		allocate (P%flags(P%ndata))
		if (nft(nf90_get_var(ncid, info%varid, P%flags, start))) then
			call rads_error (S, rads_err_nc_get, 'Error reading NetCDF array "flags" in file', P)
			return
		endif
	endif

	if (info%dataname == 'surface_type') then
		! Special setting to get surface type from bits 2, 4, 5
		! 0=ocean, 2=enclosed seas and lakes, 3=land, 4=continental ice
		! I am keeping 1 for coastal in future
		do i = 1,P%ndata
			if (btest(P%flags(i),2)) then
				data(i) = 4 ! continental ice
			else if (btest(P%flags(i),4)) then
				data(i) = 3 ! land
			else if (btest(P%flags(i),5)) then
				data(i) = 2 ! enclosed sea or lake
			else
				data(i) = 0 ! ocean
			endif
		enddo
	else
		! Just get one or more bits
		bits = (/0,1/)
		read (info%dataname, *, iostat=i) bits
		data = ibits(P%flags,bits(1),bits(2))
	endif

else if (associated(P%flags)) then
	! We already have the flagswords, so just copy
	data = P%flags

else
	! Reconstrunct flagword from various flags
	allocate (P%flags(P%ndata))
	P%flags = 0
	do j = 1,S%nvar
		k = S%var(j)%field(1)
		if (k < 2501 .or. k > 2516) cycle
		k = modulo (k - 2500, 16)
		call rads_get_var_by_var (S, P, S%var(j), data, noedit=.true.)
		if (k == 4) then ! Surface type
			do i = 1,P%ndata
				select case (nint2(data(i)))
				case (4) ! ice = land = non-ocean
					P%flags(i) = P%flags(i) + flag_masks(2) + flag_masks(4) + flag_masks(5)
				case (3) ! land = non-ocean
					P%flags(i) = P%flags(i) + flag_masks(4) + flag_masks(5)
				case (2) ! enclosed seas or lakes = non-ocean
					P%flags(i) = P%flags(i) + flag_masks(5)
				end select
			enddo
		else
			P%flags = P%flags + nint2(data) * flag_masks(k)
		endif
	enddo
	data = P%flags
endif

end subroutine rads_get_var_flags

recursive subroutine rads_get_var_math ! Get data variable from MATH statement
use rads_misc
use rads_math
type(math_ll), pointer :: top
integer(fourbyteint) :: i, i0, i1, istat

! Start with a nullified 'top'
nullify(top)

! Process the math commands left to right
i1 = 0
do
	if (.not.next_word (info%dataname, i0, i1)) exit
	if (i1 == i0) cycle
	istat = math_eval (info%dataname(i0:i1-1), P%ndata, top)
	if (istat /= 0) then  ! No command or value, likely to be a variable
		call math_push (P%ndata,top)
		call rads_get_var_by_name (S, P, info%dataname(i0:i1-1), top%data)
	endif
	if (S%error /= rads_noerr) exit
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
	write (stderr,*) i,' remaining items on stack'
	call rads_error (S, rads_noerr, 'Cleaned up')
endif
end subroutine rads_get_var_math

subroutine rads_get_var_grid ! Get data by interpolating a grid
use rads_grid
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
	S%error = rads_err_source
else if (info%datasrc == rads_src_grid_lininter) then
	forall (i = 1:P%ndata) data(i) = grid_lininter (info%grid, x(i), y(i))
else if (info%datasrc == rads_src_grid_splinter) then
	forall (i = 1:P%ndata) data(i) = grid_splinter (info%grid, x(i), y(i))
else if (info%datasrc == rads_src_grid_linphase) then
	forall (i = 1:P%ndata) data(i) = grid_lininter (info%grid, x(i), y(i), .true.)
else
	forall (i = 1:P%ndata) data(i) = grid_query (info%grid, x(i), y(i))
endif
end subroutine rads_get_var_grid

subroutine rads_get_var_constant
integer :: ios
read (info%dataname, *, iostat=ios) data(1)
if (ios == 0) then
	S%error = rads_noerr
	data = data(1)
else
	S%error = rads_err_var
	data = nan
endif
end subroutine rads_get_var_constant

subroutine rads_check_limits ! Edit the data based on limits
integer(fourbyteint) :: mask(2)

if (info%datatype == rads_type_lon) then
	! When longitude, shift the data first above the lower (west) limit
	data = info%limits(1) + modulo (data - info%limits(1), 360d0)
	! Then only check the upper (east) limit
	where (data > info%limits(2)) data = nan
else if (info%datatype == rads_type_flagmasks) then
	! Do editing of flag based on mask
	mask = nint(info%limits)
	where (info%limits /= info%limits) mask = 0 ! Because it is not guaranteed for every compiler that nint(nan) = 0
	! Reject when any of the set bits of the lower limit are set in the data
	if (mask(1) /= 0) where (iand(nint(data),mask(1)) /= 0) data = nan
	! Reject when any of the set bits of the upper limit are not set in the data
	if (mask(2) /= 0) where (iand(nint(data),mask(2)) /= mask(2)) data = nan
else
	! Do editing of any other variable
	where (data < info%limits(1) .or. data > info%limits(2)) data = nan
endif
end subroutine rads_check_limits

recursive subroutine rads_check_quality_flag
use rads_misc
real(eightbytereal) :: qual(P%ndata)
integer(fourbyteint) :: i0, i1

! Process all quality checks left to right
i1 = 0
do
	if (.not.next_word (info%quality_flag, i0, i1)) exit
	if (i1 == i0) cycle
	call rads_get_var_by_name (S, P, info%quality_flag(i0:i1-1), qual)
	if (S%error /= rads_noerr) exit
	where (isnan_(qual)) data = nan
enddo
end subroutine rads_check_quality_flag

subroutine rads_update_stat ! Update the statistics for given var
use rads_misc
real(eightbytereal) :: q, r
integer(fourbyteint) :: i

! If stats are already done for this cycle and pass, skip it
if (P%cycle == info%cycle .and. P%pass == info%pass) return
info%cycle = P%cycle
info%pass = P%pass

! Count selected and rejected measurements and update pass statistics
! Use method of West (1979)
do i = 1,P%ndata
	if (isnan_(data(i))) then
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

!****f* rads/rads_read_xml
! SUMMARY
! Read RADS4 XML file
!
! SYNOPSIS
subroutine rads_read_xml (S, filename)
use netcdf
use xmlparse
use rads_time
use rads_misc
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: filename
!
! PURPOSE
! This routine parses a RADS4 XML file and fills the <S> struct with
! information pertaining to the given satellite and all variable info
! encountered in that file.
!
! The execution terminates on any error, and also on any warning if
! fatal = .true.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! filename : XML file name
! fatal    : If .true., then all warnings are fatal.
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_xml_parse, rads_err_xml_file
!****-------------------------------------------------------------------
type(xml_parse) :: X
integer, parameter :: max_lvl = 20
character(len=rads_varl) :: tag, name, tags(max_lvl)
character(len=rads_naml) :: attr(2,max_lvl), val(max_lvl)
character(len=6) :: src
integer :: nattr, nval, i, j, ios, skip, skip_level
integer(twobyteint) :: field(2)
logical :: endtag, endskip
real(eightbytereal) :: node_rate
type(rads_varinfo), pointer :: info, info_block
type(rads_var), pointer :: var, var_block
type(rads_phase), pointer :: phase

! Initialise
S%error = rads_noerr
endskip = .true.
skip_level = 0
nullify (var_block, info_block, phase)

! Open XML file
call xml_open (X, filename, .true.)
! If failed, try with ".xml" extension (if not already)
if (X%error .and. index(filename,'.xml') == 0) call xml_open (X, trim(filename) // '.xml', .true.)
if (X%error) then
	S%error = rads_err_xml_file
	return
endif
if (rads_verbose >= 2) write (*,'(2a)') 'Parsing XML file ',trim(filename)
call xml_options (X, ignore_whitespace = .true.)

! Parse XML file, store information in S struct
do
	if (X%eof) exit
	call xml_get (X, tag, endtag, attr, nattr, val, nval)

	! Process closing tags
	if (endtag) then
		if (tag /= tags(X%level+1)) &
			call xmlparse_error ('Closing tag </'//trim(tag)//'> follows opening tag <'//trim(tags(X%level+1))//'>')
		endskip = (X%level < skip_level)
		if (endskip) skip_level = 0  ! Stop skipping when descended back below the starting level
		if (tag == 'var') nullify (var_block, info_block) ! Stop processing <var> block
		if (tag == 'phase') nullify (phase)   ! Stop processing <phase> block
		cycle  ! Ignore all other end tags
	endif

	! Set var and info to the ones for the <var> block, which means that outside of the block they will be set to null
	var => var_block
	info => info_block

	! Special actions for <else> and <elseif>
	! These will issue a 'skip' when previous <if> was not skipped
	if (tag == 'else' .or. tag == 'elseif') then
		if (tags(X%level) /= 'if' .and. tags(X%level) /= 'elseif') &
			call xmlparse_error ('Opening tag <'//trim(tag)//'> follows closing tag </'//trim(tags(X%level))//'>')
		if (.not.endskip .and. skip_level == 0) skip_level = X%level
	endif

	! Process opening tags
	tags(X%level) = tag
	if (skip_level > 0 .and. X%level >= skip_level) cycle	! Skip all data at level equal or larger than skip level

	! Error maybe?
	if (X%error) then
		call xmlparse_error ('Error parsing xml line')
		S%error = rads_err_xml_parse
		exit
	endif

	! Check if we need to skip this level
	! This level will be skipped when all three of the following are true
	! a) Tag contains attribute "sat="
	! b) The attribute value contains the satellite abbreviaton, or
	!    the attribute value starts with "!" and does not contain the satellite abbreviation
	! c) "*.r" matches all branches with extension ".r")
	! d) "6a.hr*" matches all branches that start with "6a.hr"
	! e) The satellite abbreviation is not set to "??"
	!
	! Example 1: for original TOPEX (tx)
	! sat="tx" => pass
	! sat="tx.r" => skip
	! sat="*.r" => skip
	! sat="j1" => skip
	! sat="j1 tx" => pass
	! sat="!j1" => pass
	! sat="!j1 tx" => skip
	!
	! Example 2: for TOPEX Retracked (tx.r)
	! sat="tx" => pass
	! sat="tx.r" => pass
	! sat="*.r" => pass
	! sat="!j1" => pass
	! sat="!tx" => skip
	! sat="!tx.r" => skip
	!
	! Additionally: check for var="name" option. This will temporarily overrule var and info.
	! They are reset to var_block and info_block on the next cycle of the loop.

	skip = 0
	skip_level = 0
	do i = 1,nattr
		select case (attr(1,i))
		case ('sat')
			if (skip == 0) skip = 1
			if (S%sat == '??') then
				skip = -1
			else if ((attr(2,i)(:1) == '!') .eqv. &
				(index(attr(2,i),S%sat//' ') == 0 .and. index(attr(2,i),trim(S%branch(1))//' ') == 0 &
					.and. index(attr(2,i),S%branch(1)(:5)//'* ') == 0 &
					.and. index(attr(2,i),'*'//trim(S%branch(1)(3:))//' ') == 0)) then
				skip=-1
			endif
		case ('var')
			var => rads_varptr (S, attr(2,i), null())
			info => var%info
		end select
	enddo
	if (skip == 1 .and. skip_level == 0) then
		skip_level = X%level
		cycle
	endif

	select case (tag)
	case ('global_attributes')
		allocate (S%glob_att(nval))
		S%glob_att = val(1:nval)

	case ('satellites')
		call sat_translate

	case ('satellite')
		S%satellite = val(1)(:8)

	case ('satid')
		read (val(:nval), *, iostat=ios) S%satid

	case ('phase')
		if (has_name()) then
			phase => rads_init_phase (S, name)
			phase%name = attr(2,1)(:rads_varl)
		endif

	case ('mission')
		phase%mission = val(1)(:rads_varl)

	case ('cycles')
		read (val(:nval), *, iostat=ios) phase%cycles

	case ('exclude_cycles')
		allocate (S%excl_cycles)
		S%excl_cycles%list = -1
		read (val(:nval), *, iostat=ios) S%excl_cycles%list
		S%excl_cycles%n = count(S%excl_cycles%list /= -1)

	case ('subcycles')
		allocate (phase%subcycles)
		phase%subcycles%list = 0
		phase%subcycles%i = 1
		do i = 1,nattr
			if (attr(1,i) == 'start') read (attr(2,i), *, iostat=ios) phase%subcycles%i
		enddo
		read (val(:nval), *, iostat=ios) phase%subcycles%list
		! Turn length of subcycles in passes into number of accumulative passes before subcycle
		phase%subcycles%n = count(phase%subcycles%list /= 0)
		phase%passes = maxval(phase%subcycles%list)
		do i = phase%subcycles%n, 2, -1
			phase%subcycles%list(i) = sum(phase%subcycles%list(1:i-1))
		enddo
		phase%subcycles%list(1) = 0

	case ('ref_pass') ! First part is a date string, at least 19 characters long
		i = index(val(1),' ') ! Position of first space
		phase%ref_time = strp1985f(val(1)(:i-1))
		phase%ref_pass = 1
		phase%ref_orbit = 1
		read (val(1)(i:), *, iostat=ios) phase%ref_lon, phase%ref_cycle, phase%ref_pass, phase%ref_orbit

	case ('start_time')
		phase%start_time = strp1985f(val(1))

	case ('end_time')
		phase%end_time = strp1985f(val(1))

	case ('repeat')
		phase%repeat_shift = 0d0
		read (val(:nval), *, iostat=ios) phase%repeat_days, phase%repeat_passes, phase%repeat_shift
		! Compute length of repeat in nodal days from inclination and repeat in solar days
		! This assumes 1000 km altitude to get an approximate node rate (in rad/s)
		node_rate = -1.21306d-6 * cos(S%inclination*rad)
		phase%repeat_nodal = nint(phase%repeat_days * (7.292115d-5 - node_rate) / 7.272205d-5)
		! Determine length in seconds of single pass
		phase%pass_seconds = phase%repeat_days * 86400d0 / phase%repeat_passes
		phase%passes = phase%repeat_passes

	case ('dt1hz')
		read (val(:nval), *, iostat=ios) S%dt1hz

	case ('inclination')
		read (val(:nval), *, iostat=ios) S%inclination

	case ('frequency')
		read (val(:nval), *, iostat=ios) S%frequency

	case ('xover_params')
		read (val(:nval), *, iostat=ios) S%xover_params

	case ('alias')
		if (associated(var)) then	! Within <var> block, we do not need "name" attribute
			call rads_set_alias (S, var%name, var%name // ' ' // val(1), field)
		else if (has_name (field)) then
			if (val(1)(1:1) == ',' .or. val(1)(1:1) == '/') then
				call rads_set_alias (S, name, trim(name) // val(1), field)
			else
				call rads_set_alias (S, name, val(1), field)
			endif
		endif

	case ('var')
		if (has_name (field)) then
			var_block => rads_varptr (S, name, null())
			info_block => var_block%info
			if (any(field > rads_nofield)) var_block%field = field
		endif

	case ('long_name')
		call assign_or_append (info%long_name)

	case ('standard_name')
		info%standard_name = val(1)(:rads_naml)
		select case (val(1))
		case ('time')
			info%datatype = rads_type_time
		case ('latitude')
			info%datatype = rads_type_lat
		case ('longitude')
			info%datatype = rads_type_lon
		case ('dimension')
			info%datatype = rads_type_dim
			info%standard_name = '' ! Do not use this as standard name
		case ('sea_surface_height_above_sea_level')
			info%datatype = rads_type_sla
		case default
			info%datatype = rads_type_other
		end select

	case ('source')
		info%source = val(1)(:rads_naml)

	case ('parameters')
		info%parameters = val(1)(:rads_naml)

	case ('units')
		info%units = val(1)(:rads_varl)

	case ('flag_masks')
		call assign_or_append (info%flag_meanings)
		info%datatype = rads_type_flagmasks

	case ('flag_values')
		call assign_or_append (info%flag_meanings)
		info%datatype = rads_type_flagvalues

	case ('comment')
		call assign_or_append (info%comment)

	case ('backup')
		call rads_message (trim(var%name)//': <backup> deprecated, use <default> for value, or <alias> for variable')

	case ('default')
		info%default = huge(0d0) ! Set to default first
		read (val(:nval), *, iostat=ios) info%default

	case ('quality_flag')
		call assign_or_append (info%quality_flag)

	case ('format')
		call rads_set_format (S, var%name, val(1))

	case ('compress')
		call rads_set_compress (S, var%name, val(1))

	case ('limits') ! Do not use routine rads_set_limits here!
		if (all(val == '')) info%limits = nan ! Reset limits to none if none given
		if (var%name == 'time') then
			i = index(val(1),' ') ! Position of first space
			info%limits(1) = strp1985f(val(1)(:i-1))
			info%limits(2) = strp1985f(val(1)(i+1:))
		else
			read (val(:nval), *, iostat=ios) info%limits
			! If we have an old-fashioned flagword, convert it to limits of single flags
			if (var%name == 'flags') call rads_set_limits_by_flagmask (S, info%limits)
		endif

	case ('plot_range')
		read (val(:nval), *, iostat=ios) info%plot_range

	case ('data')
		call assign_or_append (info%dataname)
		src = ''
		do i = 1, nattr
			select case (attr(1,i))
			case ('source')
				src = attr(2,i)(:6)
			case ('branch')
				! 'branch' can be like '.20hz'
				do j = 1, rads_max_branches
					if (S%branch(j) == S%sat//attr(2,i)) then
						exit
					else if (S%branch(j) == '') then
						S%branch(j) = S%sat//attr(2,i)(:rads_varl-2)
						exit
					endif
				enddo
				if (j > rads_max_branches) call rads_exit ('Too many branches: '//trim(attr(2,i)))
				info%brid = j
			end select
		enddo
		! Work out which data source
		select case (src)
		case ('constant')
			info%datasrc = rads_src_constant
		case ('grid', 'grid_l')
			info%datasrc = rads_src_grid_lininter
		case ('grid_s', 'grid_c')
			info%datasrc = rads_src_grid_splinter
		case ('grid_p')
			info%datasrc = rads_src_grid_linphase
		case ('grid_n', 'grid_q')
			info%datasrc = rads_src_grid_query
		case ('math')
			info%datasrc = rads_src_math
		case ('netcdf', 'nc_var', 'nc_att', 'nc')
			info%datasrc = rads_src_nc_var
			if (index(info%dataname,':') > 0) info%datasrc = rads_src_nc_att
		case ('flags')
			info%datasrc = rads_src_flags
		case ('tpj')
			info%datasrc = rads_src_tpj
		case default
			! Make "educated guess" of data source
			if (index(info%dataname,'.nc') > 0) then
				info%datasrc = rads_src_grid_lininter
			else if (index(trim(info%dataname),' ') > 0) then
				info%datasrc = rads_src_math
			else if (index(info%dataname,':') > 0) then
				info%datasrc = rads_src_nc_att
			else if (is_number(info%dataname)) then
				info%datasrc = rads_src_constant
			else
				info%datasrc = rads_src_nc_var
			endif
		end select
		! Additional stuff to do for grids
		if (info%datasrc / 10 * 10 == rads_src_grid_lininter) then
			info%gridx = 'lon'
			info%gridy = 'lat'
			do i = 1,nattr
				select case (attr(1,i))
				case ('x')
					info%gridx = attr(2,i)(:rads_varl)
				case ('y')
					info%gridy = attr(2,i)(:rads_varl)
				end select
			enddo
			allocate (info%grid)
			call parseenv (val(1), info%grid%filenm)
			info%grid%ntype = 0	! This signals that the grid was not loaded yet
		endif

	case ('dimensions')
		read (val(:nval), *, iostat=ios) info%ndims

	case ('if', 'elseif', 'else', '!--')
		! Dummy and comment tags

	! FOR BACKWARD COMPATIBILITY, here are tags <netcdf>, <math> and <grid>
	case ('netcdf')
		info%datasrc = rads_src_nc_var
		info%dataname = val(1)
		if (index(info%dataname, ':') > 0) info%datasrc = rads_src_nc_att

	case ('math')
		call assign_or_append (info%dataname)
		info%datasrc = rads_src_math

	case ('grid')
		info%datasrc = rads_src_grid_lininter
		info%gridx = 'lon' ; info%gridy = 'lat'
		do i = 1,nattr
			select case (attr(1,i))
			case ('x')
				info%gridx = attr(2,i)(:rads_varl)
			case ('y')
				info%gridy = attr(2,i)(:rads_varl)
			case ('inter')
				if (attr(2,i)(:1) == 'c' .or. attr(2,i)(:1) == 's') info%datasrc = rads_src_grid_splinter
				if (attr(2,i)(:1) == 'p') info%datasrc = rads_src_grid_linphase
				if (attr(2,i)(:1) == 'q') info%datasrc = rads_src_grid_query
			end select
		enddo
		allocate (info%grid)
		call parseenv (val(1), info%grid%filenm)
		info%grid%ntype = 0	! This signals that the grid was not loaded yet

	case default
		call xmlparse_error ('Unknown tag <'//trim(tag)//'>')

	end select
enddo

! Close XML file
call xml_close (X)

if (X%level > 0) call xmlparse_error ('Did not close tag <'//trim(tags(X%level))//'>')

if (S%error > rads_noerr) call rads_exit ('Fatal errors occurred while parsing XML file '//trim(filename))

contains

logical function has_name (field)
integer(twobyteint), optional :: field(2)
integer :: i, ios
name = ''
if (present(field)) field = rads_nofield
do i = 1,nattr
	if (attr(1,i) == 'name') then
		name = attr(2,i)(:rads_varl)
	else if (present(field) .and. attr(1,i) == 'field') then
		read (attr(2,i), *, iostat=ios) field
	endif
enddo
has_name = (name /= '')
if (.not.has_name) call xmlparse_error ('Tag <'//trim(tag)//'> requires ''name'' attribute')
end function has_name

subroutine assign_or_append (string)
! Assign or append concatenation of val(1:nval) to string
character(len=*), intent(inout) :: string
integer :: i, j, l
character(len=8) :: action

! Check for "action="
action = 'replace'
do i = 1,nattr
	if (attr(1,i) == 'action') action = attr(2,i)(:8)
enddo

select case (action)
case ('replace') ! Remove current content, then append new content
	string = ''
case ('append') ! Append new content
case ('delete') ! Delete from current content that what matches new content
	do i = 1,nval
		l = len_trim(val(i))
		j = index(string, val(i)(:l))
		if (j == 0) cycle
		if (string(j+l:j+l) == ' ') l = l + 1 ! Remove additional space
		string(j:) = string(j+l:)
	enddo
	return
case ('merge') ! Append new content only if it is not yet in current content
	do i = 1,nval
		l = len_trim(val(i))
		if (index(string, val(i)(:l)) > 0) cycle
		if (string == '') then
			string = val(i)
		else
			string = trim(string) // ' ' // val(i)
		endif
	enddo
	return
case default
	call xmlparse_error ('Unknown option "action='//trim(action)//'"')
end select

! Append strings val(:)
if (nval == 0) return
j = 1
if (string == '') then
	string = val(1)
	j = 2
endif
do i = j,nval
	string = trim(string) // ' ' // val(i)
enddo
end subroutine assign_or_append

subroutine xmlparse_error (string)
! Issue error message with file name and line number
character(len=*), intent(in) :: string
character(rads_naml) :: text
write (text, 1300) trim(filename), X%lineno, string
call rads_error (S, rads_err_xml_parse, text)
1300 format ('Error parsing file ',a,' at or near line ',i0,': ',a)
end subroutine xmlparse_error

subroutine sat_translate
integer :: i, j, l
! Translate the given satellite identifier into sat, tree, phase.
! Start with S%spec given on command line.
! It will be replaced by the mission phase, if any.
!
! If three characters, this may be like "e2g".
! Check if the first two characters match the list.
l = len_trim(S%spec)
if (l == 3) then
	do i = 1,nval
		if (S%spec(1:2) /= val(i)(1:2)) cycle
		S%sat = val(i)(1:2)
		S%spec = S%spec(3:3)	! Phase part
		return
	enddo
endif
! If we have a '/' or ':' or '.', then separate specification
j = scan(S%spec,'/:.')
if (j > 0) l = j - 1
! Now scan for matching strings (beginning of string only)
do i = 1,nval
	if (index(' '//val(i), ' '//strtolower(S%spec(:l))) == 0) cycle
	S%sat = val(i)(1:2)
	S%spec = S%spec(l+1:)	! Everything after <sat>
	j = scan(S%spec,'/:')
	if (j == 0) then	! No phase indication
		S%branch(1) = S%sat//S%spec(:rads_varl-2)
		S%spec = ''
	else	! With phase indication
		S%branch(1) = S%sat//S%spec(:j-1)
		S%spec = S%spec(j+1:)
	endif
	return
enddo
call rads_exit ('No satellite found based on specification "'//trim(S%spec)//'"')
end subroutine sat_translate

end subroutine rads_read_xml

!****if* rads/rads_varptr
! SUMMARY
! Returns the pointer to a given variable
!
! SYNOPSIS
function rads_varptr (S, varname, tgt) result (ptr)
use netcdf
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: varname
type(rads_var), intent(in), pointer, optional :: tgt
type(rads_var), pointer :: ptr
!
! PURPOSE
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
! ARGUMENTS
! S        : Satellite/mission dependent structure
! varname  : Name of the RADS variable
! tgt      : (optional) target pointer (see above)
! rads_varptr : Pointer to the structure for <varname>
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_var
!****-------------------------------------------------------------------
integer(fourbyteint) :: i, n
type(rads_var), pointer :: temp(:)
integer(twobyteint) :: field
character(len=rads_varl), pointer :: name
character(len=rads_naml), pointer :: long_name

S%error = rads_noerr
field = -999

! If variable name is numerical, look for field number
if (varname(1:1) >= '0' .and. varname(1:1) <= '9') then
	read (varname, *, iostat=i) field
	do i = 1,S%nvar
		if (any(S%var(i)%field == field)) exit
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

if (i <= S%nvar) then
	! Match was found: return pointer
	ptr => S%var(i)
	if (.not.present(tgt)) then
		! No association requested
		return
	else if (.not.associated(tgt)) then
		! No association requested
		return
	else if (associated(ptr%info,tgt%info)) then
		! Association is already a fact
		return
	endif
	! Reassociate existing variable
	call rads_free_var_struct (S, ptr, .true.)
else if (.not.present(tgt)) then
	! No match found, and none should be created: return null pointer and error
	nullify (ptr)
	call rads_error (S, rads_err_var, 'No variable "'//trim(varname)//'" was defined for "'//trim(S%branch(1))//'"')
	return
else
	! If we got here, we need to make a new variable. Do we also need to allocate more space?
	n = size(S%var)
	if (i > n) then
		allocate (temp(n + rads_var_chunk))
		temp(1:n) = S%var
		temp(n+1:n+rads_var_chunk) = rads_var (null(), null(), null(), null(), null(), .false., rads_nofield)
		deallocate (S%var)
		S%var => temp
		if (rads_verbose >= 3) write (*,'(a,2i5)') 'Increased S%var:',n,n+rads_var_chunk
	endif
	S%nvar = i
	ptr => S%var(i)
endif

! Assign the info struct to that of tgt or make a new one
if (associated(tgt)) then
	ptr%info => tgt%info
	allocate (name)
else
	allocate (ptr%info)
	ptr%info = rads_varinfo (varname, varname, '', '', '', '', '', '', '', '', 'f0.3', '', '', null(), &
		huge(0d0), nan, nan, 0d0, 1d0, nan, nan, 0d0, 0d0, .false., 1, 1, nf90_double, 0, 0, 0, 0, 0, 0, 0)
	name => ptr%info%name
endif
long_name => ptr%info%long_name ! This is to avoid warning in gfortran 4.8

! Assign name and long_name
ptr%name => name
ptr%long_name => long_name
if (field > rads_nofield) then ! Was given field number
	write (ptr%name, '("f",i4.4)') field
	ptr%field = field
else
	ptr%name = varname
endif

end function rads_varptr

!****f* rads/rads_set_alias
! SUMMARY
! Set alias to an already defined variable
!
! SYNOPSIS
subroutine rads_set_alias (S, alias, varname, field)
use rads_misc
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: alias, varname
integer(twobyteint), intent(in), optional :: field(2)
!
! PURPOSE
! This routine defines an alias to an existing variable, or up to three
! variables. When more than one variable is given as target, they will
! be addressed one after the other.
! If alias is already defined as an alias or variable, it will be overruled.
! The alias will need to point to an already existing variable or alias.
! Up to three variables can be specified, separated by spaces or commas.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! alias    : New alias for (an) existing variable(s)
! varname  : Existing variable name(s)
! field    : (optional) new field numbers to associate with alias
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_alias, rads_err_var
!****-------------------------------------------------------------------
type(rads_var), pointer :: tgt, src
integer :: i0, i1, n_alias, l0, l1, l2, i
character(len=rads_naml), pointer :: long_name

S%error = rads_noerr

! Check that neither alias not varname is empty
if (alias == '') then
	call rads_error (S, rads_err_alias, 'Alias is empty')
	return
else if (varname == '') then
	call rads_error (S, rads_err_var, 'Variable name is empty')
	return
endif

! Now loop through the string <varname> to find variable names
i1 = 0
n_alias = 0
nullify (src)
do
	if (.not.next_word (varname, i0, i1)) exit
	n_alias = n_alias + 1
	if (i1 == i0) cycle
	tgt => rads_varptr (S, varname(i0:i1-1))
	if (.not.associated(tgt)) then
		call rads_error (S, rads_err_var, 'Alias target "'//varname(i0:i1-1)//'" of "'//trim(alias)//'" not found')
		return
	endif
	if (.not.associated(src)) src => rads_varptr (S, alias, tgt)
	select case (n_alias)
	case (1) ! First alias
		if (S%error == rads_err_alias) then
			call rads_error (S, rads_err_alias, 'Cannot alias "'//trim(alias)//'"; it is used by other variables')
			return
		endif
		if (present(field)) then
			if (any(field /= rads_nofield)) src%field = field
		endif
		nullify (src%inf1, src%inf2)
	case (2) ! Second alias
		src%inf1 => tgt%info
	case (3) ! Third alias
		src%inf2 => tgt%info
	case default
		call rads_error (S, rads_err_alias, 'Too many aliases for "'//trim(alias)//'"')
		return
	end select
enddo

! Reset long_name to a mix of long_names
if (n_alias < 2) return
l0 = len_trim(src%info%long_name)
l1 = len_trim(src%inf1%long_name)
i0 = 0
do i = 0,min(l0,l1)-1
	if (src%info%long_name(l0-i:l0-i) /= src%inf1%long_name(l1-i:l1-i)) exit
	if (src%info%long_name(l0-i:l0-i) == ' ') i0 = i + 1
enddo
allocate (long_name)
if (n_alias == 3) then
	l2 = len_trim(src%inf2%long_name)
	long_name = src%info%long_name(1:l0-i0) // ', ' // src%inf1%long_name(1:l1-i0) // &
		', ' // src%inf2%long_name(1:l2-i0) // src%info%long_name(l0-i0+1:l0)
else
	long_name = src%info%long_name(1:l0-i0) // ', ' // src%inf1%long_name(1:l1-i0) // &
		src%info%long_name(l0-i0+1:l0)
endif
src%long_name => long_name
end subroutine rads_set_alias

!****f* rads/rads_set_limits
! SUMMARY
! Set limits on given variable
!
! SYNOPSIS
subroutine rads_set_limits (S, varname, lo, hi, string, iostat)
use rads_misc
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: varname
real(eightbytereal), intent(in), optional :: lo, hi
character(len=*), intent(in), optional :: string
integer(fourbyteint), intent(out), optional :: iostat
!
! PURPOSE
! This routine set the lower and upper limits for a given variable in
! RADS.
! The limits can either be set by giving the lower and upper limits
! as double floats <lo> and <hi> or as a character string <string> which
! contains the two numbers separated by whitespace, a comma or a slash.
! In case only one number is given, only the lower or higher bound
! (following the separator) is set, the other value is left unchanged.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! varname  : Variable name
! lo, hi   : (optional) Lower and upper limit
! string   : (optional) String of up to two values, with separating whitespace
!            or comma or slash.
! iostat   : (optional) iostat code from reading string
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_var
!****-------------------------------------------------------------------
type(rads_var), pointer :: var
var => rads_varptr (S, varname)
if (.not.associated(var)) return
call rads_set_limits_info (var%info)
call rads_set_limits_info (var%inf1)
call rads_set_limits_info (var%inf2)

contains

subroutine rads_set_limits_info (info)
type(rads_varinfo), pointer :: info
real(eightbytereal), parameter :: sec1970 = -473385600d0
if (.not.associated(info)) return
if (present(lo)) info%limits(1) = lo
if (present(hi)) info%limits(2) = hi
if (.not.present(string)) then
else if (string == '') then	! Reset limits to none
	info%limits = nan
else ! Read limits (only change the ones given)
	call read_val (string, info%limits, '/', iostat=iostat)
endif
if (info%datatype == rads_type_lat .or. info%datatype == rads_type_lon) then
	! If latitude or longitude limits are changed, recompute equator longitude limits
	call rads_traxxing (S)
else if (info%datatype == rads_type_time) then
	! If time limits are changed, also limit the cycles
	if (isan_(info%limits(1))) S%cycles(1) = max(S%cycles(1), rads_time_to_cycle (S, info%limits(1)))
	! If no upper time limit is set, use the current time
	if (isan_(info%limits(2))) then
		S%cycles(2) = min(S%cycles(2), rads_time_to_cycle (S, info%limits(2)))
	else
		S%cycles(2) = min(S%cycles(2), rads_time_to_cycle (S, time()+sec1970))
	endif
else if (var%name == 'flags') then
	call rads_set_limits_by_flagmask (S, info%limits)
endif
end subroutine rads_set_limits_info
end subroutine rads_set_limits

!****if* rads/rads_set_limits_by_flagmask
! SUMMARY
! Set limits based on flagmask
!
! SYNOPSIS
subroutine rads_set_limits_by_flagmask (S, limits)
type(rads_sat), intent(inout) :: S
real(eightbytereal), intent(inout) :: limits(2)
!****-------------------------------------------------------------------
integer :: i, ios, mask(2), bits(2)
mask = 0
where (limits == limits) mask = nint(limits) ! Because it is not guaranteed for every compiler that nint(nan)=0
! Loop through all variables to find those with field between 2501 and 2516
do i = 1,S%nvar
	if (.not.any(S%var(i)%field >= 2501 .and. S%var(i)%field <= 2516)) cycle
	bits = (/0,1/) ! Default values
	read (S%var(i)%info%dataname, *, iostat=ios) bits
	if (S%var(i)%info%dataname == 'surface_type') then
		! Special setting to get surface type from bits 2, 4, 5
		! 0=ocean, 2=enclosed seas and lakes, 3=land, 4=continental ice
		! I am keeping 1 for coastal in future
		S%var(i)%info%limits = 0
		if (btest(mask(1),5)) then
			S%var(i)%info%limits(2) = 1 ! ocean only
		else if (btest(mask(1),4)) then
			S%var(i)%info%limits(2) = 2 ! water only
		else if (btest(mask(1),2)) then
			S%var(i)%info%limits(2) = 3 ! anything but ice
		endif
		if (btest(mask(2),2)) then
			S%var(i)%info%limits(1) = 4 ! ice only
		else if (btest(mask(2),4)) then
			S%var(i)%info%limits(1) = 3 ! land or ice
		else if (btest(mask(2),5)) then
			S%var(i)%info%limits(1) = 2 ! non-ocean
		endif
	else if (S%var(i)%info%datatype == rads_type_flagmasks) then
		S%var(i)%info%limits(1) = ibits(mask(1),bits(1),bits(2))
		S%var(i)%info%limits(2) = ibits(mask(2),bits(1),bits(2))
	else ! rads_type_flagvalues
		S%var(i)%info%limits(2) = ibits(not(mask(1)),bits(1),bits(2))
		S%var(i)%info%limits(1) = ibits(mask(2),bits(1),bits(2))
	endif
enddo
end subroutine rads_set_limits_by_flagmask

!****f* rads/rads_set_region
! SUMMARY
! Set latitude/longitude limits or distance to point
!
! SYNOPSIS
subroutine rads_set_region (S, string)
use rads_misc
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: string
!
! PURPOSE
! This routine set the region for data selection (after the -R option).
! The region can either be specified as a box by four values "W/E/S/N",
! or as a circular region by three values "E/N/radius". Separators
! can be commas, slashes, or whitespace.
!
! In case of a circular region, longitude and latitude limits are set
! accordingly for a rectangular box surrounding the circle. However, when
! reading pass data, the distance to the centroid is used as well to
! edit out data.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! string   : String of three or four values with separating whitespace.
!            For rectangular region: W/E/S/N.
!            For circular region: E/N/radius (radius in degrees).
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_var
!****-------------------------------------------------------------------
real(eightbytereal) :: r(4), x
r = nan
call read_val (string, r, '/')
if (isnan_(r(4))) then ! Circular region
	r(3) = abs(r(3))
	! For longitude limits, do some trigonometry to determine furthest meridians
	x = sin(r(3)*rad) / cos(r(2)*rad)
	if (r(3) < 90d0 .and. x >= 0d0 .and. x <= 1d0) then
		x = asin(x)/rad
	else ! Circle includes either pole
		x = 180d0
	endif
	S%lon%info%limits = r(1) + (/ -x, x /)
	! For latitude limits, add and subtract radius from centroid latitude
	S%lat%info%limits(1) = max(-90d0,r(2)-r(3))
	S%lat%info%limits(2) = min( 90d0,r(2)+r(3))
	! Convert longitude, latitude and radius from degrees to radians
	S%centroid = r(1:3)*rad
else ! Rectangular region
	S%lon%info%limits = r(1:2)
	S%lat%info%limits = r(3:4)
endif
call rads_traxxing (S)
end subroutine rads_set_region

!****f* rads/rads_set_format
! SUMMARY
! Set print format for ASCII output of given variable
!
! SYNOPSIS
subroutine rads_set_format (S, varname, format)
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: varname, format
!
! PURPOSE
! This routine set the FORTRAN format specifier of output of a given
! variable in RADS.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! varname  : Variable name
! format   : FORTRAN format specifier (e.g. 'f10.3')
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_var
!****-------------------------------------------------------------------
type(rads_var), pointer :: var
var => rads_varptr(S, varname)
if (.not.associated(var)) return
var%info%format = format
var%info%boz_format = (index('bozBOZ',var%info%format(:1)) > 0)
end subroutine rads_set_format

!****if* rads/rads_set_compress
! SUMMARY
! Set format for binary storage of given variable
!
! SYNOPSIS
subroutine rads_set_compress (S, varname, format)
use netcdf
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: varname, format
!
! PURPOSE
! This routine set the format specification for binary storage of a
! variable in RADS. The specification has the form: type[,scale[,offset]]
! where
! - 'type' is one of: int1, byte, int2, short, int4, int, real, real4,
! dble, double
! - 'scale' is the scale to revert the binary back to a float (optional)
! - 'offset' is the offset to add to the result when unpacking (optional)
! Instead of commas, spaces or slashes can be used as separators.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! varname  : Variable name
! format   : Binary storage specification (e.g. int2,1e-3)
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_var
!****-------------------------------------------------------------------
character(len=len(format)) :: temp
integer :: i, ios
type(rads_var), pointer :: var

var => rads_varptr(S, varname)
var%info%scale_factor = 1d0
var%info%add_offset = 0d0

do i = 1,len(format)
	if (format(i:i) == ',' .or. format(i:i) == '/') then
		temp(i:i) = ' '
	else
		temp(i:i) = format(i:i)
	endif
enddo
i = index(temp, ' ')
select case (temp(:i-1))
case ('int1', 'byte')
	var%info%nctype = nf90_int1
case ('int2', 'short')
	var%info%nctype = nf90_int2
case ('int4', 'int')
	var%info%nctype = nf90_int4
case ('real', 'real4', 'float')
	var%info%nctype = nf90_real
case default
	var%info%nctype = nf90_double
end select
read (temp(i+1:), *, iostat=ios) var%info%scale_factor, var%info%add_offset
end subroutine rads_set_compress

!****if* rads/rads_set_quality_flag
! SUMMARY
! Set quality flag(s) for given variable
!
! SYNOPSIS
subroutine rads_set_quality_flag (S, varname, flag)
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: varname, flag
!
! PURPOSE
! This routine set appends <flag> to the variables contained in the set
! of variables to check to allow variable <varname> to pass (if <flag> is
! not already contained in that set).
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! varname  : Variable name
! flag     : Name of the variable that needs to be checked to validate
!             <varname>.
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_var
!****-------------------------------------------------------------------
type(rads_var), pointer :: var
var => rads_varptr(S, varname)
if (.not.associated(var)) then
	return
else if (var%info%quality_flag == '') then
	var%info%quality_flag = flag
else if (index(var%info%quality_flag,flag) == 0) then
	var%info%quality_flag = trim(var%info%quality_flag) // ' ' // flag
endif
end subroutine rads_set_quality_flag

!****if* rads/rads_long_name_and_units
! SUMMARY
! Get string of long_name and units (or flag_meanings)
!
! SYNOPSIS
subroutine rads_long_name_and_units (var, unit)
use rads_misc
type(rads_var), intent(in) :: var
integer(fourbyteint), optional, intent(in) :: unit
!
! PURPOSE
! This writes to the output file attached to <unit> the long_name and units (or
! flag_meanings) of a variable. The output can be one of the following:
!  long_name [units]
!  long_name [bits: flag_meanings]
!  long_name [values: flag_meanings]
!
! ARGUMENTS
! var      : Pointer to variable
! unit     : Fortran output unit (6 = stdout (default), 0 = stderr)
!****-------------------------------------------------------------------
integer :: iunit, i, i0, i1
iunit = stdout
if (present(unit)) iunit = unit
write (iunit,550,advance='no') trim(var%long_name) // ' ['
if (var%info%datatype == rads_type_flagmasks) then
	write (iunit,550,advance='no') 'bits:'
else if (var%info%datatype == rads_type_flagvalues) then
	write (iunit,550,advance='no') 'values:'
else
	write (iunit,550) trim(var%info%units) // ']'
	return
endif
i1 = 0
i = 0
do
	if (.not.next_word (var%info%flag_meanings, i0, i1)) exit
	if (i1 == i0) cycle
	write (iunit,551,advance='no') i,var%info%flag_meanings(i0:i1-1)
	i = i + 1
enddo
write (iunit,550) ']'
550 format (a)
551 format (1x,i0,'=',a)
end subroutine rads_long_name_and_units

!****if* rads/rads_stat_0d
! SUMMARY
! Print the RADS statistics for a given satellite
!
! SYNOPSIS
subroutine rads_stat_0d (S)
use rads_time
type(rads_sat), intent(in) :: S
!
! PURPOSE
! This is the scalar version of rads_stat.
!****-------------------------------------------------------------------
integer(fourbyteint) :: i
write (rads_log_unit, 700)
write (rads_log_unit, 710) trim(S%satellite), S%sat, timestamp(), trim(S%command)

write (rads_log_unit, 720) 'PASSES QUERIED', sum(S%pass_stat)
write (rads_log_unit, 724) 'REJECTED','SELECTED','LOWER','UPPER','STEP'
write (rads_log_unit, 730) 'Cycle number limits', S%pass_stat(1), sum(S%pass_stat(2:7)), S%cycles
write (rads_log_unit, 730) 'Pass number limits' , S%pass_stat(2), sum(S%pass_stat(3:7)), S%passes
write (rads_log_unit, 731) 'Time limits' , S%pass_stat(3), sum(S%pass_stat(4:7)), datestring (S%time%info%limits)
write (rads_log_unit, 732) 'Equator longitude limits (asc)' , S%pass_stat(5:7:2), S%eqlonlim(1,:)
write (rads_log_unit, 732) 'Equator longitude limits (des)' , S%pass_stat(4:6:2), S%eqlonlim(0,:)

write (rads_log_unit, 721) 'PASSES AND MEASUREMENTS READ', sum(S%pass_stat(6:7)), S%total_read
write (rads_log_unit, 724) 'REJECTED','SELECTED','LOWER','UPPER','MIN',' MAX','MEAN','STDDEV'
call rads_stat_line (S%time)
call rads_stat_line (S%lat)
call rads_stat_line (S%lon)

write (rads_log_unit, 720) 'MEASUREMENTS IN REQUESTED PERIOD AND REGION', S%total_inside
write (rads_log_unit, 724) 'REJECTED','SELECTED','LOWER','UPPER','MIN',' MAX','MEAN','STDDEV'
do i = 1,S%nvar
	if (.not.associated(S%var(i)%info)) then ! Skip undefined variables
	else if (S%var(i)%info%selected + S%var(i)%info%rejected == 0) then ! Skip "unused" variables
	else if (S%var(i)%info%datatype >= rads_type_time) then ! Skip time, lat, lon
	else if (S%var(i)%name /= S%var(i)%info%name) then ! Skip aliases
	else ! Print statistics line for whatever remains
		call rads_stat_line (S%var(i))
	endif
enddo
write (rads_log_unit, 700)

700 format (134('#'))
710 format ('# Editing statistics for ',a,' (',a,')'/'# Created: ',a,' UTC: ',a)
720 format ('#'/'# ',a,t53,i10/'#')
721 format ('#'/'# ',a,t43,2i10/'#')
724 format ('#',t43,2a10,6a12)
730 format ("# ",a,t43,2i10,3i12)
731 format ("# ",a,t43,2i10,2a)
732 format ("# ",a,t43,2i10,2f12.3)

contains

subroutine rads_stat_line (var)
type(rads_var), intent(in) :: var
real(eightbytereal), parameter :: sec2000 = 473299200d0
type(rads_varinfo), pointer :: info
info => var%info
write (rads_log_unit, '("# ",a," [",a,"]",t43,2i10)', advance='no') trim(var%long_name), trim(info%units), &
	info%rejected, info%selected
if (info%units(:13) /= 'seconds since') then
	write (rads_log_unit, '(6f12.3)') info%limits, info%xmin, info%xmax, info%mean, sqrt(info%sum2/(info%selected-1))
else if (info%units(15:18) == '1985') then
	write (rads_log_unit, '(5a,f12.0)') datestring (info%limits), datestring(info%xmin), datestring(info%xmax), &
		datestring(info%mean), sqrt(info%sum2/(info%selected-1))
else
	write (rads_log_unit, '(5a,f12.0)') datestring (info%limits+sec2000), datestring(info%xmin+sec2000), &
		datestring(info%xmax+sec2000), datestring(info%mean+sec2000), sqrt(info%sum2/(info%selected-1))
endif
end subroutine rads_stat_line

elemental function datestring (sec)
use rads_time
real(eightbytereal), intent(in) :: sec
character(len=12) :: datestring
integer(fourbyteint) :: yy, mm, dd, hh, mn
real(eightbytereal) :: ss
if (sec /= sec) then
	datestring = '         NaN'
else
	call sec85ymdhms(sec, yy, mm, dd, hh, mn, ss)
	write (datestring, '(1x,3i2.2,"/",2i2.2)') mod(yy,100), mm, dd, hh, mn
endif
end function datestring

end subroutine rads_stat_0d

subroutine rads_stat_1d (S)
type(rads_sat), intent(in) :: S(:)
integer :: i
do i = 1,size(S)
	if (S(i)%sat /= '') call rads_stat_0d (S(i))
enddo
end subroutine rads_stat_1d

!****if* rads/rads_exit
! SUMMARY
! Exit RADS with error message
!
! SYNOPSIS
subroutine rads_exit (string)
character(len=*), intent(in) :: string
!
! PURPOSE
! This routine terminates RADS after printing an error message.
!
! ARGUMENT
! string   : Error message
!****-------------------------------------------------------------------
character(len=rads_naml) :: progname
call getarg (0, progname)
call rads_message (string)
call rads_message ('Use "'//trim(progname)//' --help" for more info')
call exit (10)
end subroutine rads_exit

!****if* rads/rads_error
! SUMMARY
! Print error message and store error code
!
! SYNOPSIS
subroutine rads_error (S, ierr, string, P)
type(rads_sat), intent(inout) :: S
integer(fourbyteint), intent(in) :: ierr
character(len=*), intent(in) :: string
type(rads_pass), intent(in), optional :: P
!
! PURPOSE
! This routine prints an error message and sets the error code.
! The message is subpressed when -q is used (rads_verbose < 0)
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! ierr     : Error code
! string   : Error message
! P        : If used, it will add the pass file name at the end of <string>
!
! ERROR CODE
! S%error  : Will be set to ierr when not rads_noerr
!****-------------------------------------------------------------------
call rads_message (string, P)
if (ierr /= rads_noerr) S%error = ierr
end subroutine rads_error

!****if* rads/rads_message
! SUMMARY
! Print warning message
!
! SYNOPSIS
subroutine rads_message (string, P)
use rads_netcdf
character(len=*), intent(in) :: string
type(rads_pass), intent(in), optional :: P
!
! PURPOSE
! This routine prints a message to standard error.
! The message is subpressed when -q is used (rads_verbose < 0)
!
! ARGUMENTS
! string   : Error message
! P        : If used, it will add the pass file name at the end of <string>
!****-------------------------------------------------------------------
if (rads_verbose < 0) then
	! Remain quiet
else if (present(P)) then
	call nf90_message (trim(string) // ' ' // P%fileinfo(1)%name)
else
	call nf90_message (string)
endif
end subroutine rads_message

!****if* rads/rads_opt_error
! SUMMARY
! Print error message about failed option scanning
!
! SYNOPSIS
subroutine rads_opt_error (opt, arg)
character(len=*), intent(in) :: opt, arg
!
! PURPOSE
! This routine prints an error message to standard error in case RADS
! failed to identify or scan an option on the command line.
! The message is subpressed when -q is used (rads_verbose < 0)
!
! ARGUMENTS
! opt      : Command line option
! arg      : Argument of the command line option
!****-------------------------------------------------------------------
if (opt == ':') then ! Unknown option
	call rads_message ('Unknown option '//trim(arg)//' skipped')
else if (opt == '::') then ! Missing required argument
	call rads_exit ('Option "'//trim(arg)//'" should be followed by a required argument')
else if (len_trim(opt) == 1) then ! Short option
	call rads_exit ('Option "-'//trim(opt)//' '//trim(arg)//'" could not be parsed')
else if (len_trim(opt) == 2 .and. opt(2:2) == ':') then ! Long option
	call rads_exit ('Option "--'//trim(opt(1:1))//' '//trim(arg)//'" could not be parsed')
else ! Long option
	call rads_exit ('Option "--'//trim(opt)//' '//trim(arg)//'" could not be parsed')
endif
end subroutine rads_opt_error

!****if* rads/rads_version
! SUMMARY
! Print message about current program and version
!
! SYNOPSIS
logical function rads_version (description, unit, flag)
character(len=*), intent(in), optional :: description, flag
integer(fourbyteint), intent(in), optional :: unit
!
! PURPOSE
! This routine prints out a message in one of the following forms,
! on the first argument on the command line or the optional argument <flag>.
! 1) When first command line argument or <flag> is --version:
!    "rads_program, version <number>"
!    The program then terminates here.
! 2) When no <description> is given:
!    "rads_program (v<number>)"
!    Return value is .true.
! 3) When first command line argument or <flag> is -? or --help:
!    "rads_program (v<number>): <description>"
!    Return value is .false.
! 4) When first command line argument of <flag> is --head:
!    "rads_program (v<number>): <description>"
!    Return value is .true.
! 5) Otherwise:
!    No output
!    Return value is .true.
!
! ARGUMENTS
! description : One-line description of program
! unit        : Fortran output unit (6 = stdout (default), 0 = stderr)
! flag        : Use string in replacement of first command line argument
!
! RETURN VALUE
! rads_version: .false. if output is of type 3, otherwise .true.
!****-------------------------------------------------------------------
integer :: iunit
character(len=rads_naml) :: progname, arg
call getarg (0, progname)
call getarg (1, arg)
if (present(flag) .and. arg /= '--version' .and. arg /= '--help' .and. arg /= '-?' .and. arg /= '--head') arg = flag
rads_version = .true.
if (present(unit)) then
	iunit = unit
else
	iunit = stdout
endif
if (arg == '--version') then
	write (iunit, 1320) trim(progname), trim(rads_version_id)
	stop
else if (.not.present(description)) then
	write (iunit, 1300) trim(progname), trim(rads_version_id)
else if (arg == '--help' .or. arg == '-?') then
	write (iunit, 1310) trim(progname), trim(rads_version_id), trim(description)
	rads_version = .false.
else if (arg == '--head') then
	write (iunit, 1310) trim(progname), trim(rads_version_id), trim(description)
endif
1300 format (a,' (',a,')')
1310 format (a,' (',a,'): ',a)
1320 format (a,', version ',a)
end function rads_version

!****if* rads/rads_synopsis
! SUMMARY
! Print general usage information for all RADS programs
!
! SYNOPSIS
subroutine rads_synopsis (unit)
integer(fourbyteint), intent(in), optional :: unit
!
! PURPOSE
! This routine prints out the usage information that is common to all
! RADS programs.
!
! ARGUMENTS
! unit        : Fortran output unit (6 = stdout (default), 0 = stderr)
!****-------------------------------------------------------------------
integer :: iunit
character(len=rads_naml) :: progname
call getarg (0, progname)
if (present(unit)) then
	iunit = unit
else
	iunit = stdout
endif
write (iunit, 1300) trim(progname)
1300 format (/ &
'Usage: ',a,' [required_arguments] [rads_dataselectors] [rads_options] [program_options]' // &
'Required argument is:'/ &
'  -S, --sat SAT[/PHASE]     Specify satellite [and phase] (e.g. e1/g, tx)'// &
'Optional [rads_dataselectors] are:'/ &
'  -A, --alias VAR1=VAR2     Use variable VAR2 when VAR1 is requested'/ &
'  -C, --cycle C0[,C1[,DC]]  Specify first and last cycle and modulo'/ &
'  -F, --fmt, --format VAR=FMT'/ &
'                            Specify Fortran format used to print VAR (for ASCII output only)'/ &
'  -L, --limits VAR=MIN,MAX  Specify edit data range for variable VAR'/ &
'  --lon LON0,LON1           Specify longitude boundaries (deg)'/ &
'  --lat LAT0,LAT1           Specify latitude  boundaries (deg)'/ &
'  -P, --pass P0[,P1[,DP]]   Specify first and last pass and modulo; alternatively use -Pa'/ &
'                            (--pass asc) or -Pd (--pass des) to restrict selection to'/ &
'                            ascending or descending passes only'/ &
'  -Q, --quality-flag VAR=FLAG'/&
'                            Check variable FLAG when validating variable VAR'/ &
'  -R, --region LON0,LON1,LAT0,LAT1'/ &
'                            Specify rectangular region (deg)'/ &
'  -R, --region LON0,LAT0,RADIUS'/ &
'                            Specify circular region (deg)' / &
'  --sla SLA0,SLA1           Specify range for SLA (m)'/ &
'  --time T0,T1              Specify time selection (optionally use --ymd, --doy,'/ &
'                            or --sec for [YY]YYMMDD[HHMMSS], [YY]YYDDD, or SEC85)'/ &
'  -V, --var VAR1,...        Select variables to be read'/ &
'  -X, --xml XMLFILE         Load configuration file XMLFILE in addition to the defaults'/ &
'  -Z, --cmp, --compress VAR=TYPE[,SCALE[,OFFSET]]'/ &
'                            Specify binary output format for variable VAR (NetCDF only); TYPE'/ &
'                            is one of: int1, int2, int4, real, dble; SCALE and OFFSET are'/ &
'                            optional (def: 1,0)'// &
'Still working for backward compatibility with RADS3 are options:'/ &
'  --h H0,H1                 Specify range for SLA (m) (now --sla H0,H1)'/ &
'  --opt J                   Use field number J when J/100 requested (now -A VAR1=VAR2)'/ &
'  --opt I=J                 Make field I (range 1-99) and alias for field J (now -A VAR1=VAR2)'/ &
'  --sel VAR1,...            Select variables to read'// &
'Common [rads_options] are:'/ &
'  --args FILENAME           Get all command line arguments from FILENAME (1 argument per line)'/ &
'  --debug LEVEL             Set debug/verbosity level'/ &
'  --help                    Print this syntax massage'/ &
'  --log FILENAME            Send statistics to FILENAME (default is standard output)'/ &
'  -q, --quiet               Suppress warning messages (but keeps fatal error messages)'/ &
'  -v, --verbose             Increase verbosity level'/ &
'  --version                 Print version info only'/ &
'  --                        Terminates all options; all following command-line arguments are'/ &
'                            considered non-option arguments, even if they begin with a hyphen')
end subroutine rads_synopsis

!****if* rads/rads_init_phase
! SUMMARY
! Add new mission phase and get pointer to satellite phase info
!
! SYNOPSIS
function rads_init_phase (S, name) result (phase)
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: name
type(rads_phase), pointer :: phase
!
! PURPOSE
! Create a new phase, adding it to the <S%phases> array and create
! the array if needed. Then create pointer to the phase definitions
! for this phase.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! name     : Name of phase
!
! ERROR CODE:
! S%error  : rads_noerr, rads_err_nophase
!****-------------------------------------------------------------------
integer(fourbyteint) :: n
type(rads_phase), pointer :: temp(:)
nullify (phase)
n = 0
if (associated(S%phases)) n = size(S%phases)

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
S%phases(n) = rads_phase (name(1:1), '', (/999,0/), 0, nan, nan, nan, nan, 0, 0, 0, nan, nan, nan, 0, 0, null())
phase => S%phases(n)
end function rads_init_phase

subroutine rads_set_phase_by_name (S, name, error)
type(rads_sat), intent(inout) :: S
character(len=*), intent(in) :: name
logical, optional, intent(out) :: error
integer :: i
if (present(error)) error = .false.
S%error = rads_noerr
! Check if we are already in the right mission phase
if (name(1:1) == S%phase%name(1:1)) return
! Check all mission phases
do i = 1,size(S%phases)
	if (name(1:1) == S%phases(i)%name(1:1)) then
		S%phase => S%phases(i)
		return
	endif
enddo
if (present(error)) then
	error = .true.
	S%error = rads_err_nophase
else
	call rads_exit ('No mission phase "'//trim(name)//'" for satellite "'//S%sat//'"')
endif
end subroutine rads_set_phase_by_name

subroutine rads_set_phase_by_cycle (S, cycle, error)
type(rads_sat), intent(inout) :: S
integer(fourbyteint), intent(in) :: cycle
logical, optional, intent(out) :: error
integer :: i
character(len=3) :: name
if (present(error)) error = .false.
S%error = rads_noerr
! Check if we are already in the right mission phase
if (cycle >= S%phase%cycles(1) .and. cycle <= S%phase%cycles(2)) return
! Check all mission phases
do i = 1,size(S%phases)
	if (cycle >= S%phases(i)%cycles(1) .and. cycle <= S%phases(i)%cycles(2)) then
		S%phase => S%phases(i)
		return
	endif
enddo
if (present(error)) then
	error = .true.
	S%error = rads_err_nophase
else
	write (name, '(i3.3)') cycle
	call rads_exit ('No cycle '//name//' for any mission phase of satellite "'//S%sat//'"')
endif
end subroutine rads_set_phase_by_cycle

subroutine rads_set_phase_by_time (S, time, error)
use rads_time
type(rads_sat), intent(inout) :: S
real(eightbytereal), intent(in) :: time
logical, optional, intent(out) :: error
integer :: i
if (present(error)) error = .false.
S%error = rads_noerr
! Check if we are already in the right mission phase
if (time >= S%phase%start_time .and. time <= S%phase%end_time) return
! Check all mission phases
do i = 1,size(S%phases)
	if (time >= S%phases(i)%start_time .and. time <= S%phases(i)%end_time) then
		S%phase => S%phases(i)
		return
	endif
enddo
if (present(error)) then
	error = .true.
	S%error = rads_err_nophase
else
	call rads_exit ('Time '//strf1985f(time)//' is outside any mission phase of satellite "'//S%sat//'"')
endif
end subroutine rads_set_phase_by_time

!****f* rads/rads_predict_equator
! SUMMARY
! Predict equator crossing time and longitude
!
! SYNOPSIS
subroutine rads_predict_equator (S, P)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
!
! PURPOSE
! This routine estimates the equator time and longitude, as well
! as the start and end time of given cycle number and pass number
! stored in the pass struct <P>..
!
! The routine works for exact repeat orbits as well as drifting
! orbits.
!
! The estimated variables (equator_time, equator_lon, start_time, end_time)
! are returned in the pass struct <P>.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass structure
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_nophase
!****-------------------------------------------------------------------
integer(fourbyteint) :: cc, pp
real(eightbytereal) :: d, e
logical :: error

 cc = P%cycle
 pp = P%pass

! If the cycle is out of range for the current phase, look for a new phase
call rads_set_phase (S, cc, error)
if (error) return

! For constructions with subcycles, convert first to "real" cycle/pass number
if (associated(S%phase%subcycles)) then
	cc = cc - S%phase%subcycles%i
	pp = S%phase%subcycles%list(modulo(cc,S%phase%subcycles%n)+1) + pp
	cc = cc / S%phase%subcycles%n + 1
endif

! Now do the estimation process
d = S%phase%pass_seconds
P%equator_time = S%phase%ref_time + ((cc - S%phase%ref_cycle) * S%phase%repeat_passes + (pp - S%phase%ref_pass)) * d
P%start_time = P%equator_time - 0.5d0 * d
P%end_time = P%equator_time + 0.5d0 * d
d = -S%phase%repeat_nodal * 360d0 / S%phase%repeat_passes ! Longitude advance per pass due to precession of node and earth rotation
e = S%phase%repeat_shift * (P%equator_time - S%phase%ref_time) / S%phase%pass_seconds / S%phase%repeat_passes ! Longitude shift per cycle
P%equator_lon = modulo(S%phase%ref_lon + (pp - S%phase%ref_pass) * d + e + modulo(pp - S%phase%ref_pass,2) * 180d0, 360d0)
if (rads_verbose >= 4) write (*,'(a,3f15.3,f12.6)') 'Estimated start/end/equator time/longitude = ', &
	P%start_time, P%end_time, P%equator_time, P%equator_lon
S%error = rads_noerr
end subroutine rads_predict_equator

!****if* rads/rads_time_to_cycle_pass
! SUMMARY
! Determine cycle and pass number for given epoch
!
! SYNOPSIS
subroutine rads_time_to_cycle_pass (S, time, cycle, pass, abs_orbit)
type(rads_sat), intent(inout) :: S
real(eightbytereal), intent(in) :: time
integer(fourbyteint), intent(out) :: cycle, pass
integer(fourbyteint), intent(out), optional :: abs_orbit
!
! PURPOSE
! Given an epoch <time> in seconds since 1985, determine the <cycle>
! number and <pass> number in which that epoch falls.
! Optionally, the absolute orbit number at the equator (<abs_orbit>)
! can be computed.
!
! This routine also updates the <S%phase> pointer to point to the
! appropriate structure containing the phase information.
!
! ARGUMENTS
! S            : Satellite/mission dependent structure
! time         : Time in seconds since 1985
! cycle        : Cycle number in which <time> falls
! pass         : Pass number in which <time> falls
! abs_orbit    : Absolute orbit number (optional)
!****-------------------------------------------------------------------
integer :: i, j, n
real(eightbytereal) :: d, t0, x

! Look for the appropriate phase
i = 1
S%error = rads_noerr
do i = 1,size(S%phases)-1
	if (time < S%phases(i+1)%start_time) exit
enddo

! Link the selected phase
S%phase => S%phases(i)

! Estimate the cycle from the time
d = S%phase%pass_seconds
t0 = S%phase%ref_time - (S%phase%ref_pass - 0.5d0) * d ! Time of start of ref_cycle
x = (time - t0) / (S%phase%repeat_days * 86400d0) + S%phase%ref_cycle
cycle = floor(x)
pass = int((x - cycle) * S%phase%repeat_passes + 1)

if (present(abs_orbit)) then
	! Compute the absolute orbit number (add 2*ref_orbit first to avoid dividing negative numbers)
	abs_orbit = (S%phase%ref_orbit * 2 + (cycle - S%phase%ref_cycle) * S%phase%repeat_passes + (pass - S%phase%ref_pass)) / 2
endif

! When there are subcycles, compute the subcycle number
if (associated(S%phase%subcycles)) then
	x = (time - t0) - (cycle - S%phase%ref_cycle) * (S%phase%repeat_days * 86400d0)
	cycle = (cycle - 1) * S%phase%subcycles%n + S%phase%subcycles%i
	n = floor(x / d)
	do j = 2,S%phase%subcycles%n
		if (S%phase%subcycles%list(j) > n) exit
		cycle = cycle + 1
	enddo
	pass = n - S%phase%subcycles%list(j-1) + 1
endif
end subroutine rads_time_to_cycle_pass

!****if* rads/rads_time_to_cycle
! SUMMARY
! Determine cycle number for given epoch
!
! SYNOPSIS
function rads_time_to_cycle (S, time)
type(rads_sat), intent(inout) :: S
real(eightbytereal), intent(in) :: time
integer(fourbyteint) :: rads_time_to_cycle
!
! PURPOSE
! Given an epoch <time> in seconds since 1985, determine the cycle
! number in which that epoch falls.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! time     : Time in seconds since 1985
!
! RETURN VALUE
! rads_time_to_cycle : Cycle number in which <time> falls
!****-------------------------------------------------------------------
integer(fourbyteint) :: pass
call rads_time_to_cycle_pass (S, time, rads_time_to_cycle, pass)
end function rads_time_to_cycle

!***********************************************************************
!*rads_cycle_to_time -- Determine start or end time of cycle
!
! SYNOPSIS
function rads_cycle_to_time (S, cycle)
type(rads_sat), intent(inout) :: S
integer(fourbyteint), intent(in) :: cycle
real(eightbytereal) :: rads_cycle_to_time
!
! Given a cycle number, estimate the start time of that cycle.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! cycle    : Cycle number
!
! RETURN VALUE
! rads_cycle_to_time : Start time of <cycle> in seconds since 1985
!****-------------------------------------------------------------------
integer :: i, cc, pp
real(eightbytereal) :: d, t0
! Find the correct phase
i = 1
do i = 1,size(S%phases)-1
	if (cycle < S%phases(i+1)%cycles(1)) exit
enddo

d = S%phases(i)%pass_seconds
t0 = S%phases(i)%ref_time - (S%phases(i)%ref_pass - 0.5d0) * d ! Time of start of ref_cycle

if (associated(S%phases(i)%subcycles)) then
	! Convert subcycle numbering to "real" cycle numbering
	cc = cycle - S%phases(i)%subcycles%i
	pp = S%phases(i)%subcycles%list(modulo(cc,S%phases(i)%subcycles%n)+1)
	cc = cc / S%phases(i)%subcycles%n + 1
	rads_cycle_to_time = max(S%phases(i)%start_time, &
		t0 + (cc - S%phases(i)%ref_cycle) * S%phases(i)%repeat_days * 86400d0 + pp * d)
else
	rads_cycle_to_time = max(S%phases(i)%start_time, &
		t0 + (cycle - S%phases(i)%ref_cycle) * S%phases(i)%repeat_days * 86400d0)
endif
end function rads_cycle_to_time

!****if* rads/rads_traxxing
! SUMMARY
! Determine which tracks cross an area
!
! SYNOPSIS
subroutine rads_traxxing (S)
use rads_misc
type(rads_sat), intent(inout) :: S
!
! PURPOSE
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
! ARGUMENT
! S        : Satellite/mission dependent structure
!****-------------------------------------------------------------------
real(eightbytereal) :: u(2),l(2),latmax
integer(fourbyteint) :: i

! If either longitude limit is NaN, set eqlonlin to NaN
if (any(isnan_(S%lon%info%limits))) then
	S%eqlonlim = nan
	return
endif

! Initialize
latmax = min(S%inclination, 180d0-S%inclination)

do i = 1,2
	! Compute argument of latitude for lower and upper latitude boundary
	! This also accounts for NaNs
	if (abs(S%lat%info%limits(i)) <= latmax) then
		u(i) = asin(sin(S%lat%info%limits(i)*rad) / sin(S%inclination*rad)) / rad
		l(i) = asin(tan(S%lat%info%limits(i)*rad) / tan(S%inclination*rad)) / rad
	else
		u(i) = (2*i-3)*90d0
		l(i) = u(i)
		if (S%inclination > 90d0) l(i) = -l(i)
	endif
enddo

! Compute the longitude advance corrected for time
l = l - u * 2d0 * S%phase%repeat_nodal / S%phase%repeat_passes

! Compute the equator crossing longitudes for ascending tracks
! Add 2-degree margin
S%eqlonlim(1,1) = S%lon%info%limits(1) - max(l(1),l(2)) - 2d0
S%eqlonlim(1,2) = S%lon%info%limits(2) - min(l(1),l(2)) + 2d0

! Compute the equator crossing longitudes for descending tracks
! Add 2-degree margin
S%eqlonlim(0,1) = S%lon%info%limits(1) + min(l(1),l(2)) - 2d0
S%eqlonlim(0,2) = S%lon%info%limits(2) + max(l(1),l(2)) + 2d0

if (rads_verbose >= 3) write (*,'(a,4f12.6)') "Eqlonlim = ",S%eqlonlim

end subroutine rads_traxxing

!****if* rads/rads_progress_bar
! SUMMARY
! Print and update progress of scanning cycles/passes
!
! SYNOPSIS
subroutine rads_progress_bar (S, P, nselpass)
type(rads_sat), intent(in) :: S
type(rads_pass), intent(in) :: P
integer(fourbyteint), intent(in) :: nselpass
!
! PURPOSE
! This routine prints information on the progress of cycling through
! cycles and passes of satellite data. If used, it should preferably be
! called before every call of rads_close_pass.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass structure
! nselpass : Number of measurements selected in this pass
!****-------------------------------------------------------------------
integer :: pos_old = 50, lin_old = -1, cycle_old = -1, pos, lin, i

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

if (cycle_old < 0) write (rads_log_unit,700) trim(S%satellite),trim(S%phase%name)
pos = mod(P%pass-1,100)+1
lin = int((P%pass-1)/100)
if (P%cycle /= cycle_old) then
	write (rads_log_unit,710,advance='no')     ! Print header
	pos_old = 100
	cycle_old = P%cycle
endif
if (pos_old == 100 .or. lin_old /= lin) then
	write (rads_log_unit,720,advance='no') P%cycle,int((P%pass-1)/100)*100+1
	! Start new line, print cycle nr and pass nr
	pos_old = 0
endif
do i = pos_old+1,pos-1
	write (rads_log_unit,750,advance='no') ' ' ! Fill up to current pass with ' '
enddo
if (S%error /= rads_noerr) then
	write (rads_log_unit,750,advance='no') '-' ! Print - for non existing pass
else if (P%ndata == 0) then
	write (rads_log_unit,750,advance='no') 'x' ! Print x for empty pass
else if (nselpass == 0) then
	write (rads_log_unit,750,advance='no') 'o' ! Print o when no data is selected
else
	write (rads_log_unit,750,advance='no') '#' ! Print # when data available
endif
pos_old = pos
lin_old = lin

end subroutine rads_progress_bar

!****if* rads/rads_parse_r_option
! SUMMARY
! Parse the -r option
!
! SYNOPSIS
subroutine rads_parse_r_option (S, opt, arg, idx)
type(rads_sat), intent(in) :: S
character(len=*), intent(in) :: opt, arg
integer(fourbyteint), intent(out) :: idx
!
! PURPOSE
! This routine is used generally with the -r option, but can also be
! associated with any other command line option.
! The routine parses the argument of the option and returns the index
! value (position) of the associated variable on the -V option argument.
!
! The following cases are recognised
! Argument                Value of idx
! 'n' or 'any'            -2
! '' or '0' or 'none'     0
! 1 <= integer <= nsel    value of integer
! name of variable on -V  positional index of variable on -V option
! any other value         quit program with error message
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! opt      : Command line option name ('r' or 'reject-on-nan')
! arg      : Command line option argument
! idx      : Index value
!****-------------------------------------------------------------------
integer(fourbyteint) :: i
select case (arg)
case ('n', 'any')
	idx = -2
case ('', '0', 'none')
	idx = 0
case default
	if (arg(1:1) >= '0' .and. arg(1:1) <= '9') then
		idx = 0
		read (arg, *, iostat=i) idx
		if (idx < 1 .or. idx > S%nsel) call rads_exit ('Option -'//trim(opt)//' NR used with invalid value')
	else
		do i = 1,S%nsel
			if (arg == S%sel(i)%name .or. arg == S%sel(i)%info%name) then
				idx = i
				return
			endif
		enddo
		call rads_exit ('Option -'//trim(opt)//' <varname> does not refer to variable specified on -V option')
	endif
end select
end subroutine rads_parse_r_option

!****if* rads/rads_create_pass
! SUMMARY
! Create RADS pass (data) file
!
! SYNOPSIS
subroutine rads_create_pass (S, P, ndata, n_hz, n_wvf, name)
use netcdf
use rads_netcdf
use rads_time
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
integer(fourbyteint), intent(in), optional :: ndata, n_hz, n_wvf
character(len=*), intent(in), optional :: name
!
! PURPOSE
! This routine creates a new RADS NetCDF data file. If one of the same
! file name already exists, it is removed.
! The file is initialized with the appropriate global attributes, and
! the dimensions ('time' and optionally 'meas_ind' and 'wvf_ind') will be set up.
! The primary dimension can either be fixed (ndata > 0) or unlimited (ndata == 0).
! If an unlimited dimension is selected, then the pass-related global
! attributes, like cycle, pass, equator_time, equator_lon will not be written
! to the file.
! When the <ndata> argument is omitted, the value <P%ndata> is used.
!
! The optional argument <n_hz> gives the size of the secondary dimension for
! multi-Hertz data. When <n_hz> is omitted, the value <P%n_hz> is used.
! When n_hz == 0, no secondary dimension is created.
!
! The optional argument <name> can have one of three forms:
! - Left empty it specifies the current directory
! - Ending in '/' it specifies a named directory in which to store the pass
!   file of the form <sat>p<pass>c<cycle>.nc
! - Otherwise it specifies the file name in full.
! If <name> is omitted, the default file name is used:
! $RADSDATAROOT/<sat>/<phase>/c<cycle>/<sat>p<pass>c<cycle>.nc
! Any directory that does not yet exist will be created.
!
! Upon entry, the <P> structure needs to contain the relevant information on
! the pass (cycle, pass, equator_time, equator_lon). Upon return, the
! P%fileinfo will be updated.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass structure
! ndata    : Length of the primary ('time') dimension (use 0 for unlimited)
! n_hz     : Number of multi-Hertz data per second
! n_wvf    : Number of waveform gates
! name     : Name of directory in which to store pass file, or file name
!
! ERROR CODE
! S%error  : rads_noerr, rads_err_nc_create
!****-------------------------------------------------------------------
integer(fourbyteint) :: i, l, e
logical :: exist
real(eightbytereal), parameter :: ellipsoid_axis = 6378136.3d0, ellipsoid_flattening = 1d0/298.257d0
integer(fourbyteint) :: ncid
character(len=rads_cmdl) :: filename

! Initialise
S%error = rads_noerr
if (present(ndata)) P%ndata = ndata
if (present(n_hz)) P%n_hz = n_hz
if (present(n_wvf)) P%n_wvf = n_wvf

! Build the file name, make directory if needed
600 format (a,'/',a,'/',a,'/c',i3.3,'/',a2,'p',i4.4,'c',i3.3,'.nc')
610 format (a,'p',i4.4,'c',i3.3,'.nc')
if (.not.present(name)) then
	write (filename, 600) trim(S%dataroot), trim(S%branch(1)), trim(S%phase%name), P%cycle, S%sat, P%pass, P%cycle
	l = len_trim(filename)-15
	inquire (file = filename(:l), exist = exist)
	if (.not.exist) call system ('mkdir -p ' // filename(:l))
else if (name == '') then
	write (filename, 610) S%sat, P%pass, P%cycle
else if (name(len_trim(name):) == '/') then
	write (filename, 610) trim(name)//S%sat, P%pass, P%cycle
	inquire (file = name, exist = exist)
	if (.not.exist) call system ('mkdir -p ' // name)
else
	call rads_init_pass_struct (S, P)
	filename = name
endif
P%fileinfo(1)%name = filename

! Create the (new) data file
if (rads_verbose >= 2) write (*,'(2a,i10)') 'Creating ',trim(filename),P%ndata
if (nft(nf90_create(filename, nf90_write+nf90_nofill, ncid))) then
	call rads_error (S, rads_err_nc_create, 'Error creating file', P)
	return
endif
P%fileinfo(1)%ncid = ncid

! Define the principle dimension
if (nft(nf90_def_dim (ncid, 'time', P%ndata, l))) then
	call rads_error (S, rads_err_nc_create, 'Error creating dimension "time" in file', P)
	return
endif

! Define second and third dimension (if requested)
if (P%n_hz == 0) then
	! Do nothing
else if (nft(nf90_def_dim (ncid, 'meas_ind', P%n_hz, l))) then
	call rads_error (S, rads_err_nc_create, 'Error creating dimension "meas_ind" in file', P)
	return
endif
if (P%n_wvf == 0) then
	! Do nothing
else if (nft(nf90_def_dim (ncid, 'wvf_ind', P%n_wvf, l))) then
	call rads_error (S, rads_err_nc_create, 'Error creating dimension "wvf_ind" in file', P)
	return
endif
P%rw = .true.

! Specify the global attributes
e = 0
do i = 1,size(S%glob_att)
	l = index(S%glob_att(i),' ')
	e = e + nf90_put_att (ncid, nf90_global, S%glob_att(i)(:l-1), S%glob_att(i)(l+1:))
enddo
e = e + nf90_put_att (ncid, nf90_global, 'ellipsoid_axis', ellipsoid_axis) + &
	nf90_put_att (ncid, nf90_global, 'ellipsoid_flattening', ellipsoid_flattening)
l = index(filename, '/', .true.) + 1
e = e + nf90_put_att (ncid, nf90_global, 'filename', trim(filename(l:))) + &
	nf90_put_att (ncid, nf90_global, 'mission_name', trim(S%satellite)) + &
	nf90_put_att (ncid, nf90_global, 'mission_phase', S%phase%name(:1))
if (ndata > 0) call rads_put_passinfo (S, P)
if (P%original /= '') e = e + nf90_put_att (ncid, nf90_global, 'original', P%original)
! Temporarily also create a 'log01' record, to support RADS3
l = index(P%original, rads_linefeed) - 1
if (l < 0) l = len_trim(P%original)
e = e + nf90_put_att (ncid, nf90_global, 'log01', datestamp()//' | '//trim(S%command)//': RAW data from '//P%original(:l))

if (e /= 0) call rads_error (S, rads_err_nc_create, 'Error writing global attributes to file', P)

call rads_put_history (S, P)
end subroutine rads_create_pass

!****if* rads/rads_put_passinfo
! SUMMARY
! Write pass info to RADS data file
!
! SYNOPSIS
subroutine rads_put_passinfo (S, P)
use netcdf
use rads_time
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
!
! PURPOSE
! This routine writes cycle and pass number, equator longitude and time,
! and start and end time of the pass to the global attributes in a RADS
! pass file.
!
! Equator longitude is written in the range 0-360 and is truncated after
! six decimals.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass structure
!****-------------------------------------------------------------------
integer(fourbyteint) :: e, ncid
character(len=26) :: date(3)
ncid = P%fileinfo(1)%ncid
e = 0
date = strf1985f ((/P%equator_time,P%start_time,P%end_time/))
e = e + &
nf90_put_att (ncid, nf90_global, 'cycle_number', P%cycle) + &
nf90_put_att (ncid, nf90_global, 'pass_number', P%pass) + &
nf90_put_att (ncid, nf90_global, 'equator_longitude', 1d-6 * nint(1d6 * modulo(P%equator_lon, 360d0))) + &
nf90_put_att (ncid, nf90_global, 'equator_time', date(1)) + &
nf90_put_att (ncid, nf90_global, 'first_meas_time', date(2)) + &
nf90_put_att (ncid, nf90_global, 'last_meas_time', date(3))
if (e /= 0) call rads_error (S, rads_err_nc_create, 'Error writing global attributes to file', P)
end subroutine rads_put_passinfo

!****if* rads/rads_put_history
! SUMMARY
! Write history to RADS data file
!
! SYNOPSIS
subroutine rads_put_history (S, P)
use netcdf
use rads_time
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
!
! PURPOSE
! This routine writes a datestamp and the command line to the global
! attribute 'history', followed by any possible previous history.
!
! When we are dealing with a pass file that was previously created by RADS3,
! this routine removes any previous log?? attributes (except log01) and
! adds the attribute 'original'.
!
! This command is called from rads_create_pass, so it is not required to
! call this command following rads_create_pass. However, when updating a
! file after rads_open_pass, it is necessary to call this routine to update
! the history.
!
! Note that P%history will not be updated.
!
! ARGUMENTS
! S        : Satellite/mission dependent structure
! P        : Pass structure
!****-------------------------------------------------------------------
integer(fourbyteint) :: e, i, ncid
character(len=8) :: log

! Make sure we are in define mode and that we can write
ncid = P%fileinfo(1)%ncid
if (nf90_redef (ncid) == nf90_eperm) call rads_error (S, rads_err_nc_put, 'File not opened for writing:', P)

! Write history attribute
if (associated(P%history)) then
	e = nf90_put_att (ncid, nf90_global, 'history', trim(P%history)//rads_linefeed//timestamp()//' : '//trim(S%command))
else
	e = nf90_put_att (ncid, nf90_global, 'history', timestamp()//' : '//trim(S%command))
endif
if (e /= 0) call rads_error (S, rads_err_nc_put, 'Error writing history attribute to file', P)

! Remove "log??" entries from RADS3 and make sure P%original is written
if (P%nlogs == 0) return
do i = 2,P%nlogs
	write (log, '("log",i2.2)') i
	e = nf90_del_att (ncid, nf90_global, log)
enddo
e = nf90_put_att (ncid, nf90_global, 'original', P%original)
end subroutine rads_put_history

subroutine rads_def_var_by_var_0d (S, P, var, nctype, scale_factor, add_offset, ndims, coordinates, varid)
use netcdf
use rads_netcdf
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
integer(fourbyteint), intent(in), optional :: nctype, ndims
real(eightbytereal), intent(in), optional :: scale_factor, add_offset
logical, intent(in), optional :: coordinates
integer(fourbyteint), intent(out), optional :: varid
type(rads_varinfo), pointer :: info
integer(fourbyteint) :: e, n, xtype, ncid, varid_
integer :: j=0, j0, j1! k
character(len=5) :: hz
S%error = rads_noerr
ncid = P%fileinfo(1)%ncid

! Get some information on dimensions and scale factors
info => var%info
if (present(nctype)) info%nctype = nctype
if (present(scale_factor)) info%scale_factor = scale_factor
if (present(add_offset)) info%add_offset = add_offset
if (present(ndims)) info%ndims = ndims

! Set the range of dimensions to be referenced
j0 = 1
j1 = info%ndims
!k = index(S%time%info%dataname, ' ') - 1
!if (nft(nf90_inq_dimid(ncid,S%time%info%dataname(:k),j1))) return
!j0 = j1 - info%ndims + 1
!if (j1 .eq. 1 .and. info%ndims .eq. 2) then
!	j0 = 1
!	j1 = 2
!endif
if (info%datatype == rads_type_dim) j0 = j1 ! Single dimension that is not primary

! Make sure we are in define mode and that we can write
if (nf90_redef (ncid) == nf90_eperm) call rads_error (S, rads_err_nc_put, 'File not opened for writing:', P)

! First check if the variable already exists
if (nff(nf90_inq_varid(ncid, var%name, varid_))) then
	e = nf90_inquire_variable (ncid, varid_, xtype=xtype, ndims=n)
	if (xtype /= info%nctype .or. n /= info%ndims) then
		call rads_error (S, rads_err_nc_var, &
			'Cannot redefine variable "'//trim(var%name)//'" with different type or dimension in file', P)
		return
	endif
! Define a constant
else if (info%ndims == 0) then
	if (nft(nf90_def_var(ncid, var%name, info%nctype, varid_))) then
		call rads_error (S, rads_err_nc_var, 'Error creating variable "'//trim(var%name)//'" in file', P)
		return
	endif
! Define a 1- or 2-dimensional variable
else if (nft(nf90_def_var(ncid, var%name, info%nctype, (/(j,j=j1,j0,-1)/), varid_))) then
	call rads_error (S, rads_err_nc_var, 'Error creating variable "'//trim(var%name)//'" in file', P)
	return
endif

! (Re)set the attributes; they may have been changed since last write
e = 0
if (info%datatype == rads_type_dim) then
	! Do not write _FillValue for dimension coordinates, like meas_ind
else if (info%nctype == nf90_int1) then
	e = e + nf90_put_att (ncid, varid_, '_FillValue', huge(0_onebyteint))
else if (info%nctype == nf90_int2) then
	e = e + nf90_put_att (ncid, varid_, '_FillValue', huge(0_twobyteint))
else if (info%nctype == nf90_int4) then
	e = e + nf90_put_att (ncid, varid_, '_FillValue', huge(0_fourbyteint))
endif
e = e + nf90_put_att (ncid, varid_, 'long_name', trim(info%long_name))
if (info%standard_name /= '') e = e + nf90_put_att (ncid, varid_, 'standard_name', trim(info%standard_name))
if (info%source /= '') e = e + nf90_put_att (ncid, varid_, 'source', trim(info%source))
if (info%units /= '') e = e + nf90_put_att (ncid, varid_, 'units', trim(info%units))
if (info%datatype == rads_type_flagmasks) then
	n = count_spaces (info%flag_meanings)
	if (info%nctype == nf90_int1) then
		e = e + nf90_put_att (ncid, varid_, 'flag_masks', int(flag_masks(0:n),onebyteint))
	else
		e = e + nf90_put_att (ncid, varid_, 'flag_masks', flag_masks(0:n))
	endif
	e = e + nf90_put_att (ncid, varid_, 'flag_meanings', info%flag_meanings)
else if (info%datatype == rads_type_flagvalues) then
	n = count_spaces (info%flag_meanings)
	if (info%nctype == nf90_int1) then
		e = e + nf90_put_att (ncid, varid_, 'flag_values', flag_values(0:n))
	else
		e = e + nf90_put_att (ncid, varid_, 'flag_values', int(flag_values(0:n),twobyteint))
	endif
	e = e + nf90_put_att (ncid, varid_, 'flag_meanings', info%flag_meanings)
endif
if (info%quality_flag /= '') e = e + nf90_put_att (ncid, varid_, 'quality_flag', info%quality_flag)
if (info%scale_factor /= 1d0) e = e + nf90_put_att (ncid, varid_, 'scale_factor', info%scale_factor)
if (info%add_offset /= 0d0)  e = e + nf90_put_att (ncid, varid_, 'add_offset', info%add_offset)
if (info%datatype >= rads_type_time .or. info%dataname(:1) == ':' .or. info%ndims < 1 .or. &
	(present(coordinates) .and. .not.coordinates)) then
	! Do not add coordinate attribute for some data types
else if (info%ndims > 1 .and. S%n_hz_output .and. P%n_hz > 1) then
	! For multi-Hz data: use 'lon_#hz lat_#hz'
	write (hz, '("_",i2.2,"hz")') P%n_hz
	e = e + nf90_put_att (ncid, varid_, 'coordinates', 'lon'//hz//' lat'//hz)
else
	! All other types: use 'lon lat'
	e = e + nf90_put_att (ncid, varid_, 'coordinates', 'lon lat')
endif
if (var%field(1) /= rads_nofield) e = e + nf90_put_att (ncid, varid_, 'field', var%field(1))
if (info%comment /= '') e = e + nf90_put_att (ncid, varid_, 'comment', info%comment)
if (e /= 0) call rads_error (S, rads_err_nc_var, &
	'Error writing attributes for variable "'//trim(var%name)//'" in file', P)
info%cycle = P%cycle
info%pass = P%pass
if (present(varid)) varid = varid_

contains

pure function count_spaces (string)
character(len=*), intent(in) :: string
integer(fourbyteint) :: count_spaces, i
count_spaces = 0
do i = 2,len_trim(string)-1
	if (string(i:i) == ' ') count_spaces = count_spaces + 1
enddo
end function count_spaces

end subroutine rads_def_var_by_var_0d

subroutine rads_def_var_by_var_1d (S, P, var, nctype, scale_factor, add_offset, ndims, coordinates, varid)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var(:)
integer(fourbyteint), intent(in), optional :: nctype, ndims
real(eightbytereal), intent(in), optional :: scale_factor, add_offset
logical, intent(in), optional :: coordinates
integer(fourbyteint), intent(out), optional :: varid
integer :: i
do i = 1,size(var)
	call rads_def_var_by_var_0d (S, P, var(i), nctype, scale_factor, add_offset, ndims, coordinates, varid)
	if (S%error /= rads_noerr) return
enddo
end subroutine rads_def_var_by_var_1d

subroutine rads_def_var_by_name (S, P, varname, nctype, scale_factor, add_offset, ndims, coordinates, varid)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
integer(fourbyteint), intent(in), optional :: nctype, ndims
real(eightbytereal), intent(in), optional :: scale_factor, add_offset
logical, intent(in), optional :: coordinates
integer(fourbyteint), intent(out), optional :: varid
type(rads_var), pointer :: var
var => rads_varptr (S, varname)
if (S%error /= rads_noerr) return
call rads_def_var_by_var_0d (S, P, var, nctype, scale_factor, add_offset, ndims, coordinates, varid)
if (S%error /= rads_noerr) return
end subroutine rads_def_var_by_name

subroutine rads_put_var_by_var_0d (S, P, var, data)
use netcdf
use rads_netcdf
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
real(eightbytereal), intent(in) :: data
integer(fourbyteint) :: e, ncid, varid
varid = rads_put_var_helper (S, P, var%name)
if (varid == 0) return
ncid = P%fileinfo(1)%ncid
select case (var%info%nctype)
case (nf90_int1)
	e = nf90_put_var (ncid, varid, nint1((data - var%info%add_offset) / var%info%scale_factor))
case (nf90_int2)
	e = nf90_put_var (ncid, varid, nint2((data - var%info%add_offset) / var%info%scale_factor))
case (nf90_int4)
	e = nf90_put_var (ncid, varid, nint4((data - var%info%add_offset) / var%info%scale_factor))
case default
	e = nf90_put_var (ncid, varid, (data - var%info%add_offset) / var%info%scale_factor)
end select
if (e /= 0) call rads_error (S, rads_err_nc_put, &
	'Error writing data for variable "'//trim(var%name)//'" to file', P)
end subroutine rads_put_var_by_var_0d

subroutine rads_put_var_by_name_0d (S, P, varname, data)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
real(eightbytereal), intent(in) :: data
type(rads_var), pointer :: var
var => rads_varptr (S, varname)
if (S%error /= rads_noerr) return
call rads_put_var_by_var_0d (S, P, var, data)
end subroutine rads_put_var_by_name_0d

subroutine rads_put_var_by_var_1d (S, P, var, data)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
real(eightbytereal), intent(in) :: data(:)
call rads_put_var_by_var_1d_start (S, P, var, data, (/1/))
end subroutine rads_put_var_by_var_1d

subroutine rads_put_var_by_name_1d (S, P, varname, data)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
real(eightbytereal), intent(in) :: data(:)
type(rads_var), pointer :: var
var => rads_varptr (S, varname)
if (S%error /= rads_noerr) return
call rads_put_var_by_var_1d_start (S, P, var, data, (/1/))
end subroutine rads_put_var_by_name_1d

subroutine rads_put_var_by_var_1d_start (S, P, var, data, start)
use netcdf
use rads_netcdf
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
real(eightbytereal), intent(in) :: data(:)
integer(fourbyteint), intent(in) :: start(:)
integer(fourbyteint) :: e, ncid, varid
varid = rads_put_var_helper (S, P, var%name)
if (varid == 0) return
ncid = P%fileinfo(1)%ncid
select case (var%info%nctype)
case (nf90_int1)
	e = nf90_put_var (ncid, varid, nint1((data - var%info%add_offset) / var%info%scale_factor), start)
case (nf90_int2)
	e = nf90_put_var (ncid, varid, nint2((data - var%info%add_offset) / var%info%scale_factor), start)
case (nf90_int4)
	e = nf90_put_var (ncid, varid, nint4((data - var%info%add_offset) / var%info%scale_factor), start)
case default
	e = nf90_put_var (ncid, varid, (data - var%info%add_offset) / var%info%scale_factor, start)
end select
if (e /= 0) call rads_error (S, rads_err_nc_put, &
	'Error writing data for variable "'//trim(var%name)//'" to file', P)
end subroutine rads_put_var_by_var_1d_start

subroutine rads_put_var_by_var_2d (S, P, var, data)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
real(eightbytereal), intent(in) :: data(:,:)
call rads_put_var_by_var_2d_start (S, P, var, data, (/1,1/))
end subroutine rads_put_var_by_var_2d

subroutine rads_put_var_by_name_2d (S, P, varname, data)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
real(eightbytereal), intent(in) :: data(:,:)
type(rads_var), pointer :: var
var => rads_varptr (S, varname)
if (S%error /= rads_noerr) return
call rads_put_var_by_var_2d_start (S, P, var, data, (/1,1/))
end subroutine rads_put_var_by_name_2d

subroutine rads_put_var_by_var_2d_start (S, P, var, data, start)
use netcdf
use rads_netcdf
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
real(eightbytereal), intent(in) :: data(:,:)
integer(fourbyteint), intent(in) :: start(:)
integer(fourbyteint) :: e, ncid, varid
varid = rads_put_var_helper (S, P, var%name)
if (varid == 0) return
ncid = P%fileinfo(1)%ncid
select case (var%info%nctype)
case (nf90_int1)
	e = nf90_put_var (ncid, varid, nint1((data - var%info%add_offset) / var%info%scale_factor), start)
case (nf90_int2)
	e = nf90_put_var (ncid, varid, nint2((data - var%info%add_offset) / var%info%scale_factor), start)
case (nf90_int4)
	e = nf90_put_var (ncid, varid, nint4((data - var%info%add_offset) / var%info%scale_factor), start)
case default
	e = nf90_put_var (ncid, varid, (data - var%info%add_offset) / var%info%scale_factor, start)
end select
if (e /= 0) call rads_error (S, rads_err_nc_put, &
	'Error writing data for variable "'//trim(var%name)//'" to file', P)
end subroutine rads_put_var_by_var_2d_start

subroutine rads_put_var_by_var_3d (S, P, var, data)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
real(eightbytereal), intent(in) :: data(:,:,:)
call rads_put_var_by_var_3d_start (S, P, var, data, (/1,1,1/))
end subroutine rads_put_var_by_var_3d

subroutine rads_put_var_by_name_3d (S, P, varname, data)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
real(eightbytereal), intent(in) :: data(:,:,:)
type(rads_var), pointer :: var
var => rads_varptr (S, varname)
if (S%error /= rads_noerr) return
call rads_put_var_by_var_3d_start (S, P, var, data, (/1,1,1/))
end subroutine rads_put_var_by_name_3d

subroutine rads_put_var_by_var_3d_start (S, P, var, data, start)
use netcdf
use rads_netcdf
use rads_misc
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
type(rads_var), intent(in) :: var
real(eightbytereal), intent(in) :: data(:,:,:)
integer(fourbyteint), intent(in) :: start(:)
integer(fourbyteint) :: e, ncid, varid
varid = rads_put_var_helper (S, P, var%name)
if (varid == 0) return
ncid = P%fileinfo(1)%ncid
select case (var%info%nctype)
case (nf90_int1)
	e = nf90_put_var (ncid, varid, nint1((data - var%info%add_offset) / var%info%scale_factor), start)
case (nf90_int2)
	e = nf90_put_var (ncid, varid, nint2((data - var%info%add_offset) / var%info%scale_factor), start)
case (nf90_int4)
	e = nf90_put_var (ncid, varid, nint4((data - var%info%add_offset) / var%info%scale_factor), start)
case default
	e = nf90_put_var (ncid, varid, (data - var%info%add_offset) / var%info%scale_factor, start)
end select
if (e /= 0) call rads_error (S, rads_err_nc_put, &
	'Error writing data for variable "'//trim(var%name)//'" to file', P)
end subroutine rads_put_var_by_var_3d_start

function rads_put_var_helper (S, P, varname)
use netcdf
use rads_netcdf
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
character(len=*), intent(in) :: varname
integer(fourbyteint) :: rads_put_var_helper
integer(fourbyteint) :: e, ncid, varid
S%error = rads_noerr
ncid = P%fileinfo(1)%ncid
e = nf90_enddef (ncid) ! Make sure to get out of define mode
if (.not.P%rw) then
	call rads_error (S, rads_err_nc_put, &
	'File not opened for writing variable "'//trim(varname)//'":', P)
	rads_put_var_helper = 0
else if (nft(nf90_inq_varid (ncid, varname, varid))) then
	call rads_error (S, rads_err_nc_var, 'No variable "'//trim(varname)//'" in file', P)
	rads_put_var_helper = 0
else
	rads_put_var_helper = varid ! Use varid determined above
endif
end function rads_put_var_helper

end module rads
