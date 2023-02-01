! This is an interface to the C code of libfes, the FES ocean tide prediction software.
!-
! $Log: libfes.f90,v $
! Revision 1.2  2017/10/12 12:48:41  rads
! - Update comments and error messages
!
! Revision 1.1  2017/09/07 10:25:55  rads
! Added libfes to src/tides directory
!
!------------------------------------------------------------------------

module libfes
use iso_c_binding

enum, bind(c)			! fes_enum_access
enumerator :: &
    fes_io = 0, &		! Direct access (Grids are not loaded into memory).
    fes_mem				! Memory access (Grids are loaded into memory).
end enum

enum, bind(c)			! fes_enum_tide_type
enumerator :: &
    fes_tide = 0, &		! Ocean tide
    fes_radial			! Radial tide
end enum

enum, bind(c)			! fes_enum_error
enumerator :: &
    fes_success = 0, &	! No error reporting
    fes_no_memory, &	! Not enough memory
    fes_netcdf_error, &	! netCDF error
    fes_io_error, &		! IO error
    fes_ini_error, &	! Invalid configuration file
    fes_no_data, &		! No data available in grids for the location asked
    fes_value_error		! Function receives an argument that has the right type but an inappropriate value
end enum

type, bind(c) :: fes
	type(c_ptr) ::	info	! void* FES
end type

interface
! int fes_new(FES* info, const fes_enum_tide_type tide, const fes_enum_access mode, const char* const path);
function fes_new(info, tide, mode, path) bind(c)
import c_ptr, c_int, c_char, fes
integer(c_int) :: fes_new
type(fes) :: info
integer(c_int), value :: tide, mode
character(kind=c_char) :: path(*)
end function fes_new

! void fes_delete(FES info);
!*fes_delete -- Free up space allocated by fes_init
!+
subroutine fes_delete(info) bind(c)
import fes
type(fes), value :: info
!
! This routine frees up memory allocated by a previous call to fes_init
!
! Argument:
!  info     : info initialised by fes_init
!-----------------------------------------------------------------------
end subroutine fes_delete

! void fes_min_number(FES info);
!*fes_min_number -- Get number of points used during the interpolation
!+
function fes_min_number(info) bind(c)
import c_int, fes
integer(c_int) :: fes_min_number
type(fes), value :: info
!
! This routine returns the number of points used when interpolating the FES
! tidal models.
!
! Argument:
!  info     : info initialised by fes_init
!-----------------------------------------------------------------------
end function fes_min_number

! void fes_set_buffer_size(FES info);
!*fes_set_buffer_size -- Set size of buffer in MiB
!+
function fes_set_buffer_size(info, size) bind(c)
import c_int, fes
integer(c_int) :: fes_set_buffer_size
type(fes), value :: info
integer(c_int), value :: size
!
! This function set the size of the memory buffer in MiBytes.
! The default size is 64 MiB.
!
! Argument:
!  info     : info initialised by fes_init
!-----------------------------------------------------------------------
end function fes_set_buffer_size

!*fes_set_nodal_time -- Set the nodal time
!+
subroutine fes_set_nodal_time(info, time) bind(c)
import c_double, fes
type(fes), value :: info
real(c_double), value :: time
!
! This routine sets the nodal time to a given value.
!
! Argument:
!  info     : info initialised by fes_init
!  time     : UTC time in seconds since 1 Jan 1985
!-----------------------------------------------------------------------
end subroutine fes_set_nodal_time

! int fes_core(FES info, const double lat, const double lon, const double time, double* h, double* h_lp);
function fes_core(info, lat, lon, time, h, h_lp) bind(c)
import c_int, c_ptr, c_double, fes
integer(c_int) :: fes_core
type(fes), value :: info
real(c_double), value :: lat, lon, time
real(c_double) :: h, h_lp
end function fes_core

! fes_enum_error fes_errno(FES fes);
function fes_errno(info) bind(c)
import c_int, fes
integer(c_int) :: fes_errno
type(fes), value :: info
end function fes_errno
end interface

contains

!*fes_init -- Initialize FES2012 or FES2014 tide model
!+
function fes_init (info, tide, mode, name)
integer(c_int) :: fes_init
type(fes), intent(out) :: info
integer(c_int), intent(in) :: tide, mode
character(len=*) :: name
!
! Initialise use of FES201x tides, which sets up a structure
! <info> of type(fes). Use this structure in subsequent calls to
! fes_eval to compute the tides and fes_delete to free up memory.
!
! The argument <mode> determines whether grids are stored in
! memory or read directly from the tide grids.
! Use the value fes_io (0) to read directly from grids, but
! keep this information in a buffer for later use.
! The maximum allotted memory can be set by calling
! the routine fes_set_buffer_size.
! Use the value fes_mem (1) to allocate sufficient memory
! and read all necessary tide grids into memory.
!
! The argument <tide> determines whether only ocean tides
! (value fes_tide = 0) or radial (load) tides (value fes_radial = 1)
! are to be computed.
!
! The argument <name> relates to the *.ini file in the directory
! $(ALTIM)/data. E.g. 'FES2012/all' will be based on the content
! of $(ALTIM)/data/FES2012/all.ini
!
! Input arguments:
!  tide     : Type of tide. Either fes_tide or fes_radial
!  mode     : Type of IO method. Either fes_io or fes_mem
!  name     : Type of FES model used. Current options are
!             for example 'FES2012/all', 'FES2012/reduced',
!             'FES2014/standard', 'FES2014/extrapolated'
!
! Output argument:
!  info     : Handle initialised by fes_init
!-----------------------------------------------------------------------
character(len=256) :: path
character(len=6), parameter :: tides(0:1) = (/'ocean ', 'radial'/)
character(len=8), parameter :: modes(0:1) = (/'opening ', 'loading '/)
call getenv ('ALTIM', path)
write (0,600) modes(mode), trim(tides(tide)), trim(path)//'/data/'//trim(name)//'.ini'
600 format ('(libfes: ',a, a,' tide: ',a,')')
fes_init = fes_new (info, tide, mode, trim(path)//'/data/'//trim(name)//'.ini'//c_null_char)
call fes_perror (info)
end function fes_init

!*fes_eval -- Compute tides according to FES2012 or FES2014 model
!+
function fes_eval (info, time, lat, lon, h, h_lp)
integer(c_int) :: fes_eval
type(fes), intent(in) :: info
real(c_double), intent(in) :: time, lat, lon
real(c_double), intent(out) :: h, h_lp
!
! This routine makes the tidal predictions of ocean or load tide
! based on the FES2012 or FES2014 model.
!
! The input grids can be found in $ALTIM/data/FES2012 and
! $ALTIM/data/FES2014
!
! To initialize the computation, the function fes_init should be
! called first.
!
! Longitude and latitude are to be specified in degrees; time in UTC
! seconds since 1 Jan 1985. All predicted tides are output in meters.
! If the tide is requested in a point where it is not defined, NaN
! (Not-a-Number) is returned.
!
! The long-period tide includes both the equilibrium and non-equilibrium
! portions coming from the model grids. To these the equilibrium
! tides of an additional 89 2nd and 17 3rd order waves are added. The
! 2nd order waves exclude the 5 dynamically computed long-period tides
! (Mm, Mf, Mtm, MSqm, Ssa) and their 12 sidelines in order to avoid
! duplicating these.
!
! Input arguments:
!  info     : Handle initialised by fes_init
!  utc      : UTC time in seconds since 1 Jan 1985
!  lat      : Latitude (degrees)
!  lon      : Longitude (degrees)
!
! Output arguments:
!  h        : Predicted short-period tide (m)
!  h_lp     : Predicted long-period tide (m)
!-----------------------------------------------------------------------
real(c_double), parameter :: sec50 = -1104537600d0, day = 86400d0
fes_eval = fes_core (info, lat, lon, (time - sec50) / day, h, h_lp)
h = h * 1d-2
h_lp = h_lp * 1d-2
end function fes_eval

!*fes_perror -- Print error produced by FES routines (if any)
!+
subroutine fes_perror (info, fail)
type(fes), intent(in) :: info
logical, intent(in), optional :: fail
!
! Print error message to standard error if FES routines encountered any.
!
! Input arguments:
!  info     : Handle initialised by fes_init
!  fail     : (Optional) set to .true. if the process needs to stop on error.
!-----------------------------------------------------------------------
character(len=256) :: arg
character(len=36), parameter :: err(6) = (/ &
	'Not enough memory                  ', &
	'netCDF error                       ', &
    'IO error                           ', &
    'Configuration file contains error  ', &
    'Tide is undefined                  ', &
    'Value error                        ' /)
if (fes_errno(info) == fes_success) return
call getarg (0,arg)
write (0,600) trim(arg), trim(err(fes_errno(info)))
if (present(fail)) then
	if (fail) stop
endif
600 format (a,': libfes: ',a)
end subroutine fes_perror

end module libfes
