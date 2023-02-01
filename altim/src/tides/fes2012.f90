
module fes2012
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
    fes_success, &		! No error reporting
    fes_no_memory, &	! Not enough memory
    fes_netcdf_error, &	! netCDF error
    fes_io_error, &		! IO error
    fes_ini_error, &	! Invalid configuration file
    fes_no_data			! No data available in grids for the location asked
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
!-
end subroutine fes_delete

! void fes_set_nodal_time(FES info, double t);
!*fes_delete -- Free up space allocated by fes_init
!+
subroutine fes_set_nodal_time(info, time) bind(c)
import c_double, fes
type(fes), value :: info
real(c_double), value :: time
!
! This routine frees up memory allocated by a previous call to fes_init
!
! Argument:
!  info     : info initialised by fes_init
!-
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

! void fes_perror(FES fes);
subroutine fes_perror(info) bind(c)
import fes
type(fes), value :: info
end subroutine fes_perror
end interface

contains

!*fes_init -- Initialize FES tide model
!+
function fes_init (info, tide, mode, name)
integer(c_int) :: fes_init
type(fes), intent(out) :: info
integer(c_int), intent(in) :: tide, mode
character(len=*) :: name
!
! Initialise use of FES2012 tides.
! When mode = fes_mem: Allocate memory for FES tide modeling and
! read grids into memory. When mode = fes_io, read directly from
! grid.
! When tide = fes_tide: Compute ocean tides
! When tide = fes_radial: Compute radial (load) tides
! For consistency, however, use GOT4.8 load tides with FES2012.
!
! Input arguments:
!  tide     : Type of tide. Either fes_tide or fes_radial
!  mode     : Type of IO method. Either fes_io or fes_mem
!  name     : Type of FES model used. Current options are
!             'FES2012/all' or 'FES2012/reduced'
!
! Output argument:
!  info     : Handle initialised by fes_init
!-
character(len=256) :: path
call getenv ('ALTIM', path)
write (0,600) trim(path)//'/data/'//trim(name)//'.ini'
600 format ('(Loading FES tide: ',a,')')
fes_init = fes_new (info, tide, mode, trim(path)//'/data/'//trim(name)//'.ini'//c_null_char)
end function fes_init

!*fes_eval -- Compute tides according to FES model
!+
function fes_eval (info, time, lat, lon, h, h_lp)
integer(c_int) :: fes_eval
type(fes), intent(in) :: info
real(c_double), intent(in) :: time, lat, lon
real(c_double), intent(out) :: h, h_lp
!
! This routine makes the tidal predictions of ocean or rdial tide
! based on the FES2012 model.
!
! The input grids can be found in $ALTIM/data/<tide>.
!
! To initialize the computation, the function fes_init should be
! called first.
!
! Longitude and latitude are to be specified in degrees; time in UTC
! seconds since 1 Jan 1985. All predicted tides are output in meters.
! If the tide is requested in a point where it is not defined, NaN
! (Not-a-Number) is returned.
!
! The long-period tide already excludes the components that are captured
! in the dynamically computed short-period tides (Mm, Mf, Mtm, MSqm, Ssa)
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
!-
! $Log: fes2012.f90,v $
! Revision 1.4  2014/12/09 14:06:30  rads
! - Introduced fes_set_nodal_time
! - Set h and h_lp to NaN when interpolation fails
!
! Revision 1.3  2014/09/12 13:51:51  rads
! - Change to comment only
!
! Revision 1.2  2014/04/28 08:30:17  rads
! Manual added
!
!-----------------------------------------------------------------------
real(c_double), parameter :: sec50 = -1104537600d0, day = 86400d0
fes_eval = fes_core (info, lat, lon, (time - sec50) / day, h, h_lp)
if (fes_eval == 0) then
    h = h * 1d-2
    h_lp = h_lp * 1d-2
else
    h = 0d0
    h = h / h
    h_lp = h
endif
end function fes_eval

!*fes_error -- Print error produced by FES routines (if any)
!+
subroutine fes_error (info, fail)
type(fes), intent(in) :: info
logical, intent(in), optional :: fail
!
! Print error message to standard error if FES routines encountered any.
!
! Input arguments:
!  info     : Handle initialised by fes_init
!  fail     : (Optional) set to .true. if the process needs to stop on error.
!-
character(len=256) :: arg
if (fes_errno(info) == fes_success) return
call getarg (0,arg)
write (0,600,advance='no') trim(arg)
call fes_perror(info)
if (present(fail)) then
	if (fail) stop
endif
600 format (a,': ')
end subroutine

end module fes2012
