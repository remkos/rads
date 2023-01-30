!-----------------------------------------------------------------------
!*nf90_message -- Print message to standard error
!+
subroutine nf90_message (string, ncid)
use netcdf
character(len=*), intent(in) :: string
integer(fourbyteint), optional :: ncid
!
! This routine prints a message <string> to standard error, while
! prepending the program name.
! Optionally, it adds the file name (based on the NetCDF ID) to which
! the error pertains.
!
! Examples:
! call nf90_message ('hello world')
!     output => program: hello word
! call nf90_message ('error occurred in', ncid)
!     output => program: error occurred in filename
!
! Arguments:
!  string   : Error message string
!  ncid     : Optional: NetCDF ID providing the filename
!-----------------------------------------------------------------------
character(len=160) :: progname, filename
integer :: i, l
call getarg (0, progname)
if (present(ncid)) then
	i = nf90_inq_path (ncid, l, filename)
	write (stderr, '(a,": ",a,1x,a)') trim(progname), trim(string), filename(:l)
else
	write (stderr, '(a,": ",a)') trim(progname), trim(string)
endif
end subroutine nf90_message
