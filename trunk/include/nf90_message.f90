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
! The optional argument <ncid> is ignored.
! This is done to give some sort of backwards compatibility with
! netCDF versions that do not support nf90_inq_path.
!
! Examples:
! call nf90_message ('hello world')
!     output => program: hello word
! call nf90_message ('error occurred in file', ncid)
!     output => program: error occurred in file
!
! Arguments:
!  string   : Error message string
!  ncid     : Optional: netCDF ID (not actually used)
!-----------------------------------------------------------------------
character(len=160) :: progname
call getarg (0, progname)
write (stderr, '(a,": ",a)') trim(progname), trim(string)
end subroutine nf90_message
