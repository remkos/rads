!-----------------------------------------------------------------------
! Copyright (c) 2011-2021  Remko Scharroo
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

!*rads_add_gpd -- Add GPD along-track data to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! by adding GPD wet tropospheric correction.
!
! usage: rads_add_gpd [data-selectors] [options] < [list-of-GPD-files]
!-----------------------------------------------------------------------
program rads_add_gpd

use rads
use rads_misc
use rads_devel
use rads_netcdf
use netcdf

! Command line arguments

integer(fourbyteint) :: cyc, pass

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P
character(rads_cmdl) :: gpd_filenm
integer(fourbyteint) :: gpd_nr, gpd_done, ncid, dimid
real(eightbytereal), allocatable :: gpd_time(:), gpd_wet_tropo_cor(:), gpd_reference_height(:), gpd_source_flag(:)

! Other local variables

integer(fourbyteint) :: j, l, ios

! Initialise

call synopsis ('--head')
call rads_set_options (' all')
call rads_init (S)

! Check options

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('all')
		! Nothing to do
	end select
enddo

! Process all GPD input data files
! NOTE: This is not the most efficient way to do things in case this is run instide a SANDBOX.
!       In that case it would be futile to open and read every input file until later to discover that
!       there are no corresponding RADS pass files to process.
!       So later the opening/reading of the GPD input files should be put within the pass loop.

do
	read (*,'(a)',iostat = ios) gpd_filenm
	if (ios /= 0) exit
	l = index(gpd_filenm, '_c', .true.)
	cyc = -1
	read (gpd_filenm(l+2:l+4),*,iostat = ios) cyc
    if (cyc < S%cycles(1) .or. cyc > S%cycles(2)) cycle

	call log_string ('(' // trim(gpd_filenm) // ')')
	call nfs(nf90_open(gpd_filenm,nf90_nowrite,ncid))
	call nfs(nf90_inq_dimid(ncid, 'time_01', dimid))
	call nfs(nf90_inquire_dimension(ncid, dimid, len=gpd_nr))
	allocate (gpd_time(gpd_nr), gpd_wet_tropo_cor(gpd_nr), gpd_reference_height(gpd_nr), gpd_source_flag(gpd_nr))
	call get_var (ncid, 'time_01 473299200 ADD', gpd_time)
	call get_var (ncid, 'gpd_wet_tropo_cor_01', gpd_wet_tropo_cor)
	call get_var (ncid, 'gpd_reference_height_01', gpd_reference_height)
	call get_var (ncid, 'gpd_source_flag_01', gpd_source_flag)
	call nfs(nf90_close(ncid))
	gpd_done = 0

	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		gpd_done = gpd_done + P%ndata
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Add GPD along-track data to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  --all                     (Has no impact)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i
real(eightbytereal) :: time(n)

call log_pass (P)

! Get time and check against GPD time

call rads_get_var (S, P, 'time', time, .true.)
i = count (time /= gpd_time(gpd_done+1:gpd_done+n))
if (i > 0) then
	call log_string ('Error: time_gpd /= time')
	call log_records (0)
	return
endif

! Define the output variables

call rads_put_history (S, P)

call rads_def_var (S, P, 'gpd_wet_tropo_cor')
call rads_def_var (S, P, 'gpd_reference_height')
call rads_def_var (S, P, 'gpd_source_flag')

! Write output variable

call rads_put_var (S, P, 'gpd_wet_tropo_cor', gpd_wet_tropo_cor)
call rads_put_var (S, P, 'gpd_reference_height', gpd_reference_height)
call rads_put_var (S, P, 'gpd_source_flag', gpd_source_flag)

call log_records (n)
end subroutine process_pass

end program rads_add_gpd
