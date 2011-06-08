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

!*radsxosel -- RADS crossover extender
!+
program radsxosel

! This program adds RADS data to a file that has previously
! been created with radsxogen.
!
! The output is a new netCDF file containing all the information from
! the input file plus the data fields indicated on the sel= option.
!-
! Created by Remko Scharroo, Altimetrics LLC, based in part on previous
! programs, max and max2, developed at DEOS.
!-----------------------------------------------------------------------
use rads
use netcdf
use rads_netcdf
integer(fourbyteint), parameter :: msat = 10
real(eightbytereal) :: dt(msat) = -1d0
type(rads_sat) :: S
integer(fourbyteint) :: i, j, nsat = 0, nsel = 0, reject = -1, ios, debug, nrxo
character(len=80) :: arg, filename = 'radsxogen.nc'
type :: stat
	integer :: test, shallow, xdat, ins, gap, xdt, xout
endtype
type(stat) :: nr
real(eightbytereal), allocatable :: lat(:), lon(:), time(:,:)
integer(twobyteint), allocatable :: satid(:,:)
integer(fourbyteint), allocatable :: cycle(:,:), pass(:,:), idx(:,:)
integer(fourbyteint) :: ncid

! Initialize RADS or issue help
call synopsis
S%sat = '' ! Initialize blank
S%error = rads_noerr

! Start with this-is message
call rads_version ('$Revision$')

! Get filename
call getarg(iargc(), arg)
if (arg(:1) /= '-' .and. index(arg,'=') == 0) filename = arg

! Scan command line arguments
do i = 1,iargc()
	call getarg(i,arg)
	if (arg(:4) == 'out=') then
		outname = arg(5:)
	else if (arg(:2) == '-rn') then
		reject = -2
	else if (arg(:2) == '-r') then
		reject = 0
		read (arg(3:),*,iostat=ios) reject
	else if (arg(:2) == '-f') then
		freeform = .true.
	else if (arg(:3) == '-sp') then
		stat_mode = pass_stat
	else if (arg(:3) == '-sc') then
		stat_mode = cycle_stat
	else if (arg(:7) == 'maxrec=') then
		read (arg(8:),*) nselmax
	else if (arg(:5) == 'step=') then
		read (arg(6:),*) step
	else if (arg(:5) == 'list=') then
		listunit = getlun()
		open (listunit, file=arg(6:), status='replace')
	endif
enddo

! Read the whole file into memory
call nfs (nf90_open (filename, nf90_write, ncid))
call nfs (nf90_inquire_dimension (ncid, 2, len=nrxo))
allocate (lat(nrxo), lon(nrxo), time(2,nrxo), satid(2,nrxo), cycle(2,nrxo), pass(2,nrxo), idx(2,nrxo))
call nfs (nf90_get_var (ncid, 1, lat))
call nfs (nf90_get_var (ncid, 2, lon))
call nfs (nf90_get_var (ncid, 3, time))
call nfs (nf90_get_var (ncid, 4, satid))
call nfs (nf90_get_var (ncid, 5, cycle))
call nfs (nf90_get_var (ncid, 6, pass))
call nfs (nf90_get_var (ncid, 7, idx))

600 format (/ &
'File name             : ',a/ &
'Number of xovers read : ',i9)
write (*, 600) trim(filename), nrxo

call nfs (nf90_close (ncid))

contains

!***********************************************************************

subroutine synopsis
if (rads_version('$Revision$','Add RADS data to crossover file')) return
call rads_synopsis()
write (0,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r#               : reject xovers if data item number # on sel= specifier is NaN'/ &
'                      (default: reject if SLA field is NaN)'/ &
'  -r0, -r           : do not reject xovers with NaN values'/ &
'  -rn               : reject xovers if any value is NaN')
stop
end subroutine synopsis

!***********************************************************************

end program radsxosel
