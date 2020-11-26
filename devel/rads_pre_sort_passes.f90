!-----------------------------------------------------------------------
! Copyright (c) 2011-2020  Remko Scharroo
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

!*rads_pre_sort_passes -- Sort (combine and split) Sentinel-6 or Jason GDR-F into pass files
!
! Read Sentinel-6 standard or reduced granules or orbits and
! combine them (and split them) into pass files.
! This works also for Jason-3 GDR-F OGDR files in GDR-F format.
! The input file names are read from
! standard input. The individual pass files will be named
! <destdir>/cCCC/SSS_*_CCC_PPP_*.nc, where SSS is the satellite abbreviation
! (S6A, JA3), CCC is the cycle number and PPP the pass number.
! The directory <destdir>/cCCC will be created if needed.
!
! This program relies on the availability of S6 or JA3 ORF files.
!-----------------------------------------------------------------------
program rads_pre_sort_passes

use rads
use rads_misc
use rads_time
use rads_netcdf
use rads_devel_misc
use netcdf

! Struct for orbit info

type(orfinfo) :: orf(200000)

! Scruct to store input file information

type :: fileinfo
	integer(fourbyteint) :: ncid(0:3), rec0, rec1, nrec, type, cycle, pass
	real(eightbytereal) :: time0, time1, lat0, lat1, lon0, lon1
	character(len=rads_cmdl) :: filenm
end type
type(fileinfo) :: fin(20)

! Group names

character(len=7) :: grpnm(3) = (/ 'data_01', 'ku     ', 'c      '/)

! General variables

character(len=rads_cmdl) :: arg, filenm, dimnm, destdir, product_name, xref_orbit_data
character(len=rads_strl) :: exclude_list = ','
character(len=26) :: date(3)
character(len=15) :: newer_than = '00000000T000000'
character(len=3) :: sat
integer(fourbyteint), parameter :: mpass = 254 * 500
real(eightbytereal), parameter :: sec2000 = 473299200d0
integer(fourbyteint) :: i0, i, j, ngrps = 3, ncid1(0:3), nrec, ios, varid, in_max = huge(fourbyteint), &
	nfile = 0, orbit_type, ipass = 1, ipass0 = 0, verbose_level = 3, cycle_number, pass_number
real(eightbytereal), allocatable :: time(:), lat(:), lon(:)
real(eightbytereal) :: last_time = 0
logical :: first = .true., check_pass = .false.

! Print description, if requested

if (iargc() < 1) then
	call getarg(0,arg)
	write (*,1300) trim(arg), trim(arg)
	stop
endif
1300 format (a,' -- Sort (combine/split) Sentinel-6 or Jason GDR-F files into pass files'// &
'syntax: ',a,' [options] destdir < list'//'where'/ &
'  destdir           : Destination directory (appends c???/*.nc)'/ &
'  list              : List of input files names'// &
'where [options] are:' / &
'  -c                : Check cycle and pass number'/ &
'  -mMAXREC          : Maximum number of records allowed at input' / &
'  -pDATE            : Skip data before given processing date (yymmddThhmmss)' / &
'  -vLEVEL           : Specify verbosity level: 0 = quiet, 1 = only list created files,'/ &
'                      2 = also list unchanged files, 3 = also list input files,'/ &
'                      4 = also list skipped input files (default)'/ &
'                      5 = debug information' / &
'  -xVAR1[,VAR2,...] : Exclude variable(s) from copying')

! First determine filetype

read (*,'(a)',iostat=ios) filenm
if (ios /= 0) stop

i = index(filenm,'/',.true.) + 1
if (filenm(i:i+1) == 'JA' .and. filenm(i+10:i+10) == 'f') then
	! Jason GDR-F
else if (filenm(i:i+1) == 'S6' .and. filenm(i+4:i+8) == 'P4_2_') then
	! Sentinel-6 Level-2
	if (filenm(i+9:i+10) == 'HR') ngrps = 2
else
	call rads_exit ('Wrong filetype')
endif
sat = filenm(i:i+2)

! Get ORF file

call read_orf (sat, orf)

! Read options and destination directory

do i = 1,iargc()
	call getarg (i,arg)
	if (arg(:2) == '-c') then
		check_pass = .true.
	else if (arg(:2) == '-i') then
		! Ignore old argument from ogdrsplit
	else if (arg(:2) == '-m') then
		read (arg(3:), *, iostat=ios) in_max
	else if (arg(:2) == '-p') then
		newer_than = arg(3:)
	else if (arg(:2) == '-v') then
		read (arg(3:), *, iostat=ios) verbose_level
	else if (arg == '-x2') then
		! Ignore old argument from ogdrsplit
	else if (arg(:2) == '-x') then
		exclude_list = trim(exclude_list) // arg(3:len_trim(arg)) // ','
	else
		destdir = arg
	endif
enddo

! Cycle through all input files

do

! Read the next file name

	if (.not.first) then
		read (*,'(a)',iostat=ios) filenm
		if (ios /= 0) exit
	endif
	i = index(filenm,'/',.true.) + 1
	product_name = filenm(i:)
	first = .false.

! Open the input file

	if (nft(netcdf_open(filenm,nf90_nowrite,ncid1))) then
		call rads_message ('Error while opening file: '//filenm)
		cycle
	endif

! Read global attributes

	if (sat(:2) == 'S6') call nfs(nf90_get_att(ncid1(0),nf90_global,'product_name',product_name))
	call nfs(nf90_get_att(ncid1(0),nf90_global,'cycle_number',cycle_number))
	call nfs(nf90_get_att(ncid1(0),nf90_global,'pass_number',pass_number))

! Read the time dimension

	call nfs(nf90_inquire_dimension(ncid1(1),1,dimnm,nrec))
	if (dimnm /= 'time') call rads_exit ('Error reading time dimension')
	if (nrec > in_max) then
		call rads_message ('Too many measurements in input file, skipped: '//filenm)
		call nfs(nf90_close(ncid1(0)))
		cycle
	endif
	if (allocated(time)) deallocate (time, lat, lon)
	allocate (time(nrec),lat(nrec),lon(nrec))
	call nfs(nf90_inq_varid(ncid1(1),'time',varid))
	call nfs(nf90_get_var(ncid1(1),varid,time))

! Check processing time

	if (product_name(51:65) < newer_than) then
		call nfs(nf90_close(ncid1(0)))
		call write_line (4, '.... old', 1, nrec, nrec, time(1), time(nrec), '<', filenm)
		cycle
	endif

! Read latitude and longitude

	call nfs(nf90_inq_varid(ncid1(1),'latitude',varid))
	call nfs(nf90_get_var(ncid1(1),varid,lat))
	lat = lat * 1d-6
	call nfs(nf90_inq_varid(ncid1(1),'longitude',varid))
	call nfs(nf90_get_var(ncid1(1),varid,lon))
	lon = lon * 1d-6

! First advance to beyond the last time tag

	do i0 = 1, nrec
		if (time(i0) > last_time) exit
	enddo
	if (i0 > 1) call write_line (4, '... skip', 1, i0-1, nrec, time(1), time(i0-1), '<', filenm)

! If none of the records are after last_time, skip the whole input file

	if (i0 > nrec) then
		call nfs(nf90_close(ncid1(0)))
		cycle
	endif

! Check if the first input point is in a new pass

	call which_pass (time(i0))
	if (ipass /= ipass0) call write_output
	last_time = time(nrec)

! Two reasons to split a file into two pieces:

	do i = max(2,i0),nrec
		call which_pass (time(i))
		if (time(i) == time(i-1)) then
! - Duplicated measurements within a single file
			call fill_fin (i0, i-2)
			call write_line (4, '... skip', i-1, i-1, nrec, time(i-1), time(i-1), '<', filenm)
			i0 = i
		else if (ipass /= ipass0) then
! - Split pass when entering into a new pass
			call fill_fin (i0, i-1)
			call write_output
			i0 = i
		endif
	enddo
! - Register the remaining bit
	call fill_fin (i0, nrec)

! Deallocate time and location arrays

	deallocate (time, lat, lon)

enddo

! Dump the remainder of the input files to output

call write_output

contains

!***********************************************************************
! Fill the array of file information

subroutine fill_fin (i0, i1)
integer, intent(in) :: i0, i1
nfile = nfile + 1
if (nfile > 20) call rads_exit ('Number of granules too large (> 20)')
fin(nfile)%ncid = ncid1
fin(nfile)%nrec = nrec
fin(nfile)%filenm = filenm
fin(nfile)%type = orbit_type
fin(nfile)%cycle = cycle_number
fin(nfile)%pass = pass_number
fin(nfile)%rec0 = i0
fin(nfile)%rec1 = i1
fin(nfile)%time0 = time(i0)
fin(nfile)%time1 = time(i1)
fin(nfile)%lat0 = lat(i0)
fin(nfile)%lat1 = lat(i1)
fin(nfile)%lon0 = lon(i0)
fin(nfile)%lon1 = lon(i1)
call write_line (3, '.. input', i0, i1, nrec, time(i0), time(i1), '<', filenm)
if (verbose_level >= 5) write (*,*) "ncid1, nrec =", ncid1, nrec
end subroutine fill_fin

!***********************************************************************
! Determine the corresponding record from ORF

subroutine which_pass (time)
real(eightbytereal), intent(in) :: time
do while (time < orf(ipass)%starttime)
	ipass = ipass - 1
	if (ipass < 1) call rads_exit ('Times are beyond limits of ORF file')
enddo
do while (time > orf(ipass+1)%starttime)
	ipass = ipass + 1
	if (orf(ipass)%cycle < 0) call rads_exit ('Times are beyond limits of ORF file')
enddo
end subroutine which_pass

!***********************************************************************
! Write whatever has been buffered so far to an output file

subroutine write_output
integer(fourbyteint) :: nrec, ncid1(0:3), ncid2(0:3), varid1, varid2, i, nout, &
	cycle_number, pass_number, absolute_pass_number, absolute_rev_number
real(eightbytereal) :: equator_time, equator_longitude, x
character(len=rads_naml) :: dirnm, prdnm, outnm, varnm
logical :: exist

! Retrieve the pass variables

cycle_number = orf(ipass0)%cycle
pass_number = orf(ipass0)%pass

absolute_pass_number = orf(ipass0)%abs_pass
absolute_rev_number = orf(ipass0)%abs_rev

equator_time = orf(ipass0)%eqtime + sec2000
equator_longitude = orf(ipass0)%eqlon
ipass0 = ipass

! How many records are buffered for output?
! Skip if there is nothing left

if (nfile == 0 .or. ipass0 == 0) return
nrec = sum(fin(1:nfile)%rec1 - fin(1:nfile)%rec0 + 1)

! Open the output file. Make directory if needed.

605 format (a,'/c',i3.3)
610 format (a,'P',i3.3,'_',i3.3,'.nc') ! JA? format
612 format (a,'RED_',a,i3.3,'_',i3.3,a,'.nc') ! S6? format

write (dirnm,605) trim(destdir),cycle_number
inquire (file=dirnm,exist=exist)
if (.not.exist) call system('mkdir -p '//dirnm)
if (sat(:2) == 'JA') then
	write (prdnm,610) product_name(:11),cycle_number,pass_number
else
	write (prdnm,612) product_name(:13),product_name(92:95),cycle_number,pass_number,product_name(95:98)
endif
outnm = trim(dirnm) // '/' // trim(prdnm)
inquire (file=outnm,exist=exist)

! If exist, then keep the file if the buffer is smaller or equal in size
! If it is larger, delete the existing file

if (exist) then
	call nfs(netcdf_open(outnm,nf90_write,ncid2))
	call nfs(nf90_inquire_dimension(ncid2(1),1,len=nout))
	if (verbose_level >= 5) write (*,*) "nout =",nout
	if (sat(:2) == 'S6') then
		call nfs(nf90_get_att(ncid2(0),nf90_global,'first_measurement_time',date(1)))
		call nfs(nf90_get_att(ncid2(0),nf90_global,'last_measurement_time',date(2)))
	else
		call nfs(nf90_get_att(ncid2(0),nf90_global,'first_meas_time',date(1)))
		call nfs(nf90_get_att(ncid2(0),nf90_global,'last_meas_time',date(2)))
	endif
	call nfs(nf90_get_att(ncid2(0),nf90_global,'equator_time',date(3)))
	x = strp1985f (date(3))
	if (abs(x - equator_time) > 0.5d-3) then
		date(3) = strf1985f(equator_time)
		call nfs(nf90_put_att(ncid2(0),nf90_global,'equator_time',date(3)))
		if (verbose_level >= 1) write (*,620) 'Updating', 1, nout, nout, date(1:2), '>', trim(outnm)
		if (verbose_level >= 5) write (*,*) "ncid2, nout =", ncid2, nout
	endif
	call nfs(nf90_close(ncid2(0)))
	if (nrec <= nout) then
		do i = 1, nfile
			! Release NetCDF file when reaching end
			if (fin(i)%rec1 == fin(i)%nrec) call nfs(nf90_close(fin(i)%ncid(0)))
		enddo
		if (verbose_level >= 2) write (*,620) 'Keeping ', 1, nout, nout, date(1:2), '>', trim(outnm)
		if (verbose_level >= 5) write (*,*) "ncid2, nout =", ncid2, nout
		nfile = 0
		return
	endif
	call system('rm -f '//outnm)
endif
620 format (a,' : ',3i6,' : ',a,' - ',a,1x,a1,1x,a)

! Create a new file

call nfs(netcdf_create(outnm,nf90_write+nf90_netcdf4,ncid2))
!call nfs(nf90_set_fill(ncid2,nf90_nofill,i))

! Create the time dimension

call nfs(nf90_def_dim(ncid2(1),'time',nrec,i))

! Copy all the variable definitions and attributes

do i = 0,ngrps
	call netcdf_copy_defs(fin(1)%ncid(i),ncid2(i))
enddo

! Overwrite some attributes and product name

call write_line (1, 'Creating', 1, nrec, nrec, fin(1)%time0, fin(nfile)%time1,'>', outnm)
if (verbose_level >= 5) write (*,*) "ncid2, nrec =", ncid2, nrec
date(3) = strf1985f(equator_time)

call nfs(nf90_put_att(ncid2(0),nf90_global,'product_name',prdnm))
call nfs(nf90_put_att(ncid2(0),nf90_global,'cycle_number',cycle_number))
call nfs(nf90_put_att(ncid2(0),nf90_global,'pass_number',pass_number))
call nfs(nf90_put_att(ncid2(0),nf90_global,'absolute_rev_number',absolute_rev_number))
if (sat(:2) == 'S6') then
	call nfs(nf90_put_att(ncid2(0),nf90_global,'first_measurement_time',date(1)))
	call nfs(nf90_put_att(ncid2(0),nf90_global,'last_measurement_time',date(2)))
	call nfs(nf90_put_att(ncid2(0),nf90_global,'first_measurement_latitude',fin(1)%lat0))
	call nfs(nf90_put_att(ncid2(0),nf90_global,'last_measurement_latitude',fin(nfile)%lat1))
	call nfs(nf90_put_att(ncid2(0),nf90_global,'first_measurement_longitude',fin(1)%lon0))
	call nfs(nf90_put_att(ncid2(0),nf90_global,'last_measurement_longitude',fin(nfile)%lon1))
else
	call nfs(nf90_put_att(ncid2(0),nf90_global,'absolute_pass_number',absolute_pass_number))
	call nfs(nf90_put_att(ncid2(0),nf90_global,'first_meas_time',date(1)))
	call nfs(nf90_put_att(ncid2(0),nf90_global,'last_meas_time',date(2)))
endif
call nfs(nf90_put_att(ncid2(0),nf90_global,'equator_time',date(3)))
call nfs(nf90_put_att(ncid2(0),nf90_global,'equator_longitude',equator_longitude))
call nfs(nf90_enddef(ncid2(0)))

! If requested, check if cycles and pass numbers agree
if (check_pass) then
	do i = 1, nfile
		if (fin(i)%cycle /= cycle_number) write (*,630) 'cycle', cycle_number, fin(i)%cycle, trim(fin(i)%filenm)
		if (fin(i)%pass /= pass_number) write (*,630) ' pass', pass_number, fin(i)%pass, trim(fin(i)%filenm)
	enddo
endif
630 format ('Warning: ',a,' number; output:',i4,'; input:',i4,' from ',a)

! Copy all data elements

nout = 0
do i = 1,nfile
	do j = 0,ngrps
		call netcdf_copy_vars (fin(i)%ncid(j), fin(i)%rec0, fin(i)%rec1, ncid2(j), nout+1)
	enddo
	nout = nout + fin(i)%rec1 - fin(i)%rec0 + 1
	! Release NetCDF file when reaching end
	if (fin(i)%rec1 == fin(i)%nrec) call nfs(nf90_close(fin(i)%ncid(0)))
enddo

! Close output

call nfs(nf90_close(ncid2(0)))
nfile = 0

end subroutine write_output

!***********************************************************************
! Open NetCDF file

integer function netcdf_open (filenm, flags, ncid)
character(len=*), intent(in) :: filenm
integer(fourbyteint), intent(in) :: flags
integer(fourbyteint), intent(out) :: ncid(0:3)
integer :: i

ncid = 0
netcdf_open = nf90_open(filenm,flags,ncid(0))
if (nft(netcdf_open)) return
do i = 1,ngrps
	netcdf_open = nf90_inq_grp_ncid(ncid(i/2),grpnm(i),ncid(i))
enddo
end function netcdf_open

!***********************************************************************
! Create NetCDF file

integer function netcdf_create (filenm, flags, ncid)
character(len=*), intent(in) :: filenm
integer(fourbyteint), intent(in) :: flags
integer(fourbyteint), intent(out) :: ncid(0:3)
integer :: i

ncid = 0
netcdf_create = nf90_create(filenm,flags,ncid(0))
if (nft(netcdf_create)) return
do i = 1,ngrps
	netcdf_create = nf90_def_grp(ncid(i/2),grpnm(i),ncid(i))
enddo
end function netcdf_create

!***********************************************************************
! Copy all the variable definitions and attributes

subroutine netcdf_copy_defs(ncid1, ncid2)
integer(fourbyteint), intent(in) :: ncid1, ncid2
integer(fourbyteint) :: i, nvars, natts, ndims, xtype, dimids(2), varid1, varid2
character(len=rads_naml) :: attnm, varnm

call nfs(nf90_inquire(ncid1,nvariables=nvars,nattributes=natts))
if (verbose_level >= 5) write (*,*) "copy_defs:", ncid1, ncid2, nvars, natts
varid2 = 0
do varid1 = 0, nvars
	if (varid1 > 0) then
		call nfs(nf90_inquire_variable(ncid1,varid1,varnm,xtype,ndims,dimids,natts))
		! Skip all listed excluded variables, all 20-Hz variables
		! Skip also all uint variables because the NetCDF library doesn't work with them
		! (they occur only in enhanced_measurements.nc)
		if (excluded(varnm) .or. ndims > 1 .or. dimids(1) > 1 .or. xtype == nf90_uint) cycle
		call nfs(nf90_def_var(ncid2,varnm,xtype,dimids(1:ndims),varid2))
	endif
	do i = 1,natts
		call nfs(nf90_inq_attname(ncid1,varid1,i,attnm))
		call nfs(nf90_copy_att(ncid1,varid1,attnm,ncid2,varid2))
	enddo
enddo
end subroutine netcdf_copy_defs

!***********************************************************************
! Copy all variables

subroutine netcdf_copy_vars(ncid1, rec0, rec1, ncid2, rec2)
integer(fourbyteint), intent(in) :: ncid1, rec0, rec1, ncid2, rec2
integer(fourbyteint) :: nrec, nvars, xtype, varid1, varid2, idxin(2)=1, idxut(2)=1
real(eightbytereal), allocatable :: darr1(:)
integer(fourbyteint), allocatable :: iarr1(:)
character(len=rads_naml) :: varnm

call nfs(nf90_inquire(ncid2,nvariables=nvars))
nrec = rec1 - rec0 + 1
if (verbose_level >= 5) write (*,*) "copy_vars:", ncid1, ncid2, nvars, nrec, rec2
if (nvars == 0) return

if (nrec == 0) stop "nrec == 0"
allocate (darr1(nrec),iarr1(nrec))
idxin(2) = rec0
idxut(2) = rec2

do varid2 = 1,nvars
	call nfs(nf90_inquire_variable(ncid2,varid2,varnm,xtype))
	call nfs(nf90_inq_varid(ncid1,varnm,varid1))
	if (xtype == nf90_double) then
		call nfs(nf90_get_var(ncid1,varid1,darr1,idxin(2:2)))
		call nfs(nf90_put_var(ncid2,varid2,darr1,idxut(2:2)))
	else
		call nfs(nf90_get_var(ncid1,varid1,iarr1,idxin(2:2)))
		call nfs(nf90_put_var(ncid2,varid2,iarr1,idxut(2:2)))
	endif
enddo
deallocate (darr1,iarr1)
end subroutine netcdf_copy_vars


subroutine write_line (verbose, word, rec0, rec1, nrec, time0, time1, dir, filenm)
character(len=*), intent(in) :: word, dir, filenm
integer, intent(in) :: verbose, rec0, rec1, nrec
real(eightbytereal), intent(in) :: time0, time1
if (verbose_level < verbose) return
date(1) = strf1985f(time0+sec2000)
date(2) = strf1985f(time1+sec2000)
write (*,620) word, rec0, rec1, nrec, date(1:2), dir, trim(filenm)
620 format (a,' : ',3i6,' : ',a,' - ',a,1x,a1,1x,a)
end subroutine write_line

function excluded (varnm)
character(len=*), intent(in) :: varnm
logical :: excluded
excluded = (index(exclude_list,','//trim(varnm)//',') > 0)
end function excluded

end program rads_pre_sort_passes
