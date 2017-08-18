!-----------------------------------------------------------------------
! Copyright (c) 2011-2017  Remko Scharroo
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

!*s3combine -- Combine (and split) Sentinel-3 files into pass files
!
! Read Sentinel-3 standard_measurements.nc or reduced_measurements.nc
! granules and combine them (and split them) into pass files.
! The input file names are read from
! standard input. The individual pass files will be named
! <destdir>/cCCC/S3A_*_CCC_PPP_*.nc, where CCC is the cycle number and
! PPP the pass number. The directory <destdir>/cCCC will be created if needed.
!
! This program does not rely on ORF files, like ogdrsplit does.
!-----------------------------------------------------------------------
program s3combine

use rads
use rads_misc
use rads_netcdf
use typesizes
use netcdf

! Scruct to store input file information

type :: fileinfo
	integer(fourbyteint) :: ncid, rec0, rec1, nrec, type
	real(eightbytereal) :: time0, time1, lat0, lat1, lon0, lon1
	character(len=rads_cmdl) :: filenm
end type
type(fileinfo) :: fin(20)

! Struct to store the records from ORF file

type :: orfinfo
	integer(fourbyteint) :: cycle, pass
	real(eightbytereal) :: starttime, eqtime, eqlon
end type
type(orfinfo) :: orf(200000)

! General variables

character(len=rads_cmdl) :: arg, filenm, dimnm, destdir, product_name, xref_orbit_data
character(len=rads_strl) :: exclude_list = ','
character(len=26) :: date(3)
character(len=15) :: newer_than = '00000000T000000'
integer(fourbyteint), parameter :: mpass = 254 * 500
real(eightbytereal), parameter :: sec2000 = 473299200d0
integer(fourbyteint) :: i0, i, ncid1, nrec, ios, varid, in_max = huge(fourbyteint), nr_passes = 770, &
	nfile = 0, orbit_type, ipass = 1, ipass0 = 0
real(eightbytereal), allocatable :: time(:), lat(:), lon(:)
real(eightbytereal) :: last_time = 0
logical :: first = .true.

! Print description, if requested

if (iargc() < 1) then
	write (*,1300)
	stop
endif
1300 format ('s3combine -- Combine/split Sentinel-3 files into pass files'// &
'syntax: s3combine [options] destdir < list'//'where'/ &
'  destdir           : Destination directory (appends c???/*.nc)'/ &
'  list              : List of input files names'// &
'where [options] are:' / &
'  -mMAXREC          : Maximum number of records allowed at input' / &
'  -pDATE            : Skip data before given processing date (yymmddThhmmss)' / &
'  -xVAR1[,VAR2,...] : Exclude variable(s) from copying')

! First determine filetype

read (*,'(a)',iostat=ios) filenm
if (ios /= 0) stop
i = index(filenm,'.SEN3')
if (i == 0) stop 'Wrong filetype'

! Get ORF file

call read_orf (filenm(i-94:i-92), orf)

! Read options and destination directory

do i = 1,iargc()
	call getarg (i,arg)
	if (arg(:2) == '-m') then
		read (arg(3:), *, iostat=ios) in_max
	else if (arg(:2) == '-p') then
		newer_than = arg(3:)
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
	first = .false.

! Open the input file

	if (nft(nf90_open(filenm,nf90_nowrite,ncid1))) then
		call rads_message ('Error while opening file: '//filenm)
		cycle
	endif

! Read global attributes

	call nfs(nf90_get_att(ncid1,nf90_global,'product_name',product_name))
	call nfs(nf90_get_att(ncid1,nf90_global,'xref_orbit_data',xref_orbit_data))
	orbit_type = which_orbit_type (xref_orbit_data)

! Read the time dimension

	call nfs(nf90_inquire_dimension(ncid1,1,dimnm,nrec))
	if (dimnm /= 'time_01') stop 'Error reading time dimension'
	if (nrec > in_max) then
		call rads_message ('Too many measurements in input file, skipped: '//filenm)
		call nfs(nf90_close(ncid1))
		cycle
	endif
	if (allocated(time)) deallocate (time, lat, lon)
	allocate (time(nrec),lat(nrec),lon(nrec))
	call nfs(nf90_inq_varid(ncid1,'time_01',varid))
	call nfs(nf90_get_var(ncid1,varid,time))

! Check processing time

	if (product_name(49:63) < newer_than) then
		call nfs(nf90_close(ncid1))
		call write_line ('.... old', 1, nrec, nrec, time(1), time(nrec), '<', filenm)
		cycle
	endif

! Read latitude and longitude

	call nfs(nf90_inq_varid(ncid1,'lat_01',varid))
	call nfs(nf90_get_var(ncid1,varid,lat))
	lat = lat * 1d-6
	call nfs(nf90_inq_varid(ncid1,'lon_01',varid))
	call nfs(nf90_get_var(ncid1,varid,lon))
	lon = lon * 1d-6

! First advance to beyond the last time tag

	do i0 = 1, nrec
		if (time(i0) > last_time) exit
	enddo
	if (i0 > 1) call write_line ('... skip', 1, i0-1, nrec, time(1), time(i0-1), '<', filenm)

! If none of the records are after last_time, skip the whole input file

	if (i0 > nrec) then
		call nfs(nf90_close(ncid1))
		cycle
	endif

! Check if the first input point is in a new pass

	call which_pass (time(i0))
	if (ipass /= ipass0) call write_output
	last_time = time(nrec)

! Also duplicated measurements within a single file
! Split pass when entering into a new pass

	do i = max(2,i0),nrec
		call which_pass (time(i))
		if (time(i) == time(i-1)) then
			call fill_fin (i0, i-2)
			call write_line ('... skip', i-1, i-1, nrec, time(i-1), time(i-1), '<', filenm)
			i0 = i
		else if (ipass /= ipass0) then
			call fill_fin (i0, i-1)
			call write_output
			i0 = i
		endif
	enddo
	call fill_fin (i0, nrec)

! Deallocate time and location arrays

	deallocate (time, lat, lon)

enddo

! Dump the remainder of the input files to output

call write_output

contains

!***********************************************************************
! Read Orbital Revolution File (ORF)

subroutine read_orf (sat, orf)
character(len=3), intent(in) :: sat
type(orfinfo), intent(inout) :: orf(:)
character(len=320) :: line
integer :: hash, mjd, yy, mm, dd, hh, mn, ios, npass
real(eightbytereal) :: ss, lon, lat

! Open the equator crossing table

select case (sat)
case ('S3A')
	call parseenv ('${ALTIM}/data/ODR.SNTNL-3A/orf.txt', line)
case ('S3B')
	call parseenv ('${ALTIM}/data/ODR.SNTNL-3B/orf.txt', line)
case default
	stop 'Wrong satellite code'
end select
open (10, file=line, status='old')

! Initialise with dummy values

orf = orfinfo (-1, -1, nan, nan, nan)

! Skip until after the lines starting with #

hash = 0
do while (hash < 5)
	read (10,550,iostat=ios) line
	if (ios /= 0) stop 'Premature end of file'
	if (line(:1) == '#') hash = hash + 1
enddo

! Read the equator crossing table

npass = 1
do
	read (10,550,iostat=ios) line
	if (ios /= 0) exit
	read (line,600) yy,mm,dd,hh,mn,ss,orf(npass)%cycle,orf(npass)%pass,lon,lat
	call ymd2mjd(yy,mm,dd,mjd)
	ss = (mjd-51544)*86400d0 + hh*3600d0 + mn*60d0 + ss
	if (abs(lat) > 1) then
		orf(npass)%starttime = ss
	else
		orf(npass)%eqtime = ss
		orf(npass)%eqlon = lon
		npass=npass + 1
	endif
enddo
close (10)
550 format (a)
600 format (i4,4(1x,i2),1x,f6.3,1x,i3,1x,i5,6x,2(1x,f6.2))

end subroutine read_orf

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
fin(nfile)%rec0 = i0
fin(nfile)%rec1 = i1
fin(nfile)%time0 = time(i0)
fin(nfile)%time1 = time(i1)
fin(nfile)%lat0 = lat(i0)
fin(nfile)%lat1 = lat(i1)
fin(nfile)%lon0 = lon(i0)
fin(nfile)%lon1 = lon(i1)
call write_line ('.. input', i0, i1, nrec, time(i0), time(i1), '<', filenm)
end subroutine fill_fin

!***********************************************************************
! Determine the type of orbit used

function which_orbit_type (xref_orbit_data)
integer :: which_orbit_type
character(len=*), intent(in) :: xref_orbit_data
select case (xref_orbit_data(10:12))
case ('POE')
	which_orbit_type = 8
case ('MDO')
	which_orbit_type = 6
case ('ROE')
	which_orbit_type = 4
case ('NAV')
	which_orbit_type = 3
case ('NAT')
	which_orbit_type = 2
case ('FPO')
	which_orbit_type = 1
case ('OSF')
	which_orbit_type = 0
case default
	which_orbit_type = 127_onebyteint
end select
end function which_orbit_type

!***********************************************************************
! Determine the corresponding record from ORF

subroutine which_pass (time)
real(eightbytereal), intent(in) :: time
do while (time < orf(ipass)%starttime)
	ipass = ipass - 1
	if (ipass < 1) stop 'Times are beyond limits of ORF file'
enddo
do while (time > orf(ipass+1)%starttime)
	ipass = ipass + 1
	if (orf(ipass)%cycle < 0) stop 'Times are beyond limits of ORF file'
enddo
end subroutine which_pass

!***********************************************************************
! Write whatever has been buffered so far to an output file

subroutine write_output
integer(fourbyteint) :: nrec, ncid1, ncid2, varid1, varid2, varid3, xtype, ndims, &
	dimids(2), dimid2, natts, i, nvars, idxin(2)=1, idxut(2)=1, nout, &
	cycle_number, pass_number, absolute_pass_number, absolute_rev_number
real(eightbytereal) :: equator_time, equator_longitude
character(len=rads_naml) :: dirnm, prdnm, outnm, attnm, varnm
real(eightbytereal), allocatable :: darr1(:)
integer(fourbyteint), allocatable :: iarr1(:)
logical :: exist
integer(onebyteint), parameter :: flag_values(0:8) = int((/0,1,2,3,4,5,6,7,8/), onebyteint)

! Retrieve the pass variables

cycle_number = orf(ipass0)%cycle
pass_number = orf(ipass0)%pass
equator_time = orf(ipass0)%eqtime + sec2000
equator_longitude = orf(ipass0)%eqlon
ipass0 = ipass

! How many records are buffered for output?
! Skip if there is nothing left

if (nfile == 0 .or. ipass0 == 0) return
nrec = sum(fin(1:nfile)%rec1 - fin(1:nfile)%rec0 + 1)

! Open the output file. Make directory if needed.

605 format (a,'/c',i3.3)
610 format (a,i3.3,'_',i3.3,a,'.nc')

write (dirnm,605) trim(destdir),cycle_number
inquire (file=dirnm,exist=exist)
if (.not.exist) call system('mkdir -p '//dirnm)
write (prdnm,610) product_name(:15),cycle_number,pass_number,product_name(77:94)
outnm = trim(dirnm) // '/' // trim(prdnm)
inquire (file=outnm,exist=exist)

! If exist, then keep the file if the buffer is smaller or equal in size
! If it is larger, delete the existing file

if (exist) then
	call nfs(nf90_open(outnm,nf90_nowrite,ncid2))
	call nfs(nf90_inquire_dimension(ncid2,1,len=nout))
	call nfs(nf90_get_att(ncid2,nf90_global,'first_meas_time',date(1)))
	call nfs(nf90_get_att(ncid2,nf90_global,'last_meas_time',date(2)))
	call nfs(nf90_close(ncid2))
	if (nrec <= nout) then
		do i = 1, nfile
			! Release netCDF file when reaching end
			if (fin(i)%rec1 == fin(i)%nrec) call nfs(nf90_close(fin(i)%ncid))
		enddo
		write (*,620) 'Keeping ', 1, nout, nout, date(1:2), '>', trim(outnm)
		nfile = 0
		return
	endif
	call system('rm -f '//outnm)
endif
620 format (a,' : ',3i6,' : ',a,' - ',a,1x,a1,1x,a)

! Create a new file

call nfs(nf90_create(outnm,nf90_write+nf90_nofill,ncid2))
call nfs(nf90_set_fill(ncid2,nf90_nofill,i))

! Create the time dimension

call nfs(nf90_def_dim(ncid2,'time_01',nrec,dimid2))

! Copy all the variable definitions and attributes

ncid1 = fin(1)%ncid
call nfs(nf90_inquire(ncid1,nvariables=nvars,nattributes=natts))
varid2 = 0
do varid1 = 0, nvars
	if (varid1 > 0) then
		call nfs(nf90_inquire_variable(ncid1,varid1,varnm,xtype,ndims,dimids,natts))
		! Skip all listed excluded variables, all 20-Hz variables
		! Skip also all uint variables because the netCDF library doesn't work with them
		! (they occur only in enhanced_measurements.nc)
		if (excluded(varnm) .or. ndims > 1 .or. dimids(1) > 1 .or. xtype == nf90_uint) cycle
		call nfs(nf90_def_var(ncid2,varnm,xtype,dimids(1:ndims),varid2))
	endif
	do i = 1,natts
		call nfs(nf90_inq_attname(ncid1,varid1,i,attnm))
		call nfs(nf90_copy_att(ncid1,varid1,attnm,ncid2,varid2))
	enddo
enddo

! Add the orbit type (only for versions that do not have it already)

if (nft(nf90_inq_varid(ncid1,'orbit_type',varid1))) then
	call nfs(nf90_def_var(ncid2,'orbit_type',nf90_byte,dimid2,varid3))
	call nfs(nf90_put_att(ncid2,varid3,'long_name','Orbit type flag : 1 Hz Ku band'))
	call nfs(nf90_put_att(ncid2,varid3,'_FillValue',127_onebyteint))
	call nfs(nf90_put_att(ncid2,varid3,'flag_values',flag_values))
	call nfs(nf90_put_att(ncid2,varid3,'flag_meanings','osf fos navatt doris_nav gnss_roe pod_moe salp_moe pod_poe salp_poe'))
	call nfs(nf90_put_att(ncid2,varid3,'coordinates','lon_01 lat_01'))
else
	varid3 = 0	! This signals that there is no need to write a new orbit_type variable
endif

! Determine absolute pass, rev, and equator crossing information

absolute_pass_number = (cycle_number-1)*nr_passes+pass_number-54
absolute_rev_number = absolute_pass_number / 2

! Overwrite some attributes and product name

call write_line ('Creating', 1, nrec, nrec, fin(1)%time0, fin(nfile)%time1,'>', outnm)
call strf1985f(date(3),equator_time)

call nfs(nf90_put_att(ncid2,nf90_global,'product_name',prdnm))
call nfs(nf90_put_att(ncid2,nf90_global,'cycle_number',cycle_number))
call nfs(nf90_put_att(ncid2,nf90_global,'pass_number',pass_number))
call nfs(nf90_put_att(ncid2,nf90_global,'absolute_pass_number',absolute_pass_number))
call nfs(nf90_put_att(ncid2,nf90_global,'absolute_rev_number',absolute_rev_number))
call nfs(nf90_put_att(ncid2,nf90_global,'equator_time',date(3)))
call nfs(nf90_put_att(ncid2,nf90_global,'equator_longitude',equator_longitude))
call nfs(nf90_put_att(ncid2,nf90_global,'first_meas_time',date(1)))
call nfs(nf90_put_att(ncid2,nf90_global,'last_meas_time',date(2)))
call nfs(nf90_put_att(ncid2,nf90_global,'first_meas_lat',fin(1)%lat0))
call nfs(nf90_put_att(ncid2,nf90_global,'last_meas_lat',fin(nfile)%lat1))
call nfs(nf90_put_att(ncid2,nf90_global,'first_meas_lon',fin(1)%lon0))
call nfs(nf90_put_att(ncid2,nf90_global,'last_meas_lon',fin(nfile)%lon1))
call nfs(nf90_enddef(ncid2))

! Copy all data elements

nout = 0
do i = 1,nfile
	nrec = fin(i)%rec1 - fin(i)%rec0 + 1
	if (nrec == 0) stop "nrec == 0"
	allocate (darr1(nrec),iarr1(nrec))
	ncid1 = fin(i)%ncid
	idxin(2) = fin(i)%rec0
	idxut(2) = nout + 1
	varid2 = 0
	do varid1 = 1,nvars
		call nfs(nf90_inquire_variable(ncid1,varid1,varnm,xtype,ndims,dimids,natts))
		if (excluded(varnm) .or. ndims > 1 .or. dimids(1) > 1 .or. xtype == nf90_uint) cycle
		call nfs(nf90_inq_varid(ncid2,varnm,varid2))
		if (xtype == nf90_double) then
			call nfs(nf90_get_var(ncid1,varid1,darr1,idxin(2:2)))
			call nfs(nf90_put_var(ncid2,varid2,darr1,idxut(2:2)))
		else
			call nfs(nf90_get_var(ncid1,varid1,iarr1,idxin(2:2)))
			call nfs(nf90_put_var(ncid2,varid2,iarr1,idxut(2:2)))
		endif
	enddo
	iarr1 = fin(i)%type
	if (varid3 > 0) call nfs(nf90_put_var(ncid2,varid3,iarr1,idxut(2:2)))
	nout = nout + nrec
	deallocate (darr1,iarr1)
	! Release netCDF file when reaching end
	if (fin(i)%rec1 == fin(i)%nrec) call nfs(nf90_close(fin(i)%ncid))
enddo

! Close output

call nfs(nf90_close(ncid2))
nfile = 0

end subroutine write_output

subroutine write_line (word, rec0, rec1, nrec, time0, time1, dir, filenm)
character(len=*), intent(in) :: word, dir, filenm
integer, intent(in) :: rec0, rec1, nrec
real(eightbytereal), intent(in) :: time0, time1
call strf1985f(date(1),time0+sec2000)
call strf1985f(date(2),time1+sec2000)
write (*,620) word, rec0, rec1, nrec, date(1:2), dir, trim(filenm)
620 format (a,' : ',3i6,' : ',a,' - ',a,1x,a1,1x,a)
end subroutine write_line

function excluded (varnm)
character(len=*), intent(in) :: varnm
logical :: excluded
excluded = (index(exclude_list,','//trim(varnm)//',') > 0)
end function excluded

end program s3combine
