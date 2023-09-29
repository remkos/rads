!-----------------------------------------------------------------------
! Copyright (c) 2011-2023  Remko Scharroo
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

!*ogdrsplit -- Split Jason or SARAL OGDR files into pass files
!
! Read OSDR_SSHA or OGDR or OGDR_SSHA files from Jason-1/2/3 or SARAL and
! split them into separate pass files. The input file names are read from
! standard input. The individual pass files will be named
! <destdir>/cCCC/???_*_2P?PCCC_[P]PPP.nc, where CCC is the cycle number and
! [P]PPP the pass number. The directory <destdir>/cCCC will be created if needed.
!
! This program needs an up-to-date files $RADSROOT/ext/??/???_ORF.txt with
! equator crossing information.
!-----------------------------------------------------------------------
program ogdrsplit

use rads
use rads_misc
use rads_time
use rads_netcdf
use rads_devel_misc
use netcdf

! Struct for orbit info

type(orfinfo) :: orf(200000)

! General variables

character(len=rads_cmdl) :: arg, filenm, dimnm, filetype, destdir
character(len=rads_strl) :: exclude_list = ','
real(eightbytereal), parameter :: sec2000 = 473299200d0
integer(fourbyteint) :: l, ipass, i0, i, ncid1, &
	nrec, nhz, ios, varid, n_ignore = 0, max_dim = 2
real(eightbytereal), allocatable :: time(:)

! Print description, if requested

if (iargc() < 1) then
	write (*,1300)
	stop
endif
1300 format ('ogdrsplit -- Split OGDR or OGDR SSHA files into pass files'// &
'syntax: ogdrsplit [options] destdir < list'//'where'/ &
'  destdir           : Destination directory (appends c???/*.nc)'/ &
'  list              : List of input file names'// &
'where [options] are:' / &
'  -iNRECS           : Ignore up to NRECS record chunks to move to existing files (def: 0)' / &
'  -x2               : Exclude any multi-dimensional variables' / &
'  -xVAR1[,VAR2,...] : Exclude variable(s) from copying')

! First determine filetype

read (*,550,iostat=ios) filenm
if (ios /= 0) stop
l = index(filenm,'_2P')
if (l == 0) call rads_exit ('Wrong filetype')
do i = l-1,5,-1
	if (filenm(i:i) == '_') exit
enddo
filetype = filenm(i-3:l+3)

! Get ORF file

call read_orf (filetype, orf)

! Read options and destination directory

do i = 1,iargc()
	call getarg (i,arg)
	if (arg(:2) == '-i') then
		read (arg(3:),*) n_ignore
	else if (arg(:3) == '-x2') then
		max_dim = 1
	else if (arg(:2) == '-x') then
		exclude_list = trim(exclude_list) // arg(3:len_trim(arg)) // ','
	else
		destdir = arg
	endif
enddo

! Open the input OGDR file(s)

do
	call nfs(nf90_open(filenm,nf90_nowrite,ncid1))

! Read the time dimension

	call nfs(nf90_inquire_dimension(ncid1,1,dimnm,nrec))
	if (dimnm /= 'time') call rads_exit ('Error reading time dimension')
	if (nft(nf90_inquire_dimension(ncid1,2,dimnm,nhz))) then
		! No second dimension in OGDR-GPS or OGDR-SSH files
		nhz = 1
	else
		if (dimnm /= 'meas_ind') call rads_exit ('Error reading meas_ind dimension')
	endif
	allocate (time(nrec))
	call nfs(nf90_inq_varid(ncid1,'time',varid))
	call nfs(nf90_get_var(ncid1,varid,time))
	write (*,610) trim(filenm)

! Check time with equator crossing table

	ipass = 1
	i0 = 1
	do i = 1,nrec
		do while (time(i) > orf(ipass+1)%starttime)
			call copyfile(i0,i-1)
			ipass = ipass + 1
			if (orf(ipass)%cycle < 0) call rads_exit ('Times are beyond limits of ORF file')
			i0 = i
		enddo
	enddo
	call copyfile(i0,nrec)
	call nfs(nf90_close(ncid1))
	deallocate (time)

! Read next filename

	read (*,550,iostat=ios) filenm
	if (ios /= 0) exit
enddo

! Formats

550 format (a)
610 format ('Splitting : ',a)

contains

!***********************************************************************
! Copy the contents of the input file to a new output file

subroutine copyfile (rec0, rec1)
integer(fourbyteint), intent(inout) :: rec0
integer(fourbyteint), intent(in) :: rec1
integer(fourbyteint) :: nrec,ncid2,varid,varid2,xtype,ndims,dimids(2),natts,i,nvars,idxin(2)=1,idxut(2)=1,dimlen
character(len=rads_naml) :: outnm,attnm,varnm
real(eightbytereal) :: time2(2)
real(eightbytereal), allocatable :: darr1(:),darr2(:,:)
integer(fourbyteint), allocatable :: iarr1(:),iarr2(:,:),iarr3(:)
logical :: exist
character(len=26) :: date(3)

! Skip empty data chunks

if (rec1 < rec0) return
call nfs(nf90_inquire(ncid1,nvariables=nvars,nattributes=natts))

! Open the output file. Make directory if needed.

605 format (a,'/c',i3.3)
610 format (a,'/c',i3.3,'/',a,'P',i3.3,'_',i3.3,'.nc') ! JA? format
611 format (a,'/c',i3.3,'/',a,'P',i3.3,'_',i4.4,'.nc') ! SRL format
620 format ('... Records : ',3i6,' : ',a,1x,a,' - ',a,a)

write (outnm,605) trim(destdir),orf(ipass)%cycle
inquire (file=outnm,exist=exist)
if (.not.exist) call system('mkdir -p '//outnm)
if (filetype(:3) == 'SRL') then
	write (outnm,611) trim(destdir),orf(ipass)%cycle,trim(filetype),orf(ipass)%cycle,orf(ipass)%pass
else
	write (outnm,610) trim(destdir),orf(ipass)%cycle,trim(filetype),orf(ipass)%cycle,orf(ipass)%pass
endif
inquire (file=outnm,exist=exist)

if (exist) then

! Match up the last time in the existing file with times in the input file

	call nfs(nf90_open(outnm,nf90_write,ncid2))
	call nfs(nf90_set_fill(ncid2,nf90_nofill,i))
	call nfs(nf90_inquire_dimension(ncid2,1,len=dimlen))
	idxin(2)=dimlen
	call nfs(nf90_get_var(ncid2,1,time2(1),idxin(1:1)))
	call nfs(nf90_get_var(ncid2,1,time2(2),idxin(2:2)))
	do while (time(rec0) < time2(2)+0.5d0 .and. rec0 <= rec1)
		rec0 = rec0 + 1
	enddo
	nrec = rec1 - rec0 + 1
	if (nrec <= n_ignore) then
		if (nrec > 0) then
			date(1) = strf1985f(time(rec0)+sec2000)
			date(2) = strf1985f(time(rec1)+sec2000)
			write (*,620) rec0,rec1,nrec,trim(outnm),date(1:2),' (skipped)'
		endif
		call nfs(nf90_close(ncid2))
		return
	endif
	call nfs(nf90_redef(ncid2))

else

	call nfs(nf90_create(outnm,nf90_write+nf90_nofill,ncid2))
	call nfs(nf90_set_fill(ncid2,nf90_nofill,i))
	dimlen = 0

! Create the dimensions

	call nfs(nf90_def_dim(ncid2,'time',nf90_unlimited,i))
	if (nhz > 1 .and. max_dim > 1) &
		call nfs(nf90_def_dim(ncid2,'meas_ind',nhz,i)) ! Create only if there is a 2nd dimension
	time2(1)=time(rec0)

! Copy all the variable definitions and attributes

	varid2 = 0
	do varid = 0,nvars
		if (varid > 0) then
			call nfs(nf90_inquire_variable(ncid1,varid,varnm,xtype,ndims,dimids,natts))
			if (excluded(varnm) .or. ndims > max_dim .or. dimids(1) > max_dim) cycle
			call nfs(nf90_def_var(ncid2,varnm,xtype,dimids(1:ndims),varid2))
		endif
		do i = 1,natts
			call nfs(nf90_inq_attname(ncid1,varid,i,attnm))
			call nfs(nf90_copy_att(ncid1,varid,attnm,ncid2,varid2))
		enddo
	enddo
endif

! Initialize

nrec = rec1 - rec0 + 1
idxin(2) = rec0
idxut(2) = dimlen + 1
time2(2) = time(rec1)
date(1) = strf1985f(time(rec0)+sec2000)
date(2) = strf1985f(time(rec1)+sec2000)
write (*,620) rec0,rec1,nrec,trim(outnm),date(1:2)
date(1) = strf1985f(time2(1)+sec2000)
date(2) = strf1985f(time2(2)+sec2000)
date(3) = strf1985f(orf(ipass)%eqtime+sec2000)

allocate (darr1(nrec),darr2(nhz,nrec),iarr1(nrec),iarr2(nhz,nrec),iarr3(nhz))

! Overwrite cycle/pass attributes

call nfs(nf90_put_att(ncid2,nf90_global,'cycle_number',orf(ipass)%cycle))
call nfs(nf90_put_att(ncid2,nf90_global,'pass_number',orf(ipass)%pass))
call nfs(nf90_put_att(ncid2,nf90_global,'absolute_pass_number',orf(ipass)%abs_pass))
call nfs(nf90_put_att(ncid2,nf90_global,'absolute_rev_number',orf(ipass)%abs_rev))
call nfs(nf90_put_att(ncid2,nf90_global,'equator_time',date(3)))
call nfs(nf90_put_att(ncid2,nf90_global,'equator_longitude',orf(ipass)%eqlon))
call nfs(nf90_put_att(ncid2,nf90_global,'first_meas_time',date(1)))
call nfs(nf90_put_att(ncid2,nf90_global,'last_meas_time',date(2)))
call nfs(nf90_put_att(ncid2,nf90_global,'original',trim(filenm)))
call nfs(nf90_enddef(ncid2))

! Copy all data elements

varid2 = 0
do varid = 1,nvars
	call nfs(nf90_inquire_variable(ncid1,varid,varnm,xtype,ndims,dimids,natts))
	if (excluded(varnm) .or. ndims > max_dim .or. dimids(1) > max_dim) cycle
	varid2 = varid2 + 1
	if (xtype == nf90_double) then
		if (ndims == 2) then
			call nfs(nf90_get_var(ncid1,varid ,darr2,idxin))
			call nfs(nf90_put_var(ncid2,varid2,darr2,idxut))
		else
			call nfs(nf90_get_var(ncid1,varid ,darr1,idxin(2:2)))
			call nfs(nf90_put_var(ncid2,varid2,darr1,idxut(2:2)))
		endif
	else
		if (ndims == 2) then
			call nfs(nf90_get_var(ncid1,varid ,iarr2,idxin))
			call nfs(nf90_put_var(ncid2,varid2,iarr2,idxut))
		else if (dimids(1) == 2) then
			call nfs(nf90_get_var(ncid1,varid ,iarr3))
			call nfs(nf90_put_var(ncid2,varid2,iarr3))
		else
			call nfs(nf90_get_var(ncid1,varid ,iarr1,idxin(2:2)))
			call nfs(nf90_put_var(ncid2,varid2,iarr1,idxut(2:2)))
		endif
	endif
enddo

deallocate (darr1,iarr1,iarr2,iarr3)

call nfs(nf90_close(ncid2))

end subroutine copyfile

function excluded (varnm)
character(len=*), intent(in) :: varnm
logical :: excluded
excluded = (index(exclude_list,','//trim(varnm)//',') > 0)
end function excluded

end program ogdrsplit
