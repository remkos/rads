!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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
! Read OSDR_SSHA or OGDR or OGDR_SSHA files from Jason-1, Jason-2 or SARAL and
! split them into separate pass files. The OGDR file names are read from
! standard input. The individual pass files will be named
! <destdir>/cCCC/???_*_2P?PCCC_[P]PPP.nc, where CCC is the cycle number and
! [P]PPP the pass number. The directory <destdir>/cCCC will be created if needed.
!
! This program needs an up-to-date files $RADSROOT/ext/??/???_ORF.txt with
! equator crossing information.
!-----------------------------------------------------------------------
program ogdrsplit

use rads
use rads_netcdf
use typesizes
use netcdf

character(80) :: radsroot,arg,orf,dimnm,filetype,destdir
character(26) :: date(2)
integer(fourbyteint), parameter :: mpass = 254 * 500
real(eightbytereal), parameter :: sec2000 = 473299200d0
integer(fourbyteint) :: l,yy,mm,dd,hh,mn,cycle(mpass),pass(mpass),npass,ipass,i0,i,ncid1, &
	nrec,nhz,hash,mjd,ios,varid
real(eightbytereal) :: ss,lon,lat,eqtime(mpass),eqlon(mpass),starttime(mpass)
real(eightbytereal), allocatable :: time(:)

! Print description, if requested

if (iargc() /= 1) then
	write (*,1300)
	stop
endif
1300 format ('ogdrsplit -- Split OGDR or OGDR SSHA file into pass files'// &
'syntax: ogdrsplit destdir < list'//'where'/ &
'  destdir : Destination directory (appends c???/*.nc)'/ &
'  list    : Input OGDR or OGDR SSHA files')

! First determine filetype

read (*,550,iostat=ios) arg
if (ios /= 0) stop
l = index(arg,'_2P')
if (l == 0) stop 'Wrong filetype'
do i = l-1,5,-1
	if (arg(i:i) == '_') exit
enddo
filetype = arg(i-3:l+3)

! Open the equator crossing table

call checkenv('RADSROOT',radsroot,l)
select case (filetype(:3))
case ('JA1')
	orf = radsroot(:l)//'/ext/j1/JA1_ORF.txt'
case ('JA2')
	orf = radsroot(:l)//'/ext/j2/JA2_ORF.txt'
case ('SRL')
	orf = radsroot(:l)//'/ext/sa/SRL_ORF.txt'
case default
	stop 'Wrong filetype'
end select
open (10, file=orf, status='old')

! Skip until after the lines starting with #

hash = 0
do while (hash < 5)
	read (10,550,iostat=ios) orf
	if (ios /= 0) stop 'Premature end of file'
	if (orf(:1) == '#') hash = hash + 1
enddo

! Read the equator crossing table

npass = 1
do
	read (10,550,iostat=ios) orf
	if (ios /= 0) exit
	read (orf,600) yy,mm,dd,hh,mn,ss,cycle(npass),pass(npass),lon,lat
	if (abs(lat) > 1) then
		call ymd2mjd(yy,mm,dd,mjd)
		starttime(npass) = (mjd-51544)*86400d0 + hh*3600d0 + mn*60d0 + ss
	else
		eqtime(npass) = (mjd-51544)*86400d0 + hh*3600d0 + mn*60d0 + ss
		eqlon(npass) = lon
		npass=npass + 1
	endif
enddo
close (10)
npass = npass - 1
550 format (a)
600 format (i4,4(1x,i2),1x,f6.3,1x,i3,1x,i5,6x,2(1x,f6.2))
610 format ('Splitting : ',a)

! Read destination directory

call getarg (1,destdir)

! Open the input OGDR file(s)

do
	call nfs(nf90_open(arg,nf90_nowrite,ncid1))

! Read the time dimension

	call nfs(nf90_inquire_dimension(ncid1,1,dimnm,nrec))
	if (dimnm /= 'time') stop 'Error reading time dimension'
	if (nft(nf90_inquire_dimension(ncid1,2,dimnm,nhz))) then
		! No second dimension in OGDR-GPS or OGDR-SSH files
		nhz = 1
	else
		if (dimnm /= 'meas_ind') stop 'Error reading meas_ind dimension'
	endif
	allocate (time(nrec))
	call nfs(nf90_inq_varid(ncid1,'time',varid))
	call nfs(nf90_get_var(ncid1,varid,time))
	call nfs(nf90_get_att(ncid1,nf90_global,'first_meas_time',date(1)))
	call nfs(nf90_get_att(ncid1,nf90_global,'last_meas_time',date(2)))
	write (*,610) trim(arg)

! Check time with equator crossing table

	if (time(nrec) > starttime(npass)) stop 'Times are beyond equator crossing table'

	ipass = 1
	i0 = 1
	do i = 1,nrec
		do while (time(i) > starttime(ipass+1))
			call copyfile(i0,i-1)
			ipass = ipass + 1
			i0 = i
		enddo
	enddo
	call copyfile(i0,nrec)
	call nfs(nf90_close(ncid1))
	deallocate (time)

! Read next filename

	read (*,550,iostat=ios) arg
	if (ios /= 0) exit
enddo

contains

!***********************************************************************
! Copy the contents of the input file to a new output file

subroutine copyfile (rec0, rec1)
integer(fourbyteint), intent(inout) :: rec0
integer(fourbyteint), intent(in) :: rec1
integer(fourbyteint) :: nrec,ncid2,varid,xtype,ndims,dimids(2),natts,i,nvars,idxin(2)=1,idxut(2)=1,dimlen
character(80) :: outnm,attnm,varnm
real(eightbytereal) :: time2(2)
real(eightbytereal), allocatable :: darr1(:),darr2(:,:)
integer(fourbyteint), allocatable :: iarr1(:),iarr2(:,:),iarr3(:)
logical :: exist

! Skip empty data chunks

if (rec1 < rec0) return
call nfs(nf90_inquire(ncid1,nvariables=nvars,nattributes=natts))

! Open the output file. Make directory if needed.

605 format (a,'/c',i3.3)
610 format (a,'/c',i3.3,'/',a,'P',i3.3,'_',i3.3,'.nc') ! JA1/JA2 format
611 format (a,'/c',i3.3,'/',a,'P',i3.3,'_',i4.4,'.nc') ! SRL format
620 format ('... Records : ',3i6,' : ',a,1x,a,' - ',a)

write (outnm,605) trim(destdir),cycle(ipass)
inquire (file=outnm,exist=exist)
if (.not.exist) call system('mkdir -p '//outnm)
if (filetype(:3) == 'SRL') then
	write (outnm,611) trim(destdir),cycle(ipass),trim(filetype),cycle(ipass),pass(ipass)
else
	write (outnm,610) trim(destdir),cycle(ipass),trim(filetype),cycle(ipass),pass(ipass)
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
	if (rec1 < rec0) then
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
	if (nhz > 1) call nfs(nf90_def_dim(ncid2,'meas_ind',nhz,i)) ! Create only if there is a 2nd dimension
	time2(1)=time(rec0)

! Copy all the variable definitions and attributes

	do varid = 0,nvars
		if (varid > 0) then
			call nfs(nf90_inquire_variable(ncid1,varid,varnm,xtype,ndims,dimids,natts))
			call nfs(nf90_def_var(ncid2,varnm,xtype,dimids(1:ndims),varid))
		endif
		do i = 1,natts
			call nfs(nf90_inq_attname(ncid1,varid,i,attnm))
			call nfs(nf90_copy_att(ncid1,varid,attnm,ncid2,varid))
		enddo
	enddo
endif

! Overwrite cycle/pass attributes

call nfs(nf90_put_att(ncid2,nf90_global,'cycle_number',cycle(ipass)))
call nfs(nf90_put_att(ncid2,nf90_global,'pass_number',pass(ipass)))
call nfs(nf90_put_att(ncid2,nf90_global,'absolute_pass_number',(cycle(ipass)-1)*254+pass(ipass)))
call strf1985f(date(1),eqtime(ipass)+sec2000)
call nfs(nf90_put_att(ncid2,nf90_global,'equator_time',date(1)))
call nfs(nf90_put_att(ncid2,nf90_global,'equator_longitude',eqlon(ipass)))
call strf1985f(date(1),time2(1)+sec2000)
call nfs(nf90_put_att(ncid2,nf90_global,'first_meas_time',date(1)))
call strf1985f(date(2),time(rec1)+sec2000)
call nfs(nf90_put_att(ncid2,nf90_global,'last_meas_time',date(2)))
call strf1985f(date(1),time(rec0)+sec2000)
call nfs(nf90_enddef(ncid2))

! Initialize

nrec = rec1 - rec0 + 1
idxin(2)=rec0
idxut(2)=dimlen + 1
write (*,620) rec0,rec1,rec1-rec0+1,trim(outnm),date

allocate (darr1(nrec),darr2(nhz,nrec),iarr1(nrec),iarr2(nhz,nrec),iarr3(nhz))

! Copy all data elements

do varid = 1,nvars
	call nfs(nf90_inquire_variable(ncid1,varid,varnm,xtype,ndims,dimids,natts))
	if (xtype == nf90_double) then
		if (ndims == 2) then
			call nfs(nf90_get_var(ncid1,varid,darr2,idxin))
			call nfs(nf90_put_var(ncid2,varid,darr2,idxut))
		else
			call nfs(nf90_get_var(ncid1,varid,darr1,idxin(2:2)))
			call nfs(nf90_put_var(ncid2,varid,darr1,idxut(2:2)))
		endif
	else
		if (ndims == 2) then
			call nfs(nf90_get_var(ncid1,varid,iarr2,idxin))
			call nfs(nf90_put_var(ncid2,varid,iarr2,idxut))
		else if (dimids(1) == 2) then
			call nfs(nf90_get_var(ncid1,varid,iarr3))
			call nfs(nf90_put_var(ncid2,varid,iarr3))
		else
			call nfs(nf90_get_var(ncid1,varid,iarr1,idxin(2:2)))
			call nfs(nf90_put_var(ncid2,varid,iarr1,idxut(2:2)))
		endif
	endif
enddo

deallocate (darr1,iarr1,iarr2)

call nfs(nf90_close(ncid2))

end subroutine copyfile

end program ogdrsplit
