!-----------------------------------------------------------------------
! $Id$
!
! Copyright (C) 2012  Remko Scharroo (Altimetrics LLC)
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

!*radsxolist -- RADS crossover lister
!+
program radsxolist

! This program lists the contents of RADS netCDF crossover files and
! computes statistics.
!-----------------------------------------------------------------------
use netcdf
use rads
use rads_misc
use rads_netcdf
use rads_time
integer(fourbyteint) :: ncid, i, icol, ios, nvar, nxo_in, nxo_out, ntrk, id_satid, id_track, id_sla, id_new, id_old
real(eightbytereal), allocatable :: var(:,:,:), stat(:,:,:)
integer(fourbyteint), allocatable :: track(:,:)
character(len=rads_naml), allocatable :: long_name(:)
logical, allocatable :: mask(:)
type :: trk_
	real(eightbytereal) :: equator_lon, equator_time, start_time, end_time
	integer(twobyteint) :: nr_alt, nr_xover, satid, cycle, pass
end type
type(trk_), allocatable :: trk(:)
character(len=640) :: fmt_string
character(len=rads_naml) :: optopt, optarg
character(len=rads_cmdl) :: command
character(len=1) :: order = ''
character(len=4) :: statname(5) = (/ 'min ', 'max ', 'mean', 'rms ', 'std ' /)
integer(fourbyteint), parameter :: msat = 20, maxtrk = 32768
type :: sat_
	character(len=4) :: name
	real(eightbytereal) :: period, altsig, orberr, inclination
end type sat_
type(sat_) :: sat(msat)
character(len=3*msat) :: satlist = ''
type(rads_sat) :: S
character(len=*), parameter :: optlist = &
	'e::dslnotr: lon: lat: dt: edit:: dual single nolist both-legs order both-times time: ymd: doy: sec:'

integer(fourbyteint) :: var0 = 0
logical :: diff = .true., stat_only = .false., singles = .true., duals = .true., xostat = .false.
real(eightbytereal) :: t0 = 0d0, t1 = 0d0, lon0 = 0d0, lon1 = 0d0, lat0 = 0d0, lat1 = 0d0, dt0 = 0d0, dt1 = 0d0, edit = -1d0

! Check operation mode
call getarg (0, command)
xostat = (index(command, 'radsxostat') > 0)
stat_only = xostat

! Initialize RADS or issue help
call synopsis

! Write header
call get_command (command, status = i)
if (xostat) then
	write (*, 600) 'Statistics', timestamp(), trim(command)
else
	write (*, 600) 'List', timestamp(), trim(command)
endif
600 format ('# ',a,' of RADS crossovers'/'# Created: ',a,' UTC: ',a)

! Scan command line arguments
do
	call getopt (optlist, optopt, optarg)
	select case (optopt)
	case ('!') ! End of arguments
		exit
	case (':', 'r') ! Ignore -r for the time being, ignore unknown option
	case ('lon')
		read (optarg, *, iostat=ios) lon0,lon1
	case ('lat')
		read (optarg, *, iostat=ios) lat0,lat1
	case ('dt')
		read (optarg, *, iostat=ios) dt0,dt1
		if (dt0 > dt1) then
			dt1 = dt0
			dt0 = 0d0
		endif
		dt0 = dt0 * 86400d0
		dt1 = dt1 * 86400d0
	case ('e', 'edit')
		edit = 3.5d0
		read (optarg, *, iostat=ios) edit
	case ('d', 'dual')
		singles = .false.
	case ('s', 'single')
		duals = .false.
	case ('l', 'both-legs')
		diff = .false.
		var0 = 1
	case ('n', 'nolist')
		stat_only = .true.
	case ('o', 'order')
		order = optarg(1:1)
	case ('t', 'both-times')
		var0 = 1
	case (' ')
		call process (optarg) ! Process each file
	case default
		if (dateopt(optopt, optarg, t0, t1)) cycle
	end select
enddo

contains

!***********************************************************************

subroutine process (filename)
character(len=*) :: filename
integer(fourbyteint) :: i, j, k
character(len=rads_naml) :: legs
real(eightbytereal) :: mean, sigma

! Open netCDF file
write (*, 600) trim(filename)
600 format ('#'/'# File name     = ',a)
call nfs (nf90_open (filename, nf90_write, ncid))

! Get the number of xovers and number of tracks
call nfs (nf90_inquire_dimension (ncid, 2, len=nxo_in))
call nfs (nf90_inquire_dimension (ncid, 3, len=ntrk))
if (ntrk >= maxtrk) then
	write (*,'(a)') 'Output format does not allow track numbers exceeding 32767. We will modulo them.'
endif

write (*, 605) nxo_in, ntrk
605 format ('# Xovers,tracks = ',2i9)
if (edit > 0d0) write (*, 606) edit
606 format ('# Edit sigma    = ',f9.3)
fmt_string = '('

! Read all the "base variables" into memory
id_track = get_varid('track')
id_satid = get_varid('satid')
nvar = (id_satid - id_track - 1)
allocate (track(2,nxo_in), trk(ntrk), var(2,nxo_in,-1:nvar), stat(2,-1:nvar,5), mask(nxo_in), long_name(-2:nvar))
mask = .true.
call get_var_1d (get_varid('lat'), var(1,:,-1), long_name(-2))
call get_var_1d (get_varid('lon'), var(2,:,-1), long_name(-1))
call get_var_2d (get_varid('time'), var(:,:,0), long_name(0))
call nfs (nf90_get_var (ncid, id_track, track))
call nfs (nf90_get_var (ncid, id_satid, trk%satid))
call nfs (nf90_get_var (ncid, get_varid('cycle'), trk%cycle))
call nfs (nf90_get_var (ncid, get_varid('pass'), trk%pass))
call nfs (nf90_get_var (ncid, get_varid('equator_lon'), trk%equator_lon))
call nfs (nf90_get_var (ncid, get_varid('equator_time'), trk%equator_time))
call nfs (nf90_get_var (ncid, get_varid('start_time'), trk%start_time))
call nfs (nf90_get_var (ncid, get_varid('end_time'), trk%end_time))
call nfs (nf90_get_var (ncid, get_varid('nr_xover'), trk%nr_xover))
call nfs (nf90_get_var (ncid, get_varid('nr_alt'), trk%nr_alt))
if (nft (nf90_get_att (ncid, nf90_global, 'legs', legs))) legs = 'undetermined undetermined'

! Get essential satellite information
! Older files only have IDs, newer have the satellite abbreviations
if (nf90_get_att (ncid, get_varid('satid'), 'flag_meanings', satlist) == nf90_noerr) then
	do i = 1,len_trim(satlist)/3+1
		call rads_init (S, satlist(i*3-2:i*3-1))
		sat(S%satid) = sat_ (S%sat, 2*S%phase%pass_seconds, &
			S%xover_params(1), S%xover_params(2), S%inclination)
	enddo
else
	satlist = 'g3 ss gs e1 tx pn e2 g1 j1 n1 j2 c2'
	do i = 1,12
		if (.not.any(trk%satid == i)) cycle
		call rads_init (S, satlist(i*3-2:i*3-1))
		sat(S%satid) = sat_ (S%sat, 2*S%phase%pass_seconds, &
			S%xover_params(1), S%xover_params(2), S%inclination)
	enddo
endif

! Now load all the "data variables"
do i = 1,nvar
	call get_var_2d (id_track + i, var(:,:,i), long_name(i))
enddo
fmt_string(len_trim(fmt_string)-3:) = ')'

! Exchange any correction if requested
do
	call getopt (optlist, optopt, optarg)
	if (optopt == '!') exit
	if (optopt == 'r') then
		j = index(optarg,'=')
		id_old = get_varid(optarg(:j-1)) - id_track
		id_new = get_varid(optarg(j+1:)) - id_track
		id_sla = get_varid('sla') - id_track
		if (optarg(:3) == 'alt') then	! Add change to altitude
			var(:,:,id_sla) = var(:,:,id_sla) + (var(:,:,id_new) - var(:,:,id_old))
		else	! Subtract change to corrections
			var(:,:,id_sla) = var(:,:,id_sla) - (var(:,:,id_new) - var(:,:,id_old))
		endif
	endif
enddo

! Close netCDF file
call nfs (nf90_close (ncid))

! Mask out data not in specified range
if (lat1 > lat0) where (var(1,:,-1) < lat0 .or. var(1,:,-1) > lat1) mask = .false.
if (lon1 > lon0) where (var(2,:,-1) < lon0 .or. var(2,:,-1) > lon1) mask = .false.
if (t1   > t0  ) where (var(1,:, 0) < t0   .or. var(1,:, 0) > t1   .or. &
						var(2,:, 0) < t0   .or. var(2,:, 0) > t1  ) mask = .false.
if (dt1  > dt0 ) where (abs(var(1,:,0)-var(2,:,0)) < dt0 .or. abs(var(1,:,0)-var(2,:,0)) > dt1) mask = .false.

! Mask out singles or duals
if (.not.singles) then
	mask = (trk(track(1,:))%satid /= trk(track(2,:))%satid)
else if (.not.duals) then
	mask = (trk(track(1,:))%satid == trk(track(2,:))%satid)
endif

! If editing is requested, mask based on sigma-editing of first variable only
! - Compute the standard deviation of the first
! - Edit based on outliers in first variable
if (edit > 0d0) then
	do j = 1,1	! For time being: first variable only
		call mean_variance (pack(var(1,:,j)-var(2,:,j),mask), mean, sigma)
		sigma = sqrt(sigma)*edit
		mask = mask .and. (abs(var(1,:,j)-var(2,:,j)-mean) <= sigma)
	enddo
endif

! Print number of xovers selected
nxo_out = count(mask)
write (*,610) nxo_out
610 format ('# Xovers sel''d  = ',i9)

! If no xovers skip the rest
if (nxo_out == 0) then
	deallocate (track, trk, var, stat, mask, long_name)
	return
endif

! Write column info
if (.not.xostat) then
	icol = 0
	do i = -2,-1
		icol = icol + 1
		write (*,611) icol,trim(long_name(i))
	enddo
	do i = 0,var0-1
		icol = icol + 2
		write (*,612) icol-1,icol,trim(long_name(i))
	enddo
	if (diff) then
		do i = var0,nvar
			icol = icol + 1
			write (*,611) icol,trim(long_name(i))
		enddo
	else
		do i = var0,nvar
			icol = icol + 2
			write (*,612) icol-1,icol,trim(long_name(i))
		enddo
	endif
endif
611 format ('# Column  ',i2,'    = ',a)
612 format ('# Columns ',i2,'-',i2,' = ',a)

! Sort first and last depending on order
select case (order)
case ('A')
	call reorder (mod(trk(track(1,:))%pass,2) < mod(trk(track(2,:))%pass,2))
	legs = 'ascending_pass descending_pass'
case ('a')
	call reorder (mod(trk(track(2,:))%pass,2) < mod(trk(track(1,:))%pass,2))
	legs = 'descending_pass ascending_pass'
case ('T')
	call reorder (var(1,:,0) < var(2,:,0))
	legs = 'later_measurement earlier_measurement'
case ('t')
	call reorder (var(2,:,0) < var(1,:,0))
	legs = 'earlier_measurement later_measurement'
case ('S')
	call reorder (trk(track(1,:))%satid < trk(track(2,:))%satid)
	legs = 'higher_satellite_id lower_satellite_id'
case ('s')
	call reorder (trk(track(2,:))%satid < trk(track(1,:))%satid)
	legs = 'lower_satellite_id higher_satellite_id'
case ('H')
	call reorder (sat(trk(track(1,:))%satid)%period < sat(trk(track(2,:))%satid)%period)
	legs = 'higher_satellite lower_satellite'
case ('h')
	call reorder (sat(trk(track(2,:))%satid)%period < sat(trk(track(1,:))%satid)%period)
	legs = 'lower_satellite higher_satellite'
case default ! Native order
	if (.not.duals .or. len_trim(satlist) <= 2) legs = 'ascending_pass descending_pass'
	if (.not.singles) legs = satlist
end select

! Specify order of values
if (diff) then
	write (*,615) 'Difference',trim(legs)
else
	write (*,615) 'Order     ',trim(legs)
endif
615 format ('# ',a,'    = ',a)

! Now collapse the data if needing difference
if (diff) var(1,:,var0:nvar) = var(1,:,var0:nvar) - var(2,:,var0:nvar)

! Print out data
if (stat_only) then
	! Skip
else if (diff) then
	do i = 1,nxo_in
		if (mask(i)) write (*,fmt_string) var(:,i,-1:var0-1),var(1,i,var0:nvar)
	enddo
else
	do i = 1,nxo_in
		if (mask(i)) write (*,fmt_string) var(:,i,:)
	enddo
endif

! Do statistics
do j = -1,nvar
	do k = 1,2
		stat(k,j,1) = minval(var(k,:,j),mask)
		stat(k,j,2) = maxval(var(k,:,j),mask)
		call mean_variance (pack(var(k,:,j),mask), stat(k,j,3), stat(k,j,5))
		stat(k,j,4) = sqrt(stat(k,j,3)**2+stat(k,j,5)*(nxo_out-1)/nxo_out)
		stat(k,j,5) = sqrt(stat(k,j,5))
	enddo
enddo

! Print statistics
if (xostat) then
	write (*,640) 'MIN','MAX','MEAN','RMS','STDDEV'
	write (*,645) trim(long_name(-2)),stat(1,-1,:)
	write (*,645) trim(long_name(-1)),stat(2,-1,:)
	do i = 0,nvar
		write (*,645) trim(long_name(i)),stat(1,i,:)
		if (.not.diff .or. i < var0) write (*,645) trim(long_name(i)),stat(2,i,:)
	enddo
else
	do i = 1,5
		write (*,650,advance='no') statname(i)
		if (diff) then
			write (*,fmt_string) stat(:,-1:var0-1,i),stat(1,var0:nvar,i)
		else
			write (*,fmt_string) stat(:,:,i)
		endif
	enddo
endif
640 format ('# ',t43,5a16)
645 format ('# ',a,t43,5f16.4)
650 format ('# ',a,' : ')

deallocate (track, trk, var, stat, mask, long_name)
end subroutine process

!***********************************************************************

subroutine synopsis
if (rads_version ('$Revision$','RADS crossover file lister')) return
if (xostat) then
	write (stderr,1300) 'stat'
else
	write (stderr,1300) 'list'
	write (stderr,1301)
endif
1300 format (/ &
'usage: radsxo',a,' [options] FILENAME ...' // &
'Required argument:' / &
'  FILENAME                  Name of input netCDF xover file'// &
'Optional arguments [options] are:'/ &
'  --lon=LON0,LON1           Specify longitude boundaries (deg)'/ &
'  --lat=LAT0,LAT1           Specify latitude boundaries (deg)'/ &
'  --t=T0,T1                 Specify time selection (optionally use --ymd=, --doy=,'/ &
'                            or --sec= for [YY]YYMMDD[HHMMSS], YYDDD, or SEC85)'/ &
'  --dt=[DTMIN,]DTMAX        Use only xovers with [DTMIN <] dt < DTMAX (days)'/ &
'  -d, --dual                Select duals satellite crossovers only'/ &
'  -s, --single              Select single satellite crossovers only'/ &
'  -l, --both-legs           Write out both xover values (default is differences)'/ &
'  -t, --both-times          Write out both times (default is difference)'/ &
'  -oTYPE, --order=TYPE      Order of the xover values (or difference), where TYPE is one of:'/ &
'                              A|a = ascending-descending or vv'/ &
'                              H|h = higher-lower satellite or vv'/ &
'                              S|s = higher-lower satellite ID or vv'/ &
'                              T|t = later-earlier measurement or vv')
1301 format (/ &
'  -n, --no-list             Do not print listing (print statistics only)')
stop
end subroutine synopsis

!***********************************************************************

function get_varid (name)
character(len=*), intent(in) :: name
integer(fourbyteint) :: get_varid
call nfs (nf90_inq_varid(ncid,name,get_varid))
end function get_varid

!***********************************************************************

subroutine get_var_1d (varid, array, long_name)
integer(fourbyteint), intent(in) :: varid
real(eightbytereal), intent(out) :: array(:)
character(len=*), intent(out) :: long_name
real(eightbytereal) :: val
character(len=rads_varl) :: fmt, units
call nfs (nf90_get_var (ncid, varid, array))
if (nf90_get_att (ncid, varid, 'scale_factor', val) == nf90_noerr) array = array * val
if (nf90_get_att (ncid, varid, 'add_offset', val) == nf90_noerr) array = array + val
if (nf90_get_att (ncid, varid, 'long_name', long_name) /= nf90_noerr) long_name = ''
if (nf90_get_att (ncid, varid, 'units', units) /= nf90_noerr) units = ''
if (nf90_get_att (ncid, varid, 'format', fmt) /= nf90_noerr) fmt = 'f0.3'
fmt_string = trim(fmt_string) // trim(fmt) // ',1x,'
long_name = trim(long_name)//' ['//trim(units)//']'
end subroutine get_var_1d

!***********************************************************************

subroutine get_var_2d (varid, array, long_name)
integer(fourbyteint), intent(in) :: varid
real(eightbytereal), intent(out) :: array(:,:)
character(len=*), intent(out) :: long_name
real(eightbytereal) :: val
character(len=rads_varl) :: fmt, units
call nfs (nf90_get_var (ncid, varid, array))
if (nf90_get_att (ncid, varid, 'scale_factor', val) == nf90_noerr) array = array * val
if (nf90_get_att (ncid, varid, 'add_offset', val) == nf90_noerr) array = array + val
if (nf90_get_att (ncid, varid, 'long_name', long_name) /= nf90_noerr) long_name = ''
if (nf90_get_att (ncid, varid, 'units', units) /= nf90_noerr) units = ''
if (nf90_get_att (ncid, varid, 'format', fmt) /= nf90_noerr) fmt = 'f0.3'
if (.not.diff .or. (var0 == 1 .and. long_name == 'time')) then
	fmt_string = trim(fmt_string) // trim(fmt) // ',1x,' // trim(fmt) // ',1x,'
else
	fmt_string = trim(fmt_string) // trim(fmt) // ',1x,'
endif
long_name = trim(long_name)//' ['//trim(units)//']'
end subroutine get_var_2d

!***********************************************************************

subroutine reorder (order)
logical, intent(in) :: order(:)
integer(fourbyteint) :: i, j
do i = 1,nxo_in
	if (order(i)) then
		do j = 0,nvar
			call flip(var(:,i,j))
		enddo
	endif
enddo
end subroutine reorder

subroutine flip (a)
real(eightbytereal), intent(inout) :: a(2)
real(eightbytereal) :: b
b = a(1)
a(1) = a(2)
a(2) = b
end subroutine flip

end program radsxolist
