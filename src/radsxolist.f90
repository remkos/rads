!-----------------------------------------------------------------------
! Copyright (c) 2011-2022  Remko Scharroo
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

!*radsxolist -- RADS crossover lister and statistics generator
!+
program radsxolist

! This program lists the contents of RADS NetCDF crossover files and
! computes statistics.
!-----------------------------------------------------------------------
use netcdf
use rads
use rads_misc
use rads_netcdf
use rads_time
integer(fourbyteint) :: ncid, i, icol, ios, nvar, nxo_in, nxo_out, ntrk, nxml = 0, &
	id_flag, id_satid, id_track, id_sla, id_new, id_old, id_offset, id_alt_rate, nvar_replace = 0
real(eightbytereal), allocatable :: var(:,:,:), stat(:,:,:), binstat(:,:,:), tbiasstat(:,:)
integer(fourbyteint), allocatable :: track(:,:), binnr(:,:), bincount(:)
character(len=rads_naml), allocatable :: long_name(:)
logical, allocatable :: is_alt(:), boz(:), mask(:), binmask(:)
type :: trk_
	real(eightbytereal) :: equator_lon, equator_time, start_time, end_time
	integer(twobyteint) :: nr_alt, nr_xover, satid, cycle, pass
end type
type(trk_), allocatable :: trk(:)
character(len=640) :: fmt_str

character(len=rads_varl) :: optopt
character(len=rads_naml) :: optarg, str_replace(10)
character(len=rads_cmdl) :: command, xml(10)
character(len=1) :: order = ''
character(len=4) :: statname(5) = (/ 'min ', 'max ', 'mean', 'rms ', 'std ' /)
integer(fourbyteint), parameter :: msat = 20
type :: sat_
	character(len=2) :: sat
	real(eightbytereal) :: period, altsig, orberr, inclination
	logical :: mask
end type sat_
type(sat_) :: sat(msat)
character(len=3*msat) :: satlist = '', use_sats = ''
type(rads_sat) :: S
character(len=*), parameter :: optlist = &
	'S:X:e::b::o:r:dsnltp sat: xml: lon: lat: dt: edit:: bin:: order: replace: dual single nolist both-legs both-times ' // &
	' time: ymd: doy: sec: check-flag: full-year tbias add-tbias: sub-tbias: pass-info dual-asc dual-des'

integer(fourbyteint) :: var0 = 0, check_flag = -1
logical :: diff = .true., stat_only, singles = .true., duals = .true., xostat, fullyear = .false., &
	ltbias = .false., pass_info = .false., dual_asc = .false., dual_des = .false.
real(eightbytereal) :: t0 = 0d0, t1 = 0d0, lon0 = 0d0, lon1 = 0d0, lat0 = 0d0, lat1 = 0d0, &
	dt0 = 0d0, dt1 = 0d0, edit = -1d0, bin = 0d0, tbias = 0d0

! Check operation mode
call getarg (0, command)
xostat = (index(command, 'radsxostat') > 0)
stat_only = xostat

! Initialize RADS or issue help
call synopsis
xml = ''

! Scan command line arguments
do
	call getopt (optlist, optopt, optarg)
	select case (optopt)
	case ('!') ! End of arguments
		exit
	case (':', '::')
		call rads_opt_error (optopt, optarg)
	case ('S', 'sat')
		use_sats = optarg(:3*msat)
	case ('X', 'xml')
		nxml = nxml + 1
		xml(nxml) = optarg
	case ('lon')
		read (optarg, *, iostat=ios) lon0,lon1
		if (ios /= 0) call rads_opt_error (optopt, optarg)
	case ('lat')
		read (optarg, *, iostat=ios) lat0,lat1
		if (ios /= 0) call rads_opt_error (optopt, optarg)
	case ('dt')
		read (optarg, *, iostat=ios) dt0,dt1
		if (ios > 0) call rads_opt_error (optopt, optarg)
		if (dt0 > dt1) then
			dt1 = dt0
			dt0 = 0d0
		endif
		dt0 = dt0 * 86400d0
		dt1 = dt1 * 86400d0
	case ('e', 'edit')
		edit = 3.5d0
		read (optarg, *, iostat=ios) edit
		if (ios > 0) call rads_opt_error (optopt, optarg)
	case ('dual-asc')
		singles = .false.
		dual_asc = .true.
	case ('dual-des')
		singles = .false.
		dual_des = .true.
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
	case ('b', 'bin')
		bin = 1d0
		read (optarg, *, iostat=ios) bin
		if (ios > 0) call rads_opt_error (optopt, optarg)
	case ('r', 'replace')
		nvar_replace = nvar_replace + 1
		if (nvar_replace > 10) call rads_exit ('Too many -r options')
		str_replace(nvar_replace) = optarg
	case ('full-year')
		fullyear = .true.
	case ('tbias')
		ltbias = .true.
	case ('add-tbias', 'sub-tbias')
		read (optarg, *, iostat=ios) tbias
		if (ios > 0) call rads_opt_error (optopt, optarg)
		tbias = tbias * 1d-3
		if (optopt == 'sub-tbias') tbias = -tbias
	case ('p', 'pass-info')
		pass_info = .true.
	case ('check-flag')
		read (optarg, *, iostat=ios) check_flag
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case (' ')
		! Skip filenames for now
	case ('time', 't:', 'sec', 'mjd', 'doy', 'ymd') ! Finally try date arguments
		call dateopt (optopt, optarg, t0, t1, iostat=ios)
		if (ios > 0) call rads_opt_error (optopt, optarg)
	end select
enddo

! Write header
call get_command (command, status=i)
if (xostat) then
	write (*, 600) 'Statistics', timestamp(), trim(command)
else
	write (*, 600) 'List', timestamp(), trim(command)
endif
600 format ('# ',a,' of RADS crossovers'/'# Created: ',a,' UTC: ',a)

! Some options exclude or imply others
if (bin > 0d0) then
	diff = .false.
	var0 = 1
	stat_only = .false.
	xostat = .false.
endif

! Now process all files
call getopt_reset ()
do
	call getopt (optlist, optopt, optarg)
	select case (optopt)
	case ('!') ! End of arguments
		exit
	case (' ')
		call process (optarg)
	end select
enddo

contains

!***********************************************************************

subroutine process (filename)
character(len=*), intent(in) :: filename
integer(fourbyteint) :: i, j, k, k1, k2, yy, mm, dd, trkid(2)
character(len=rads_naml) :: legs, varnm
real(eightbytereal) :: mean, sigma
real(eightbytereal), parameter :: day = 86400d0
type(rads_var), pointer :: ptr
logical :: old

! Open NetCDF file
600 format (/'# File name     = ',a)
if (nft (nf90_open (filename, nf90_nowrite, ncid))) then
	call rads_message ('error opening file '//trim(filename))
	return
endif
write (*, 600) trim(filename)

! Get the number of xovers and number of tracks
call nfs (nf90_inq_dimid (ncid, 'xover', i))
call nfs (nf90_inquire_dimension (ncid, i, len=nxo_in))
call nfs (nf90_inq_dimid (ncid, 'track', i))
call nfs (nf90_inquire_dimension (ncid, i, len=ntrk))

write (*, 605) nxo_in, ntrk
605 format ('# Xovers,tracks = ',2i9)
if (edit > 0d0) write (*, 606) edit
606 format ('# Edit sigma    = ',f9.3)

! Read all the "base variables" into memory
id_track = get_varid('track')
id_satid = get_varid('satid')
if (id_track == 4) then ! Old style xover file
	nvar = id_satid - 5
	id_offset = 4
else ! New style xover file
	nvar = id_track - 4
	id_offset = 3
endif
allocate (track(2,nxo_in), trk(ntrk), var(2,nxo_in,-1:nvar), stat(2,-1:nvar,5), mask(nxo_in), long_name(-2:nvar), &
		is_alt(nvar), boz(nvar))
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
old = nft (nf90_get_att (ncid, get_varid('satid'), 'flag_meanings', satlist))
if (old) satlist = 'g3 ss gs e1 tx pn e2 g1 j1 n1 j2 c2 sa'
do i = 1,len_trim(satlist),3
	if (old) then
		 if (.not.any(trk%satid == i/3+1)) cycle
	endif
	call rads_init (S, satlist(i:i+1), xml(1:nxml))
	sat(S%satid) = sat_ (S%sat, 2*S%phase%pass_seconds, &
		S%xover_params(1), S%xover_params(2), S%inclination, .true.)
	if (use_sats /= '') sat(S%satid)%mask = (index(use_sats, S%sat) > 0)
enddo

! Start format string with lat, lon, time
fmt_str = '(' // trim(S%lat%info%format) // ',1x,' // trim(S%lon%info%format) // ',1x,' // &
	trim(S%time%info%format) // ',1x,'
if (var0 == 1) fmt_str = trim(fmt_str) // trim(S%time%info%format) // ',1x,'

! Load all the "data variables" and complete the format string
id_sla = 0
id_alt_rate = 0
do i = 1,nvar
	call get_var_2d (i + id_offset, var(:,:,i), long_name(i))
	j = nf90_inquire_variable (ncid, i + id_offset, varnm)
	ptr => rads_varptr (S, varnm)
	if (diff) then
		fmt_str = trim(fmt_str) // trim(ptr%info%format) // ',1x,'
	else
		fmt_str = trim(fmt_str) // trim(ptr%info%format) // ',1x,' // trim(ptr%info%format) // ',1x,'
	endif
	boz(i) = ptr%info%boz_format
	is_alt(i) = ptr%info%standard_name == 'height_above_reference_ellipsoid' .or. &
		ptr%info%standard_name == 'sea_surface_height_above_sea_level'
	if (varnm == 'sla') id_sla = i
	if (varnm == 'alt_rate') id_alt_rate = i
enddo
fmt_str(len_trim(fmt_str)-3:) = ')'

! Try if we can generate a time tag bias
if (ltbias .and. (id_sla == 0 .or. id_alt_rate == 0)) &
	call rads_exit ('Cannot compute time tag bias without both SLA and ALT_RATE present')
if (tbias /= 0 .and. id_alt_rate == 0) &
	call rads_exit ('Cannot adjust for time tag bias without ALT_RATE present')

! Exchange any correction if requested
do i = 1,nvar_replace
	if (id_sla == 0) call rads_exit ('Cannot replace fields without SLA present')
	j = index(optarg,'=')
	id_old = get_varid(optarg(:j-1)) - id_offset
	id_new = get_varid(optarg(j+1:)) - id_offset
	if (optarg(:3) == 'alt') then	! Add change to altitude
		var(:,:,id_sla) = var(:,:,id_sla) + (var(:,:,id_new) - var(:,:,id_old))
	else	! Subtract change to corrections
		var(:,:,id_sla) = var(:,:,id_sla) - (var(:,:,id_new) - var(:,:,id_old))
	endif
enddo
id_flag = 0
if (check_flag >= 0) id_flag = get_varid('flag_alt_oper_mode') - id_offset

! Adjust for a time tag bias (tbias in ms) on TIME, SLA and all ALT fields, if requested
if (tbias /= 0d0) then
	var(:,:,0) = var(:,:,0) + tbias
	do i = 1,nvar
		if (is_alt(i)) var(:,:,i) = var(:,:,i) + tbias * var(:,:,id_alt_rate)
	enddo
endif

! Close NetCDF file
call nfs (nf90_close (ncid))

! Mask out singles or duals
if (.not.singles) then
	mask = (trk(track(1,:))%satid /= trk(track(2,:))%satid)
else if (.not.duals) then
	mask = (trk(track(1,:))%satid == trk(track(2,:))%satid)
else
	mask = .true.
endif

! Mask out ascending or descending passes on first satellite
if (dual_asc) where (modulo(trk(track(1,:))%pass,2) == 0) mask = .false.
if (dual_des) where (modulo(trk(track(1,:))%pass,2) == 1) mask = .false.

! Mask out data not in specified range
if (lat1 > lat0) where (var(1,:,-1) < lat0 .or. var(1,:,-1) > lat1) mask = .false.
if (lon1 > lon0) where (var(2,:,-1) < lon0 .or. var(2,:,-1) > lon1) mask = .false.
if (t1   > t0  ) where (var(1,:, 0) < t0   .or. var(1,:, 0) > t1   .or. &
						var(2,:, 0) < t0   .or. var(2,:, 0) > t1  ) mask = .false.
if (dt1  > dt0 ) where (abs(var(1,:,0)-var(2,:,0)) < dt0 .or. abs(var(1,:,0)-var(2,:,0)) > dt1) mask = .false.

! Select only specified satellites
if (use_sats /= '') &
	where (.not.sat(trk(track(1,:))%satid)%mask .or. .not.sat(trk(track(2,:))%satid)%mask) mask = .false.

! Mask out based on data flag (for use with CryoSat)
if (id_flag > 0) where (var(1,:,id_flag) /= check_flag) mask = .false.

! If editing is requested, mask based on sigma-editing of first variable only
! - Compute the standard deviation of the first
! - Edit based on outliers in first variable
if (edit > 0d0) then
	do j = 1,1	! For time being: first variable only
		call mean_variance (pack(var(1,:,j)-var(2,:,j),mask), mean, sigma)
		sigma = sqrt(sigma)*edit
		where (abs(var(1,:,j)-var(2,:,j)-mean) > sigma) mask = .false.
	enddo
endif

! Print number of xovers selected
nxo_out = count(mask)
write (*,610) nxo_out
610 format ('# Xovers sel''d  = ',i9)

! If no xovers: skip the rest
if (nxo_out == 0) then
	deallocate (track, trk, var, stat, mask, long_name, is_alt, boz)
	return
endif

! If stats are to be binned, assign bins
if (bin > 0d0) then
	allocate (binnr(2,nxo_in), binmask(nxo_in))
	binnr = floor(var(:,:,0)/bin/day)
endif

! Write column info
if (.not.xostat) then
	icol = 0
	if (bin > 0d0) then
		write (*,611) 1, 'date [YYMMDD]'
		write (*,611) 2, 'number of crossovers'
		icol = 2
		write (*,613) 'Mean of the following variables:'
	endif
	do i = -2,-1
		icol = icol + 1
		write (*,611) icol,trim(long_name(i))
	enddo
	do i = 0,var0-1
		icol = icol + 2
		write (*,612) icol-1,icol,trim(long_name(i))
	enddo
	if (bin > 0d0) write (*,613) 'Mean and std dev of the xover differences of variables:'
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
	if (bin > 0d0 .and. ltbias) then
		icol = icol + 2
		write (*,612) icol-1,icol,'time tag bias [ms]'
	endif
	if (pass_info) then
		icol = icol + 2
		write (*,612) icol-1,icol,'satellite id'
		icol = icol + 2
		write (*,612) icol-1,icol,'cycle number'
		icol = icol + 2
		write (*,612) icol-1,icol,'pass number'
	endif
endif
611 format ('# Column  ',i2,'    = ',a)
612 format ('# Columns ',i2,'-',i2,' = ',a)
613 format ('# ',a)

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
! Collapse data if needed
if (diff .or. bin > 0d0) then
	write (*,615) 'Difference',trim(legs)
	var(1,:,var0:nvar) = var(1,:,var0:nvar) - var(2,:,var0:nvar)
else
	write (*,615) 'Order     ',trim(legs)
endif
615 format ('# ',a,'    = ',a)

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

! Print out data (warning, this alters variables with boz format)
if (.not.stat_only .and. bin == 0d0) then
	do i = 1,nxo_in
		if (.not.mask(i)) cycle
		call print_data (var(:,i,:), .not.pass_info)
		if (.not.pass_info) cycle
		trkid = track(:,i)
		write (*,632) trk(trkid)%satid,trk(trkid)%cycle,trk(trkid)%pass
	enddo
endif

! Print statistics
if (bin > 0d0) then
	! Statistics per bin
	k1 = minval(binnr)
	k2 = maxval(binnr)
	allocate (bincount(k1:k2), binstat(2,-1:nvar,k1:k2), tbiasstat(2,k1:k2))
	bincount = 0
	binstat = 0d0
	tbiasstat = 0d0
	do i = 1,nxo_in
		if (.not.mask(i)) cycle
		call update_stat (binnr(1,i),var(:,i,:))
		if (binnr(1,i) /= binnr(2,i)) call update_stat (binnr(2,i),var(:,i,:))
	enddo
	do k = k1,k2
		if (bincount(k) == 0) cycle
		binstat(2,var0:nvar,k) = sqrt((binstat(2,var0:nvar,k) - binstat(1,var0:nvar,k)**2/bincount(k))/(bincount(k)-1))
		binstat(1,:,k) = binstat(1,:,k) / bincount(k)
		binstat(2,-1:var0-1,k) = binstat(2,-1:var0-1,k) / bincount(k)
		call mjd2ymd (nint(k*bin)+46066,yy,mm,dd)
		if (.not.fullyear) yy = modulo(yy,100)
		write (*,630,advance='no') yy,mm,dd,bincount(k)
		call print_data (binstat(:,:,k), .not.ltbias)
		if (ltbias) write(*,631) tbiasstat(2,k)/tbiasstat(1,k)*1d3, sqrt(1d0/tbiasstat(1,k))*1d3
	enddo
	deallocate (bincount, binstat, tbiasstat)
else if (xostat) then
	write (*,640) 'MIN','MAX','MEAN','RMS','STDDEV'
	write (*,645) trim(long_name(-2)),stat(1,-1,:)
	write (*,645) trim(long_name(-1)),stat(2,-1,:)
	do i = 0,nvar
		write (*,645) trim(long_name(i)),stat(1,i,:)
		if (.not.diff .or. i < var0) write (*,645) trim(long_name(i)),stat(2,i,:)
	enddo
	if (ltbias .and. diff) then
		call solve_tbias (pack(var(1,:,id_sla),mask),pack(var(1,:,id_alt_rate),mask),mean,sigma)
		write (*,646) 'time tag bias [ms]', mean*1d3, sigma*1d3
	endif
else
	do i = 1,5
		write (*,650,advance='no') statname(i)
		call print_data (stat(:,:,i), .true.)
	enddo
endif
630 format (3i0.2,i9,1x)
631 format (1x,2f8.4)
632 format (2i3,2i4,2i5)
640 format ('# ',t43,5a16)
645 format ('# ',a,t43,5f16.4)
646 format ('# ',a,t43,32x,f16.4,16x,f16.4)
650 format ('# ',a,' : ')

deallocate (track, trk, var, stat, mask, long_name, is_alt, boz)
if (bin > 0d0) deallocate (binnr, binmask)
end subroutine process

!***********************************************************************

subroutine synopsis
if (rads_version ('RADS crossover file lister')) return
if (xostat) then
	write (*,1300) 'stat'
else
	write (*,1300) 'list'
	write (*,1301)
endif
1300 format (/ &
'usage: radsxo',a,' [options] [--] FILENAME ...' // &
'Required argument:' / &
'  FILENAME                  Name of input NetCDF xover file'// &
'Optional arguments [options] are:'/ &
'  -S, --sat SAT1[,SAT2,..]  Comma-separated list of satellites (default: all)'/ &
'  -X, --xml XMLFILE         Load XMLFILE in addition to RADS defaults'/ &
'  --lon LON0,LON1           Specify longitude boundaries (deg)'/ &
'  --lat LAT0,LAT1           Specify latitude boundaries (deg)'/ &
'  --time T0,T1              Specify time selection (optionally use --ymd, --doy,'/ &
'                            or --sec for [YY]YYMMDD[HHMMSS], YYDDD, or SEC85)'/ &
'  --dt [DTMIN,]DTMAX        Use only xovers with [DTMIN <] dt < DTMAX (days)'/ &
'  -d, --dual                Select duals satellite crossovers only'/ &
'  -s, --single              Select single satellite crossovers only'/ &
'  -l, --both-legs           Write out both xover values (default is differences)'/ &
'  -t, --both-times          Write out both times (default is difference)'/ &
'  --tbias                   Determine time tag bias (requires variables sla and alt_rate)'/ &
'  --add-tbias VAL           Add (apply) effect of time tag bias (VAL in ms)'/ &
'                            (requires variable alt_rate)'/ &
'  --sub-tbias VAL           Subtract (remove) effect of time tag bias (VAL in ms)'/ &
'                            (requires variable alt_rate)'/ &
'  -o, --order TYPE          Order of the xover values (or difference), where TYPE is one of:'/ &
'                              A|a = ascending-descending or vv'/ &
'                              H|h = higher-lower satellite or vv'/ &
'                              S|s = higher-lower satellite ID or vv'/ &
'                              T|t = later-earlier measurement or vv'/ &
'  --dual-asc                Select only ascending passes on first satellite'/ &
'  --dual-des                Select only descending passes on first satellite'/ &
'  -e, --edit [VAL]          Edit the data beyond VAL [3.5] times the std dev of the 1st variable')
1301 format ( &
'  -b, --bin [DAYS]          Bin data by number of DAYS [1] and print mean and std dev.'/ &
'  --full-year               Print data as YYYYMMDD instead of YYMMDD'/ &
'  -p, --pass-info           Write out satellite ID, cycle and pass number for every record'/ &
'  -n, --no-list             Do not print listing (print overall statistics only)')
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
character(len=rads_varl) :: units
call nfs (nf90_get_var (ncid, varid, array))
if (nf90_get_att (ncid, varid, 'scale_factor', val) == nf90_noerr) array = array * val
if (nf90_get_att (ncid, varid, 'add_offset', val) == nf90_noerr) array = array + val
if (nf90_get_att (ncid, varid, 'long_name', long_name) /= nf90_noerr) long_name = ''
if (nf90_get_att (ncid, varid, 'units', units) /= nf90_noerr) units = ''
long_name = trim(long_name)//' ['//trim(units)//']'
end subroutine get_var_1d

!***********************************************************************

subroutine get_var_2d (varid, array, long_name)
integer(fourbyteint), intent(in) :: varid
real(eightbytereal), intent(out) :: array(:,:)
character(len=*), intent(out) :: long_name
real(eightbytereal) :: val
character(len=rads_varl) :: units
call nfs (nf90_get_var (ncid, varid, array))
if (nf90_get_att (ncid, varid, 'scale_factor', val) == nf90_noerr) array = array * val
if (nf90_get_att (ncid, varid, 'add_offset', val) == nf90_noerr) array = array + val
if (nf90_get_att (ncid, varid, 'long_name', long_name) /= nf90_noerr) long_name = ''
if (nf90_get_att (ncid, varid, 'units', units) /= nf90_noerr) units = ''
long_name = trim(long_name)//' ['//trim(units)//']'
end subroutine get_var_2d

!***********************************************************************
! Print a single line of data (take care of flags)

subroutine print_data (x, advance)
real(eightbytereal), intent(inout) :: x(2,-1:nvar)
logical, intent(in) :: advance
integer :: j
character(len=3) :: yes_no
! In the following we would have liked to use nint8 in the transfer
! function, but an 8-byte integer is not guaranteed to work, so we use
! padded 4-byte integers instead.
if (advance) then
	yes_no = 'yes'
else
	yes_no = 'no'
endif
do j = 1,nvar
	if (boz(j)) call bit_transfer (x(:,j))
enddo
if (diff) then
	write (*,fmt_str,advance=yes_no) x(:,-1:var0-1), x(1,var0:nvar)
else
	write (*,fmt_str,advance=yes_no) x(:,:)
endif
end subroutine print_data

!***********************************************************************

subroutine reorder (order)
logical, intent(in) :: order(:)
integer(fourbyteint) :: i, j, k
real(eightbytereal) :: a
do i = 1,nxo_in
	if (order(i)) then
		do j = 0,nvar
			a = var(1,i,j)
			var(1,i,j) = var(2,i,j)
			var(2,i,j) = a
		enddo
		k = track(1,i)
		track(1,i) = track(2,i)
		track(2,i) = k
	endif
enddo
end subroutine reorder

!***********************************************************************

subroutine solve_tbias (sla, alt_rate, tbias, sig_tbias)
real(eightbytereal), intent(in) :: sla(:), alt_rate(:)
real(eightbytereal), intent(out) :: tbias, sig_tbias
real(eightbytereal) :: atwa, atwb
real(eightbytereal), parameter :: w = 1d2 ! Measurement sigma = 10 cm
atwa = sum(alt_rate * w * alt_rate)
atwb = sum(sla * w * alt_rate)
tbias = atwb / atwa
sig_tbias = sqrt(1d0/atwa)
end subroutine solve_tbias

!***********************************************************************

subroutine update_stat (k, val)
integer(fourbyteint), intent(in) :: k
real(eightbytereal), intent(in) :: val(2,-1:nvar)
real(eightbytereal), parameter :: w = 1d2 ! Measurement sigma = 10 cm
binstat(1,:,k) = binstat(1,:,k) + val(1,:)
binstat(2,-1:var0-1,k) = binstat(2,-1:var0-1,k) + val(2,-1:var0-1)
binstat(2,var0:nvar,k) = binstat(2,var0:nvar,k) + val(1,var0:nvar) ** 2
bincount(k) = bincount(k) + 1
if (ltbias) tbiasstat(:,k) = tbiasstat(:,k) + w * val(1,id_alt_rate) * (/ val(1,id_alt_rate), val(1,id_sla) /)
end subroutine update_stat
end program radsxolist
