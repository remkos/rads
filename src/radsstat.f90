!-----------------------------------------------------------------------
! Copyright (c) 2011-2025  Remko Scharroo
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

program radsstat

! This program reads the RADS data base and computes statistics
! by pass, cycle or day of a number of RADS data variables.
!
! Output can be either as an ASCII table or as a NetCDF file.
! The output will always include both the mean and standard deviation of
! each of the selected variables per the selected period (N cycles, N
! passes or N days). Optionally also minimum and maximum values can be
! produced.
!
! Statistics are produced in one of various ways:
! - Weighting boxes of given size by their respective area
! - Equal weight to each measurement
! - Weight measurements by cosine of latitude
! - Use inclination-dependent weight
!
! When choosing statistics spanning N days, the N-day intervals are
! at fixed locations, starting at 1.0 Jan 1985. For example, 7-day
! periods will thus always start on Tuesday (as is 1 Jan 1985).
!
! usage: radsstat [RADS_options] [options]
!-----------------------------------------------------------------------
use netcdf
use rads
use rads_time
use rads_misc
use rads_netcdf
integer(fourbyteint) ::  lstat = 2, reject = -1
character(len=*), parameter :: wtype(0:3)=(/ &
	'box weight                  ', 'constant weight             ', &
	'area weighted               ', 'inclination-dependent weight'/)
character(len=rads_strl) :: format_string
character(len=rads_cmdl) :: filename = ''
integer(fourbyteint), parameter :: period_day=1, period_pass=2, period_cycle=3
integer(fourbyteint) :: nr, minnr=2, cycle, pass, i, output_format=0, &
	period=period_day, wmode=0, nx, ny, kx, ky, ios, sizes(2), ncid, varid(2), grpid(4) = 0
real(eightbytereal), allocatable :: lat_w(:)
real(eightbytereal) :: sini, step=1d0, x0, y0, res(2)=(/3d0,1d0/), start_time, dt = 0d0
type :: stat
	integer(fourbyteint) :: nr
	real(eightbytereal) :: wgt, z(4)
end type
type(stat), allocatable :: box(:,:,:), tot(:)
type(rads_sat) :: S(2)
type(rads_pass) :: P(2), Pout
logical :: ascii = .true., fullyear = .false., boz_format = .false., force = .false.
logical :: echofilepaths = .false. , groups = .false., mask(4) = .false.

! Initialize RADS or issue help
call synopsis
call rads_set_options ('ab::c::d::flmo::p::q::r::s ' // &
	'dt echo-file-paths force format-cycle format-day format-pass groups full-year min: mean-only minmax no-stddev ' // &
	'output:: reject-on-nan: res:')
call rads_init (S)
if (any(S%error /= rads_noerr)) call rads_exit ('Fatal error')

! Determine how many satellites and cycles.
! Also check that the same number of variables have been selected for each satellite
do i = 1,2
	if (S(i)%sat == '') exit
	if (S(i)%nsel == 0) call rads_parse_varlist (S(i), 'sla')
	if (S(i)%nsel /= S(1)%nsel) call rads_exit ('Unequal amount of variables requested for different missions')
	if (S(i)%cycles(3) /= S(1)%cycles(3)) call rads_exit ('Cycle step size should be the same for all missions')
	if (S(i)%inclination /= S(1)%inclination) &
		call rads_exit ('Missions need to have the same inclination to be considered collinear')
	dt = max(dt,S(i)%dt1hz)
enddo

! Scan command line arguments
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('d')
		period = period_day
		read (rads_opt(i)%arg, *, iostat=ios) step
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('p')
		period = period_pass
		read (rads_opt(i)%arg, *, iostat=ios) step
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('c')
		if (rads_opt(i)%arg(1:1) == '/') then ! When stepping by fraction of cycles, redefine it as a number of passes.
			period = period_pass
			read (rads_opt(i)%arg(2:), *, iostat=ios) step
			if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
			step = S(1)%passes(2) / step
		else
			period = period_cycle
			read (rads_opt(i)%arg, *, iostat=ios) step
			if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
		endif
	case ('b')
		wmode = 0
		call read_val (rads_opt(i)%arg, res, '/-+x', iostat=ios)
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('m')
		wmode = 1
	case ('a')
		wmode = 2
	case ('s')
		wmode = 3
	case ('no-stddev')
		if (lstat == 2 .or. lstat == 4) lstat = lstat - 1
	case ('mean-only')
		lstat = 0
	case ('l', 'minmax')
		if (lstat == 1 .or. lstat == 2) lstat = lstat + 2
	case ('format-cycle')
		output_format = period_cycle
	case ('format-day')
		output_format = period_day
	case ('format-pass')
		output_format = period_pass
	case ('groups')
		groups = .true.
	case ('full-year')
		fullyear = .true.
	case ('echo-file-paths')
		echofilepaths = .true.
	case ('min')
		read (rads_opt(i)%arg, *, iostat=ios) minnr
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('res')
		call read_val (rads_opt(i)%arg, res, '/-+x', iostat=ios)
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('f', 'force')
		force = .true.
	case ('o', 'output')
		filename = rads_opt(i)%arg
		if (filename == '') filename = 'radsstat.nc'
	case ('r', 'reject-on-nan')
		call parse_r_option (S(1), rads_opt(i)%opt, rads_opt(i)%arg)
	case ('dt')
		read (rads_opt(i)%arg, *, iostat=ios) dt
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	end select
enddo

! At this point, lstat indicates the number of output columns, being either
! 0 = mean only in 1-D (netCDF) array
! 1 = mean only
! 2 = mean. srddev
! 3 = mean, min, amx
! 4 = mean, stddev, min, max
mask(1) = .true.
if (lstat == 2 .or. lstat == 4) mask(2) = .true.
if (lstat >= 3) mask(3:4) = .true.

! If no filename is given, we have ASCII output
ascii = (filename == '')

if (period == period_day) step = step * 86400d0 ! Convert from days to seconds

! Determine output format if undefined
if (output_format == 0) output_format = period

! If SLA is among the selected variables, remember index
! Also check if we have boz-formats
do i = 1,S(1)%nsel
	if (reject == -1 .and. S(1)%sel(i)%info%datatype == rads_type_sla) reject = i
	if (S(1)%sel(i)%info%boz_format) boz_format = .true.
enddo

! Set up the boxes
x0 = S(1)%lon%info%limits(1)
y0 = S(1)%lat%info%limits(1)
nx = int((S(1)%lon%info%limits(2)-x0)/res(1)+0.999d0)
ny = int((S(1)%lat%info%limits(2)-y0)/res(2)+0.999d0)
allocate (box(0:S(1)%nsel,nx,ny),tot(0:S(1)%nsel),lat_w(ny))

! Set up the weights
sini = sin(S(1)%inclination*rad)
if (wmode == 1) then
	lat_w = 1d0
else if (wmode == 3) then
	forall (ky=1:ny) lat_w(ky) = sqrt(max(1d-2,1d0-(sin((y0+(ky-0.5d0)*res(2))*rad)/sini)**2))
else
	forall (ky=1:ny) lat_w(ky) = cos((y0+(ky-0.5d0)*res(2))*rad)
endif

! Initialize statistics
call init_stat
if (ascii) then
	call ascii_header
else
	call netcdf_header
endif

! Start looping through cycles and passes
do cycle = S(1)%cycles(1), S(1)%cycles(2), S(1)%cycles(3)
	! Process passes one-by-one
	do pass = S(1)%passes(1), S(1)%passes(2), S(1)%passes(3)
		call rads_open_pass (S(1), P(1), cycle, pass, echofilepaths=echofilepaths)
		! Process the pass data
		if (P(1)%ndata > 0) call process_pass (P(1)%ndata, S(1)%nsel)
		! Print the statistics at the end of the data pass (if requested)
		if (period == period_pass .and. floor((pass-1)/step) /= floor(pass/step)) call output_stat
		! Close the pass file
		call rads_close_pass (S(1), P(1))
	enddo

	! Print the statistics at the end of the cycle (if requested)
	if (period == period_cycle .and. nint(modulo(dble(cycle),step)) == 0) call output_stat
enddo

! Flush the statistics and close RADS4
cycle = S(1)%cycles(2)
pass = S(1)%passes(2)
call output_stat
if (rads_verbose >= 1) call rads_stat (S)
call rads_end (S)
deallocate (box,tot,lat_w)

! Close NetCDF file
if (.not.ascii) call nfs (nf90_close (ncid))

contains

!***********************************************************************

subroutine parse_r_option (S, opt, arg)
type(rads_sat), intent(in) :: S
character(len=*), intent(in) :: opt, arg
integer(fourbyteint) :: i
select case (arg)
case ('n', 'any')
	reject = -2
case ('', '0', 'none')
	reject = 0
case default
	if (arg(1:1) >= '0' .and. arg(1:1) <= '9') then
		reject = 0
		read (arg, *, iostat=i) reject
	else
		do i = 1,S%nsel
			if (arg == S%sel(i)%name .or. arg == S%sel(i)%info%name) then
				reject = i
				return
			endif
		enddo
		call rads_exit ('Option -'//trim(opt)//' <varname> does not refer to variable specified on -V option')
	endif
end select
end subroutine parse_r_option

!***********************************************************************

subroutine synopsis
if (rads_version ('Print RADS statistics per cycle, pass or day(s)')) return
call rads_synopsis
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r, --reject-on-nan VAR   Reject records if variable VAR is NaN (default: reject if SLA is NaN)'/ &
'  -r NR                     Reject records if data item number NR on -V specifier is NaN'/ &
'  -r 0, -r none, -r         Do not reject measurement records with NaN values'/ &
'                            In all cases above, all non-NaN values per variable are averaged'/ &
'  -r n, -r any              Reject measurement records if any value is NaN'/ &
'                      Note: If no -r option is given, -r sla is assumed'/ &
'  -c [N]                    Statistics per cycle or N cycles or use /N to specify a fraction of cycles'/ &
'  -p [N]                    Statistics per pass or N passes'/ &
'  -d [N]                    Statistics per day (default) or N days'/ &
'  -b [DX,DY]                Average by boxes with size (default on: 3x1 degrees)'/ &
'  -m                        Give all measurements equal weight'/ &
'  -a                        Weight measurements by cosine of latitude'/ &
'  -s                        Use inclination-dependent weight'/ &
'  -l, --minmax              Output min and max in addition to mean and stddev'/ &
'  --no-stddev               Output mean values only' / &
'  --mean-only               Output mean values only in 1-D array (applies to NetCDF output only)' / &
'  --dt DT                   Set minimum bin size in seconds (default is determined by satellite)'/ &
'  --full-year               Write date as YYYYMMDD instead of the default YYMMDD'/ &
'  --min MINNR               Minimum number of measurements per statistics record (default: 2)'/ &
'  --res DX,DY               Size of averaging boxes in degrees (default: 3,1)'/ &
'  --echo-file-paths         Write to STDOUT the paths of the files before being read'/ &
'  --format-cycle            Output format starts with CYCLE, [YY]YYMMDD (default with -c)'/ &
'  --format-day              Output format starts with [YY]YYMMDD (default with -d)'/ &
'  --format-pass             Output format starts with CYCLE, PASS (default with -p)'/ &
'  --groups                  Create groups (applies to NetCDF output only)' / &
'  -f, --force               Force comparison, even when missions are not considered collinear'/ &
'  -o, --output [OUTNAME]    Create NetCDF output instead of ASCII (default output'/ &
'                            filename is "radsstat.nc")')
stop
end subroutine synopsis

!***********************************************************************
! Process data for a single pass

subroutine process_pass (ndata, nsel)
integer(fourbyteint), intent(in) :: ndata, nsel
real(eightbytereal) :: z(ndata,0:nsel)
real(eightbytereal), allocatable :: a(:)
integer, allocatable :: idx(:)
integer :: i, j
real(eightbytereal) :: equator_time_2

! Read the data for this pass
z(:,0) = P(1)%tll(:,1)	! Store time
do j = 1,nsel
	call rads_get_var (S(1), P(1), S(1)%sel(j), z(:,j))
enddo

if (S(2)%sat /= '') then
	call rads_open_pass (S(2), P(2), S(2)%cycles(1) - S(1)%cycles(1) + cycle, pass, echofilepaths=echofilepaths)
	if (force) then
		! Skip the next check, do collinear anyway
	else if (S(2)%phase%passes /= S(1)%phase%passes .or. S(2)%phase%ref_lon /= S(1)%phase%ref_lon) then
		! Pass ranges should be the same for all satellites, otherwise we do not have collinear tracks
		call rads_exit ('Satellite missions '//S(2)%sat//'/'//trim(S(2)%phase%name)// &
			' and '//S(1)%sat//'/'//trim(S(1)%phase%name)//' are not collinear')
	endif
	if (P(2)%ndata > 0) then
		allocate (a(0:P(2)%ndata), idx(P(1)%ndata))
		a(0) = nan
		! If the equator times are almost the same, use the same
		equator_time_2 = P(2)%equator_time
		if (abs(P(1)%equator_time - P(2)%equator_time) < 1d0) equator_time_2 = P(1)%equator_time
		! The little delta of 1 microsecond is to avoid bin jump at equator when it is xx.5000000
		call make_idx (P(1)%ndata, nint((P(1)%tll(:,1) - P(1)%equator_time + 1d6)/dt), &
					   P(2)%ndata, nint((P(2)%tll(:,1) - equator_time_2    + 1d6)/dt), idx)
		do j = 1,nsel
			call rads_get_var (S(2), P(2), S(2)%sel(j), a(1:))
			z(:,j) = z(:,j) - a(idx(:))
		enddo
		deallocate (a, idx)
	else
		z(:,1:) = nan
	endif
	call rads_close_pass (S(2), P(2))
endif

! Update the statistics with data in this pass
do i = 1,ndata
	if (reject > 0) then
		if (isnan_(z(i,reject))) cycle ! Reject if selected variable is NaN
	else if (reject == -2) then
		if (any(isnan_(z(i,:)))) cycle ! Reject if any variable is NaN
	endif

	if (period == period_day) then	! If "daily" statistics are requested
		if (P(1)%tll(i,1) >= start_time + step) call output_stat	! Output stat when beyond end of "day"
		! First call (after start or statistics reset) sets the start time to integer multiple of "step"
		if (isnan_(start_time)) start_time = floor(P(1)%tll(i,1)/step) * step
	else
		! First call (after start or statistics reset) saves the start time
		if (isnan_(start_time)) start_time = P(1)%tll(i,1)
	endif

	! Update the box statistics
	kx = floor((P(1)%tll(i,3)-x0)/res(1) + 1d0)
	ky = floor((P(1)%tll(i,2)-y0)/res(2) + 1d0)
	kx = max(1,min(kx,nx))
	ky = max(1,min(ky,ny))
	do j = 0,nsel
		call update_stat (box(j,kx,ky), z(i,j))
	enddo
	nr = nr + 1
enddo
end subroutine process_pass

!***********************************************************************
! Create index from pass A to pass B by linking the common times

subroutine make_idx (na, ta, nb, tb, idx)
integer(fourbyteint), intent(in) :: na, ta(:), nb, tb(:)
integer(fourbyteint), intent(out) :: idx(:)
integer(fourbyteint) :: ia, ib
idx = 0
ib = 1
do ia = 1,na
	do while (tb(ib) < ta(ia) .and. ib < nb)
		ib = ib + 1
	enddo
	if (tb(ib) == ta(ia)) idx(ia) = ib
	if (ib < nb .and. tb(ib+1) == tb(ib)) ib = ib + 1 ! Encourage indexing to the next point if two points are in the same bin
!	write (*,*) ia, idx(ia), ia - idx(ia), ta(ia), tb(ia), P(1)%tll(ia,1), P(2)%tll(ia,1)
enddo
! write (*,*) P(1)%equator_time, P(2)%equator_time
end subroutine make_idx

!***********************************************************************
! Update statistics inside a given box

subroutine update_stat (box, z)
type(stat), intent(inout) :: box
real(eightbytereal), intent(in) :: z
if (isnan_(z)) return
box%nr   = box%nr   + 1
box%wgt  = box%wgt  + 1d0
box%z(1) = box%z(1) + z
box%z(2) = box%z(2) + z*z
box%z(3) = min(box%z(3), z)
box%z(4) = max(box%z(4), z)
end subroutine update_stat

!***********************************************************************
! Output statistics for one batch of data

subroutine output_stat
integer(fourbyteint) :: j, k, yy, mm, dd, start(2) = 1
real(eightbytereal) :: w
type(rads_var), pointer :: var

! If not enough data available, simply reset the statistics
if (nr < minnr) then
	call init_stat
	return
endif

! Init global stats
tot = stat(0, 0d0, (/ 0d0, 0d0, huge(1d0), -huge(1d0) /))

! Cycle trough all boxes and determine overall weighted mean
do ky=1,ny
	do kx=1,nx
		if (box(0,kx,ky)%nr == 0) cycle ! Skip empty boxes

		! Determine the weight
		if (wmode == 0) then
			w = lat_w(ky) / box(0,kx,ky)%wgt
		else
			w = lat_w(ky)
		endif

		! Update overall statistics
		tot(:)%nr   = tot(:)%nr   + box(:,kx,ky)%nr
		tot(:)%wgt  = tot(:)%wgt  + w * box(:,kx,ky)%wgt
		tot(:)%z(1) = tot(:)%z(1) + w * box(:,kx,ky)%z(1)
		tot(:)%z(2) = tot(:)%z(2) + w * box(:,kx,ky)%z(2)
		tot(:)%z(3) = min(tot(:)%z(3), box(:,kx,ky)%z(3))
		tot(:)%z(4) = max(tot(:)%z(4), box(:,kx,ky)%z(4))
	enddo
enddo

! Divide by total weight to compute overall weighted mean and standard deviation

tot%z(1) = tot%z(1) / tot%wgt
where (tot%nr > 1)
	tot%z(2) = sqrt(max(tot%z(2) / tot%wgt - tot%z(1)*tot%z(1), 0d0) * tot%nr / (tot%nr - 1d0)) ! max() avoids tiny negatives
else where
	tot%z(2) = nan
end where
where (tot%nr == 0)
	tot%z(3) = nan
	tot%z(4) = nan
end where

! Write out statistics in ASCII
if (ascii) then
	call mjd2ymd(floor(start_time/86400d0)+46066,yy,mm,dd)
	if (.not.fullyear) yy = modulo(yy,100)
	select case (output_format)
	case (period_day)
		write (*,600,advance='no') yy,mm,dd
	case (period_pass)
		write (*,601,advance='no') cycle,pass
	case default
		write (*,602,advance='no') cycle,yy,mm,dd
	end select
	if (boz_format) then
		do j=1,S(1)%nsel
			if (.not.S(1)%sel(j)%info%boz_format) cycle
			do k = 1,4
				call bit_transfer (tot(j)%z(k))
			enddo
		enddo
	endif
	write (*,format_string) nr, tot(0)%z(1), (pack(tot(j)%z(:),mask),j=1,S(1)%nsel)

! Write output to NetCDF if requested
else
	call nfs (nf90_put_var (ncid, varid(1), nr, start(2:2)))
	Pout%fileinfo(1) = rads_file (ncid, filename)
	if (output_format /= period_day) then
		var => rads_varptr (S(1), 'cycle')
		call rads_put_var (S(1), Pout, var, (/ dble(cycle) /), start(2:2))
	endif
	if (output_format == period_pass) then
		var => rads_varptr (S(1), 'pass')
		call rads_put_var (S(1), Pout, var, (/ dble(pass) /), start(2:2))
	endif
	call rads_put_var (S(1), Pout, S(1)%time, (/ tot(0)%z(1) /), start(2:2))
	if (groups .or. lstat == 0) then
		do k = 1,4
			Pout%fileinfo(1) = rads_file (grpid(k), filename)
			if (mask(k)) then
				do j = 1,S(1)%nsel
					call rads_put_var (S(1), Pout, S(1)%sel(j), (/ tot(j)%z(k) /), start(2:2))
				enddo
			endif
		enddo
	else
		do j = 1,S(1)%nsel
			call rads_put_var (S(1), Pout, S(1)%sel(j), pack(tot(j)%z(:),mask), start)
		enddo
	endif
	start(2) = start(2) + 1
endif

! Reset statistics
call init_stat

600 format (3i0.2)
601 format (i3,i5)
602 format (i3,1x,3i0.2)
end subroutine output_stat

!***********************************************************************
! Initialise statistics

subroutine init_stat
box = stat(0, 0d0, (/ 0d0, 0d0, huge(1d0), -huge(1d0) /))
nr = 0
start_time = nan
end subroutine init_stat

!***********************************************************************
! Write out ASCII the header

subroutine ascii_header
integer :: j0, j, l

600 format ('# Statistics of RADS variables (',a,')'/ &
'# Created: ',a,' UTC: ',a/ &
'#'/'# Satellite : ',a/'# Cycles    :',i5,' -',i5/ &
'# Passes    :',i5,' -',i5/'#'/'# Output columns:')
610 format ('#    ( 1) date [YYMMDD]')
611 format ('# ( 1, 2) cycle and pass')
612 format ('#    ( 1) cycle'/'#    ( 2) date at beginning of cycle [YYMMDD]')
620 format ('#    (',i2,') nr of measurements'/'#    (',i2,') mean time [',a,']')
621 format ('#    (',i2,') mean of ')
622 format ('# (',i2,'-',i2,') mean and stddev of ')
623 format ('# (',i2,'-',i2,') mean, min and max of ')
624 format ('# (',i2,'-',i2,') mean, stddev, min and max of ')

write (*,600) trim(wtype(wmode)), timestamp(), trim(S(1)%command), &
	trim(S(1)%branch(1)), S(1)%cycles(1:2), S(1)%passes(1:2)
select case (output_format)
case (period_day)
	write (*,610)
	j0 = 1
case (period_pass)
	write (*,611)
	j0 = 2
case default
	write (*,612)
	j0 = 2
end select
write (*,620) j0+1,j0+2,trim(S(1)%time%info%units)

format_string = '(i9,f12.0'
do j = 1,S(1)%nsel
	! Write description of variables
	select case (lstat)
	case (0, 1)
		write (*,621,advance='no') (j-1)+j0+3
	case (2)
		write (*,622,advance='no') 2*(j-1)+j0+3,2*j+j0+2
	case (3)
		write (*,623,advance='no') 3*(j-1)+j0+3,3*j+j0+2
	case default
		write (*,624,advance='no') 4*(j-1)+j0+3,4*j+j0+2
	end select
	call rads_long_name_and_units(S(1)%sel(j))
	! Assemble format for statistics lines
	l = len_trim(format_string)
	! Add one decimal to an f-format, or copy boz-format
	if (S(1)%sel(j)%info%boz_format) then
		write (format_string(l+1:),'(",",i0,"(1x,",a,")")') max(lstat,1),trim(S(1)%sel(j)%info%format)
	else
		call read_val (S(1)%sel(j)%info%format(2:), sizes, '.')
		if (sizes(1) > 0) sizes(1) = sizes(1) + 1
		sizes(2) = sizes(2) + 1
		write (format_string(l+1:),'(",",i0,"(1x,f",i0,".",i0,")")') max(lstat,1),sizes
	endif
enddo

l = len_trim(format_string)
format_string(l+1:) = ')'
end subroutine ascii_header

!***********************************************************************
! Write out the NetCDF header

subroutine netcdf_header
integer :: j, k, dimid(2)
integer(onebyteint) :: istat(4) = (/ 1_onebyteint, 2_onebyteint, 3_onebyteint, 4_onebyteint /)

! Open output NetCDF file
call nfs (nf90_create (filename, nf90_write + nf90_netcdf4, ncid))
call nfs (nf90_def_dim (ncid, 'time', nf90_unlimited, dimid(1)))
if (lstat > 0) call nfs (nf90_def_dim (ncid, 'stat', lstat, dimid(2)))

! To use general NetCDF creation machinary, we trick the library a bit here
Pout%fileinfo(1) = rads_file (ncid, filename)
Pout%rw = .true.
grpid(1) = ncid

! Define "time" variables
call nfs (nf90_def_var (ncid, 'nr', nf90_int4, dimid(1:1), varid(1)))
call nfs (nf90_put_att (ncid, varid(1), 'long_name', 'number of measurements'))
if (output_format /= period_day) call rads_def_var (S(1), Pout, 'cycle')
if (output_format == period_pass) call rads_def_var (S(1), Pout, 'pass')
call rads_def_var (S(1), Pout, S(1)%time)

! Create groups when requested

if (groups) then
	call nfs (nf90_def_grp (ncid, 'mean', grpid(1)))
	call nfs (nf90_put_att (grpid(1), nf90_global, 'comment', 'Data variables in this group contain mean values'))
	if (mask(2)) then
		call nfs (nf90_def_grp (ncid, 'stddev', grpid(2)))
		call nfs (nf90_put_att (grpid(2), nf90_global, 'comment', 'Data variables in this group contain standard deviations'))
	endif
	if (mask(3)) then
		call nfs (nf90_def_grp (ncid, 'min', grpid(3)))
		call nfs (nf90_put_att (grpid(3), nf90_global, 'comment', 'Data variables in this group contain minimum values'))
	endif
	if (mask(4)) then
		call nfs (nf90_def_grp (ncid, 'max', grpid(4)))
		call nfs (nf90_put_att (grpid(4), nf90_global, 'comment', 'Data variables in this group contain maximum values'))
	endif

	! Define selected variables
	do k = 1,4
		if (mask(k)) then
			Pout%fileinfo(1) = rads_file (grpid(k), filename) ! Set to right group
			do j = 1,S(1)%nsel
				call rads_def_var (S(1), Pout, S(1)%sel(j), coordinates=.false.)
			enddo
		endif
	enddo
	Pout%fileinfo(1) = rads_file (ncid, filename)

! When not using groups

else
	! Define "stat" variable
	if (lstat > 0) then
		call nfs (nf90_def_var (ncid, 'stat', nf90_int1, dimid(2:2), varid(2)))
		call nfs (nf90_put_att (ncid, varid(2), 'long_name', 'statistics type'))
		call nfs (nf90_put_att (ncid, varid(2), 'flag_values', istat))
		call nfs (nf90_put_att (ncid, varid(2), 'flag_meanings', 'mean standard_deviation minimum maximum'))
	endif
	select case (lstat)
	case (0)
		call nfs (nf90_put_att (ncid, nf90_global, 'comment', 'Data variables contain mean values'))
	case (1)
		call nfs (nf90_put_att (ncid, varid(2), 'comment', 'Data columns are mean values'))
	case (2)
		call nfs (nf90_put_att (ncid, varid(2), 'comment', 'Data columns are mean and standard deviation'))
	case (3)
		call nfs (nf90_put_att (ncid, varid(2), 'comment', 'Data columns are mean, minimum, and maximum'))
	case (4)
		call nfs (nf90_put_att (ncid, varid(2), 'comment', 'Data columns are mean, standard deviation, minimum, and maximum'))
	end select

	! Define selected variables
	do j = 1,S(1)%nsel
		if (lstat > 0) S(1)%sel(j)%info%ndims = 2
		call rads_def_var (S(1), Pout, S(1)%sel(j), coordinates=.false.)
	enddo
endif

! Define global attibutes
call nfs (nf90_put_att (ncid, nf90_global, 'Conventions', 'CF-1.7'))
call nfs (nf90_put_att (ncid, nf90_global, 'title', 'RADS 4 statistics file'))
call nfs (nf90_put_att (ncid, nf90_global, 'institution', 'EUMETSAT / NOAA / TU Delft'))
call nfs (nf90_put_att (ncid, nf90_global, 'references', 'RADS Data Manual, Version ' // trim(rads_version_id)))
call nfs (nf90_put_att (ncid, nf90_global, 'satellite', trim(S(1)%branch(1))))
call nfs (nf90_put_att (ncid, nf90_global, 'cycles', S(1)%cycles(1:2)))
call nfs (nf90_put_att (ncid, nf90_global, 'passes', S(1)%passes(1:2)))
call nfs (nf90_put_att (ncid, nf90_global, 'weights', trim(wtype(wmode))))
call nfs (nf90_put_att (ncid, nf90_global, 'box_size', res))
call nfs (nf90_put_att (ncid, nf90_global, 'history', timestamp()//' UTC: '//trim(S(1)%command)))
call nfs (nf90_enddef (ncid))

! Write "stat" coordinate
if (.not.groups .and. lstat > 0) call nfs (nf90_put_var (ncid, varid(2), istat(:lstat)))

end subroutine netcdf_header

!***********************************************************************

end program radsstat
