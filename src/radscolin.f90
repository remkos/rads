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

!*radscolin -- Make collinear data sets from RADS
!+
program radscolin
!
! This program provides a quick and crude way to make collinear data
! sets (in ASCII).
!
! Usage: radscolin [RADS_options] [options]
!-----------------------------------------------------------------------
use rads
use rads_misc
use rads_time

integer(fourbyteint), parameter :: msat = 20
type(rads_sat) :: S(msat)
integer(fourbyteint) :: nsel = 0, reject = -2, cycle, pass, i, j, ios, &
	nbins, nsat = 0, ntrx = 0, ntrx1, ntrx2, ntrx3, ntrx4, type_sla = 1, out_diff = 0, step = 1, ncols
real(eightbytereal) :: dt = 0d0
character(len=rads_naml) :: prefix = 'radscolin_p', suffix = '.nc', satlist
logical :: ascii = .true., out_data = .true., out_mean = .false., out_sdev = .false., out_extr = .false., &
	out_track = .true., out_cumul = .false., keep = .false., force = .false., boz_format = .false.
real(eightbytereal), allocatable :: data(:,:,:)
logical, allocatable :: mask(:,:)
integer(fourbyteint), allocatable :: nr_in_bin(:), idx(:)
type :: stat_
	integer(fourbyteint) :: nr
	real(eightbytereal) :: mean, var
end type
type(stat_), allocatable :: stat(:,:), cumul_stat(:,:)
type :: info_
	character(len=6) :: sat
	integer(twobyteint) :: satid, cycle
	integer(fourbyteint) :: ndata
end type
type(info_), allocatable :: info(:), info_backup(:)
character(len=rads_strl) :: format_string

! Initialize RADS or issue help
call synopsis
call rads_set_options ('acdefklnNo::r::st ' // &
	'cumul diff diff-no-coord diff1 dt: reject-on-nan:: extremes force keep mean minmax no-pass no-track ' // &
	'output:: stddev step:')
call rads_init (S)
if (any(S%error /= rads_noerr)) call rads_exit ('Fatal error')

! Determine how many satellites and cycles.
! Also check that the same number of variables have been selected for each satellite
do i = 1,msat
	if (S(i)%sat == '') exit
	if (S(i)%nsel == 0) call rads_parse_varlist (S(i), 'sla')
	if (S(i)%nsel /= S(1)%nsel) call rads_exit ('Unequal amount of variables requested for different missions')
	if (S(i)%cycles(3) /= S(1)%cycles(3)) call rads_exit ('Cycle step size should be the same for all missions')
	ntrx = ntrx + (S(i)%cycles(2) - S(i)%cycles(1)) / S(i)%cycles(3) + 1
	dt = max(dt,S(i)%dt1hz)
	nsat = i
enddo
nsel = S(1)%nsel

! String of sat names
satlist = S(1)%sat
do i = 2,nsat
	satlist = trim(satlist) // ' ' // S(i)%sat
enddo

! Determine which column has the SLA (if any)
! Also check for boz_format
do j = 1,nsel
	if (S(nsat)%sel(j)%info%datatype == rads_type_sla) type_sla = j
	if (S(nsat)%sel(j)%info%boz_format) boz_format = .true.
enddo

! Scan command line arguments
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('r', 'reject-on-nan')
		call parse_r_option (S(1), rads_opt(i)%opt, rads_opt(i)%arg)
	case ('step')
		read (rads_opt(i)%arg, *, iostat=ios) step
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('dt')
		read (rads_opt(i)%arg, *, iostat=ios) dt
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('a', 'mean')
		out_mean = .true.
	case ('s', 'stddev')
		out_sdev = .true.
	case ('e', 'extremes')
		call rads_message ('warning: option -e|--extremes has been replaced by -l|--minmax')
		out_extr = .true.
	case ('l', 'minmax')
		out_extr = .true.
	case ('diff')
		out_diff = 1
		keep = .true.
	case ('diff-no-coord')
		out_diff = 2
		keep = .true.
	case ('diff1')
		out_diff = 3
		keep = .true.
	case ('d', 'no-pass')
		out_data = .false.
	case ('t', 'no-track')
		out_track = .false.
	case ('c', 'cumul')
		out_cumul = .true.
		keep = .true.
	case ('k', 'keep')
		keep = .true.
	case ('f', 'force')
		force = .true.
	case ('o', 'output')
		ascii = .false.
		if (rads_opt(i)%arg == '') cycle
		j = index (rads_opt(i)%arg, '#')
		if (j == 0) call rads_exit ('Output file name needs to include at least one "#"')
		prefix = rads_opt(i)%arg(:j-1)
		j = index (rads_opt(i)%arg, '#', .true.)
		suffix = rads_opt(i)%arg(j+1:)
	end select
enddo

! --cumul requires --diff or --r0 and not allowed together with --output
if (.not.out_cumul) then
	! Continue
else if (.not.ascii) then
	call rads_exit ('--cumul cannot be used together with -o or --output')
else if  (out_diff == 0 .and. reject /= 0) then
	call rads_exit ('--cumul requires --diff or -r')
endif

! --diff1 not allowed together with --output, -a, -s, -l
if (out_diff /= 3) then
	! Continue
else if (.not.ascii) then
	call rads_exit ('--diff1 cannot be used together with -o or --output')
else if (out_mean .or. out_sdev .or. out_extr) then
	call rads_exit ('--diff1 cannot be used together with statistics (-a, -s, or -l)')
endif

! Allocate data arrays
nbins = nint(S(1)%phase%pass_seconds/dt * 0.6d0) ! Number of bins on either side of equator (20% margin)
allocate (data(ntrx+4,nsel,-nbins:nbins), mask(ntrx+4,-nbins:nbins), nr_in_bin(-nbins:nbins), &
	idx(-nbins:nbins), stat(ntrx+4,nsel), cumul_stat(ntrx+4,nsel), info(ntrx+4), info_backup(ntrx+4))

forall (i=-nbins:nbins)	idx(i) = i
cumul_stat = stat_ (0, 0d0, 0d0)

if (out_diff == 3) then
	! Read one pass for collinear cycles across missions
	do i = 0, S(1)%cycles(2)-S(1)%cycles(1), S(1)%cycles(3)
		do pass = S(1)%passes(1), S(1)%passes(2), S(1)%passes(3)
			call process_pass (S(:)%cycles(1)+i, S(:)%cycles(1)+i, 1)
		enddo
	enddo
else
	! Read one pass for all cycles and missions at a time
	do pass = S(1)%passes(1), S(1)%passes(2), S(1)%passes(3)
		call process_pass (S(:)%cycles(1), S(:)%cycles(2), S(1)%cycles(3))
	enddo
endif

! Write out cumulative stats, if requested
if (out_cumul) then
	cumul_stat%mean = cumul_stat%mean / cumul_stat%nr
	cumul_stat%var = sqrt((cumul_stat%var - cumul_stat%nr * cumul_stat%mean**2) / (cumul_stat%nr - 1))
	call write_pass_stat (cumul_stat(1:ncols,:), &
		S(1)%cycles(2)-S(1)%cycles(1)+1, S(1)%passes(2)-S(1)%passes(1)+1, ' # cumul_')
endif

! End RADS
call rads_end (S)

! Deallocate data arrays
deallocate (data, mask, nr_in_bin, idx, stat, info, info_backup)

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
case ('all')
	type_sla = 0
case default
	if (arg(1:1) >= '0' .and. arg(1:1) <= '9') then
		reject = 0
		read (arg, *, iostat=i) reject
	else
		do i = 1,S%nsel
			if (arg == S%sel(i)%name .or. arg == S%sel(i)%info%name) then
				type_sla = i
				return
			endif
		enddo
		call rads_exit ('Option -'//trim(opt)//' <varname> does not refer to variable specified on -V option')
	endif
end select
end subroutine parse_r_option

!***********************************************************************

subroutine synopsis
if (rads_version ('Make collinear data sets from RADS')) return
call rads_synopsis
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  --dt DT                   Set minimum bin size in seconds (default is determined by satellite)'/ &
'  --step N                  Write out only one out of N bins along track'/ &
'  -r, --reject-on-nan VAR   Base rejection criteria (below) on VAR'/ &
'  -r all                    Base rejection criteria (below) on all selected variables'/ &
'                      Note: By default rejection of data is based on ''sla'' (if selected) or 1st on -V'/ &
'                      Note: The two options above can be combined with any of the three below'/ &
'  -r NR                     Reject stacked data when fewer than NR tracks with valid values'/ &
'  -r 0, -r none, -r         Keep all stacked data points, even NaN'/ &
'  -r n, -r any              Reject stacked data when data on any track is NaN (default)'/ &
'                      Note: If no -r option is given, -r any is assumed'/ &
'  -k, --keep                Keep all passes, even those that do not have data in the selected area'/ &
'  -a, --mean                Output mean in addition to pass data'/ &
'  -s, --stddev              Output standard deviation in addition to pass data'/ &
'  -l, --minmax              Output minimum and maximum in addition to pass data'/ &
'  -d, --no-pass             Do not output pass data'/ &
'  -t, --no-track            Do not print along-track data (ascii output only)'/ &
'  -c, --cumul               Output cumulative statistics (ascii output only) (implies --keep)'/ &
'      --diff                Compute difference between first and second half of selected passes'/ &
'                            (implies --keep)'/ &
'      --diff-no-coord       Same as --diff, but excluding coordinates from computing difference'/ &
'      --diff1               Same as --diff-no-coord, writing one cycle at a time, chronologically'/ &
'                            (cannot be combined with -o, -a, -s, or -l)' / &
'  -f, --force               Force comparison, even when missions are not considered collinear'/ &
'  -o, --output [FILENAME]   Create NetCDF output by pass (default is ascii output to stdout).'/ &
'                            Optionally specify FILENAME including "#", to be replaced by the pass'/ &
'                            number. Default is "radscolin_p#.nc"')
stop
end subroutine synopsis

!***********************************************************************
! Process the data for a single pass

subroutine process_pass (cycle_start, cycle_stop, cycle_step)
integer(fourbyteint) :: cycle_start(:), cycle_stop(:), cycle_step
real(eightbytereal), allocatable :: temp(:)
integer, allocatable :: bin(:)
integer :: i, j, j0, j1, k, m, ndata
type(rads_pass) :: P

! Initialize
data = nan
ntrx = 0
nr_in_bin = 0
mask = .false.
stat = stat_ (0, 0d0, 0d0)
ndata = 0

! Read in data
do m = 1,nsat
	do cycle = cycle_start(m), cycle_stop(m), cycle_step
		call rads_open_pass (S(m), P, cycle, pass)
		if (force) then
			! Skip the next check, do collinear anyway
		else if (S(m)%phase%passes /= S(1)%phase%passes .or. S(m)%phase%ref_lon /= S(1)%phase%ref_lon) then
			! Pass ranges should be the same for all satellites, otherwise we do not have collinear tracks
			call rads_exit ('Satellite missions '//S(m)%sat//'/'//trim(S(m)%phase%name)// &
				' and '//S(1)%sat//'/'//trim(S(1)%phase%name)//' are not collinear')
		endif
		if (P%ndata > 0 .or. keep) then
			ntrx = ntrx + 1 ! track counter
			ndata = ndata + P%ndata ! data counter
			info(ntrx) = info_ ('    '//S(m)%sat, S(m)%satid, int(cycle,twobyteint), P%ndata)
		endif
		if (P%ndata > 0) then
			allocate (temp(P%ndata),bin(P%ndata))
			bin = nint((P%tll(:,1) - P%equator_time) / dt) ! Store bin nr associated with measurement
			! Guard against rogue timings
			if (minval(bin) < -nbins .or. maxval(bin) > nbins) then
				call rads_message ('Skipping pass because of rogue timings (time(:) or equator_time):', P)
			else
				do j = 1,nsel
					call rads_get_var (S(m), P, S(m)%sel(j), temp)
					data(ntrx,j,bin(:)) = temp(:)
				enddo
				! For time being, set to "true" ANY incoming data point, even if NaN
				mask(ntrx,bin(:)) = .true.
			endif
			deallocate (temp,bin)
		endif
		call rads_close_pass (S(m), P)
	enddo
enddo

! Skip this pass if no data is found
if (ndata == 0) return

! If requested, do difference of first half and second half of tracks
if (out_diff > 0) then
	ntrx = ntrx / 2
	do i = 1,ntrx
		do j = 1,nsel
			if (out_diff < 2 .or. S(nsat)%sel(j)%info%datatype < rads_type_time) then
				! Coordinates skipped if requested
				data(i,j,:) = data(i,j,:) - data(i+ntrx,j,:)
				if (S(nsat)%sel(j)%info%datatype == rads_type_lon) then
					! For longitudes, bring differences to within ±180º
					where (data(i,j,:) < -180d0) data(i,j,:) = data(i,j,:) + 360d0
					where (data(i,j,:) >  180d0) data(i,j,:) = data(i,j,:) - 360d0
				endif
			endif
		enddo
		mask(i,:) = mask(i,:) .or. mask(i+ntrx,:)
	enddo
	do i = 2*ntrx,ntrx+1,-1
		info(i+2) = info(i)
	enddo
endif

! Skip if there are no tracks selected
if (ntrx == 0) return

! Specify the columns for statistics
info_backup = info
ntrx1 = ntrx + 1
ntrx2 = ntrx + 2
ntrx3 = ntrx + 3
ntrx4 = ntrx + 4
info(ntrx1) = info_ ('  mean', 0, 9001, 9999)
info(ntrx2) = info_ ('stddev', 0, 9002, 9999)
info(ntrx3) = info_ ('   min', 0, 9003, 9999)
info(ntrx4) = info_ ('   max', 0, 9004, 9999)

! Determine which range of selected variables to scan for NaNs
if (type_sla == 0) then ! This signals: check all variables
	j0 = 1; j1 = nsel
else
	j0 = type_sla; j1 = type_sla
endif

! If reject == 0, count number of SLA measurements per bin, also the NaNs
! Else, count the number of non-NaN SLA measurements per bin
! In both cases, set mask to non-NaN SLA measurements only
if (reject == 0) then ! This signals: ignore NaNs
	forall (k=-nbins:nbins)
		nr_in_bin(k) = count(mask(1:ntrx,k))
		forall (i=1:ntrx) mask(i,k) = all(isan_(data(i,j0:j1,k)))
	end forall
else
	forall (k=-nbins:nbins)
		forall (i=1:ntrx) mask(i,k) = all(isan_(data(i,j0:j1,k)))
		nr_in_bin(k) = count(mask(1:ntrx,k))
	end forall
	! Set to zero the bins that not reach the threshold number
	k = reject
	if (reject == -2) k = ntrx
	where (nr_in_bin < k) nr_in_bin = 0
endif

! If no valid measurements at all, return
if (sum(nr_in_bin) == 0) return

! Mask out bins with zero measurements
do k = -nbins,nbins
	if (nr_in_bin(k) == 0) mask(1:ntrx,k) = .false.
enddo

! Compute per-bin statistics (horizonally)
do k = -nbins,nbins
	do j = 1,nsel
		if (S(nsat)%sel(j)%info%datatype == rads_type_lon) then
			! For longitudes, use different routine to create proper mean
			call mean_variance (pack(data(1:ntrx,j,k),mask(1:ntrx,k)), data(ntrx1,j,k), data(ntrx2,j,k), &
				data(ntrx3,j,k), data(ntrx4,j,k), S(nsat)%sel(j)%info%limits)
		else
			! For all other variables
			call mean_variance (pack(data(1:ntrx,j,k),mask(1:ntrx,k)), data(ntrx1,j,k), data(ntrx2,j,k), &
				data(ntrx3,j,k), data(ntrx4,j,k))
		endif
	enddo
enddo
data(ntrx2,:,:) = sqrt(data(ntrx2,:,:)) ! Variance to std dev
! Mask out NaN statistics
forall (i=ntrx1:ntrx4,k=-nbins:nbins) mask(i,k) = all(isan_(data(i,j0:j1,k)))

! Compute per-track statistics (vertically)
do i = 1,ntrx4
	do j = 1,nsel
		call mean_variance (pack(data(i,j,:),mask(i,:)), stat(i,j)%mean, stat(i,j)%var)
	enddo
	stat(i,:)%nr = count(mask(i,:))
enddo
stat%var = sqrt(stat%var) ! Variance to std dev

! Do cumulative statistics, if requested
if (out_cumul) then
	where (stat%nr >= 2)
		cumul_stat%mean = cumul_stat%mean + stat%nr * stat%mean
		cumul_stat%var = cumul_stat%var + stat%nr * stat%mean**2 + (stat%nr - 1) * stat%var**2
		cumul_stat%nr = cumul_stat%nr + stat%nr
	endwhere
endif

! Determine column ranges for output
ncols = 0
if (out_data) ncols = ntrx
if (out_mean) ncols = ncols + 1
if (out_sdev) ncols = ncols + 1
if (out_extr) ncols = ncols + 2

! If not printing standard dev, move the rest of the stats forward
if (.not.out_sdev) then
	data(ntrx2:ntrx3,:,:) = data(ntrx3:ntrx4,:,:)
	info(ntrx2:ntrx3) = info(ntrx3:ntrx4)
	stat(ntrx2:ntrx3,:) = stat(ntrx3:ntrx4,:)
endif

! If not printing mean, move the rest of the stats forward
if (.not.out_mean) then
	data(ntrx1:ntrx3,:,:) = data(ntrx2:ntrx4,:,:)
	info(ntrx1:ntrx3) = info(ntrx2:ntrx4)
	stat(ntrx1:ntrx3,:) = stat(ntrx2:ntrx4,:)
endif

! If not printing data, move statistics forward
if (.not.out_data) then
	data (1:ncols,:,:) = data(ntrx1:ntrx+ncols,:,:)
	info(1:ncols) = info(ntrx1:ntrx+ncols)
	stat(1:ncols,:) = stat(ntrx1:ntrx+ncols,:)
endif

! Print out the pass
if (ascii) then
	call write_pass_ascii
else
	call write_pass_netcdf
endif

end subroutine process_pass

!***********************************************************************
! Write the pass in NetCDF

subroutine write_pass_netcdf
use netcdf
use rads_netcdf
character(len=rads_naml) :: filename
integer(fourbyteint) :: ncid, dimid(2), j, k, n, start(2) = 1, varid(4)
type(rads_pass) :: P
real(eightbytereal), allocatable :: tmp(:,:)

! Count number of bins
n = count (nr_in_bin(-nbins:nbins:step) > 0)

! Construct filename
write (filename,'(a,i4.4,a)') trim(prefix),pass,trim(suffix)

! Open output NetCDF file
call nfs (nf90_create (filename, nf90_write, ncid))
call nfs (nf90_def_dim (ncid, 'bin', n, dimid(1)))
call nfs (nf90_def_dim (ncid, 'track', ncols, dimid(2)))

! Define track info
call nfs (nf90_def_var (ncid, 'satid', nf90_int1, dimid(2:2), varid(1)))
call nfs (nf90_put_att (ncid, varid(1), 'long_name', 'satellite ID'))
call nfs (nf90_put_att (ncid, varid(1), 'flag_values', int(S(1:nsat)%satid, onebyteint)))
call nfs (nf90_put_att (ncid, varid(1), 'flag_meanings', trim(satlist)))
call nfs (nf90_put_att (ncid, varid(1), 'comment', 'Satellite IDs relate to the different missions'))
call nfs (nf90_def_var (ncid, 'cycle', nf90_int2, dimid(2:2), varid(2)))
call nfs (nf90_put_att (ncid, varid(2), 'long_name', 'cycle number'))
call nfs (nf90_put_att (ncid, varid(2), &
	'comment', 'Cycle number; 9001 = mean; 9002 = standard deviation; 9003 = minimum; 9004 = maximum'))

! Define bin info
call nfs (nf90_def_var (ncid, 'bin', nf90_int2, dimid(1:1), varid(3)))
call nfs (nf90_put_att (ncid, varid(3), 'long_name', 'bin number'))
call nfs (nf90_put_att (ncid, varid(3), 'comment', 'Bin number is 0 at equator, adding/subtracting 1 for each 1-Hz time step'))
call nfs (nf90_def_var (ncid, 'nr', nf90_int2, dimid(1:1), varid(4)))
call nfs (nf90_put_att (ncid, varid(4), 'long_name', 'number of collinear measurements in bin'))

! Define global attibutes
call nfs (nf90_put_att (ncid, nf90_global, 'Conventions', 'CF-1.5'))
call nfs (nf90_put_att (ncid, nf90_global, 'title', 'RADS 4.0 collinear tracks file'))
call nfs (nf90_put_att (ncid, nf90_global, 'institution', 'EUMETSAT / NOAA / TU Delft'))
call nfs (nf90_put_att (ncid, nf90_global, 'references', 'RADS Data Manual, Version ' // trim(rads_version_id)))
call nfs (nf90_put_att (ncid, nf90_global, 'pass_number', pass))
call nfs (nf90_put_att (ncid, nf90_global, 'history', timestamp()//' UTC: '//trim(S(1)%command)))

! To use general NetCDF creation machinary, we trick the library a bit here
P%fileinfo(1) = rads_file (ncid, filename)
P%rw = .true.
S(1)%time%info%ndims = 2

! Define selected variables
do j = 1,nsel
	call rads_def_var (S(1), P, S(1)%sel(j), ndims = 2)
enddo

call nfs (nf90_enddef (ncid))

! Write track info
call nfs (nf90_put_var (ncid, varid(1), info(1:ncols)%satid))
call nfs (nf90_put_var (ncid, varid(2), info(1:ncols)%cycle))

! Write bin info
call nfs (nf90_put_var (ncid, varid(3), pack(idx(-nbins:nbins:step), nr_in_bin(-nbins:nbins:step) > 0)))
call nfs (nf90_put_var (ncid, varid(4), pack(nr_in_bin(-nbins:nbins:step), nr_in_bin(-nbins:nbins:step) > 0)))

! Write data
allocate (tmp(ncols,n))

do j = 1,nsel
	i = 0
	do k = -nbins,nbins,step
		if (nr_in_bin(k) == 0) cycle
		i = i + 1
		tmp(:,i) = data(1:ncols,j,k)
	enddo
	call rads_put_var (S(1), P, S(1)%sel(j), tmp, start)
enddo

deallocate (tmp)

call nfs (nf90_close (ncid))

end subroutine write_pass_netcdf

!***********************************************************************
! Write the pass in ASCII

subroutine write_pass_ascii
logical :: first = .true.
integer :: i, j, k

600 format ('# RADS collinear track file'/'# Created: ',a,' UTC: ',a/'#')
601 format ('# Pass      = ',i4.4)
610 format ('# Satellite =',999(1x,a6))
615 format ('# Cycles    =',999(4x,i3.3))
618 format ('# ... minus ...')
620 format ('#'/'# Column ranges for each variable:')
622 format ('# ',i4,' -',i4,' : ')
625 format ('# ',i4,7x,': ',a)
630 format ('#')

if (first .or. out_track) then
	if (.not.first) write (*,*) ! Skip line between passes

	! Describe data set per variable
	write (*,600) timestamp(), trim(S(1)%command)
	if (out_track) write (*,601) pass
	if (out_data .or. out_mean .or. out_sdev) then
		write (*,610) info(1:ncols)%sat
		write (*,615) info(1:ncols)%cycle
	else
		write (*,610) info(1:ntrx)%sat
		write (*,615) info(1:ntrx)%cycle
	endif
	if (out_diff > 0 .and. out_data) then
		write (*,618)
		write (*,610) info_backup(ntrx+3:2*ntrx+2)%sat
		write (*,615) info_backup(ntrx+3:2*ntrx+2)%cycle
	endif

	! Describe variables
	write (*,620)
	i = 1
	if (ncols > 0) then
		do j = 1,nsel
			write (*,622,advance='no') i,i+ncols-1
			call rads_long_name_and_units(S(nsat)%sel(j))
			i = i + ncols
		enddo
	endif
	write (*,625) i,'number of measurements'
	i = i + 1
	write (*,625) i,'record number'
	write (*,630)
endif

! Build format string
if (ncols == 0) then
	format_string = '(a1,'
else if (ncols == 1) then
	write (format_string,'("(a1,",a,",")') trim(S(nsat)%sel(1)%info%format)
else
	write (format_string,'("(a1,",a,",",i0,"(1x,",a,"),")') &
		trim(S(nsat)%sel(1)%info%format),ncols-1,trim(S(nsat)%sel(1)%info%format)
endif
do i = 2,nsel
	write (format_string(len_trim(format_string)+1:),'(i0,"(1x,",a,"),")') &
		ncols,trim(S(nsat)%sel(i)%info%format)
enddo
format_string(len_trim(format_string)+1:) = '2i9,a)'

! Do a transfer of bit patterns if needed
if (boz_format) then
	do i = 1,nsel
		if (.not.S(nsat)%sel(i)%info%boz_format) cycle
		call bit_transfer (data(:,i,:))
		call bit_transfer (stat(:,i)%mean)
		call bit_transfer (stat(:,i)%var)
	enddo
endif

! Print out data that are common to some passes
if (out_track) then
	do k = -nbins,nbins,step
		if (nr_in_bin(k) == 0) cycle
		write (*,format_string) ' ', data(1:ncols,:,k), nr_in_bin(k), k
	enddo
endif

! Write per-pass stats
call write_pass_stat (stat(1:ncols,:), S(1)%cycles(1), pass, ' # ')

first = .false.

end subroutine write_pass_ascii

!***********************************************************************
! Write per-pass stats

subroutine write_pass_stat (stat, cyc, pass, string)
type(stat_) :: stat(:,:)
integer(fourbyteint) :: cyc, pass
character(len=*) :: string
650 format('#',i8,999i9)
write (*,format_string) '#', stat%mean, cyc, pass, string//'avg'
write (*,format_string) '#', stat%var, cyc, pass, string//'std'
write (*,650,advance='no') stat%nr, cyc, pass
write (*,'(a)') string//'nr'
end subroutine write_pass_stat

!***********************************************************************

end program radscolin
