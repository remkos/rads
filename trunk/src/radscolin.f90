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

integer(fourbyteint), parameter :: msat = 5
type(rads_sat) :: S(msat)
integer(fourbyteint) :: nsel = 0, reject = 9999, cycle, pass, i, j, ios, &
	nbins, nsat = 0, ntrx = 0, ntrx1, ntrx2, type_sla = 1, step = 1, ncols
real(eightbytereal) :: dt = 0.97d0
character(len=rads_naml) :: prefix = 'radscolin_p', suffix = '.nc', satlist
logical :: ascii = .true., out_data = .true., out_mean = .false., out_sdev = .false.
real(eightbytereal), allocatable :: data(:,:,:)
logical, allocatable :: mask(:,:)
integer(fourbyteint), allocatable :: nr_in_bin(:), bin(:)
type :: stat_
	integer(fourbyteint) :: nr
	real(eightbytereal) :: mean,sum2
end type
type(stat_), allocatable :: stat(:,:)
type :: info_
	character(len=6) :: sat
	integer(twobyteint) :: satid, cycle
	integer(fourbyteint) :: ndata
end type
type(info_), allocatable :: info(:)

! Initialize RADS or issue help
call synopsis
call rads_set_options ('adsnNo::r:: step: dt: output:')
call rads_init (S)
if (any(S%error /= rads_noerr)) call rads_exit ('Fatal error')

! Determine how many satellites and cycles.
! Also check that the same number of variables have been selected for each satellite
do i = 1,msat
	if (S(i)%sat == '') exit
	if (S(i)%nsel == 0) call rads_parse_varlist (S(i), 'sla')
	if (S(i)%nsel /= S(1)%nsel) call rads_exit ('Unequal amount of variables on sel= for different satellites')
	if (any(S(i)%passes /= S(1)%passes)) call rads_exit ('Unequal number of passes per cycle for different satellites')
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

! Determine which column to check for NaN (default: 1st)
do j = 1,nsel
	if (S(nsat)%sel(j)%info%datatype == rads_type_sla) type_sla = j
enddo

! Scan command line arguments
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('r')
		if (rads_opt(i)%arg /= 'n') then
			reject = 0
			read (rads_opt(i)%arg, *, iostat=ios) reject
		endif
	case ('step')
		read (rads_opt(i)%arg, *, iostat=ios) step
	case ('dt')
		read (rads_opt(i)%arg, *, iostat=ios) dt
	case ('a')
		out_mean = .true.
	case ('s')
		out_sdev = .true.
	case ('d')
		out_data = .false.
	case ('o', 'out')
		ascii = .false.
		if (rads_opt(i)%arg == '') cycle
		j = index (rads_opt(i)%arg, '#')
		if (j == 0) call rads_exit ('Output file name needs to include at least one "#"')
		prefix = rads_opt(i)%arg(:j-1)
		j = index (rads_opt(i)%arg, '#', .true.)
		suffix = rads_opt(i)%arg(j+1:)
	end select
enddo

! Allocate data arrays
nbins = nint((S(1)%phase%pass_seconds + 60d0)/dt/2d0) ! Number of bins on either side of equator
allocate (data(ntrx+2,nsel,-nbins:nbins), mask(ntrx+2,-nbins:nbins), nr_in_bin(-nbins:nbins), &
	bin(-nbins:nbins), stat(ntrx+2,nsel), info(ntrx+2))

forall (i=-nbins:nbins) bin(i) = i

! Read one pass for each satellites at a time
do pass = S(1)%passes(1), S(1)%passes(2), S(1)%passes(3)
	call process_pass
enddo

! End RADS
call rads_end (S)

! Deallocate data arrays
deallocate (data, mask, nr_in_bin, bin, stat, info)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('$Revision$','Make collinear data sets from RADS')) return
call rads_synopsis ()
write (stderr,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r#                       Reject data when there are fewer than # tracks with valid SLA'/ &
'                            (default: # = number of selected cycles)'/ &
'  -r0, -r                   Keep all stacked data points, even NaN'/ &
'  -rn                       Reject data when any track is NaN (default)'/ &
'  --dt=DT                   Set minimum bin size in seconds (default is determined by satellite)'/ &
'  --step=N                  Write out only every N points'/ &
'  -a                        Output mean in addition to pass data'/ &
'  -s                        Output standard deviation in addition to pass data'/ &
'  -d                        Do not output pass data'/ &
'  -o, --out[=OUTNAME]       Create netCDF output by pass. Optionally specify filename including "#", which'/ &
'                            is to be replaced by the psss number. Default is "radscolin_p#.nc"')
stop
end subroutine synopsis

!***********************************************************************
! Process the data for a single pass

subroutine process_pass
real(eightbytereal), allocatable :: temp(:)
integer, allocatable :: bin(:)
integer :: i, j, k, m
type(rads_pass) :: P

! Initialize
data = S(1)%nan
i = 0
ntrx = 0
nr_in_bin = 0
mask = .false.
stat = stat_ (0, 0d0, 0d0)

! Read in data
do m = 1,nsat
	do cycle = S(m)%cycles(1), S(m)%cycles(2), S(m)%cycles(3)
		call rads_open_pass (S(m), P, cycle, pass)
		if (P%ndata > 0) then
			ntrx = ntrx + 1 ! track counter
			info(ntrx) = info_ ('    '//S(m)%sat, S(m)%satid, int(cycle,twobyteint), P%ndata)
			allocate (temp(P%ndata), bin(P%ndata))
			bin = nint((P%tll(:,1) - P%equator_time) / dt) ! Store bin nr associated with measurement
			do j = 1,nsel
				call rads_get_var (S(m), P, S(m)%sel(j), temp)
				data(ntrx,j,bin(:)) = temp(:)
			enddo
			! For time being, set to "true" ANY incoming data point, even if NaN
			mask(ntrx,bin(:)) = .true.
			deallocate (temp,bin)
		endif
		call rads_close_pass (S(m), P)
	enddo
enddo

! Specify the columns for statistics
ntrx1 = ntrx + 1
ntrx2 = ntrx + 2
info(ntrx1) = info_ ('  mean', 0, 9001, 9999)
info(ntrx2) = info_ ('stddev', 0, 9002, 9999)

! If reject == 0, count number of SLA measurements per bin, also the NaNs
! Else, count the number of non-NaN SLA measurements per bin
! In both cases, set mask to non-NaN SLA measurements only
if (reject == 0) then
	forall (k=-nbins:nbins) nr_in_bin(k) = count(mask(1:ntrx,k))
	mask = .not.isnan(data(:,type_sla,:))
else
	mask = .not.isnan(data(:,type_sla,:))
	forall (k=-nbins:nbins) nr_in_bin(k) = count(mask(1:ntrx,k))
	! Set to zero the bins that not reach the threshold number
	where (nr_in_bin < min(ntrx,reject)) nr_in_bin = 0
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
		call mean_variance (pack(data(1:ntrx,j,k),mask(1:ntrx,k)), data(ntrx1,j,k), data(ntrx2,j,k))
	enddo
enddo
data(ntrx2,:,:) = sqrt(data(ntrx2,:,:)) ! Variance to std dev
! Mask out NaN statistics
mask(ntrx1:ntrx2,:) = .not.isnan(data(ntrx1:ntrx2,type_sla,:))

! Compute per-track statistics (vertically)
do i = 1,ntrx2
	do j = 1,nsel
		call mean_variance (pack(data(i,j,:),mask(i,:)), stat(i,j)%mean, stat(i,j)%sum2)
	enddo
	stat(i,:)%nr = count(mask(i,:))
enddo
stat%sum2 = sqrt(stat%sum2) ! Variance to std dev

! Determine column ranges for output
ncols = 0
if (out_data) ncols = ntrx
if (out_mean) ncols = ncols + 1
if (out_sdev) ncols = ncols + 1

! If printing standard deviation but not mean, put std dev in place of mean
if (out_sdev.and..not.out_mean) then
	data (ntrx1,:,:) = data(ntrx2,:,:)
	info(ntrx1) = info(ntrx2)
	stat(ntrx1,:) = stat(ntrx2,:)
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
! Write the pass in netCDF

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

! Open output netCDF file
call nfs (nf90_create (filename, nf90_write, ncid))
call nfs (nf90_def_dim (ncid, 'bin', n, dimid(1)))
call nfs (nf90_def_dim (ncid, 'track', ncols, dimid(2)))

! To use general netCDF creation machinary, we trick the library a bit here
P%ncid = ncid
P%filename = filename
S(1)%time%info%ndims = 2

! Define selected variables
do j = 1,nsel
	S(1)%sel(j)%info%ndims = 2
	call rads_def_var (S(1), P, S(1)%sel(j))
enddo

! Define track info
call nfs (nf90_def_var (ncid, 'satid', nf90_int1, dimid(2:2), varid(1)))
call nfs (nf90_put_att (ncid, varid(1), 'long_name', 'satellite ID'))
call nfs (nf90_put_att (ncid, varid(1), 'flag_values', int(S(1:nsat)%satid, onebyteint)))
call nfs (nf90_put_att (ncid, varid(1), 'flag_meanings', trim(satlist)))
call nfs (nf90_put_att (ncid, varid(1), 'comment', 'Satellite IDs relate to the different missions'))
call nfs (nf90_def_var (ncid, 'cycle', nf90_int2, dimid(2:2), varid(2)))
call nfs (nf90_put_att (ncid, varid(2), 'long_name', 'cycle number'))
call nfs (nf90_put_att (ncid, varid(2), 'comment', 'Cycle number 9001 denotes "mean", 9002 denotes "standard deviation"'))

! Define bin info
call nfs (nf90_def_var (ncid, 'nr', nf90_int2, dimid(1:1), varid(3)))
call nfs (nf90_put_att (ncid, varid(3), 'long_name', 'number of collinear measurements in bin'))
call nfs (nf90_def_var (ncid, 'bin', nf90_int2, dimid(1:1), varid(4)))
call nfs (nf90_put_att (ncid, varid(4), 'long_name', 'bin number'))
call nfs (nf90_put_att (ncid, varid(4), 'comment', 'Bin number is 0 at equator, adding/subtracting 1 for each 1-Hz time step'))

! Define global attibutes
call nfs (nf90_put_att (ncid, nf90_global, 'Conventions', 'CF-1.5'))
call nfs (nf90_put_att (ncid, nf90_global, 'title', 'RADS 4.0 collinear tracks file'))
call nfs (nf90_put_att (ncid, nf90_global, 'institution', 'Altimetrics / NOAA / TU Delft'))
call nfs (nf90_put_att (ncid, nf90_global, 'references', 'RADS Data Manual, Issue 4.0'))
call nfs (nf90_put_att (ncid, nf90_global, 'pass_number', pass))
call nfs (nf90_put_att (ncid, nf90_global, 'history', timestamp()//' UTC: '//trim(S(1)%command)))
call nfs (nf90_enddef (ncid))

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

! Write track info
call nfs (nf90_put_var (ncid, varid(1), info(1:ncols)%satid))
call nfs (nf90_put_var (ncid, varid(2), info(1:ncols)%cycle))

! Write bin info
call nfs (nf90_put_var (ncid, varid(3), pack(nr_in_bin, nr_in_bin > 0)))
call nfs (nf90_put_var (ncid, varid(4), pack(bin, nr_in_bin > 0)))

call nfs (nf90_close (ncid))
end subroutine write_pass_netcdf

!***********************************************************************
! Write the pass in ASCII

subroutine write_pass_ascii
logical :: first = .true.
integer :: i, j, k
character(len=640) :: format_string

600 format('# RADS collinear track file'/'# Created: ',a,' UTC: ',a)
610 format('#'/'# Pass      = ',i4.4/'# Satellite =',999(1x,a6))
615 format('# Cycles    =',999(4x,i3.3))
620 format('#'/'# Column ranges for each variable:')
622 format('# ',i4,' -',i4,' : ',a,' [',a,']')
625 format('# ',i4,7x,': ',a)
630 format('#')
635 format(2(1x,i5))
640 format('# ',a,': ')
645 format(i4,i5)
650 format('# nr : ',999i6)

if (.not.first) write (*,*) ! Skip line between passes
first = .false.

! Describe data set per variable
write (*,600) timestamp(), trim(S(1)%command)
write (*,610) pass, info(1:ncols)%sat
write (*,615) info(1:ncols)%cycle

! Describe variables
write (*,620)
i = 1
do j = 1,nsel
	write (*,622) i,i+ncols-1,trim(S(nsat)%sel(j)%info%long_name),trim(S(nsat)%sel(j)%info%units)
	i = i + ncols - 1
enddo
write (*,625) i,'number of measurements'
i = i + 1
write (*,625) i,'record number'
write (*,630)

! Build format string
if (ncols == 1) then
	write (format_string,'("(",a,",")') trim(S(nsat)%sel(1)%info%format)
else
	write (format_string,'("(",a,",",i3,"(1x,",a,"),")') &
		trim(S(nsat)%sel(1)%info%format),ncols-1,trim(S(nsat)%sel(1)%info%format)
endif
do i = 2,nsel
	write (format_string(len_trim(format_string)+1:),'(i3,"(1x,",a,"),")') &
		ncols,trim(S(nsat)%sel(i)%info%format)
enddo
format_string(len_trim(format_string):) = ')'

! Print out data that are common to some passes
do k = -nbins,nbins,step
	if (nr_in_bin(k) == 0) cycle
	write (*,format_string,advance='no') data(1:ncols,:,k)
	write (*,635,advance='no') nr_in_bin(k), k
	write (*,*)
enddo

! Write per-pass stats
write (*,640,advance='no') 'avg'
write (*,format_string,advance='no') stat(1:ncols,:)%mean
write (*,645) S(1)%cycles(1),pass
write (*,640,advance='no') 'std'
write (*,format_string,advance='no') stat(1:ncols,:)%sum2
write (*,645) S(1)%cycles(1),pass
write (*,650) stat(1:ncols,:)%nr,S(1)%cycles(1),pass

end subroutine write_pass_ascii

!***********************************************************************

end program radscolin
