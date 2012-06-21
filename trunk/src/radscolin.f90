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
	nbins, nsat = 0, ntrx = 0, ntrx1, ntrx2, ptrx0, ptrx1, type_sla = 1, step = 1
real(eightbytereal) :: dt = 0.97d0
character(len=rads_naml) :: arg, opt, optarg
character(len=640) :: format_string
logical :: numbered = .false., counter = .false.
real(eightbytereal), allocatable :: data(:,:,:)
logical, allocatable :: mask(:,:)
integer(fourbyteint), allocatable :: nr_in_bin(:)
type :: stat_
	integer(fourbyteint) :: nr
	real(eightbytereal) :: mean,sum2
end type
type(stat_), allocatable :: stat(:,:)

! Initialize RADS or issue help
call synopsis
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

! Determine which column to check for NaN (default: 1st)
do j = 1,nsel
	if (S(nsat)%sel(j)%info%datatype == rads_type_sla) type_sla = j
enddo

! Set default column ranges
ntrx1 = ntrx + 1
ntrx2 = ntrx + 2
ptrx0 = 1
ptrx1 = ntrx
reject = ntrx

! Scan command line arguments
do i = 1,iargc()
	call getarg (i, arg)
	call splitarg (arg, opt, optarg)
	select case (opt)
	case ('-r')
		if (optarg == 'n') then
			reject = ntrx
		else
			reject = 0
			read (optarg, *, iostat=ios) reject
		endif
	case ('--step')
		read (optarg, *, iostat=ios) step
	case ('--dt')
		read (optarg, *, iostat=ios) dt
	case ('-s') ! for backward compatibility only
		ptrx1 = ntrx2
	case ('-a')
		ptrx1 = ntrx1
		if (optarg == 's') ptrx1 = ntrx2
	case ('-A')
		ptrx0 = ntrx1
		ptrx1 = ntrx1
		if (optarg == 's') ptrx1 = ntrx2
	case ('-n')
		numbered = .true.
	case ('-N')
		counter = .true.
	end select
enddo
reject = max(0,min(reject,ntrx))

! Build format string
if (ptrx1 > ptrx0) then
	write (format_string,'("(",a,",",i3,"(1x,",a,"),")') &
		trim(S(nsat)%sel(1)%info%format),(ptrx1-ptrx0),trim(S(nsat)%sel(1)%info%format)
else
	write (format_string,'("(",a,",")') trim(S(nsat)%sel(1)%info%format)
endif
do i = 2,nsel
	write (format_string(len_trim(format_string)+1:),'(i3,"(1x,",a,"),")') &
		(ptrx1-ptrx0+1),trim(S(nsat)%sel(i)%info%format)
enddo
format_string(len_trim(format_string):) = ')'

! Allocate data arrays
nbins = nint((S(1)%phase%pass_seconds + 60d0)/dt/2d0) ! Number of bins on either side of equator
allocate (data(ntrx2,nsel,-nbins:nbins), mask(ntrx2,-nbins:nbins), nr_in_bin(-nbins:nbins), stat(ptrx0:ptrx1,nsel))

! Read one pass for each satellites at a time
do pass = S(1)%passes(1), S(1)%passes(2), S(1)%passes(3)
	call process_pass
enddo

! Deallocate data arrays
deallocate (data, mask, nr_in_bin, stat)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('$Revision$','Make collinear data sets from RADS')) return
call rads_synopsis ()
write (stderr,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r#                 reject data when there are fewer than # tracks with valid SLA'/ &
'                      (default: # = number of selected cycles)'/ &
'  -r0, -r             keep all stacked data points, even NaN'/ &
'  -rn                 reject data when any track is NaN (default)'/ &
'  --dt=DT             set minimum bin size in seconds (default is determined by satellite)'/ &
'  --step=N            write out only every N points'/ &
'  -a                  print mean in addition to pass data'/ &
'  -as                 print mean and standard deviation in addition to pass data'/ &
'  -A                  print only mean (no pass data)'/ &
'  -As                 print only mean and standard deviation (no pass data)'/ &
'  -n                  add record number' / &
'  -N                  add number of measurements in each bin')
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
nr_in_bin = 0
mask = .false.
stat = stat_ (0, 0d0, 0d0)

! Read in data
do m = 1,nsat
	do cycle = S(m)%cycles(1), S(m)%cycles(2), S(m)%cycles(3)
		i = i + 1 ! track counter
		call rads_open_pass (S(m), P, cycle, pass)
		if (P%ndata > 0) then
			allocate (temp(P%ndata), bin(P%ndata))
			bin = nint((P%tll(:,1) - P%equator_time) / dt) ! Store bin nr associated with measurement
			do j = 1,nsel
				call rads_get_var (S(m), P, S(m)%sel(j), temp)
				data(i,j,bin(:)) = temp(:)
			enddo
			! For time being, set to "true" ANY incoming data point, even if NaN
			mask(i,bin(:)) = .true.
			deallocate (temp,bin)
		endif
		call rads_close_pass (S(m), P)
	enddo
enddo

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
	where (nr_in_bin < reject) nr_in_bin = 0
endif

! If no valid measurements at all, return
if (sum(nr_in_bin) == 0) return

! Mask out bins with zero measurements
do k = -nbins,nbins
	if (nr_in_bin(k) == 0) mask(:,k) = .false.
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

! Print out the pass
call write_pass_ascii
! call write_pass_nc

call rads_end (S)

end subroutine process_pass

!***********************************************************************
! Write the pass in netCDF

subroutine write_pass_nc
use netcdf
use rads_netcdf
character(len=rads_naml) :: filename = 'radscolin.nc'
integer(fourbyteint) :: ncid, dimid(2), j, k, n, start(2) = 1
type(rads_pass) :: P
real(eightbytereal), allocatable :: tmp(:,:)

! Count number of bins
n = count (nr_in_bin > 0)

! Open output netCDF file
call nfs (nf90_create (filename, nf90_write, ncid))
call nfs (nf90_def_dim (ncid, 'time', n, dimid(1)))
call nfs (nf90_def_dim (ncid, 'track', ntrx, dimid(2)))

! To use general netCDF creation machinary, we trick the library a bit here
P%ncid = ncid
P%filename = filename
S(1)%time%info%ndims = 2

! Define selected variables
do j = 1,nsel
	S(1)%sel(j)%info%ndims = 2
	call rads_def_var (S(1), P, S(1)%sel(j))
enddo

! Define other stuff
call nfs (nf90_put_att (ncid, nf90_global, 'Conventions', 'CF-1.5'))
call nfs (nf90_put_att (ncid, nf90_global, 'title', 'RADS 4.0 colinear tracks file'))
call nfs (nf90_put_att (ncid, nf90_global, 'institution', 'Altimetrics / NOAA / TU Delft'))
call nfs (nf90_put_att (ncid, nf90_global, 'references', 'RADS Data Manual, Issue 4.0'))
call nfs (nf90_put_att (ncid, nf90_global, 'history', timestamp()//' UTC: '//trim(S(1)%command)))
call nfs (nf90_enddef (ncid))

allocate (tmp(ntrx,n))

do j = 1,nsel
	i = 0
	do k = -nbins,nbins
		if (nr_in_bin(k) == 0) cycle
		i = i + 1
		tmp(:,i) = data(1:ntrx,j,k)
	enddo
	call rads_put_var (S(1), P, S(1)%sel(j), tmp, start)
enddo

deallocate (tmp)

call nfs (nf90_close (ncid))
end subroutine write_pass_nc

!***********************************************************************
! Write the pass in ASCII

subroutine write_pass_ascii
logical :: first = .true.
integer :: i, j, k

600 format('# RADS collinear track file'/'# Created: ',a,' UTC: ',a/'#'/'# Satellite data selections:')
610 format('#   sat=',a,1x,'cycle=',i3.3,'-',i3.3,' pass=',i4.4)
611 format(' (',i3,'-',i3,')')
615 format('#   ',a,' (',i3,')')
620 format('#'/'# Column ranges for each variable:')
622 format('# ',i4,' -',i4,' : ',a,' [',a,']')
625 format('# ',i4,7x,': ',a)
630 format('#')
635 format(1x,i5)
640 format('# ',a,': ')
645 format(i4,i5)
650 format('# nr : ',400i6)

if (.not.first) write (*,*) ! Skip line between passes
first = .false.

! Describe data set per variable
write (*,600) timestamp(), trim(S(1)%command)
i = 1
do j = 1,nsat
	if (ptrx0 == 1) then !
		write (*,610,advance='no') S(j)%sat,S(j)%cycles(1:2),pass
		write (*,611) i,i + (S(j)%cycles(2) - S(j)%cycles(1)) / S(j)%cycles(3)
		i = i + S(j)%cycles(2) - S(j)%cycles(1) + 1
	else
		write (*,610) S(j)%sat,S(j)%cycles(1:2),pass
	endif
enddo
if (ptrx1 > ntrx) then
	write (*,615) 'average of all cycles',i
	i = i + 1
endif
if (ptrx1 > ntrx1) write (*,615) 'standard deviation of all cycles',i

! Describe variables
write (*,620)
i = 1
do j = 1,nsel
	write (*,622) i,i+ptrx1-ptrx0,trim(S(nsat)%sel(j)%info%long_name),trim(S(nsat)%sel(j)%info%units)
	i = i + ptrx1 - ptrx0 + 1
enddo
if (counter) then
	write (*,625) i,'number of measurements'
	i = i + 1
endif
if (numbered) write (*,625) i,'record number'
write (*,630)

! Print out data that are common to some passes
do k = -nbins,nbins,step
	if (nr_in_bin(k) == 0) cycle
	write (*,format_string,advance='no') data(ptrx0:ptrx1,:,k)
	if (counter) write (*,635,advance='no') nr_in_bin(k)
	if (numbered) write (*,635,advance='no') k
	write (*,*)
enddo

! Write per-pass stats
write (*,640,advance='no') 'avg'
write (*,format_string,advance='no') stat%mean
write (*,645) S(1)%cycles(1),pass
write (*,640,advance='no') 'std'
write (*,format_string,advance='no') stat%sum2
write (*,645) S(1)%cycles(1),pass
write (*,650) stat%nr,S(1)%cycles(1),pass

end subroutine write_pass_ascii

!***********************************************************************

end program radscolin
