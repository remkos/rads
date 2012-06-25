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

program radsstat

! This program reads the RADS data base and computes statistics
! by pass, cycle or day of a number of RADS data variables.
!
! usage: radsstat [RADS_options] [options]
!-----------------------------------------------------------------------
use rads
use rads_time
use rads_misc
integer(fourbyteint) :: lstat=1
character(len=32) :: wtype(0:3)=(/ &
	'box weight                  ', 'constant weight             ', &
	'area weighted               ', 'inclination-dependent weight'/)
character(len=640) :: format_string
integer(fourbyteint), parameter :: period_day=0, period_pass=1, period_cycle=2, day_init=-999999
integer(fourbyteint) :: nr=0, day=day_init, day_old=day_init, cycle, pass, i, l, &
	period=period_day, wmode=0, nx, ny, kx, ky, ios, sizes(2)
real(eightbytereal), allocatable :: z(:,:), lat_w(:)
real(eightbytereal) :: sini, step=1d0, x0, y0, res(2)=(/3d0,1d0/)
type :: stat
	real(eightbytereal) :: wgt, mean, sum2, xmin, xmax
end type
type(stat), allocatable :: box(:,:,:), tot(:)
type(rads_sat) :: S
type(rads_pass) :: P

! Initialize RADS or issue help
call synopsis
call rads_set_options ('c::d::p::b::masl res:')
call rads_init (S)
if (S%error /= rads_noerr) call rads_exit ('Fatal error')

! If no sel= is given, use sla
if (S%nsel == 0)  call rads_parse_varlist (S, 'sla')

! Scan command line arguments
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('d')
		period = period_day
		read (rads_opt(i)%arg, *, iostat=ios) step
	case ('p')
		period = period_pass
		read (rads_opt(i)%arg, *, iostat=ios) step
		step = dble(S%passes(2))/nint(S%passes(2)/step)
	case ('c')
		period = period_cycle
		read (rads_opt(i)%arg, *, iostat=ios) step
	case ('b')
		wmode = 0
		call read_val (rads_opt(i)%arg, res, '/-+x')
	case ('m')
		wmode = 1
	case ('a')
		wmode = 2
	case ('s')
		wmode = 3
	case ('l')
		lstat = 2
	case ('res')
		call read_val (rads_opt(i)%arg, res, '/-+x')
	end select
enddo

! Set up the boxes
x0 = S%lon%info%limits(1)
y0 = S%lat%info%limits(1)
nx = int((S%lon%info%limits(2)-x0)/res(1)+0.999d0)
ny = int((S%lat%info%limits(2)-y0)/res(2)+0.999d0)
allocate (box(0:S%nsel,nx,ny),tot(0:S%nsel),lat_w(ny))

! Set up the weights
sini = sin(S%inclination*rad)
if (wmode == 1) then
	lat_w = 1d0
else if (wmode == 3) then
	forall (ky=1:ny) lat_w(ky) = sqrt(max(1d-2,1d0-(sin((y0+(ky-0.5d0)*res(2))*rad)/sini)**2))
else
	forall (ky=1:ny) lat_w(ky) = cos((y0+(ky-0.5d0)*res(2))*rad)
endif

! Initialize statistics
box = stat(0d0, 0d0, 0d0, S%nan, S%nan)
call write_header

! Determine format for statistics
format_string = '(i9,f12.0'
do i = 1,S%nsel
	l = len_trim(format_string)
	! Add one decimal to the format
	call read_val (S%sel(i)%info%format(2:), sizes, '.')
	write (format_string(l+1:),'(",",i0,"(1x,f",i0,".",i0,")")') 2*lstat,sizes+1
enddo
l = len_trim(format_string)
format_string(l+1:) = ')'

! Start looping through cycles and passes
do cycle = S%cycles(1), S%cycles(2), S%cycles(3)
	! Process passes one-by-one
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cycle, pass)
		if (P%ndata > 0) call process_pass

		! Print the statistics at the end of the data pass (if requested)
		if (period == period_pass .and. nint(modulo(dble(pass),step)) == 0) call print_stat
		call rads_close_pass (S, P)
	enddo

	! Print the statistics at the end of the cycle (if requested)
	if (period == period_cycle .and. nint(modulo(dble(cycle),step)) == 0) call print_stat
enddo

! Flush the statistics and close RADS4
cycle = S%cycles(2)
pass = S%passes(2)
call print_stat
if (S%debug >= 1) call rads_stat (S)
call rads_end (S)
deallocate (box,tot,lat_w)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('$Revision$','Print RADS statistics per cycle, pass or day(s)')) return
call rads_synopsis ()
write (stderr,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -c[N]                     Statistics per cycle or N cycles'/ &
'  -d[N]                     Statistics per day (default) or N days'/ &
'  -p[N]                     Statistics per pass or N passes'/ &
'  -b[DX,DY]                 Average by boxes with size (default on: 3x1 degrees)'/ &
'  -m                        Give all measurements equal weight'/ &
'  -a                        Weight measurements by cosine of latitude'/ &
'  -s                        Use inclination-dependent weight'/ &
'  -l                        Print min and max in addition to mean and stddev'/ &
'  --res=DX,DY               Size of averaging boxes (default = 3x1 degrees)')
stop
end subroutine synopsis

!***********************************************************************
! Process data for a single pass

subroutine process_pass
integer :: i, j

! Read the data for this pass
allocate (z(P%ndata,0:S%nsel))
z(:,0) = P%tll(:,1)	! Store time
do j = 1,S%nsel
	call rads_get_var (S, P, S%sel(j), z(:,j))
enddo

! Initialise day_old if not done before
if (day_old == day_init) day_old = floor(P%tll(1,1)/86400d0)

! Update the statistics with data in this pass
do i = 1,P%ndata
	! Print the statistics (if in "daily" mode)
	day = floor(P%tll(i,1)/86400d0)
	if (period == period_day .and. day >= nint(day_old+step)) call print_stat
   	if (any(isnan(z(i,:)))) cycle ! Reject outliers

	! Update the box statistics
   	kx = floor((P%tll(i,3)-x0)/res(1) + 1d0)
   	ky = floor((P%tll(i,2)-y0)/res(2) + 1d0)
   	kx = max(1,min(kx,nx))
   	ky = max(1,min(ky,ny))
   	box(:,kx,ky)%wgt  = box(:,kx,ky)%wgt  + 1d0
   	box(:,kx,ky)%mean = box(:,kx,ky)%mean + z(i,:)
   	box(:,kx,ky)%sum2 = box(:,kx,ky)%sum2 + z(i,:)*z(i,:)
   	box(:,kx,ky)%xmin = min(box(:,kx,ky)%xmin, z(i,:))
   	box(:,kx,ky)%xmax = max(box(:,kx,ky)%xmax, z(i,:))
   	nr = nr + 1
enddo

deallocate (z)
end subroutine process_pass

!***********************************************************************
! Print statistics for one batch of data

subroutine print_stat
integer(fourbyteint) :: j,yy,mm,dd
real(eightbytereal) :: w

if (nr == 0) then
	day_old = day
	return
endif

! Init global stats
tot = stat(0d0, 0d0, 0d0, S%nan, S%nan)

! Cycle trough all boxes and determine overall weighted mean
do ky=1,ny
	do kx=1,nx
		if (box(1,kx,ky)%wgt == 0d0) cycle ! Skip empty boxes

		! Determine the weight
		if (wmode == 0) then
	    	w = lat_w(ky) / box(1,kx,ky)%wgt
		else
	    	w = lat_w(ky)
		endif

		! Update overall statistics
		tot(:)%wgt  = tot(:)%wgt  + w * box(:,kx,ky)%wgt
		tot(:)%mean = tot(:)%mean + w * box(:,kx,ky)%mean
		tot(:)%sum2 = tot(:)%sum2 + w * box(:,kx,ky)%sum2
		tot(:)%xmin = min(tot(:)%xmin, box(:,kx,ky)%xmin)
		tot(:)%xmax = max(tot(:)%xmax, box(:,kx,ky)%xmax)
	enddo
enddo
	
! Divide by total weight to compute overall weighted mean and standard deviation

tot%mean = tot%mean / tot%wgt
if (nr > 1) then
	tot%sum2 = sqrt((tot%sum2 / tot%wgt - tot%mean*tot%mean) * nr / (nr - 1d0))
else
	tot%sum2 = S%nan
endif

! Print results
call mjd2ymd(day_old+46066,yy,mm,dd)
select case (period)
case (period_day)
	write (*,600,advance='no') modulo(yy,100),mm,dd
case (period_pass)
	write (*,601,advance='no') cycle,pass
case default
	write (*,602,advance='no') cycle,modulo(yy,100),mm,dd
endselect
if (lstat == 1) then
	write (*,format_string) nr,tot(0)%mean,(tot(j)%mean,tot(j)%sum2,j=1,S%nsel)
else
	write (*,format_string) nr,tot(0)%mean,(tot(j)%mean,tot(j)%sum2,tot(j)%xmin,tot(j)%xmax,j=1,S%nsel)
endif

! Reset statistics
box = stat(0d0, 0d0, 0d0, S%nan, S%nan)
nr  = 0
day_old = day

600 format (3i2.2)
601 format (i3,i5)
602 format (i3,1x,3i2.2)
end subroutine print_stat

!***********************************************************************
! Write out the header

subroutine write_header
integer :: j0, j

600 format ('# Statistics of RADS variables (',a,')'/ &
'# Created: ',a,' UTC: ',a/ &
'#'/'# Satellite : ',a,'/',a/'# Cycles    :',i5,' -',i5/ &
'# Passes    :',i5,' -',i5/'#'/'# Output columns:')
610 format ('#    ( 1) date [YYMMDD]')
611 format ('# ( 1, 2) cycle and pass')
612 format ('#    ( 1) cycle'/'#    ( 2) date at beginning of cycle [YYMMDD]')
620 format ('#    (',i2,') nr of measurements'/'#    (',i2,') mean time [',a,']')
621 format ('# (',i2,'-',i2,') mean and stddev of ',a,' [',a,']')
622 format ('# (',i2,'-',i2,') mean, stddev, min and max of ',a,' [',a,']')

write (*,600) trim(wtype(wmode)),timestamp(),trim(S%command), &
	trim(S%sat),trim(S%phase%name),S%cycles(1:2),S%passes(1:2)
select case (period)
case (period_day)
	write (*,610)
	j0 = 1
case (period_pass)
	write (*,611)
	j0 = 2
case default
	write (*,612)
	j0 = 2
endselect
write (*,620) j0+1,j0+2,trim(S%time%info%units)
do j = 1,S%nsel
	if (lstat == 1) then
		write (*,621) 2*j+j0+1,2*j+j0+2,trim(S%sel(j)%info%long_name),trim(S%sel(j)%info%units)
	else
		write (*,622) 4*j+j0+1,4*j+j0+2,trim(S%sel(j)%info%long_name),trim(S%sel(j)%info%units)
	endif
enddo
end subroutine write_header

!***********************************************************************

end program radsstat
