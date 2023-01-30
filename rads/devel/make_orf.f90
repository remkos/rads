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

!*make_orf -- Create an ORF file from orbits
!
! This program creates an Orbital Revolution File (ORF), similar to
! those released by CNES for the Jason and SARAL missions.
! Output is to stdout.
!
! usage: make_orf [data-selectors] [options]
!-----------------------------------------------------------------------
program make_orf

use rads
use rads_misc
use rads_devel
use rads_time

! Define orbit type

type :: orbit
	real(eightbytereal) :: time, lat, lon
end type

! Data variables

integer(fourbyteint), parameter :: morf = 500000
integer(fourbyteint) :: norf = 0
type(rads_sat) :: S
type(orbit) :: info(-1:1), orf(morf), diff

! Command line arguments

integer(fourbyteint) :: i
real(eightbytereal) :: dt = 1d0
character(len=rads_cmdl) :: dir = ''

! Other variables

integer(fourbyteint) :: ios, cycle, pass, abs_orbit
real(eightbytereal) :: min_lat, max_lat, time, time_step = 0d0
character(len=26) :: date

! Scan command line for options

call synopsis
call rads_set_options (' dir: dt: ext:')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('dt')
		read (rads_opt(i)%arg, *, iostat=ios) dt
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('dir')
		dir = rads_opt(i)%arg
	end select
enddo

! If cycles given, and no time limits, set the time limits based on the cycles

if (isnan_(S%time%info%limits(1))) S%time%info%limits(1) = max(S%phases(1)%start_time, rads_cycle_to_time (S, S%cycles(1)))
if (isnan_(S%time%info%limits(2))) S%time%info%limits(2) = rads_cycle_to_time (S, S%cycles(2)+1)

! If dir is not given, figure out directory from variable selection

if (dir /= '') then
	! skip
else if (S%nsel == 0) then
	call rads_exit ('Specify --dir <directory> or --var <alt_variable>')
else
	dir = S%sel(1)%info%parameters
endif

! Add prefix ${ALTIM}/data/ODR.<satellite>/ to directory name

if (dir(:1) == '/' .or. dir(:2) == './') then
else
	call parseenv ('${ALTIM}/data/ODR.' // trim(S%satellite) // '/' // dir, dir)
endif

! Determine threshold of latitudes

min_lat = dt * 0.1d0	! Keep significant margin with 0.1 degrees per second
max_lat = min(S%inclination, 180d0 - S%inclination) - min_lat

! Start with start time

time = floor(S%time%info%limits(1))

! Run until no orbits left, in search of equator crossings and polar crossovers

do
	time = time + time_step
	if (time > S%time%info%limits(2)) exit
	if (time_step == dt) then	! Normal step by dt
		info(-1:0) = info(0:1)
		info(1)%time = time + dt
	else	! Upon initialisation (time_step = 0) or after large jump (time_step > dt)
		info(-1)%time = time - dt
		call get_orbit (-1)
		info( 0)%time = time
		call get_orbit ( 0)
		time_step = dt
	endif
	info( 1)%time = time + dt
	call get_orbit ( 1)
	if (ios > 0) exit
	if (abs(info(0)%lat) < min_lat) then
		call find_equator
	else if (abs(info(0)%lat) > max_lat) then
		call find_pole
	endif
enddo

! Extend the list until the upper time limit
! For computing steps into future, start with an equator crossing; it is more accurate

if (abs(orf(norf)%lat) > min_lat) norf = norf - 1
if (norf >= 5) then
	diff%time = orf(norf)%time - orf(norf-4)%time
	diff%lon = orf(norf)%lon - orf(norf-4)%lon
	do
		if (orf(norf-3)%time + diff%time > S%time%info%limits(2)) exit
		norf = norf + 1
		orf(norf)%time = orf(norf-4)%time + diff%time
		orf(norf)%lon = orf(norf-4)%lon + diff%lon
		orf(norf)%lat = orf(norf-4)%lat
	enddo
endif

! Print file header

write (*,600)
600 format (199('#')/ &
'# TABLE OF ORBIT REFERENCE PARAM: DATEUTC, NCYCLE, NPAS, NREV, LONG, LAT'/ &
'# DATEUTC: Date of the event in UTC; (NCYCLE, NPAS): Cycle and pass number at time event'/ &
'# (LONG, LAT): Geodetic longitude and latitude of the event.  '/ &
199('#'))

! Print out all the equator crossings and polar crossovers

do i = 1,norf
	! To allow for some sloppiness in the orbit scenario specified in the rads.xml file, add 1/8
	! of an orbital revolution time to the time when retreiving the cycle and pass number
	call rads_time_to_cycle_pass (S, orf(i)%time + 0.25d0 * S%phase%pass_seconds, cycle, pass, abs_orbit)
	date = strf1985f (orf(i)%time)
	date(5:5) = '/'
	date(8:8) = '/'
	! Subtract 1 from absolute orbit number for the start of a pass (large negative latitude)
	if (orf(i)%lat < -min_lat) abs_orbit = abs_orbit - 1
	write (*,610) date(1:23),rads_tab,cycle,rads_tab,pass,rads_tab,abs_orbit,rads_tab, &
		modulo(orf(i)%lon,360d0),rads_tab,orf(i)%lat
enddo
610 format(a,a1,i3.3,a1,i5.5,a1,i5.5,a1,f6.2,a1,f6.2)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Create ORF file', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:'/ &
'  --dir STR                 Specify orbit directory path'/ &
'  --dt DT                   Specify time stepping interval (s), default: 1')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Get orbit information
!-----------------------------------------------------------------------

subroutine get_orbit (i)
integer(fourbyteint), intent(in) :: i
integer(fourbyteint) :: getorb
real(eightbytereal) :: alt
ios = getorb (info(i)%time, info(i)%lat, info(i)%lon, alt, dir, rads_verbose > 0)
end subroutine get_orbit

!-----------------------------------------------------------------------
! Find equator crossing
!-----------------------------------------------------------------------

subroutine find_equator
real(eightbytereal) :: time, lat, lon, alt
integer(fourbyteint) :: getorb
if (sign(1d0,info(0)%lat) == sign(1d0,info(1)%lat)) return
! Find the root with two-step secant method
time = info(0)%time - info(1)%lat / (info(1)%lat - info(0)%lat) * dt
ios = getorb (time, lat, lon, alt, dir, rads_verbose > 0)
time = time - lat / (lat - info(0)%lat) * (time - info(0)%time)
call store_orf (time)
end subroutine find_equator

!-----------------------------------------------------------------------
! Find rollover point near pole
!-----------------------------------------------------------------------

subroutine find_pole
real(eightbytereal) :: a2, b, time
b = (info(1)%lat - info(-1)%lat) / 2d0
a2 = info(1)%lat + info(-1)%lat - 2d0 * info(0)%lat
time = - b / a2
if (time < -0.5d0 .or. time > 0.5d0) return
time = info(0)%time + time * dt
call store_orf (time)
end subroutine find_pole

!-----------------------------------------------------------------------
! Store ORF record
!-----------------------------------------------------------------------

subroutine store_orf (time)
use rads_time
real(eightbytereal), intent(in) :: time
real(eightbytereal) :: alt
integer(fourbyteint) :: getorb
norf = norf + 1
if (norf > morf) call rads_exit ('Too many ORF records')
orf(norf)%time = time
ios = getorb (time, orf(norf)%lat, orf(norf)%lon, alt, dir, rads_verbose > 0)
time_step = 1200d0	! Jump 20 minutes (quarter revolution ahead)
end subroutine store_orf

end program make_orf
