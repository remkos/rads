!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2015  Remko Scharroo
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

!*rads_add_iono -- Add ionosphere models to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! by adding one or more ionosphere models.
!
! usage: rads_add_iono [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_iono

use rads
use rads_misc
use rads_devel
use gimsubs
use nicsubs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P
type(rads_var), pointer :: var
type(giminfo) :: info

! Command line arguments

character(rads_cmdl) :: path
integer(fourbyteint) :: j, cyc, pass
real(eightbytereal) :: f, f_scaled
integer(fourbyteint), parameter :: nmod = 3
logical :: model(nmod) = .false.

! Initialise

call synopsis ('--head')
call rads_set_options ('gin gim iri2007 nic09 all')
call rads_init (S)

! Determine conversion factor from TEC units to ionospheric delay in metres
! and mean altitude.

f = -0.4028d0/S%frequency(1)**2	! freq in GHz from rads.xml; f = factor from TEC to meters
var => rads_varptr (S, 'iono_gim')
read (var%info%parameters, *) f_scaled ! Get the ionosphere scaling parameter
f_scaled = f_scaled * f

! Get ${ALTIM}/data directory

call parseenv ('${ALTIM}/data/', path)

! Check options

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('g', 'gim')
		model(1) = .true.
	case ('i', 'iri2007')
		model(2) = .true.
	case ('n', 'nic09')
		model(3) = .true.
	case ('all')
		model = .true.
	end select
enddo

if (model(1)) info = giminit(trim(path)//'gim/jplg_',1)
if (model(3)) call nicinit(trim(path)//'nic09/nic09_clim.nc',trim(path)//'nic09/nic09_gtec.nc')

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Add ionosphere models to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -g, --gim                 Add JPL GIM model data' / &
'  -i  --iri2007             Add IRI2007 ionosphere model' / &
'  -n, --nic09               Add NIC09 ionosphere model' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, j, ii, iold
real(eightbytereal) :: time(n), lat(n), lon(n), alt(n), tec1(n), z(n,nmod), tec2, dtime, d
logical :: ok(nmod)

call log_pass (P)

! Get time and position

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)
call rads_get_var (S, P, 'alt', alt, .true.)

! Now do all models

if (model(1)) then ! JPL GIM
	do i = 1,n
		z(i,1) = f_scaled * gimtec(time(i),lat(i),lon(i),info)
	enddo
endif

if (model(2)) then ! IRI2007
	! Compute IRI every 10 seconds, then interpolate linearly in time
	iold = 1
	do i = 1,n
		dtime = time(i) - time(iold)
		if (dtime < 10d0 .and. i > 1 .and. i < n) cycle
		call iri2007tec(0,time(i),lat(i),lon(i),alt(i),alt(i),tec1(i),tec2)
		do ii = iold+1,i-1
			d = (time(ii)-time(iold)) / dtime
			tec1(ii) = (1d0-d) * tec1(iold) + d * tec1(i)
		enddo
		iold = i
	enddo
	z(:,2) = f * tec1
endif

if (model(3)) then ! NIC09
	do i = 1,n
		z(i,3) = f_scaled * nictec(time(i),lat(i),lon(i))
	enddo
endif

! Store all (non-NaN) data fields

do j = 1,nmod
	ok(j) = model(j) .and. .not.all(isnan_(z(:,j)))
enddo
if (.not.any(ok)) then
	call log_records (0)
	return
endif

call rads_put_history (S, P)

if (ok(1)) call rads_def_var (S, P, 'iono_gim')
if (ok(2)) call rads_def_var (S, P, 'iono_iri2007')
if (ok(3)) call rads_def_var (S, P, 'iono_nic09')

if (ok(1)) call rads_put_var (S, P, 'iono_gim', z(:,1))
if (ok(2)) call rads_put_var (S, P, 'iono_iri2007', z(:,2))
if (ok(3)) call rads_put_var (S, P, 'iono_nic09', z(:,3))

call log_records (n)
end subroutine process_pass

end program rads_add_iono
