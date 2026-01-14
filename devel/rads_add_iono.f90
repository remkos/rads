!-----------------------------------------------------------------------
! Copyright (c) 2011-2026  Remko Scharroo
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

! Command line arguments

character(rads_cmdl) :: path
integer(fourbyteint) :: j, cyc, pass, ios
real(eightbytereal) :: f, f_scaled
integer(fourbyteint), parameter :: nmod = 3
logical :: model(nmod) = .false.

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P
type(rads_var), pointer :: var
type(giminfo) :: info(2)

! Other parameters

real(eightbytereal), parameter :: jplg_start = 283996800d0, codg_start = 418438800d0 ! 1994-01-01T00 and 1998-04-06T01

! Initialise

call synopsis ('--head')
call rads_set_options ('cgjn code jpl gim nic09 all scale:')
call rads_init (S)

! Determine conversion factor from TEC units to ionospheric delay in metres
! and mean altitude.

f = -0.4028d0/S%frequency(1)**2	! freq in GHz from rads.xml; f = factor from TEC to meters
var => rads_varptr (S, 'iono_gim')
read (var%info%parameters, *) f_scaled ! Get the ionosphere scaling parameter

! Get ${ALTIM}/data directory

call parseenv ('${ALTIM}/data/', path)

! Check options

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('g', 'gim', 'j', 'jpl')
		model(1) = .true.
	case ('c', 'code')
		model(2) = .true.
	case ('n', 'nic09')
		model(3) = .true.
	case ('all')
		model(1) = .true.
		model(3) = .true.
	case ('scale')
		read (rads_opt(j)%arg, *, iostat=ios) f_scaled ! Overrule config file
	end select
enddo
f_scaled = f_scaled * f

if (model(1)) info(1) = giminit(trim(path)//'gim/jplg_',max(1,rads_verbose))
if (model(2)) info(2) = giminit(trim(path)//'gim/codg_',max(1,rads_verbose))
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
if (rads_version ('Add ionosphere models to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -g, --gim, -j, --jpl      Add JPL GIM model data (only for 1994-01-01 and later)' / &
'  -n, --nic09               Add NIC09 ionosphere model' / &
'  --all                     All of the above' / &
'  -c, --code                Add CODE GIM model data (only for 1998-01-01 and later)' / &
'  --scale=SCALE             Set scale factor (default is from rads.xml)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, j
real(eightbytereal) :: time(n), lat(n), lon(n), z(n,nmod)
logical :: ok(nmod)

call log_pass (P)

! Get time and position

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)

! Now do all models

z = nan

if (model(1) .and. time(1) > jplg_start) then ! JPL GIM
	do i = 1,n
		z(i,1) = f_scaled * gimtec(time(i),lat(i),lon(i),info(1))
	enddo
endif

if (model(2) .and. time(1) > codg_start) then ! CODE GIM
	do i = 1,n
		z(i,2) = f_scaled * gimtec(time(i),lat(i),lon(i),info(2))
	enddo
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
if (ok(2)) call rads_def_var (S, P, 'iono_gim')	! Temporarily uses the same variable
if (ok(3)) call rads_def_var (S, P, 'iono_nic09')

if (ok(1)) call rads_put_var (S, P, 'iono_gim', z(:,1))
if (ok(2)) call rads_put_var (S, P, 'iono_gim', z(:,2)) ! Temporarily uses the same variable
if (ok(3)) call rads_put_var (S, P, 'iono_nic09', z(:,3))

call log_records (n)
end subroutine process_pass

end program rads_add_iono
