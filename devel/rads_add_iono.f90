!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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
type(giminfo) :: info

! Command line arguments

character(rads_naml) :: dir
integer(fourbyteint) :: j, cyc, pass
real(eightbytereal) :: f, f_scaled, hgt

! Initialise

call synopsis
call rads_init (S)

! Determine conversion factor from TEC units to ionospheric delay in metres
! and mean altitude.

f = -0.40250d0/S%frequency(1)**2	! freq in GHz from rads.xml; f = factor from TEC to meters
select case (S%sat)
case ('c2')
	f_scaled = 0.844 * f	! CryoSat
	hgt = 750d3
case ('tx', 'pn', 'j1', 'j2', 'j3', 'j4')
	f_scaled = 0.925 * f	! TOPEX, Jason
	hgt = 1350d3
case default
	f_scaled = 0.856 * f	! Other LEOs
	hgt = 800d3
end select

! Get $ALTIM directory

call getenv ('ALTIM', dir)

! Initialise the selected models

do j = 1,S%nsel
	select case (S%sel(j)%name)
	case ('iono_gim')
		info = giminit(trim(dir)//'/data/gim/jplg_',1)
	case ('iono_nic09')
		call nicinit(trim(dir)//'/data/nic09/nic09_clim.nc',trim(dir)//'/data/nic09/nic09_gtec.nc')
	end select
enddo

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata,S%nsel)
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
'  -V, --var=NAME[,...]      Select variable name(s) of ionospere models (required)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n, nmod)
integer(fourbyteint), intent(in) :: n, nmod
integer(fourbyteint) :: i, j, ii, iold
real(eightbytereal) :: time(n), lat(n), lon(n), tec1(n), z(n, nmod), tec2, dtime, d
logical :: allnan(nmod)

! Formats

551  format (a,' ...',$)
552  format (i5,' records changed')

write (*,551) trim(P%filename)

! Get time and position

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)

! Now do all models

do j = 1,nmod
	select case (S%sel(j)%name)
	case ('iono_gim')
		do i = 1,n
			z(i,j) = f_scaled * gimtec(time(i),lat(i),lon(i),info)
		enddo
	case ('iono_iri2007')
		! Compute IRI every 10 seconds, then interpolate linearly in time
		iold = 1
		do i = 1,n
			dtime = time(i) - time(iold)
			if (dtime < 10d0 .and. i > 1 .and. i < n) cycle
			call iri2007tec(0,time(i),lat(i),lon(i),hgt,hgt,tec1(i),tec2)
			do ii = iold+1,i-1
				d = (time(ii)-time(iold)) / dtime
				tec1(ii) = (1d0-d) * tec1(iold) + d * tec1(i)
			enddo
			iold = i
		enddo
		z(:,j) = f * tec1
	case ('iono_nic09')
		do i = 1,n
			z(i,j) = f_scaled * nictec(time(i),lat(i),lon(i))
		enddo
	end select
enddo

! Store all (non-NaN) data fields

do j = 1,nmod
	allnan(j) = all(isnan_(z(:,j)))
enddo
if (all(allnan)) then
	write (*,552) 0
	return
endif
call rads_put_history (S, P)
do j = 1,nmod
	if (.not.allnan(j)) call rads_def_var (S, P, S%sel(j))
enddo
do j = 1,nmod
	if (.not.allnan(j)) call rads_put_var (S, P, S%sel(j), z(:,j))
enddo

write (*,552) n
end subroutine process_pass

end program rads_add_iono
