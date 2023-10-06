!-----------------------------------------------------------------------
! Copyright (c) 2011-2023  Remko Scharroo
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

!*rads_add_refframe -- Add reference frame correction to RADS data
!+
! This program adds a new field to the RADS data, which is meant
! to model a constant reference frame offset between the satellite
! and the others. A total of 5 spherical harmonic coefficients can
! be used to model the offset: C00, C10, C11, S11 and C20.
!
! In practice we will use a zero offset for TOPEX.
!
! Special provisions are made for the following missions:
! J1: The reference frame offset for Phase A and B is given in the
!     configuration. For Phase C, 5 mm is added to the C00 term.
!
! usage: rads_add_refframe [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_refframe

use rads
use rads_grid
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P
type(grid) :: info

! Command line arguments

integer(fourbyteint) :: cyc, pass
real(eightbytereal) :: a(5)

! Other variables

integer(fourbyteint) :: i, k, ios
type :: var_
	real(eightbytereal) :: coef(5)
	logical :: constant
end type
type(var_), allocatable :: var(:)

! Initialise

call synopsis ('--head')
call rads_set_options ('x:: ext:: coef:')
call rads_init (S)

! Check for -x or --ext options

do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('x', 'ext')
		if (rads_opt(i)%arg == '') then
			call rads_parse_varlist (S, 'ref_frame_offset')
		else
			call rads_parse_varlist (S, 'ref_frame_offset_' // rads_opt(i)%arg)
		endif
	end select
enddo
! Default to adding 'ref_frame_offset' only
if (S%nsel == 0) call rads_parse_varlist (S, 'ref_frame_offset')
allocate(var(S%nsel))

! Get default coefficients from rads.xml file

do k = 1,S%nsel
	var(k)%coef = 0d0
	read (S%sel(k)%info%parameters, *, iostat=ios) var(k)%coef
enddo

! Check for --coef option

do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('coef')
		var(1)%coef = 0d0
		read (rads_opt(i)%arg, *, iostat=ios) var(1)%coef
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
		do k = 2,S%nsel
			var(k)%coef = var(1)%coef
		enddo
	end select
enddo

! Convert coefficients from mm to m

do k = 1,S%nsel
	var(k)%coef = var(k)%coef*1d-3
	var(k)%constant = all(var(k)%coef(2:5) == 0d0)
enddo

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata, S%nsel)
		call rads_close_pass (S, P)
	enddo
enddo

! Free the allocated grid

call grid_free(info)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Add reference frame correction to RADS data', flag=flag)) return
call synopsis_devel ('')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -V, --var=VAR1[,VAR2,...] Add specified variables (e.g. ref_frame_offset,ref_frame_offset_mle3)'/ &
'  -x, --ext EXT             Produce field ssha_EXT (e.g. "-x -x mle3" for "-Vref_frame_offset,ref_frame_offset_mle3")' / &
' --coef C00,C10,C11,S11,C20 Set reference frame offset coefficients (mm); applies to all variables' / &
'                            (default values come from rads.xml)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n, m)
integer(fourbyteint), intent(in) :: n, m
real(eightbytereal) :: lon(n), lat(n), cor(n,m)
real(eightbytereal), parameter :: rad=atan(1d0)/45d0
integer(fourbyteint) :: i, k

call log_pass (P)

! Get lat, lon

if (all(var(:)%constant)) then
	do k = 1,m
		cor(:,k) = var(k)%coef(1)
	enddo
else
	call rads_get_var (S, P, 'lon', lon, .true.)
	call rads_get_var (S, P, 'lat', lat, .true.)
	do i = 1,n
		call get_spharm(lat(i)*rad,lon(i)*rad,a)
		do k = 1,m
			cor(i,k) = dot_product(var(k)%coef,a)
		enddo
	enddo
endif

! If Jason-1 phase C, add 5 mm

if (S%sat == 'j1' .and. S%phase%name == 'c') cor = cor + 5d-3

! Store all data fields.

call rads_put_history (S, P)
do k = 1,m
	call rads_def_var (S, P, S%sel(k))
enddo
do k = 1,m
	call rads_put_var (S, P, S%sel(k), cor(:,k))
enddo

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Determine first five spherical harmonics
!-----------------------------------------------------------------------

subroutine get_spharm (lat, lon, a)
use typesizes
real(eightbytereal), intent(in) :: lat, lon
real(eightbytereal), intent(out) :: a(5)
real(eightbytereal) :: plm(0:2)

a(1) = 1d0					! C00
call p_lm (1, lat, plm)
a(2) = plm(0)				! C10
a(3) = plm(1) * cos(lon)	! C11
a(4) = plm(1) * sin(lon)	! S11
call p_lm (2, lat, plm)
a(5) = plm(0)				! C20
end subroutine get_spharm

end program rads_add_refframe
