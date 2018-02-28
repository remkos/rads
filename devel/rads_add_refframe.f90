!-----------------------------------------------------------------------
! Copyright (c) 2011-2018  Remko Scharroo
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
! C2: In order to keep on properly patching baseline B data, there is
!     a special provision made to update reference frame offsets that
!     were previously < -500 mm (C00 term only). If so, then we use a value
!     673 mm less than the configured value.
! J1: The reference frame offset for Phase A and B is given in the
!     configuration. For Phase C, 5 mm is added to the C00 term.
! J2: The reference frame offset for the SLA from MLE4 measurements is
!     given in the configuration. For MLE3 measurements 28.5 mm is added
!     to the C00 term.
! 3A: The reference frame offset is for the latest product baseline.
!     For data before PB2.19 (SM-2 6.10) increase the reference frame offset
!     by 23 mm. See notes of 2018-01-26.
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
real(eightbytereal) :: coef(5) = 0d0, a(5)
logical :: constant
type(rads_var), pointer :: var

! Other variables

integer(fourbyteint) :: i, ios

! Initialise

call synopsis ('--head')
call rads_set_options (' coef:')
call rads_init (S)

! Get default coefficients from rads.xml file

var => rads_varptr (S, 'ref_frame_offset')
read (var%info%parameters, *, iostat=ios) coef

! Check for options

do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('coef')
		coef = 0d0
		read (rads_opt(i)%arg, *, iostat=ios) coef
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	end select
enddo

coef = coef*1d-3	! Convert mm to m
constant = all(coef(2:5) == 0d0)

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
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
' --coef C00,C10,C11,S11,C20 Set reference frame offset coefficients (mm)' / &
'                            (default values come from rads.xml)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: lon(n), lat(n), cor(n), ref_frame_offset(n)
real(eightbytereal), parameter :: rad=atan(1d0)/45d0
integer(fourbyteint) :: i

call log_pass (P)

! Get lat, lon

call rads_get_var (S, P, 'lon', lon, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)

! Compute correction

if (constant) then
	cor = coef(1)
else
	do i = 1,n
		call get_spharm(lat(i)*rad,lon(i)*rad,a)
		cor(i) = dot_product(coef,a)
	enddo
endif

! If CryoSat-2 and we had values less than -500 mm, then use C00 - 673 mm.

if (S%sat == 'c2' .and. cyc < 63) then
	call rads_get_var (S, P, var, ref_frame_offset, .true.)
	if (ref_frame_offset(1) < -500d-3) cor = cor - 673d-3
endif

! If Jason-1 phase C, add 5 mm

if (S%sat == 'j1' .and. S%phase%name == 'c') cor = cor + 5d-3

! If Sentinel-3 prior to SM-2 version 6.10, add 23 mm

if (S%sat == '3a') then
	i = index(P%original, 'IPF-SM-2')
	if (i > 0 .and. P%original(i+9:i+13) < '06.10') cor = cor + 25d-3
endif

! Store all data fields.

call rads_put_history (S, P)
call rads_def_var (S, P, var)
if (S%sat == 'j2') call rads_def_var (S, P, 'ref_frame_offset_mle3')
call rads_put_var (S, P, var, cor)
if (S%sat == 'j2') call rads_put_var (S, P, 'ref_frame_offset_mle3', cor+28.5d-3)

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
