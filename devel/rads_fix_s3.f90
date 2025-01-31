!-----------------------------------------------------------------------
! Copyright (c) 2011-2025  Remko Scharroo
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

!*rads_fix_s3 -- Patch RADS altimeter files of Sentinel-3 for various anomalies
!
! This program makes numerous patches to the Sentinel-3 RADS data processed
! by rads_gen_s3. These patches include:
!
!  --uso                     Apply USO correction (S3B only)
!  --all                     All the above
!
! usage: rads_fix_s3 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_s3

use rads
use rads_misc
use rads_devel
use meteo_subs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! To store USO table

integer(fourbyteint), parameter :: m_uso = 10000
integer(fourbyteint) :: n_uso = 0, i_uso = 1, t_uso(m_uso)
real(eightbytereal) :: f_uso(m_uso)

! Other local variables

integer(fourbyteint) :: i, cyc, pass, ios, yyyy, mm, dd
real(eightbytereal) :: uso_scale
logical :: luso = .false.
character(len=320) :: uso_filenm

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' uso all')
call rads_init (S)
if (S%sat /= "3b") stop

do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('uso')
		luso = .true.
	case ('all')
		luso = .true.
	end select
enddo

if (.not.luso) stop

! Load the USO table

call parseenv ('${ALTIM}/data/tables/S3B.cor_uso_freq.daily.csv', uso_filenm)
call log_string ('(' // trim(uso_filenm) // ')')
i = getlun()
open (unit=i, file=uso_filenm, status='old', iostat=ios)
if (ios /= 0) call rads_exit ('Error reading USO file ' // trim(uso_filenm))
read (i,600)
do
	read (i,600,iostat=ios) yyyy, mm, dd, uso_scale
	if (ios /= 0) exit
	n_uso = n_uso + 1
	call ymd2mjd (yyyy, mm, dd, t_uso(n_uso))
	f_uso(n_uso) = uso_scale
enddo
close (i)
t_uso(:n_uso) = t_uso(:n_uso) - 46066	! Convert to days since 1985

600 format (i4,1x,i2,1x,i2,1x,f18.16)

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	do pass = S%passes(1),S%passes(2),S%passes(3)
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
if (rads_version ('Patch Sentinel-3B data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --uso                     Apply USO correction to range' / &
'  --all                     All the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
use rads_time
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: range_ku(n), range_ku_plrm(n), range_c(n), uso_scale

call log_pass (P)

! Do this routine only for Baseline < 005

if (P%original(38:40) >= '005') then
	call log_records(0)
	return
endif

! Get the USO scale factor

call get_uso (P%equator_time / 86400d0, uso_scale)

! Update ranges

call rads_get_var (S, P, 'range_ku', range_ku, .true.)
call rads_get_var (S, P, 'range_ku_plrm', range_ku_plrm, .true.)
call rads_get_var (S, P, 'range_c', range_c, .true.)

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Write out all the data

call rads_put_var (S, P, 'range_ku', range_ku * uso_scale)
call rads_put_var (S, P, 'range_ku_plrm', range_ku_plrm * uso_scale)
call rads_put_var (S, P, 'range_c', range_c * uso_scale)

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Find and interpolate USO scale factor at given date
!-----------------------------------------------------------------------

subroutine get_uso (t, uso_scale)
real(eightbytereal), intent(in) :: t
real(eightbytereal), intent(out) :: uso_scale
real(eightbytereal) :: x
do while (i_uso < n_uso - 1 .and. t_uso(i_uso) < t)
	i_uso = i_uso + 1
enddo
x = t - t_uso(i_uso)
uso_scale = f_uso(i_uso) * (1-x) + f_uso(i_uso+1) * x
end subroutine get_uso

end program rads_fix_s3
