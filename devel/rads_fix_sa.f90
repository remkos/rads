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

!*rads_fix_sa -- Patch RADS altimeter files of SARAL for various anomalies
!
! This program makes numerous patches to the SARAL RADS data processed
! by rads_gen_saral. These patches include:
!
! sig0:
! - Add attenuation correction from Patch2 version of NN algorithm
! ssb:
! - Add hybrid SSB to the RADS data
! wet:
! - Shift MWR wet prior to 2013-10-22: subtract 6.4 mm
! wind:
! - Update the wind speed
!
! usage: rads_fix_sa [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_sa

use rads
use rads_misc
use rads_grid
use rads_devel
use meteo_subs
use netcdf
use rads_netcdf

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

character(len=rads_cmdl) :: path
integer(fourbyteint) :: i, cyc, pass, ncid, varid
integer(fourbyteint), allocatable :: rad_pass(:), rad_dsig0(:)
logical, allocatable :: mask(:)
logical :: lsig0 = .false., usig0, lssb = .false., lwet = .false., lwind = .false.
type(grid) :: issb_hyb
character(len=9) :: str = '123456789'

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' sig0 ssb wet wind all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('sig0')
		lsig0 = .true.
	case ('ssb')
		lssb = .true.
	case ('wet')
		lwet = .true.
	case ('wind')
		lwind = .true.
	case ('all')
		lssb = .true.
		lwet = .true.
		lwind = .true.
	end select
enddo

! Load SSB model

if (lssb) then
	call parseenv ('${ALTIM}/data/models/sa_ssb_hyb.nc?ssb_hyb', path)
	if (grid_load(path,issb_hyb) /= 0) call rads_exit ('Error loading '//trim(path))
endif

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	! Load radiometer patch file (if requested and if available)
	if (lsig0) then
		call parseenv ('${RADSROOT}/ext/sa/mwr/AL_RadiometerL2_CLS_Patch2draft_c'//str(cyc:cyc)//'.nc', path)
		usig0 = nff(nf90_open(path, nf90_nowrite, ncid))
	else
		usig0 = .false.
	endif
	if (usig0) then
		call nfs(nf90_inquire_dimension (ncid, 1, len=i))
		allocate (rad_pass(i), rad_dsig0(i), mask(i))
		call nfs(nf90_inq_varid (ncid, 'pass', varid))
		call nfs(nf90_get_var (ncid, varid, rad_pass))
		call nfs(nf90_inq_varid (ncid, 'atmos_corr_sig0', varid))
		call nfs(nf90_get_var (ncid, varid, rad_dsig0))
		call nfs(nf90_close (ncid))
	endif

	! Process passes
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo

	! Deallocate patch file info
	if (usig0) deallocate (rad_pass, rad_dsig0, mask)
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Patch SARAL data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --ssb                     Add hybrid SSB' / &
'  --wet                     Shift MWR wet prior to 2013-10-22' / &
'  --wind                    Compute wind speed' / &
'  --all                     All of the above' / &
'  --sig0                    Update sigma0 with new (Patch2 draft) attenuation')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n),wet(n),swh(n),ssb(n),sig0(n),dsig0(n),dsig0_old(n),wind(n),t
real(eightbytereal), parameter :: time0 = 909058527d0 ! 2013-10-22 12:15:27

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename(len_trim(S%dataroot)+2:))

! Process data records

call rads_get_var (S, P, 'sig0_ka', sig0, .true.)

! Update sigma0

if (usig0) then
	call rads_get_var (S, P, 'dsig0_atmos_ka', dsig0_old, .true.)
	mask = (rad_pass == pass)
	if (count(mask) == n) then
		dsig0 = pack(rad_dsig0, mask) * 1d-4
		sig0 = sig0 + dsig0_old - dsig0
	else
		dsig0 = dsig0_old
	endif
endif

! Compute SSB

if (lssb) then
	call rads_get_var (S, P, 'swh_ka', swh, .true.)
	do i = 1,n
		t = swh(i)
		if (t < 0d0) t = 0d0
		if (t > 12d0) t = 12d0
		ssb(i) = grid_lininter (issb_hyb,sig0(i),t)
	enddo
endif

! Shift MWR wet tropo prior to 2013-10-22 12:15:27

if (P%start_time > time0) lwet = .false.
if (lwet) then
	call rads_get_var (S, P, 'time', time, .true.)
	call rads_get_var (S, P, 'wet_tropo_rad', wet, .true.)
	where (time < time0) wet = wet - 6.4d-3
endif

! Compute wind speed

if (lwind) wind = wind_ecmwf (sig0, .true.)

! If nothing changed, stop here

if (.not.(usig0 .or. lssb .or. lwet .or. lwind)) then
	write (*,552) 0
	return
endif

! Write out all the data

call rads_put_history (S, P)
! if (usig0) call rads_put_var (S, P, 'sig0', sig0)
if (usig0) call rads_put_var (S, P, 'dsig0_atmos_nn_ka', dsig0)
if (lssb) call rads_def_var (S, P, 'ssb_hyb')
if (lssb) call rads_put_var (S, P, 'ssb_hyb', ssb)
if (lwet) call rads_put_var (S, P, 'wet_tropo_rad', wet)
if (lwind) call rads_put_var (S, P, 'wind_speed_alt', wind)

write (*,552) n
end subroutine process_pass

end program rads_fix_sa
