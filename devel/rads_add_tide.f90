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

!*rads_add_tide -- Add tide models to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! with values computed from tide models. The tide models provide
! both ocean and load tide. Additionally, solid earth, long-period and
! pole tide can be computed.
!
! usage: rads_add_tide [data-selectors] [options]
!
! [options] indicate the tide models to be used.
!-----------------------------------------------------------------------
program rads_add_tide

use rads
use rads_grid
use rads_misc
use rads_devel
use tides

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint), parameter :: mfes = 2, mgot = 3
integer(fourbyteint) :: fes_mode = fes_mem
type(gottideinfo) :: gotinfo(mgot)
character(len=6), parameter :: nfes(mfes) = (/'fes14 ', 'fes22 '/)
character(len=6), parameter :: ngot(mgot) = (/'got48 ', 'got410', 'got51 '/)
type(grid) :: sininfo, cosinfo
type(fes) :: fes_long(mfes), fes_shrt(mfes), fes_load(mfes)
type(hrettideinfo) :: hretinfo
type(rads_var), pointer :: var

! Command line arguments

integer(fourbyteint) :: cyc, pass
character(len=rads_cmdl) :: models = '', path
logical :: do_ptide=.false., do_stide=.false., do_lptide=.false., do_hret=.false., do_annual=.false., &
	do_fes(mfes)=.false., do_got(mgot)=.false.

! Other variables

integer(fourbyteint) :: j, jdum, i0, i1

! Initialise

call synopsis ('--head')
call rads_set_options ('m: models: fes-io')
call rads_init (S)

! Check for options

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('m', 'models')
		models = rads_opt(j)%arg
	case ('fes-io')
		fes_mode = fes_io
	end select
enddo

! Which models are to be used?

i1 = 0
do
	if (.not.next_word (models, i0, i1)) exit
	select case (models(i0:i1-1))
	case ('ptide')
		do_ptide = .true.
		call poletideinit
	case ('stide')
		do_stide = .true.
	case ('lptide')
		do_lptide = .true.
	case ('hret')
		do_hret = .true.
		var => rads_varptr (S, 'tide_internal')
		call parseenv ('${ALTIM}/data/' // var%info%parameters, path)
		call hrettideinit(path, hretinfo)
	case ('annual')
		do_annual = .true.
		var => rads_varptr (S, 'mss_annual')
		jdum = index(var%info%parameters, ' ')
		call parseenv ('${ALTIM}/data/' // var%info%parameters(:jdum-1), path)
		if (grid_load (path,cosinfo) /= 0) call rads_exit ('Error loading grid')
		call parseenv ('${ALTIM}/data/' // var%info%parameters(jdum+1:), path)
		if (grid_load (path,sininfo) /= 0) call rads_exit ('Error loading grid')
	case ('fes14', 'fes2014')
		do_fes(1) = .true.
		var => rads_varptr (S, 'tide_ocean_fes14')
		jdum = fes_init(fes_long(1),fes_tide,  fes_mode,'FES2014/long_period_ocean_tide_extrapolated')
		jdum = fes_init(fes_shrt(1),fes_tide,  fes_mode,'FES2014/short_period_ocean_tide_extrapolated')
		jdum = fes_init(fes_load(1),fes_radial,fes_mode,'FES2014/load_tide')
	case ('fes22', 'fes2022')
		do_fes(2) = .true.
		var => rads_varptr (S, 'tide_ocean_fes22')
		jdum = fes_init(fes_long(2),fes_tide,  fes_mode,'FES2022b/long_period_ocean_tide_extrapolated')
		jdum = fes_init(fes_shrt(2),fes_tide,  fes_mode,'FES2022b/short_period_ocean_tide_extrapolated')
		jdum = fes_init(fes_load(2),fes_radial,fes_mode,'FES2022b/load_tide')
	case ('got48')
		do_got(1) = .true.
		call gottideinit('GOT4.8',.true.,gotinfo(1))
	case ('got410')
		do_got(2) = .true.
		call gottideinit('GOT4.10c_extrapolated',.true.,gotinfo(2))
	case ('got51')
		do_got(3) = .true.
		call gottideinit('GOT5.1',.true.,gotinfo(3))
	end select
enddo

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

! Free the allocated grids

if (do_ptide) call poletidefree
do j = 1,mfes
	if (.not.do_fes(j)) cycle
	call fes_delete(fes_long(j))
	call fes_delete(fes_shrt(j))
	call fes_delete(fes_load(j))
enddo
do j = 1,mgot
	if (do_got(j)) call gottidefree(gotinfo(j))
enddo
if (do_hret) call hrettidefree(hretinfo)
if (do_annual) then
	call grid_free(cosinfo)
	call grid_free(sininfo)
endif

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Add tide models to RADS data', flag=flag)) return
call synopsis_devel ('')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  --fes-io                  Do not load entire FES tide models into memory' / &
'  -m, --models MODEL[,...]  Select tide models' // &
'Currently available MODELs are:'/ &
'  fes14  : FES2014 ocean and load tide'/ &
'  fes22  : FES2022b ocean and load tide'/ &
'  got48  : GOT4.8 ocean and load tide'/ &
'  got49  : GOT4.9 ocean and load tide'/ &
'  got410 : GOT4.10 ocean and load tide'// &
'In addition, several of the following MODEL indicators can be used:'/ &
'  ptide  : Pole tide'/ &
'  stide  : Solid earth tide'/ &
'  lptide : Long-period tides'/ &
'  hret   : HRET internal tide'/ &
'  annual : Annual sea level variation')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: phase, co, si
real(eightbytereal), parameter :: pi = 4d0 * atan(1d0), t_2000 = 473299200d0, t_year = 365.25d0 * 86400d0
real(eightbytereal) :: time(n), lon(n), lat(n), &
	otide_sp(n), otide_lp(n), ltide_sp(n), ltide_lp(n), lptide_eq(n), lptide_mf(n), itide(n), itide_comp(6)
integer(fourbyteint) :: i, j

call log_pass (P)

! Get time and location

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)

! Long-period tide for GOT models when there is no FES model
if (any(do_got) .and. .not.any(do_fes)) then
	call rads_get_var (S, P, 'tide_equil', lptide_eq, .true.) ! Read existing field
	if (S%error /= 0) then ! If not existing, compute it with the analytical formula
		do i = 1,n
			call lpetide (time(i), lat(i), 1, lptide_eq(i), lptide_mf(i))
		enddo
	endif
endif

! Reset time reference at the start of each pass.
! This makes sure that the nodal arguments are always recomputed per pass, so it does not
! matter if the job run for one pass only or several.

gotinfo(:)%t_nodal = 1d30
hretinfo%t_nodal = 1d30
do j = 1,mfes
	if (.not.do_fes(j)) cycle
	call fes_set_nodal_time (fes_long(j), 1d30)
	call fes_set_nodal_time (fes_shrt(j), 1d30)
	call fes_set_nodal_time (fes_load(j), 1d30)
enddo

! Define output variables

call rads_put_history (S, P)

if (do_ptide) call rads_def_var (S, P, 'tide_pole')
if (do_stide) call rads_def_var (S, P, 'tide_solid')

do j = 1,mfes
	if (.not.do_fes(j)) cycle
	call rads_def_var (S, P, 'tide_ocean_' // nfes(j))
	call rads_def_var (S, P, 'tide_load_' // nfes(j))
enddo

do j = 1,mgot
	if (.not.do_got(j)) cycle
	call rads_def_var (S, P, 'tide_ocean_' // ngot(j))
	call rads_def_var (S, P, 'tide_load_' // ngot(j))
enddo

if (do_lptide) then
	call rads_def_var (S, P, 'tide_equil')
	if (any(do_fes)) call rads_def_var (S, P, 'tide_non_equil')
endif

if (do_hret) call rads_def_var (S, P, 'tide_internal')

if (do_annual) call rads_def_var (S, P, 'mss_annual')

! Process data records

! Pole tide
if (do_ptide) then
	do i = 1,n
		otide_sp(i) = poletide (time(i)/86400d0+46066d0, lat(i), lon(i))
	enddo
	call rads_put_var (S, P, 'tide_pole', otide_sp)
endif

! Solid earth tide
if (do_stide) then
	do i = 1,n
		otide_sp(i) = etide_ce(time(i)/86400d0+46066d0,lat(i),lon(i))
	enddo
	call rads_put_var (S, P, 'tide_solid', otide_sp)
endif

! FES2014 and later models
do j = 1,mfes
	if (.not.do_fes(j)) cycle
	! otide_lp already includes both non-equilibrium and equilibrium long-period tides
	! otide_sp is ignored here
	do i = 1,n
		jdum = fes_eval(fes_long(j), time(i), lat(i), lon(i), otide_sp(i), otide_lp(i))
	enddo
	! otide_sp is computed here. The long_period component is the LPE tide.
	do i = 1,n
		jdum = fes_eval(fes_shrt(j), time(i), lat(i), lon(i), otide_sp(i), lptide_eq(i))
	enddo
	do i = 1,n
		jdum = fes_eval(fes_load(j), time(i), lat(i), lon(i), ltide_sp(i), ltide_lp(i))
	enddo
	call rads_put_var (S, P, 'tide_ocean_'//nfes(j), otide_sp + otide_lp)
	call rads_put_var (S, P, 'tide_load_'//nfes(j), ltide_sp + ltide_lp)
enddo

! Write equilibrium and non-equilibrium long-period tides (from FES solution)
if (do_lptide) then
	call rads_put_var (S, P, 'tide_equil', lptide_eq)
	if (any(do_fes)) call rads_put_var (S, P, 'tide_non_equil', otide_lp + ltide_lp - lptide_eq)
endif

! GOT models
do j = 1,mgot
	if (.not.do_got(j)) cycle
	do i = 1,n
		call gottide(gotinfo(j), time(i), lat(i), lon(i), otide_sp(i), ltide_sp(i))
	enddo
	! Add equilibrium long-period tide to ocean tide
	call rads_put_var (S, P, 'tide_ocean_'//ngot(j), otide_sp + lptide_eq)
	call rads_put_var (S, P, 'tide_load_'//ngot(j), ltide_sp)
enddo

! Internal tide
if (do_hret) then
	do i = 1,n
		call hrettide (hretinfo, time(i), lat(i), lon(i), itide(i), itide_comp)
		itide(i) = sum(itide_comp(1:4))
!		write (*,'(f16.2,2f16.6,7f10.5)') time(i), lat(i), lon(i), itide(i), itide_comp
	enddo
	call rads_put_var (S, P, 'tide_internal', itide)
endif

! Annual mean sea surface variations
if (do_annual) then
	phase = (P%equator_time - t_2000) / t_year * 2d0 * pi
	co = cos(phase)
	si = sin(phase)
	do i = 1,n
		otide_lp(i) = co * grid_lininter(cosinfo,lon(i),lat(i)) + si * grid_lininter(sininfo,lon(i),lat(i))
	enddo
	call rads_put_var (S, P, 'mss_annual', otide_lp)
endif

call log_records (n)

end subroutine process_pass

end program rads_add_tide
