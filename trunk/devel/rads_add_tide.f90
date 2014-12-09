!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2014  Remko Scharroo (Altimetrics LLC)
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
integer(fourbyteint), parameter :: mfes = 2, mgot = 5
type(festideinfo) :: fesinfo(mfes)
type(gottideinfo) :: gotinfo(mgot)
character(len=6), parameter :: nfes(mfes) = (/'fes04 ','fes12 '/)
character(len=6), parameter :: ngot(mgot) = (/'got00 ','got47 ','got48 ','got49 ','got410'/)
type(grid) :: sininfo, cosinfo
type(fes) :: fesinfo1

! Command line arguments

integer(fourbyteint) :: cyc, pass
character(len=rads_cmdl) :: models = '', path
logical :: do_ptide=.false., do_stide=.false., do_lptide=.false., do_annual=.false., do_fes(mfes)=.false., do_got(mgot)=.false.

! Other variables

integer(fourbyteint) :: j, jdum, i0, i1

! Initialise

call synopsis ('--head')
call rads_set_options ('m: models:')
call rads_init (S)

! Check for options

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('m', 'models')
		models = rads_opt(j)%arg
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
	case ('annual')
		do_annual = .true.
		call parseenv ('${ALTIM}/data/DTU10/DTU10ANN_', path)
		if (grid_load (trim(path)//'cos.nc',cosinfo) /= 0 .or. &
			grid_load (trim(path)//'sin.nc',sininfo) /= 0) call rads_exit ('Error loading grid')
	case ('fes04', 'fes2004')
		do_fes(1) = .true.
		call festideinit('FES2004',.true.,fesinfo(1))
	case ('fes12', 'fes2012')
		do_fes(2) = .true.
		jdum = fes_init(fesinfo1,fes_tide,fes_mem,'FES2012/all')
		fesinfo(2)%haveload = .false. ! No load tide in FES2012
	case ('got00')
		do_got(1) = .true.
		call gottideinit('GOT00.2',.true.,gotinfo(1))
	case ('got47')
		do_got(2) = .true.
		call gottideinit('GOT4.7',.true.,gotinfo(2))
	case ('got48')
		do_got(3) = .true.
		call gottideinit('GOT4.8',.true.,gotinfo(3))
	case ('got49')
		do_got(4) = .true.
		call gottideinit('GOT4.9',.true.,gotinfo(4))
	case ('got410')
		do_got(5) = .true.
		call gottideinit('GOT4.10c',.true.,gotinfo(5))
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

if (do_annual) then
	call grid_free(cosinfo)
	call grid_free(sininfo)
endif
do j = 1,mfes
	if (do_fes(j)) call festidefree(fesinfo(j))
enddo
do j = 1,mgot
	if (do_got(j)) call gottidefree(gotinfo(j))
enddo
if (do_ptide) call poletidefree

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Add tide models to RADS data', flag=flag)) return
call synopsis_devel ('')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -m, --models=MODEL[,...]  Select tide models' // &
'Currently available MODELs are:'/ &
'  fes04  : FES2004 ocean and load tide'/ &
'  fes12  : FES2012 ocean tide'/ &
'  got00  : GOT00.2 ocean and load tide'/ &
'  got47  : GOT4.7 ocean and load tide'/ &
'  got48  : GOT4.8 ocean and load tide'/ &
'  got49  : GOT4.9 ocean and load tide'/ &
'  got410 : GOT4.10 ocean and load tide'// &
'In addition, several of the following MODEL indicators can be used:'/ &
'  ptide  : Pole tide'/ &
'  stide  : Solid earth tide'/ &
'  lptide : Long-period tides'/ &
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
real(eightbytereal), parameter :: k2 = 0.302, h2 = 0.609, h2k2 = h2 / (1 + k2)
real(eightbytereal) :: time(n), lon(n), lat(n), surface_type(n), &
	otide_sp(n), otide_lp(n), ltide_sp(n), ltide_lp(n), lptide_eq(n), lptide_mf(n)
integer(fourbyteint) :: i, j

call log_pass (P)

! Get time, location, and surface_type

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'surface_type', surface_type, .true.)

! Reset time reference at the start of each pass

fesinfo(:)%t_nodal = 1d30
gotinfo(:)%t_nodal = 1d30
if (do_fes(2)) call fes_set_nodal_time (fesinfo1, 1d30)

! Define output variables

call rads_put_history (S, P)

if (do_ptide) call rads_def_var (S, P, 'tide_pole')
if (do_stide) call rads_def_var (S, P, 'tide_solid')

do j = 1,mfes
	if (do_fes(j)) then
		call rads_def_var (S, P, 'tide_ocean_' // nfes(j))
		if (fesinfo(j)%haveload) call rads_def_var (S, P, 'tide_load_' // nfes(j))
	endif
enddo

do j = 1,mgot
	if (do_got(j)) then
		call rads_def_var (S, P, 'tide_ocean_' // ngot(j))
		if (gotinfo(j)%haveload) call rads_def_var (S, P, 'tide_load_' // ngot(j))
	endif
enddo

if (do_lptide) then
	call rads_def_var (S, P, 'tide_equil')
	if (any(do_fes)) call rads_def_var (S, P, 'tide_non_equil')
endif

if (do_annual) call rads_def_var (S, P, 'mss_annual')

! Process data records

if (do_lptide .or. any(do_fes) .or. any(do_got)) then
	do i = 1,n
		call lpetide (time(i), lat(i), 1, lptide_eq(i), lptide_mf(i))
	enddo
else
	lptide_eq = 0d0
	lptide_mf = 0d0
endif

if (do_ptide) then
	do i = 1,n
		otide_sp(i) = poletide (time(i)/86400d0+46066d0, lat(i), lon(i))
	enddo
	! Over land/lake use Love number h2, not 1+k2
	where (surface_type > 1.5d0) otide_sp = otide_sp * h2k2
	call rads_put_var (S, P, 'tide_pole', otide_sp)
endif

if (do_stide) then
	do i = 1,n
		otide_sp(i) = etide_ce(time(i)/86400d0+46066d0,lat(i),lon(i))
	enddo
	call rads_put_var (S, P, 'tide_solid', otide_sp)
endif

! FES models
do j = 1,mfes
	if (.not.do_fes(j)) cycle
	if (j == 1) then ! FES2014
! In order to allow parallelisation we do one measurement first (which initialises), then do the next n-1 in parallel
		call festide(fesinfo(j),time(1),lat(1),lon(1),otide_sp(1),otide_lp(1),ltide_sp(1),ltide_lp(1))
!$omp parallel do shared(fesinfo,time,lat,lon,otide_sp,otide_lp,ltide_sp,ltide_lp,n) private(i)
		do i = 2,n
			call festide(fesinfo(j),time(i),lat(i),lon(i),otide_sp(i),otide_lp(i),ltide_sp(i),ltide_lp(i))
		enddo
!$omp end parallel do
		! FES2004: First remove equilibrium part of Mm,Mf,Mtm,MSqm from long-period tides, then add equilibrium part
		otide_lp = otide_lp - lptide_mf + lptide_eq
		if (fesinfo(j)%haveload) call rads_put_var (S, P, 'tide_load_'//nfes(j), ltide_sp + ltide_lp)
	else ! Other
		do i = 1,n
			jdum = fes_eval(fesinfo1,time(i),lat(i),lon(i),otide_sp(i),otide_lp(i))
		enddo
		! FES2012: otide_lp already includes both non-equilibrium and equilibrium tide
	endif
	call rads_put_var (S, P, 'tide_ocean_'//nfes(j), otide_sp + otide_lp)
enddo

! Write equilibrium and non-equilibrium long-period tides (from FES solution)
if (do_lptide) then
	call rads_put_var (S, P, 'tide_equil', lptide_eq)
	if (any(do_fes)) call rads_put_var (S, P, 'tide_non_equil', otide_lp - lptide_eq)
endif

! GOT models
do j = 1,mgot
	if (.not.do_got(j)) cycle
! In order to allow parallelisation we do one measurement first (which initialises), then do the next n-1 in parallel
	call gottide(gotinfo(j),time(1),lat(1),lon(1),otide_sp(1),ltide_sp(1))
!$omp parallel do shared(gotinfo,time,lat,lon,otide_sp,ltide_sp,n) private(i)
	do i = 2,n
		call gottide(gotinfo(j),time(i),lat(i),lon(i),otide_sp(i),ltide_sp(i))
	enddo
!$omp end parallel do
	if (j == 2 .and. any(do_fes)) then ! GOT4.7
		! GOT4.7: Take long-period ocean and load tide from FES2004 (for consistency with old practices only)
		call rads_put_var (S, P, 'tide_ocean_'//ngot(j), otide_sp + otide_lp)
		if (gotinfo(j)%haveload) call rads_put_var (S, P, 'tide_load_'//ngot(j), ltide_sp + ltide_lp)
	else
		! Other: just add equilibrium tide
		call rads_put_var (S, P, 'tide_ocean_'//ngot(j), otide_sp + lptide_eq)
		if (gotinfo(j)%haveload) call rads_put_var (S, P, 'tide_load_'//ngot(j), ltide_sp)
	endif
enddo

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
