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

!*rads_add_orbit -- Add new orbit to RADS data
!
! This program adjusts the contents of RADS altimeter data files
! with orbital altitudes inpterpolated from external files,
! applies time tag bias, and/or stores orbital altitude rate.
!
! usage: rads_add_orbit [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_orbit

use rads
use rads_grid
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P
type(grid) :: info

! Command line arguments

integer(fourbyteint) :: i,cyc,pass
logical :: equator=.false.,range=.false.,doppler=.false.,range2=.false.,rate=.false.
character(160) :: name='',dir='',gridnm='',var_old='',arg
real(eightbytereal) :: tbias=0d0,maxrms=1d30,loc=0d0,dt=1d0,chirp=0d0

! Other variables

integer(fourbyteint) :: ellipse=0,ios

! Formats

550  format (a)
551  format (a,' ...',$)

! Scan command line for options

info%ntype = 0
call synopsis
call rads_set_options (' model: name: dir: tbias:: rate loc-6 loc-7 equator range range2 doppler grid:: chirp: dt: all flag:')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('loc-7')
		loc = 1d-7
	case ('loc-6')
		loc = 1d-6
	case ('rate')
		rate = .true.
	case ('dt')
		read (rads_opt(i)%arg,*,iostat=ios) dt
	case ('doppler')
		doppler = .true.
	case ('grid')
		gridnm = rads_opt(i)%arg
		if (gridnm == '') gridnm='EGM96/egm96.nc'
	case ('chirp')
		read (rads_opt(i)%arg,*,iostat=ios) chirp
	case ('equator')
		equator = .true.
	case ('tbias')
		tbias = 1d30
		read (rads_opt(i)%arg,*,iostat=ios) tbias
	case ('flag')
		read (rads_opt(i)%arg,*,iostat=ios) var_old,maxrms
		maxrms = maxrms/1d2
	case ('range2')
		range = .true.
		range2 = .true.
	case ('range')
		range = .true.
	case ('all')
		loc = 1d-6
		rate = .true.
		tbias = 1d30
		equator = .true.
	case ('name')
		name = rads_opt(i)%arg
	case ('dir')
		dir = rads_opt(i)%arg
	end select
enddo

! Overrule long_name if name is specified

if (name /= '') S%sel(1)%info%long_name = trim(name) // ' orbital altitude'

! If dir is not given, figure out directory from variable selection

if (dir == '') dir = S%sel(1)%info%source_file

! Add prefix $ALTIM/data/ODR.<satellite>/ to directory name

if (dir(:1) == '/' .or. dir(:2) == './') then
else
	call getenv ('ALTIM', arg)
	dir = trim(arg) // '/data/ODR.' // trim(S%satellite) // '/' // dir
endif

! Figure out time bias

if (tbias < 1d20) then
	tbias = tbias * 1d-3
else if (S%sat == 'e1') then
	tbias = 1.5d-3
else if (S%sat == 'e2') then
	tbias = 1.3d-3
else if (S%sat == 'n1') then
	tbias = 0d0
else if (S%sat == 'gs') then
	tbias = 2.382d-3
else if (S%sat == 'g1') then
	tbias = 0.945d-3
else if (S%sat == 'c2') then
	tbias = -8.2d-3
else
	call rads_exit ('Timing bias not known for this satellite')
endif

! Determine base ellipsoid on orbit files (Geosat ellipsoid = 3, TOPEX = 0)

select case (S%sat)
case ('e1', 'e2', 'n1')
	if (S%sel(1)%name == 'alt_jgm3' .or. S%sel(1)%name == 'alt_dgme04') ellipse = 3
end select

! Load grid when requested.
! Prefix $ALTIM/data/ if grid filenname does not start in / or ./
! Store limits. Set x-limits to mid-point +/- 180.

if (gridnm /= '') then
	write (*,551) 'Loading grid '//trim(gridnm)
	if (gridnm(:1) == '/' .or. gridnm(:2) == './') then
		if (grid_load(gridnm,info) /= 0) call rads_exit ('Error loading grid.')
	else
		call getenv ('ALTIM', arg)
		if (grid_load(trim(arg)//'/data/'//gridnm,info) /= 0) call rads_exit ('Error loading grid.')
	endif
	write (*,550) 'done'
endif

! Determine factor for Doppler correction

if (.not.doppler .and. info%ntype == 0) then
	chirp = 0d0		! No Doppler correction requested
else if (chirp /= 0d0) then
	chirp = chirp * 1d-3
else if (S%sat == 'gs') then
	chirp = 4.32d-3
else if (S%sat == 'e1' .or. S%sat == 'e2') then
	chirp = 0.836d-3
else
	call rads_exit ('Chirp parameter not known for this satellite')
endif

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
if (rads_version ('$Revision$', 'Add orbit to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:'/ &
'  --var=NAME                Set variable name for output orbit (required)'/ &
'  --name=STR                Specify name of orbit'/ &
'  --dir=STR                 Specify orbit directory path'/ &
'  --tbias[=VAL]             Add timing bias (ms) (defaults for various satellites)'/ &
'  --rate                    Compute and store orbital altitude rate'/ &
'  --loc-6                   Change latitude/longtitude (10^-6 deg)'/ &
'  --loc-7                   Change latitude/longtitude (10^-7 deg)'/ &
'  --equator                 Change equator crossing time and longitude'/ &
'  --range                   Change field 601 from SSH to range'/ &
'  --doppler                 Correct range for Doppler, based on altitude rate and chirp parameter'/ &
'  --grid                    Use EGM96 grid to add slope-induced Doppler'/ &
'  --grid=NAME               Use named grid to add slope-induced Doppler'/ &
'  --chirp=F                 Specify chirp parameter other than default (ms)'/ &
'  --dt=DT                   Specify time interval for rate, etc (s), default: 1'/ &
'  --all                     Same as: --tbias --rate --loc-6 --equator'// &
'The orbit error flag can be set using the --flag= option:'/ &
'  --flag=VAR,RMS            Specify variable of reference orbit and maximum rms with new orbit (cm)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: getorb,kstep=0,k,i,i0,i1
integer(twobyteint) :: flag
real(eightbytereal) :: utc(n),lat(n),lon(n),alt(n),alt_old(n),alt_rate(n),flags(n),alt_cnes(n),range_ku(n), &
	range_c(n),drange_fm(n),dalt(-1:1),dlat(-1:1),dlon(-1:1),t,f,rms,xx,yy,zz,xx0,yy0,dhellips
logical :: asc, cryofix

! Formats

550  format (a)
551  format (a,' ...',$)
552  format (i5,' records changed')

write (*,551) trim(P%filename)

! Determine if we need to fix Cryosat data
! Only if original data is LRM L2 and this is the first radsp_orbit run

cryofix = (S%sat == 'c2' .and. index(P%history,'SIR_LRM_2') > 0 .and. index(P%history,'radsp_orbit') == 0 &
	.and. index(P%history,'rads_add_orbit') == 0)

! Figure out step size for orbit rate or Doppler correction

if (rate .or. doppler .or. info%ntype > 0 .or. cryofix) kstep = 1

! Load time and determine ascending/descending

call rads_get_var (S, P, 'time', utc, .true.)
asc = (mod(pass,2) == 1)

! If comparison with old orbit is requested, load it together with the flags.

if (var_old /= '') then
	call rads_get_var (S, P, var_old, alt_old, .true.)
	call rads_get_var (S, P, 'flags', flags, .true.)
endif
rms = 0d0

! If Doppler correction or conversion of SSH to range is requested, load field 601.

if (cryofix) call rads_get_var (S, P, 'alt_cnes', alt_cnes, .true.)
if (range .or. doppler .or. info%ntype > 0 .or. cryofix) call rads_get_var (S, P, 'range_ku', range_ku, .true.)
if (range2) call rads_get_var (S, P, 'range_c' ,range_c, .true.)

! Try to load variable drange_fm also when slope-induced Doppler is requested

if (info%ntype > 0) then
	call rads_get_var (S, P, 'drange_fm', drange_fm, .true.)
	if (S%error /= rads_noerr) drange_fm = 0d0
endif

! Process data records

do i = 1,n

! Convert raw SSH to range using the orbit at the original time tag.
! No ellipsoid correction since both SSH and orbit are over the same ellipsoid.

	if (range) then
		if (getorb(utc(i),yy,xx,zz,dir,.true.) > 0) then
			write (*,552) 0
			return
		endif
		range_ku(i) = zz - range_ku(i)
		if (range2) range_c(i) = zz - range_c(i)
	endif

! Correct time tag for timing bias

	utc(i) = utc(i) + tbias

! Compute locations and orbital altitudes

	do k = -kstep,kstep
		if (getorb(utc(i)+k*dt/2,dlat(k),dlon(k),dalt(k),dir,.true.) > 0) then
			write (*,552) 0
			return
		endif
	enddo

! Compute orbit rate

	if (kstep > 0) alt_rate(i) = (dalt(1) - dalt(-1)) / dt

! Compute Doppler correction based on orbit rate

	if (doppler) range_ku(i) = range_ku(i) + chirp * alt_rate(i)

! Do fix for Cryosat: add 'half a second' of orbit rate to original orbit and range
! The orbit is, in addition, corrected for timing bias.

	if (cryofix) then
		alt_cnes(i) = alt_cnes(i) + (0.448130d0 + tbias) * alt_rate(i)
		range_ku(i) = range_ku(i) + (0.448130d0) * alt_rate(i)
	endif

! Compute grid slope and corresponding Doppler correction

	if (info%ntype > 0) then
		do k = -1,1,2
			dalt(k) = grid_splinter(info,dlon(k),dlat(k))
		enddo
		xx = chirp * (dalt(-1) - dalt(1))
		range_ku(i) = range_ku(i) - drange_fm(i) + xx
		drange_fm(i) = xx
	endif

! Compute orbital altitude, latitude and longitude

	alt(i) = dalt(0) + dhellips(ellipse,dlat(0))
	lat(i) = dlat(0)
	if (loc > 5d-7 .and. dlon(0) < 0d0) dlon(0) = dlon(0) + 360d0
	lon(i) = dlon(0)

! Keep record of this pass' orbit difference RMS

	rms = rms + (alt_old(i) - alt(i))**2

enddo

! Adjust start and end time

P%start_time = utc(1)
P%end_time   = utc(n)

! Scan for the equator crossing.
! If there is already a value for equator_time, assume it is an educated guess.
! Otherwise, it always has to be at most 1/4 rev before the end or
! 1/4 rev (approx 1509) after the start of the track.

if (equator) then
	xx0 = 0
	yy0 = 0
	if (P%equator_time > 0d0) then
		i0 = nint(P%equator_time) - 120
		i1 = nint(P%equator_time) + 120
	else
		i0 = nint(P%end_time) - 1510
		i1 = nint(P%start_time) + 1510
	endif
	do i = i0,i1
		t = i
		if (getorb(t,yy,xx,zz,dir,.true.) > 0) then
			write (*,552) 0
			return
		endif
		if (yy*yy0 >= 0) then
			! no equator crossing (dlat and yy0 same sign)
		else if (yy > yy0 .eqv. asc) then
			! crossing is in the right sense
			f = yy/(yy-yy0)
			P%equator_time = t - f
			P%equator_lon = xx - f*(xx-xx0)
			exit
		else
			write (*,*) yy0,yy,asc
		endif
		xx0 = xx
		yy0 = yy
		if (i == i1) write (*,550) 'error scanning for equator crossing'
	enddo
endif
if (S%debug >= 1) write (*,*) 'eq:', P%equator_time, P%equator_lon

! Update history and (re)define variables that may be new and will be written

call rads_put_passinfo (S, P)
call rads_put_history (S, P)
if (loc /= 0d0) then
	call rads_def_var (S, P, 'lat', scale_factor=loc)
	call rads_def_var (S, P, 'lon', scale_factor=loc)
endif
call rads_def_var (S, P, S%sel(1))
if (rate) call rads_def_var (S, P, 'alt_rate')
if (range .or. range2) call rads_def_var (S, P, 'range_ku')
if (range2) call rads_def_var (S, P, 'range_c')
if (info%ntype > 0) call rads_def_var (S, P, 'drange_fm')

! Now write out data

if (tbias /= 0d0) call rads_put_var (S, P, 'time', utc)
if (loc /= 0d0) then
	call rads_put_var (S, P, 'lat', lat)
	call rads_put_var (S, P, 'lon', lon)
endif
if (cryofix) call rads_put_var (S, P, 'alt_cnes', alt_cnes)
call rads_put_var (S, P, S%sel(1), alt)
if (rate) call rads_put_var (S, P, 'alt_rate', alt_rate)
if (range2) then
	call rads_put_var (S, P, 'range_ku', range_ku)
	call rads_put_var (S, P, 'range_c', range_c)
else if (range .or. doppler .or. info%ntype > 0 .or. cryofix) then
	call rads_put_var (S, P, 'range_ku', range_ku)
endif
if (info%ntype > 0) call rads_put_var (S, P, 'drange_fm', drange_fm)

! If orbit difference RMS exceeds maximum, set all orbit error flags;
! otherwise: clear them.

rms = sqrt(rms/n)
if (var_old /= '') then
	do i = 1,n
		flag = nint(flags(i),twobyteint)
		if (rms > maxrms) then
			flag = ibset(flag,15)
		else
			flag = ibclr(flag,15)
		endif
		flags(i) = flag
	enddo
	call rads_put_var (S, P, 'flags', flags)
endif
if (S%debug >= 1) write (*,*) 'rms:', pass, rms

write (*,552) n
end subroutine process_pass

!-----------------------------------------------------------------------
! Interpolate a table using 8th order Legendre interpolation
!-----------------------------------------------------------------------

subroutine intab8 (ndim, nvec, table, t1, tn, t, vec)
integer(fourbyteint) :: ndim, nvec, itrel
real(eightbytereal) ::  table(ndim,nvec), t1, tn, t, vec(ndim), trel
trel = (t-t1)/(tn-t1)*(nvec-1) + 1
itrel = max(1,min(nint(trel)-4,nvec-8))
trel = nint((trel-itrel)*60d6)/60d6   ! Fix time resolution to 1 microsecond
call inter8(ndim,trel,table(1,itrel),vec)
end subroutine intab8

end program rads_add_orbit
