!-----------------------------------------------------------------------
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

!*radsxogen -- RADS crossover generator
!+
program radsxogen

! This program "generates" crossovers from RADS data, i.e.,
! it searches crossovers within a batch of passes of altimeter
! data (single satellite crossovers) or between different
! batches of altimeter data (dual satellite crossovers).
!
! The output is a netCDF file containing the crossover times
! and locations of the crossovers, with the information of
! the respective passes. With a consecutive run of radsxosel
! more information can be added to the crossover file.
!
! This program is based in part on previous programs, max and max2,
! developed by Marc Naeije and Remko Scharroo at TU Delft/DEOS.
!-----------------------------------------------------------------------
use rads
use netcdf
use rads_netcdf
use rads_misc
use rads_time

integer(fourbyteint), parameter :: msat = 20, vbase = 13, mtrk = 500000
real(eightbytereal) :: dt(msat)
type(rads_sat) :: S(msat)
type(rads_pass) :: P
integer(fourbyteint) :: i, j, nsat = 0, nsel = 0, reject = -1, ios, ncid, dimid(3), start(2), varid(vbase)
logical :: duals = .true., singles = .true., batches = .false., l
character(len=rads_cmdl) :: satlist, filename = 'radsxogen.nc'
character(len=rads_naml) :: legs = 'undetermined - undetermined'
character(len=1) :: interpolant = 'q'
integer(fourbyteint) :: inter_half = 3, inter_points, max_gap = 3

type :: nr_
	integer :: test		! Number potential xovers tested
	integer :: shallow	! Number xovers with shallow angle crossing
	integer :: xdat		! Number of potential xovers with too few data on pass
	integer :: ins		! Number of potential xovers with no intersection
	integer :: gap		! Number xovers with too few data around xover
	integer :: xdt		! Number xovers with time difference too large
	integer :: xout		! Number of xovers written
	integer :: trk		! Number of tracks
endtype
type(nr_) :: nr

type :: trk_
	real(eightbytereal) :: equator_lon, equator_time, start_time, end_time
	integer(twobyteint) :: nr_alt, nr_xover, satid, cycle, pass
endtype
type(trk_) :: trk(mtrk)

integer(fourbyteint), allocatable :: key(:), idx(:), trkid(:,:)

! Initialize RADS or issue help
call synopsis
call rads_set_options ('dsbi:g:r:: dual single batch interpolant: gap: reject-on-nan:: dt:')
call rads_init (S)
if (any(S%error /= rads_noerr)) call rads_exit ('Fatal error')
dt = nan
dt(1) = -0.5d0

! Start with this-is message
l = rads_version ()

! If no sat= is given, exit
if (S(1)%sat == '') call rads_exit ('Need at least one sat= option')

! Determine how many satellites and cycles.
! Also check that the same number of variables have been selected for each satellite
do i = 1,msat
	if (S(i)%sat == '') exit
	if (S(i)%nsel == 0) call rads_parse_varlist (S(i), 'sla')
	if (S(i)%nsel /= S(1)%nsel) call rads_exit ('Unequal amount of variables on -V for different satellites')
	nsat = i
enddo
nsel = S(1)%nsel

! Scan command line arguments
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('r', 'reject-on-nan')
		call rads_parse_r_option (S(1), rads_opt(i)%opt, rads_opt(i)%arg, reject)
	case ('s', 'single')
		duals = .false.
	case ('d', 'dual')
		singles = .false.
	case ('b', 'batch')
		batches = .true.
	case ('i', 'interpolant')
		select case (rads_opt(i)%arg)
		case ('a':'z')
			interpolant = rads_opt(i)%arg(1:1)
			read (rads_opt(i)%arg(2:), *, iostat=ios) inter_half
		case default
			read (rads_opt(i)%arg, *, iostat=ios) inter_half
		end select
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('g', 'gap')
		read (rads_opt(i)%arg, *, iostat=ios) max_gap
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
		if (max_gap < 1) call rads_exit ('Option -g <gap> : needs <gap> > 0')
	case ('dt')
		read (rads_opt(i)%arg, *, iostat=ios) dt
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case (' ')
		filename = rads_opt(i)%arg
	end select
enddo

inter_points = inter_half * 2
select case (interpolant)
case ('q', 'c', 's')
	if (inter_half < 2) call rads_exit ('Option -i[c|q|s]<n> : needs <n> > 1')
case default
	if (inter_half < 1) call rads_exit ('Option -i[a|l|n]<n> : needs <n> > 0')
end select

! Set all unspecified dt equal to the previous and create a string of sat names
satlist = S(1)%sat
do i = 2,nsat
	if (isnan_(dt(i))) dt(i) = dt(i-1)
	satlist = trim(satlist) // ' ' // S(i)%sat
enddo

! If SLA is among the results, remember which index that is
do i = 1,nsel
	if (S(1)%sel(i)%info%datatype == rads_type_sla) then
		if (reject == -1) reject = i
	endif
enddo

! Write info about this job
if (nsat == 1) duals = .false.
write (*,600) timestamp(), trim(S(1)%command), trim(filename)
if (singles .and. duals) then
	write (*,550) 'Processing single and dual satellite crossovers'
	if (nsat == 2) legs = 'ascending_pass_or_' // S(1)%sat // ' descending_pass_or_' // S(2)%sat
else if (singles) then
	write (*,550) 'Processing single satellite crossovers'
	legs = 'ascending_pass descending_pass'
else if (duals) then
	write (*,550) 'Processing dual satellite crossovers'
	if (nsat == 2) legs = S(1)%sat // ' ' // S(2)%sat
else
	write (*,550) 'No single or duals satellite crossovers are processed. Aborting'
	stop
endif
write (*,610)
do i = 1,nsat
	write (*,620) S(i)%satellite, S(i)%cycles(1:2), S(i)%passes(1:2), dt_secs(S(i),dt(i))/86400d0
enddo
550 format (a)
600 format (/'Created: ',a,' UTC: ',a/'Output file name: ',a)
610 format (/'Satellite  Cycles    Passes  Delta-time cutoff (days)')
620 format (a,i5,i4,2i5,f8.3)

! Do further initialisations
nr = nr_ (0, 0, 0, 0, 0, 0, 0, 0)

! Open output netCDF file
call nfs (nf90_create (filename, nf90_write, ncid))
call nfs (nf90_def_dim (ncid, 'xover', nf90_unlimited, dimid(1)))
call nfs (nf90_def_dim (ncid, 'leg', 2, dimid(2)))

! To use general netCDF creation machinary, we trick the library a bit here
P%ncid = ncid
P%filename = filename
P%rw = .true.
S(1)%time%info%ndims = 2

! Define lat and lon (1D) and time (2D)
call rads_def_var (S(1), P, S(1)%lat)
call rads_def_var (S(1), P, S(1)%lon)
call rads_def_var (S(1), P, S(1)%time)
varid(1) = S(1)%lat%info%varid
varid(2) = S(1)%lon%info%varid
varid(3) = S(1)%time%info%varid

! Define selected variables
do i = 1,nsel
	S(1)%sel(i)%info%ndims = 2
	call rads_def_var (S(1), P, S(1)%sel(i))
enddo

! Define other stuff
call def_var (ncid, 'track', 'track number', '', nf90_int4, dimid(2:1:-1), varid(4))
call nfs (nf90_put_att (ncid, varid(4), 'comment', 'Internal track number relating to satid/cycle/pass below'))
call nfs (nf90_put_att (ncid, nf90_global, 'Conventions', 'CF-1.5'))
call nfs (nf90_put_att (ncid, nf90_global, 'title', 'RADS 4.0 crossover file'))
call nfs (nf90_put_att (ncid, nf90_global, 'institution', 'Altimetrics / NOAA / TU Delft'))
call nfs (nf90_put_att (ncid, nf90_global, 'references', 'RADS Data Manual, Issue 4.0'))
call nfs (nf90_put_att (ncid, nf90_global, 'history', timestamp()//' UTC: '//trim(S(1)%command)))
call nfs (nf90_put_att (ncid, nf90_global, 'legs', trim(legs)))
call nfs (nf90_enddef (ncid))

! We are now ready to compare the different batches of passes
do i = 1,nsat
	do j = i,nsat
		if (batches .and. i == j) cycle
		call xogen_batch (S(i), S(j), dt(i), dt(j))
	enddo
enddo

! Print statistics
write (*, 500) nr
500 format (/ &
'Number of crossings tested  :',i9/ &
'- Shallow crossings         :',i9/ &
'- Too few data on pass      :',i9/ &
'- No intersection           :',i9/ &
'- Too few data around xover :',i9/ &
'- Too large time difference :',i9/ &
'= Number crossovers written :',i9/ &
'Number of tracks            :',i9/)

call rads_stat (S)

if (nr%trk == 0) then
	call nfs (nf90_close (ncid))
	stop
endif

! Sort the track information in order satid/cycle/pass
allocate (key(nr%trk), idx(nr%trk), trkid(2,nr%xout))
forall (i = 1:nr%trk)
	idx(i) = i
	key(i) = trk(i)%satid * 10000000 + trk(i)%cycle * 10000 + trk(i)%pass
end forall
call iqsort (idx, key)
if (idx(1) == 0) call rads_exit ('stack size for iqsort is too small')

! Add the track info to the netCDF file
call nfs (nf90_redef (ncid))
call nfs (nf90_def_dim (ncid, 'track', nr%trk, dimid(3)))
call def_var (ncid, 'satid', 'satellite ID', '', nf90_int1, dimid(3:3), varid(5))
call nfs (nf90_put_att (ncid, varid(5), 'flag_values', int(S(1:nsat)%satid, onebyteint)))
call nfs (nf90_put_att (ncid, varid(5), 'flag_meanings', trim(satlist)))
call nfs (nf90_put_att (ncid, varid(5), 'comment', 'Satellite IDs relate to the different missions'))
call def_var (ncid, 'cycle', 'cycle number', '', nf90_int2, dimid(3:3), varid(6))
call def_var (ncid, 'pass', 'pass number', '', nf90_int2, dimid(3:3), varid(7))
call def_var (ncid, 'equator_lon', 'longitude of equator crossing', 'degrees_east', nf90_double, dimid(3:3), varid(8))
call def_var (ncid, 'equator_time', 'time of equator crossing', S(1)%time%info%units, nf90_double, dimid(3:3), varid(9))
call def_var (ncid, 'start_time', 'start time of track', S(1)%time%info%units, nf90_double, dimid(3:3), varid(10))
call def_var (ncid, 'end_time', 'end time of track', S(1)%time%info%units, nf90_double, dimid(3:3), varid(11))
call def_var (ncid, 'nr_xover', 'number of crossovers along track', '', nf90_int2, dimid(3:3), varid(12))
call def_var (ncid, 'nr_alt', 'number of measurements along track', '', nf90_int2, dimid(3:3), varid(13))
call nfs (nf90_enddef (ncid))

call nfs (nf90_put_var (ncid, varid( 5), trk(idx)%satid))
call nfs (nf90_put_var (ncid, varid( 6), trk(idx)%cycle))
call nfs (nf90_put_var (ncid, varid( 7), trk(idx)%pass))
call nfs (nf90_put_var (ncid, varid( 8), trk(idx)%equator_lon))
call nfs (nf90_put_var (ncid, varid( 9), trk(idx)%equator_time))
call nfs (nf90_put_var (ncid, varid(10), trk(idx)%start_time))
call nfs (nf90_put_var (ncid, varid(11), trk(idx)%end_time))
call nfs (nf90_put_var (ncid, varid(12), trk(idx)%nr_xover))
call nfs (nf90_put_var (ncid, varid(13), trk(idx)%nr_alt))

! Renumber the track number for the xover data
forall (i = 1:nr%trk)
	key(idx(i)) = i	! key becomes the inverse of idx
end forall
call nfs (nf90_get_var (ncid, varid(4), trkid))
trkid(1,:) = key(trkid(1,:))
trkid(2,:) = key(trkid(2,:))
call nfs (nf90_put_var (ncid, varid(4), trkid))
deallocate (key, idx, trkid)

call nfs (nf90_close (ncid))

call rads_end (S)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('Generate altimeter crossovers from RADS')) return
call rads_synopsis
write (stderr,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -d, --dual                Do dual satellite crossovers only'/ &
'  -s, --single              Do single satellite crossovers only'/ &
'  -b, --batch               Do crossovers between different batches only'/ &
'  -i, --interpolant [n|s|a|l|q|c][N]'/ &
'                            Interpolate 2N along-track values to crossover by picking (n)earest neighbor,'/ &
'                            or by cubic (s)pline, or by (a)veraging, or by (l)inear, (q)uadratic or (c)ubic'/ &
'                            polynomial fit; optionally add number of points N required on BOTH sides of the'/ &
'                            crossover for interpolation. Default: q3'/ &
'  -g, --gap GAP             Specify the maximum gap between two nearest points to crossover, in 1-Hz intervals;' / &
'                            also sets maximum gaps between 1st and last point of interpolation window to (2N+1)*GAP' / &
'  -r, --reject-on-nan VAR   Reject xovers if variable VAR on -V specifier is NaN'/ &
'  -r #                      Reject xovers if data item number # on -V specifier is NaN'/ &
'  -r 0, -r none, -r         Do not reject xovers with NaN values'/ &
'  -r n, -r any              Reject xovers if any value is NaN'/ &
'                      Note: If no -r option is given -r sla is assumed'/ &
'  --dt DT                   Limit crossover time interval to number of days'/ &
'                            use negative number to specify interval in fraction of cycles (default = -0.5)'/ &
'  FILENAME                  Specify output FILENAME (default: radsxogen.nc)')
stop
end subroutine synopsis

!***********************************************************************
! Function to convert requested time interval to seconds
! Note that we DO NOT use repeat_days because that would give a terribly
! long time for geodetic orbits (like c2 and j1/c). Thus we use the length
! of the subcycle instead

function dt_secs (S, dt)
type(rads_sat), intent(in) :: S
real(eightbytereal), intent(in) :: dt
real(eightbytereal) :: dt_secs
if (dt > 0d0) then
	dt_secs = dt * 86400d0
else
	dt_secs = -dt * S%phase%pass_seconds * S%phase%passes
endif
end function dt_secs

!***********************************************************************
! Define variables for the output file

subroutine def_var (ncid, var, long_name, units, nctype, dimid, varid)
use netcdf
use rads_netcdf
integer(fourbyteint), intent(in) :: ncid, nctype, dimid(:)
character(len=*), intent(in) :: var, long_name, units
integer(fourbyteint), intent(out) :: varid
call nfs (nf90_def_var (ncid, trim(var), nctype, dimid, varid))
call nfs (nf90_put_att (ncid, varid, 'long_name', trim(long_name)))
if (units /= '') call nfs (nf90_put_att (ncid, varid, 'units', trim(units)))
end subroutine def_var

!***********************************************************************
! Batch process passes selected for satellites S1 and S2 (which could be the same).
! Passes of satellite S1 are cycled through (opened and closed) in order.
! Passes of satellite S2 are held in memory until no longer needed, but still processed in order.

subroutine xogen_batch (S1, S2, dt1, dt2)
type(rads_sat), intent(inout) :: S1, S2
real(eightbytereal), intent(in) :: dt1, dt2
integer(fourbyteint) :: cycle1, pass1, cycle2, pass2, step
type(rads_pass) :: P1
type(rads_pass), pointer :: P2, top, prev
real(eightbytereal) :: t0, t1, t2, t3, dt

! Skip singles or duals when they are not wanted
! If singles, match only ascending and descending, hence step = 2
if (S1%sat == S2%sat) then
	if (.not.singles) return
	step = 2
else
	if (.not.duals) return
	step = 1
endif

! See if there is any chance of a temporal overlap between S1 and S2
dt = min (dt_secs(S1,dt1), dt_secs(S2,dt2))
t0 = rads_cycle_to_time (S1, S1%cycles(1)) - dt
t1 = rads_cycle_to_time (S1, S1%cycles(2)) + S1%phase%repeat_days * 86400d0 + dt
t2 = rads_cycle_to_time (S2, S2%cycles(1))
t3 = rads_cycle_to_time (S2, S2%cycles(2)) + S2%phase%repeat_days * 86400d0
if (t0 > t3 .or. t1 < t2) return

! Initialize pass and cycle numbers of S2
! Base the starting cycle number for S2 on the start time of the first cycle of S1
pass2 = S2%passes(1) - 1
cycle2 = max(S2%cycles(1), rads_time_to_cycle (S2, t0))

! Nullify all pointers
nullify (P2, top, prev)

600 format (a,3i6)

! Cycle through cycles and passes for S1
do cycle1 = S1%cycles(1), S1%cycles(2)
	do pass1 = S1%passes(1), S1%passes(2), step

		! Open a pass for S1, to be crossed with any pass of S2
		if (rads_verbose > 2) write (*,600) 'open S1', cycle1, pass1
		call rads_open_pass (S1, P1, cycle1, pass1)
		if (P1%ndata <= 0) then
			if (rads_verbose > 2) write (*,600) 'empty S1', P1%cycle, P1%pass
			call rads_close_pass (S1, P1)
			cycle
		endif
		call load_data (S1, P1)
		if (rads_verbose > 2) write (*,600) 'close S1', P1%cycle, P1%pass, P1%ncid
		call rads_close_pass (S1, P1, .true.) ! Close the pass file, but keep all its info

		! Limit the time selection for P2 based on the time limits of P1
		t0 = P1%start_time - dt
		t1 = P1%end_time   + dt

		! Walk through a linked list of already loaded passes
		P2 => top
		do while (associated(P2))
			if (P2%end_time < t0) then ! Release passes that are far in the past
				top => P2%next ! Reassign the top of the list to the next pass (if any)
				nullify (prev)
				if (rads_verbose > 2) write (*,600) 'release S2', P2%cycle, P2%pass
				call rads_close_pass (S2, P2)
				deallocate (P2)
				P2 => top
			else if (P2%start_time > t1) then ! Do not consider (yet) passes in the far future
				exit
			else ! Cross P1 and P2 if P2 holds any data
				if (P2%ndata > 0) call xogen_passes (S1, P1, S2, P2, dt)
				pass2 = P2%pass
				cycle2 = P2%cycle
				prev => P2
				P2 => P2%next
			endif
		enddo

		! Now load new passes within the time range
		do while (.not.associated(P2))
			pass2 = pass2 + step
			if (pass2 > S2%passes(2) .or. pass2 > S2%phase%passes) then ! Wrap to next cycle
				cycle2 = cycle2 + 1
				if (cycle2 > S2%cycles(2)) exit
				pass2 = S2%passes(1) - 1 + step ! Start with pass "1" or "2"
			endif
			allocate (P2)
			if (rads_verbose > 2) write (*,600) 'open S2', cycle2, pass2
			call rads_open_pass (S2, P2, cycle2, pass2)
			if (P2%end_time < t0) then ! There seems to be a reason why P2%ndata == 0 is kept
				if (rads_verbose > 2) write (*,600) 'release S2', P2%cycle, P2%pass
				call rads_close_pass (S2, P2)
				deallocate (P2)
				cycle
			endif
			if (associated(prev)) then ! Append to linked list
				prev%next => P2
			else ! Start new linked list
				top => P2
			endif
			call load_data (S2, P2)
			if (rads_verbose > 2) write (*,600) 'close S2', P2%cycle, P2%pass, P2%ncid
			call rads_close_pass (S2, P2, .true.) ! Close the pass file, but keep all its info
			if (P2%start_time > t1) exit
			if (P2%ndata > 0) call xogen_passes (S1, P1, S2, P2, dt)
			! Point to next (empty) slot
			prev => P2
			P2 => P2%next
		enddo

		! Clear any memory of pass P1
		if (rads_verbose > 2) write (*,600) 'release S1', P1%cycle, P1%pass
		call rads_close_pass (S1, P1)
	enddo
enddo

! Dump remaining open passes

do while (associated(top))
	P2 => top
	if (rads_verbose > 2) write (*,600) 'release S2', P2%cycle, P2%pass
	call rads_close_pass (S2, P2)
	top => P2%next
	deallocate (P2)
enddo
end subroutine xogen_batch

!***********************************************************************
! Load the data for a pass, remove invalids, store it back into an
! expanded P%tll.

subroutine load_data (S, P)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
real(eightbytereal), pointer :: data(:,:), tll(:,:)
real(eightbytereal) :: x
integer(fourbyteint) :: i
logical, allocatable :: keep(:)

if (P%ndata == 0) return

! Allocate enough memory for the data and get all of it
allocate (data(P%ndata,nsel),keep(P%ndata))
do i = 1,nsel
	call rads_get_var (S, P, S%sel(i), data(:,i))
enddo

! Determine which records are kept
if (reject < 0) then	! Reject if any is NaN
	forall (i=1:P%ndata)
		keep(i) = .not.any(isnan_(data(i,:)))
	end forall
else if (reject > 0) then
	keep = isan_(data(:,reject))
else
	keep = .true.
endif

! Create new tll array with only kept records, expanded with all data fields
P%ndata = count(keep)
if (P%ndata == 0) return
allocate (tll(P%ndata,nsel+3))
forall (i=1:3)
	tll(:,i) = pack (P%tll(:,i),keep)
end forall
forall (i=1:nsel)
	tll(:,i+3) = pack (data(:,i),keep)
end forall

! Group all longitudes to within 180 degrees of the equator crossing longitude
x = P%equator_lon - 180d0
tll(:,3) = modulo (tll(:,3) - x, 360d0) + x

! Get rid of old tll array and replace it with the new one
deallocate (P%tll,data,keep)
P%tll => tll

! Store track info
nr%trk = nr%trk + 1
if (nr%trk > mtrk) call rads_exit ('Too many tracks')
trk(nr%trk) = trk_ (P%equator_lon, P%equator_time, tll(1,1), tll(P%ndata,1), &
	int(P%ndata,twobyteint), 0_twobyteint, S%satid, int(P%cycle,twobyteint), int(P%pass,twobyteint))
P%trkid = nr%trk
end subroutine load_data

!***********************************************************************
! Compute the intersection of pass P1 of satellite S1 and pass P2 of satellite S2.

subroutine xogen_passes (S1, P1, S2, P2, dt)
type(rads_sat), intent(inout) :: S1, S2
type(rads_pass), intent(inout) :: P1, P2
real(eightbytereal), intent(in) :: dt
real(eightbytereal) :: shiftlon, x, y, t1, t2
real(eightbytereal), allocatable :: data(:,:)
integer :: i, i1, i2, j1, j2

! Count number of calls
nr%test = nr%test + 1

! If inclinations of S1 and S2 differ less than 1 degree, match only ascending and descending passes
if (abs(S1%inclination - S2%inclination) < 1d0 .and. modulo(P1%pass - P2%pass, 2) == 0) then
	nr%shallow = nr%shallow + 1
	return
endif

! If there are fewer than 2 points on either pass, then there is no crossing
if (P1%ndata < 2 .or. P2%ndata < 2) then
	nr%xdat = nr%xdat + 1
	return
endif

! Determine if the tracks have common latitude
if (min(P1%tll(1,2), P1%tll(P1%ndata,2)) > max(P2%tll(1,2), P2%tll(P2%ndata,2)) .or. &
	max(P1%tll(1,2), P1%tll(P1%ndata,2)) < min(P2%tll(1,2), P2%tll(P2%ndata,2))) then
	nr%ins = nr%ins + 1
	return
endif

! Determine a longitude shift for the second pass to make it closer than 180 degrees from the first
if (P2%equator_lon - P1%equator_lon > 180d0) then
	shiftlon = -360d0
else if (P2%equator_lon - P1%equator_lon < -180d0) then
	shiftlon = 360d0
else
	shiftlon = 0d0
endif

600 format (a,6f12.6)

if (rads_verbose > 2) then
	write (*,'(a,4i6)') 'processing',P1%cycle,P1%pass,P2%cycle,P2%pass
	write (*,600) 'equator',P1%equator_lon,P2%equator_lon,shiftlon
	write (*,600) 'x-ranges',P1%tll(1,3),P1%tll(P1%ndata,3),P2%tll(1,3)+shiftlon,P2%tll(P2%ndata,3)+shiftlon
endif

! Determine if the passes have common longitudes
if (min(P1%tll(1,3), P1%tll(P1%ndata,3)) > max(P2%tll(1,3), P2%tll(P2%ndata,3))+shiftlon .or. &
	max(P1%tll(1,3), P1%tll(P1%ndata,3)) < min(P2%tll(1,3), P2%tll(P2%ndata,3))+shiftlon) then
	nr%ins = nr%ins + 1
	return
endif

! Use a modified sweep-line algorithm to look for intersections
call sweep_line (P1%tll(:,3), P1%tll(:,2), P1%tll(:,1), P2%tll(:,3)+shiftlon, P2%tll(:,2), P2%tll(:,1), &
	i1, j1, i2, j2, x, y, t1, t2)
if (i1 == 0) then ! No crossing
	nr%ins = nr%ins + 1
	return
endif

! Move x back to proper range
if (x < S1%lon%info%limits(1)) x = x + 360d0
if (x > S1%lon%info%limits(2)) x = x - 360d0

if (rads_verbose > 2) then
	write (*,'(a,6f15.3)') 'time',P1%tll(i1,1),t1,P1%tll(j1,1),P2%tll(i2,1),t2,P2%tll(j2,1)
	write (*,600) 'lat ',P1%tll(i1,2),y ,P1%tll(j1,2),P2%tll(i2,2),y ,P2%tll(j2,2)
	write (*,600) 'lon ',P1%tll(i1,3),x ,P1%tll(j1,3),P2%tll(i2,3),x ,P2%tll(j2,3)
endif

! See if time interval exceeds limits
if (abs(t1-t2) <= dt) then
	! Continue only for small time interval, not if NaN
else
	if (rads_verbose > 2) write (*,'(a,4f15.3)') 'rejected: dt: ',t1,t2,abs(t1-t2),dt
	nr%xdt = nr%xdt + 1
	return
endif

! For the later interpolation of the variables along-track we will use <inter_points> points,
! <inter_half> on each side of the crossover.
i1 = max(i1, j1) - inter_half ! Indices of first of <inter_points> points for interpolation
i2 = max(i2, j2) - inter_half
j1 = i1 + inter_points - 1    ! Indices of last of <inter_points> points for interpolation
j2 = i2 + inter_points - 1

! Check that interval (i1:j1) and (i2:j2) or still entirely within the pass.
! The central two points of the interval have to be maximum <max_gap> 1-Hz intervals apart.
! The furthest two points of the interval have to be maximum <gap>*<inter_half+1> 1-Hz intervals apart.
if (large_gap (S1, P1, i1, j1) .or. large_gap (S2, P2, i2, j2)) then
	nr%gap = nr%gap + 1
	return
endif

! Interpolate data along each track
allocate (data(2,nsel))
call interpolate (t1, P1%tll(i1:j1,:), data(1,:))
call interpolate (t2, P2%tll(i2:j2,:), data(2,:))

! Write the data to file
nr%xout = nr%xout + 1
trk(P1%trkid)%nr_xover = trk(P1%trkid)%nr_xover + 1_twobyteint
trk(P2%trkid)%nr_xover = trk(P2%trkid)%nr_xover + 1_twobyteint
start = (/ 1, nr%xout /)
call nfs (nf90_put_var (ncid, varid(1), nint4(y / S1%lat%info%scale_factor), start(2:2)))
call nfs (nf90_put_var (ncid, varid(2), nint4(x / S1%lon%info%scale_factor), start(2:2)))
call nfs (nf90_put_var (ncid, varid(3), (/ t1, t2 /), start))
call nfs (nf90_put_var (ncid, varid(4), (/ P1%trkid, P2%trkid /), start))
do i = 1,nsel
	call rads_put_var (S1, P, S1%sel(i), data(:,i), start)
enddo
deallocate (data)

end subroutine xogen_passes

!***********************************************************************
! Interpolate the data to the crossover along the track

subroutine interpolate (t, tll, xoval)
use spline
real(eightbytereal), intent(in) :: t, tll(1-inter_half:,:)	! This makes indices of central points: 0 and 1
real(eightbytereal), intent(out) :: xoval(:)
integer(fourbyteint) :: i
integer(fourbyteint), parameter :: n = 4
real(eightbytereal) :: xc, x(1-inter_half:inter_half), y(1-inter_half:inter_half,nsel), aa(1-inter_half:inter_half)
real(eightbytereal) :: w(inter_points), b(n), c(n), d(n)

! x = time, y = data value
! Reduce time and values to numbers relative to the left (early) side of the central interval
x = tll(:,1) - tll(0,1)
xc = t - tll(0,1)
forall (i = 1:nsel)
	y(:,i) = tll(:,3+i) - tll(0,3+i)
end forall

select case (interpolant)
case ('a')   ! Average
	do i = 1,nsel
		xoval(i) = sum(y(:,i)) / inter_points
	enddo
case ('l')   ! Linear fit
	w = 1d0
	do i = 1,nsel
		call least_set (inter_points,x,y(:,i),w,2,b,c,d)
		call least_val (2,b,c,d,xc,xoval(i))
	enddo
case ('q')   ! Quadratic polynomial fit
	w = 1d0
	do i = 1,nsel
		call least_set (inter_points,x,y(:,i),w,3,b,c,d)
		call least_val (3,b,c,d,xc,xoval(i))
	enddo
case ('c')   ! Cubic polynomial fit
	w = 1d0
	do i = 1,nsel
		call least_set (inter_points,x,y(:,i),w,4,b,c,d)
		call least_val (4,b,c,d,xc,xoval(i))
	enddo
case ('s')   ! Cubic spline interpolation
	do i = 1,nsel
		call spline_cubic_set (inter_points,x,y(:,i),0,0d0,0,0d0,aa)
		xoval(i) = xc * (y(1,i)/x(1) - (aa(1)/6d0+aa(0)/3d0) * x(1)  &
			+ xc * (0.5d0 * aa(0) + xc * (aa(1)-aa(0)) / (6d0 * x(1))))
	enddo
case default ! Nearest neighbour
	if (xc <= 0.5d0 * x(1)) then
		xoval = 0d0 ! Technically, value at left of central interval
	else
		xoval = y(1,:) ! Value at right of centre point
	endif
end select

! Restore value at left of central interval
xoval = xoval + tll(0,4:)
end subroutine interpolate

!***********************************************************************
! Determine if gap between nearest points to the crossover or furthest of six points is too large

function large_gap (S, P, i, j)
type(rads_sat), intent(in) :: S
type(rads_pass), intent(in) :: P
integer(fourbyteint), intent(in) :: i, j
logical :: large_gap
large_gap = (i < 1 .or. j > P%ndata .or. &
	nint((P%tll(j,1) - P%tll(i,1)) / S%dt1hz) > (inter_half+1)*max_gap .or. &
	nint((P%tll(i+inter_half,1) - P%tll(j-inter_half,1)) / S%dt1hz) > max_gap)
end function large_gap

!***********************************************************************
! This is a modified type of sweep-line algorithm to determine whether two
! line segments cross, and if so, at what indices.
! x1 and x2 need to be continuously increasing or decreasing.
! Upon return i1,j1 and i2,j2 are the indices of the two sides of the intervals
! at which the lines intersect, or 0 when there is no intersection.
! The returned values xc,yc are the coordinates of the crossover computed
! by crossing two spherical arcs through the end points.
! The times of the crossovers for the two passes are tc1 and tc2.

subroutine sweep_line (x1, y1, t1, x2, y2, t2, i1, j1, i2, j2, xc, yc, tc1, tc2)
real(eightbytereal), intent(in) :: x1(:), y1(:), t1(:), x2(:), y2(:), t2(:)
integer, intent(out) :: i1, j1, i2, j2
real(eightbytereal), intent(out) :: xc, yc, tc1, tc2
integer(fourbyteint) :: n1, n2, d1, d2
n1 = size(x1)
n2 = size(x2)
if (x1(2) > x1(1)) then
	i1 = 1
	d1 = 1
else
	i1 = n1
	d1 = -1
endif
if (x2(2) > x2(1)) then
	i2 = 1
	d2 = 1
else
	i2 = n2
	d2 = -1
endif
! i1,i2 are the indices of the left (west) sides of the segments
! j1,j2 are the indices of the right (east) sides of the segments
! d1,d2 are the directions of advance (-1 or +1) making x increase
do
	j1 = i1 + d1
	j2 = i2 + d2
	if (j1 < 1 .or. j1 > n1 .or. j2 < 1 .or. j2 > n2) exit
	! The lines certainly do not cross when the x-ranges do not overlap, so make sure of that first
	if (x1(j1) < x2(i2)) then
		i1 = j1
		cycle
	else if (x2(j2) < x1(i1)) then
		i2 = j2
		cycle
	endif
	if (intersect(x1(i1),y1(i1),t1(i1),x1(j1),y1(j1),t1(j1),x2(i2),y2(i2),t2(i2),x2(j2),y2(j2),t2(j2),xc,yc,tc1,tc2)) then
		if (rads_verbose > 2) write (*,'(a,4i6,2f12.6,2f15.3)') 'crossing',i1,j1,i2,j2,xc,yc,tc1,tc2
		return
	endif

	! Move the leftmost interval further to the right
	if (x1(j1) < x2(j2)) then
		i1 = j1
	else
		i2 = j2
	endif
enddo
if (rads_verbose > 2) write (*,'(a)') 'no crossing'
i1 = 0
i2 = 0
j1 = 0
j2 = 0
xc = 0d0
xc = xc / xc
yc = xc
end subroutine sweep_line

!***********************************************************************
! Does the line from (x1,y1) via (x2,y2) to (x3,y3) go counter-clockwise?

function ccw (x1, y1, x2, y2, x3, y3)
real(eightbytereal), intent(in) :: x1, y1, x2, y2, x3, y3
logical :: ccw
ccw = ((x2-x1) * (y3-y1) > (y2-y1) * (x3-x1))
end function ccw

!***********************************************************************
! The lines between two consecutive points on line 1 and line 2 intersect when neither combination
! of 3 points goes both clockwise or both counter-clockwise

function intersect (x11, y11, t11, x12, y12, t12, x21, y21, t21, x22, y22, t22, xc, yc, tc1, tc2)
use rads_misc
real(eightbytereal), intent(in) :: x11, y11, t11, x12, y12, t12, x21, y21, t21, x22, y22, t22
real(eightbytereal), intent(out) :: xc, yc, tc1, tc2
real(eightbytereal) :: v11(3), v12(3), v21(3), v22(3), vc(3)
logical :: intersect
if (ccw(x11,y11,x21,y21,x22,y22) .eqv. ccw(x12,y12,x21,y21,x22,y22)) then
	intersect = .false.
else if (ccw(x21,y21,x11,y11,x12,y12) .eqv. ccw(x22,y22,x11,y11,x12,y12)) then
	intersect = .false.
else
	intersect = .true.
	! Convert longitude and latitude of the end points to vectors
	v11 = xy2v (x11,y11)
	v12 = xy2v (x12,y12)
	v21 = xy2v (x21,y21)
	v22 = xy2v (x22,y22)
	! Let v1 and v2 be vectors perpendicular to planes through arcs 1 and 2
	! The cross product of v1 and v2 points to the crossing point
	vc = cross_product (cross_product (v11, v12), cross_product (v21, v22))
	! The result may be pointing the opposite direction to what is intended
	if (dot_product(vc,v11) < 0d0) vc = -vc
	! Now convert back to longitude and latitude
	call v2xy (vc, xc, yc)
	! Compute the times of the crossing points using proportionalities
	tc1 = acos(dot_product(vc,v11)) / acos(dot_product(v11,v12)) * (t12-t11) + t11
	tc2 = acos(dot_product(vc,v21)) / acos(dot_product(v21,v22)) * (t22-t21) + t21
endif
end function intersect

!***********************************************************************

function xy2v (x, y) result (v)
real(eightbytereal), intent(in) :: x, y
real(eightbytereal) :: v(3)
real(eightbytereal) :: xx,yy
xx = x * rad
yy = y * rad
v(1) = cos(xx)
v(2) = sin(xx)
v(1:2) = v(1:2) * cos(yy)
v(3) = sin(yy)
end function xy2v

!***********************************************************************

subroutine v2xy (v, x, y)
real(eightbytereal), intent(inout) :: v(3)
real(eightbytereal), intent(out) :: x, y
real(eightbytereal) :: c
c = sqrt(sum(v*v))
v = v / c
x = atan2(v(2),v(1)) / rad
y = asin(v(3)) / rad
end subroutine v2xy

!***********************************************************************

end program radsxogen
