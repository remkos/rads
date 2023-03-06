!-----------------------------------------------------------------------
! Copyright (c) 2011-2022  Remko Scharroo
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

!*rads_add_mfwam -- Add variables from Meteo-France's WAM model to RADS data
!+
! This program adds up to 3 variables from the WAM model as run by Meteo-France
! These variables are:
! VHM0 = swh_mfwam = Significant wave height in m
! VMDR = mean_wave_direction = Mean wave direction in deg
! VTM02 = mean_wave_period = Mean wave period in s
!
! The variables are stored in daily files with 5'x5' grids at 3-hour intervals.
! All three variables and all 8 daily temporal intervals can be found
! in a single daily NetCDF file. These files can be found in
! ${ALTIM}/data/mfwam/%Y/%m/mfwamglocep_%Y%m%d.nc.
!
! The grids contain invalid values, whose mask changes with time, because
! of ice cover. The grids are linearly interpolated in space and time
! and require a minimum weight of 0.5 for the value to be used.
!
! The UTC times in the files are 03, 06, 09, 12, 15, 18, 21, 24 (00 of next day)
!
! usage: rads_add_mfwam [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_mfwam

use rads
use rads_misc
use rads_devel
use rads_netcdf
use netcdf

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cyc, pass
logical :: update = .false.

! Data elements

character(rads_cmdl) :: path
integer(fourbyteint), parameter :: nvar=3, i2min=-32767
integer(fourbyteint) :: mjd, mjdold=-99999, j, mvar = 0, nx = 0, ny = 0, nt = 0, ios, smear=0
real(eightbytereal) :: xmin, xmax, ymin, ymax, wmin = 0.5d0

type :: var_
	character(len=5) :: ncname
	character(len=20) :: radsname
	integer(twobyteint), allocatable :: grid(:,:,:)
	integer(twobyteint) :: fillvalue
	real(eightbytereal) :: add_offset, scale_factor
endtype
type(var_) :: var(nvar)

! Initialise

call synopsis ('--head')
call rads_set_options ('sdpu swh direction period all wmin: smear: update')
call rads_init (S)

! Get template for path name

call parseenv ('${ALTIM}/data/mfwam/%Y/%m/mfwamglocep_%Y%m%d.nc', path)

! Check all options

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('s', 'swh')
		call add_var ('VHM0', 'swh_mfwam')
	case ('d', 'direction')
		call add_var ('VMDR', 'mean_wave_direction')
	case ('p', 'period')
		call add_var ('VTM02', 'mean_wave_period')
	case ('all')
		call add_var ('VHM0', 'swh_mfwam')
		call add_var ('VMDR', 'mean_wave_direction')
		call add_var ('VTM02', 'mean_wave_period')
	case ('wmin')
		read (rads_opt(j)%arg, *, iostat=ios) wmin
	case ('smear')
		read (rads_opt(j)%arg, *, iostat=ios) smear
	case ('u', 'update')
		update = .true.
	end select
enddo
if (mvar == 0) call rads_exit ('No variables selected. Use the appropriate options.')

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata, mvar)
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Add variables from Meteo-France WAM model to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -s, --swh                 Add SWH (swh_mfwam)' / &
'  -d, --direction           Add mean wave direction (mean_wave_direction)' / &
'  -p, --period              Add mean wave peiod (mean_wave_period)' / &
'  --all                     All of the above' / &
'  --wmin=WMIN               Minumum total weight for interpolation (default: 0.5)' / &
'  --smear=SEC               Smear values into NaN areas by up to SEC seconds' / &
'  -u, --update              Update files only when there are changes')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n, mvar)
integer(fourbyteint), intent(in) :: n, mvar
real(eightbytereal) :: time(n), lat(n), lon(n), wave(n,mvar), tmp(n), x, y, t, w(2,2,2), z(2,2,2), wsum, dt, dtmin
integer(fourbyteint) :: i, j, k, ik, ix, iy, it, idx(n)
logical :: err

call log_pass (P)

! Get time and location

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)

! Process data records

do i = 1,n

! Determine the current day

	t = time(i)/86400
	mjd = int(t)

! Load new grids when entering new day

	if (mjd /= mjdold) then
		if (mjd == mjdold+1) then
			! Replace first grid with last grid and load new set of grids
			do j = 1,mvar
				var(j)%grid(:,:,1) = var(j)%grid(:,:,nt+1)
			enddo
			err = get_mfwam(mjd, 2)
		else
			! Replace both sets of grids
			err = get_mfwam(mjd-1, 1) .or. get_mfwam(mjd, 2)
		endif
		if (err) then
			call log_string ('Warning: no MFWAM data for current time')
			call log_records (0)
			stop
		endif
		mjdold = mjd
	endif

! Skip data beyond latitude limits

	if (lat(i) < ymin .or. lat(i) > ymax) then
		wave(i,:) = nan
		cycle
	endif

! Linearly interpolate in space and time
! First set the 3D grid indices

	x = modulo(lon(i) - xmin, 360d0) / 360d0 * nx + 1
	y = (lat(i) - ymin) / (ymax - ymin) * (ny - 1) + 1
	t = (t - mjd) * nt + 1
	ix = int(x)
	iy = int(y)
	it = int(t)
	x = x - ix
	y = y - iy
	t = t - it

! Compute the weights

	w(1,1,:) = (1d0 - x) * (1d0 - y)
	w(1,2,:) = (1d0 - x) * (      y)
	w(2,1,:) = (      x) * (1d0 - y)
	w(2,2,:) = (      x) * (      y)
	w(:,:,1) =  w(:,:,1) * (1d0 - t)
	w(:,:,2) =  w(:,:,2) * (      t)

! Set weights to zero for defaulted points, then sum the weights

	where (var(1)%grid(ix:ix+1,iy:iy+1,it:it+1) == var(1)%fillvalue) w = 0
	wsum = sum(w)

! If total weight is less than 0.5 set to nan

	if (wsum == 0 .or. wsum < wmin) then
		wave(i,:) = nan
		cycle
	endif

! Compute weighted averages;

	do j = 1,mvar
		z = var(j)%grid(ix:ix+1,iy:iy+1,it:it+1) * var(j)%scale_factor + var(j)%add_offset
		if (var(j)%ncname == 'VMDR') then ! Treat direction through vector sum
			z = z * rad
			wave(i,j) = atan2 (sum(sin(z) * w), sum(cos(z) * w)) / rad
		else
			wave(i,j) = sum(z * w) / wsum
		endif
	enddo
enddo

! "Smear" points into areas with NaNs

if (smear > 0) then
	do i = 1,n
		idx(i) = i
		if (isan_(wave(i,1))) cycle
		dtmin = (smear + 0.5d0) * S%dt1hz
		do k = 1, smear
			ik = i-k
			if (ik < 1) exit
			dt = abs(time(ik)-time(i))
			if (isan_(wave(ik,1)) .and. dt < dtmin) then
				dtmin = dt
				idx(i) = ik
				exit
			endif
		enddo
		do k = 1, smear
			ik = i+k
			if (ik > n) exit
			dt = abs(time(ik)-time(i))
			if (isan_(wave(ik,1)) .and. dt < dtmin) then
				dtmin = dt
				idx(i) = ik
				exit
			endif
		enddo
	enddo
	wave(:,:) = wave(idx(:),:)
endif

! If requested, check for changes first

if (update) then
	i = rads_verbose; rads_verbose = -1 ! Temporarily suspend warning
	call rads_get_var (S, P, var(1)%radsname, tmp, .true.)
	rads_verbose = i
	do i = 1,n
		if (isnan_(tmp(i)) .and. isnan_(wave(i,1))) cycle
		if (isnan_(tmp(i)) .or. isnan_(wave(i,1))) exit
		if (nint(tmp(i)/var(1)%scale_factor) /= nint(wave(i,1)/var(1)%scale_factor)) exit
	enddo
	if (i > n) then	! No changes
		call log_records (0)
		return
	endif
endif

! Define data fields

call rads_put_history (S, P)
do i = 1,mvar
	call rads_def_var (S, P, var(i)%radsname)
enddo

! Store data fields

do i = 1,mvar
	call rads_put_var (S, P, var(i)%radsname, wave(:,i))
enddo

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Define new variable to use
!-----------------------------------------------------------------------

subroutine add_var (ncname, radsname)
character(len=*) :: ncname, radsname
mvar = mvar + 1
var(mvar)%ncname = ncname
var(mvar)%radsname = radsname
end subroutine add_var

!-----------------------------------------------------------------------
! Get WaveWatch3 data
!-----------------------------------------------------------------------

logical function get_mfwam (mjd, istart)
integer(fourbyteint), intent(in) :: mjd, istart
integer(fourbyteint) ::	ncid, dimid, varid, l, strf1985, j
real(eightbytereal), allocatable :: x(:), y(:), t(:)
character(len=rads_naml) :: filenm

600 format ('(',a,')')
1300 format (a,': ',a)

! Determine file name

get_mfwam = .true.

l = strf1985(filenm, path, mjd*86400)

! Open input file

write (*,600,advance='no') trim(basename(filenm))
if (nft(nf90_open(filenm,nf90_nowrite,ncid))) then
	write (*,1300) 'Error opening file',filenm(:l)
	return
endif

! For the very first grid, get dimentsions and ranges

if (nx == 0) then
	call nfs(nf90_inq_dimid(ncid,'longitude',dimid))
	call nfs(nf90_inquire_dimension(ncid,dimid,len=nx))
	call nfs(nf90_inq_dimid(ncid,'latitude',dimid))
	call nfs(nf90_inquire_dimension(ncid,dimid,len=ny))
	call nfs(nf90_inq_dimid(ncid,'time',dimid))
	call nfs(nf90_inquire_dimension(ncid,dimid,len=nt))
	allocate (x(nx),y(ny),t(nt))
	call nfs(nf90_inq_varid(ncid,'longitude',varid))
	call nfs(nf90_get_var(ncid,varid,x))
	xmin = x(1) ; xmax = x(nx)
	call nfs(nf90_inq_varid(ncid,'latitude',varid))
	call nfs(nf90_get_var(ncid,varid,y))
	ymin = y(1) ; ymax = y(ny)
	call nfs(nf90_inq_varid(ncid,'time',varid))
	call nfs(nf90_get_var(ncid,varid,t))
	deallocate (x,y,t)
	do j = 1,mvar
		allocate (var(j)%grid(nx+1,ny,nt+1))
	enddo
endif

! Load from NetCDF file all selected variables, with scale_factor, add_offset and fillvalue
! If istart = 1, read the last temporal grid into the first slot
! If istart = 2, read the whole grid into the second and following slots
! Also copy the first column to the last to simplefy interpolation

do j = 1,mvar
	if (nft(nf90_inq_varid(ncid,var(j)%ncname,varid))) call rads_exit ('Error finding variable')
	if (nft(nf90_get_att(ncid,varid,'scale_factor',var(j)%scale_factor))) var(j)%scale_factor = 1d0
	if (nft(nf90_get_att(ncid,varid,'add_offset',var(j)%add_offset))) var(j)%add_offset = 0d0
	if (nft(nf90_get_att(ncid,varid,'_FillValue',var(j)%fillvalue))) var(j)%fillvalue = 32767
	if (istart == 1) then
		if (nft(nf90_get_var(ncid,varid,var(j)%grid(1:nx,1:ny,1:1),(/1,1,nt/)))) &
			call rads_exit ('Error reading data grid')
		var(j)%grid(nx+1,:,1) = var(j)%grid(1,:,1)
	else
		if (nft(nf90_get_var(ncid,varid,var(j)%grid(1:nx,1:ny,2:nt+1)))) call rads_exit ('Error reading data grid')
		var(j)%grid(nx+1,:,2:nt+1) = var(j)%grid(1,:,2:nt+1)
	endif
enddo

l = nf90_close(ncid)
get_mfwam = .false.
end function get_mfwam

end program rads_add_mfwam
