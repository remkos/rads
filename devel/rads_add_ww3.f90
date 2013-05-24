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

!*rads_add_ww3 -- Add variables from WaveWatch3 model to RADS data
!+
! This program adds 5 to 7 variables from the WaveWatch 3 model to
! the RADS data base. These variables are:
! m0 = wave height variance in m^2
! m1 = first order moment in m^2/s
! m2 = second order moment (velocity variance) in m^2/s^2
! m4 = fourth order moment (slope variance) (unitless)
! shs = swell height in m
! u = zonal wind speed in m/s
! v = meridionalwind speed in m/s
!
! The variables are stored in 1x1 degree grids at 6 hourly intervals.
! All five variables and all 4 daily temporal intervals can be found
! in a single netCDF file.
!
! These files can be found in $ALTIM/data/wwatch3. They are compressed
! by bzip2. The required files are unpacked and placed in /tmp before
! being used in this program.
!
! The grids contain invalid values, whose mask changes with time, because
! of ice cover. The grids are linearly interpolated in space and time
! and require a minimum weight of 0.5 for the value to be used.
!
! usage: rads_add_ww3 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_ww3

use rads
use rads_misc
use rads_devel
use rads_netcdf
use netcdf

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cyc, pass, var0=9, var1=0
logical :: update = .false.

! Data elements

character(rads_cmdl) :: path
integer(fourbyteint), parameter :: nx=361, ny=156, nhex=8, nvar=7, i2min=-32768
integer(fourbyteint) :: mjd, mjdold=-99999, j
real(eightbytereal) :: x0=-180d0, x1=180d0, y0=-77.5d0, y1=77.5d0, dx=1d0, dy=1d0, dz(nvar)
integer(twobyteint) :: grids(nx,ny,nhex,nvar), subgrid(2,2,2,nvar)
character(4) :: varnm(nvar) = (/ 'm0  ', 'm1  ', 'm2  ', 'm4  ', 'shs ', 'wu  ', 'wv  ' /)

! Initialise

call synopsis ('--head')
call rads_set_options ('wvu ww3 wind all update')
call rads_init (S)

! Get template for path name

call parseenv ('${ALTIM}/data/wwatch3/%Y/ww3_%Y%m%d.nc', path)

! Check all options
do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('w', 'ww3')
		var0 = 1 ; var1 = max(var1,5)
	case ('v', 'wind')
		var0 = min(var0,6) ; var1 = nvar
	case ('all')
		var0 = 1 ; var1 = nvar
	case ('u', 'update')
		update = .true.
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

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Add variables from WaveWatch3 model to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -w, --ww3                 Add WaveWatch3 variables m0, m1, m2, m4 and shs' / &
'  -v, --wind                Add U and V wind speed components' / &
'  --all                     All of the above' / &
'  -u, --update              Update files only when there are changes')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n), lat(n), lon(n), tmp(n), z(nvar), w, x, y, &
	f1, f2, f(2,2,2), ww3(n,nvar)
integer(fourbyteint) :: i, ix, iy, hex, var
logical :: err

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename)

! Get time and location

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)

! Process data records
do i = 1,n

! Skip data beyond latitude limits

	if (lat(i) < y0 .or. lat(i) > y1) then
		ww3(i,:) = nan
		cycle
	endif
	if (lon(i) > x1) lon(i) = lon(i) - 360d0

! Determine number of 6-hourly blobks since 1985

	f2 = time(i)/21600d0
	hex = int(f2)
	mjd = hex/4

! Load new grids when entering new day

	if (mjd == mjdold) then
		! Nothing
	else if (mjd == mjdold+1) then
		! Replace first set of grids with second set grids and load new set
		grids(:,:,1:4,:) = grids(:,:,5:8,:)
		err = get_ww3(mjd+1,grids(:,:,5:8,:))
		mjdold = mjd
	else
		! Replace both sets of grids
		err = get_ww3(mjd  ,grids(:,:,1:4,:)) .or. get_ww3(mjd+1,grids(:,:,5:8,:))
		mjdold = mjd
	endif
	if (err) then
		write (*,'(a,$)') 'No WW3 data for current time ...'
		write (*,552) 0
		stop
	endif

! Linearly interpolate in space and time

	f2 = f2 - hex
	f1 = 1d0 - f2

	x = (lon(i)-x0)/dx + 1
	y = (lat(i)-y0)/dy + 1
	ix = int(x)
	iy = int(y)
	x = x - ix
	y = y - iy

	f(1,1,:) = (1d0-x)*(1d0-y)
	f(1,2,:) = (1d0-x)*(    y)
	f(2,1,:) = (    x)*(1d0-y)
	f(2,2,:) = (    x)*(    y)
	f(:,:,1) = f(:,:,1)*f1
	f(:,:,2) = f(:,:,2)*f2

	hex = modulo(hex,4)+1
	subgrid = grids(ix:ix+1,iy:iy+1,hex:hex+1,:)

	w = 0d0
	z = 0d0
	do hex = 1,2
		do iy = 1,2
			do ix = 1,2
				if (subgrid(ix,iy,hex,1) /= i2min) then
					w = w + f(ix,iy,hex)
					z = z + f(ix,iy,hex) * subgrid(ix,iy,hex,:)
				endif
			enddo
		enddo
	enddo

	if (w > 0.5d0) then
		ww3(i,:) = z/w
	else
		ww3(i,:) = nan
	endif
enddo

! If requested, check for changes first

if (update) then
	call rads_get_var (S, P, 'wave_'//varnm(var0), tmp)
	do i = 1,n
		if (isnan(tmp(i)) .and. isnan(ww3(i,var0))) cycle
		if (isnan(tmp(i))) exit
		if (nint(tmp(i)/dz(var0)) .ne. nint(ww3(i,var0))) exit
	enddo
	if (i > n) then	! No changes
		write (*,552) 0
		return
	endif
endif

! Define data fields

call rads_put_history (S, P)
do var = var0,var1
	call rads_def_var (S, P, 'wave_'//varnm(var))
enddo

! Store data fields

do var = var0,var1
	call rads_put_var (S, P, 'wave_'//varnm(var), ww3(:,var)*dz(var))
enddo

write (*,552) n
end subroutine process_pass

!-----------------------------------------------------------------------
! Get WaveWatch3 data
!-----------------------------------------------------------------------

function get_ww3 (mjd, grid)
logical :: get_ww3
integer(fourbyteint), intent(in) :: mjd
integer(twobyteint), intent(out) :: grid(:,:,:,:)
integer(fourbyteint) ::	ncid, v_id, l, strf1985, var
character(len=rads_naml) :: filenm

600 format ('(',a,')')
1300 format (a,': ',a)

! Determine file name

get_ww3 = .true.

l = strf1985(filenm, path, mjd*86400)

! Open input file

write (*,600,advance='no') filenm(l-14:l)
if (nft(nf90_open(filenm,nf90_nowrite,ncid))) then
	write (*,1300) 'Error opening file',filenm(:l)
	return
endif

! Check netCDF file for all variables, long_name, scale_factor, units, and grids

do var = 1,nvar
	if (nft(nf90_inq_varid(ncid,varnm(var),v_id))) call rads_exit ('Error finding variable')
	if (nft(nf90_get_att(ncid,v_id,'scale_factor',dz(var)))) call rads_exit ('Missing attribute scale_factor')
	if (nft(nf90_get_var(ncid,v_id,grid(:,:,:,var)))) call rads_exit ('Error reading data grid')
enddo

l = nf90_close(ncid)
get_ww3 = .false.
end function get_ww3

end program rads_add_ww3
