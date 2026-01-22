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

!*rads_add_seaice -- Add OSIAF sea ice concentration to RADS data
!+
! This program adds the Ocean and Sea Ice Satellite Application
! Facility (OSI SAF) sea ice concentration to the RADS data. This field
! is only provided as reference or can be used for data screening or
! validation purposed.
!
! The OSI SAF sea ice concentration grids are produced daily. The OSI SAF grids
! for the Northern and Southern Hemispheres in a polar stereographic projection
! are read and interpolated along track.
! The OSI SAF NetCDF files are stored in ${ALTIM}/data/OSIAF_conc/[yyyy]/[mm]/
!
! usage: rads_add_seaice [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_seaice

use rads
use rads_misc
use rads_devel
use rads_netcdf
use rads_time
use netcdf

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cyc, pass
logical :: update=.false.

! Data variables

integer(fourbyteint) :: day1985old=-99999, j
real(eightbytereal), parameter :: lat_nh = 31d0, lat_sh = -39d0

type :: proj
	real(eightbytereal) :: e, akm1, lat_ts, lon_0
end type
type(proj) :: q(2)

type :: hemi
	character(len=rads_cmdl) :: path
	integer(twobyteint), allocatable :: grid(:,:,:)
	integer(fourbyteint) :: nx, ny
	real(eightbytereal) :: x0, dx, y0, dy
end type
type(hemi) :: h(2)

! Initialise

call synopsis ('--head')
call rads_set_options ('u update')
call rads_init (S)

! Check all options
do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('u', 'update')
		update = .true.
	end select
enddo

! Set the projections

call polar_stereo_init (q(1), 6378273d0, 6356889.44891d0, 70d0, -45d0) ! North Pole
call polar_stereo_init (q(2), 6378273d0, 6356889.44891d0, -70d0, 0d0)  ! South Pole

! Define the grids

h(1)%nx = 760
h(1)%ny = 1120
h(1)%x0 = -3845d3
h(1)%y0 = 5845d3

h(2)%nx = 790
h(2)%ny = 830
h(2)%x0 = -3950d3
h(2)%y0 = 4340d3

h(:)%dx = 10d3
h(:)%dy = -10d3

do j = 1,2
	allocate (h(j)%grid(h(j)%nx,h(j)%ny,2))
enddo

! Get template for path name

call parseenv ('${ALTIM}/data/OSIAF_conc/%Y/%m/ice_conc_nh_polstere-100_multi_%Y%m%d1200.nc',h(1)%path)
call parseenv ('${ALTIM}/data/OSIAF_conc/%Y/%m/ice_conc_sh_polstere-100_multi_%Y%m%d1200.nc',h(2)%path)

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
if (rads_version ('Add  ice concentration to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -u, --update              Update files only when there are changes')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, j, ix, iy, day1985
real(eightbytereal) :: time(n), lat(n), lon(n), surf(n), ice(n), tmp(n)
real(eightbytereal) :: wx, wy, wt, f(2,2,2), f2, fsum
integer(fourbyteint) :: t1, t2
integer(twobyteint) :: g(2,2,2)
real(eightbytereal), parameter :: dz = 1d-2
logical :: err

call log_pass (P)

! Get time and location

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)
call rads_get_var (S, P, 'surface_class', surf, .true.)

! Ensure that longitude is between 0 and 360

where (lon < 0d0) lon = lon + 360d0

! Process data records

do i = 1,n

! If over land, set sea ice concentration to NaN
! Else if outside of the grids, set it to 0

	if (surf(i) == 1 .or. surf(i) > 2) then
		ice(i) = nan
		cycle
	else if (lat(i) > lat_sh .and. lat(i) < lat_nh) then
		ice(i) = 0d0
		cycle
	endif

! Today and tomorrow

	t1 = floor(time(i)/86400)
	t2 = t1 + 1
	f2 = time(i)/86400d0 ! Number of days since 1985
	day1985 = int(f2)

! Load new grids when entering new day

	if (day1985 /= day1985old) then
		if (day1985 == day1985old+1) then
			! Copy the newest grids (2) to the oldest position (1)
			h(1)%grid(:,:,1) = h(1)%grid(:,:,2)
			h(2)%grid(:,:,1) = h(2)%grid(:,:,2)
			! Replace newest grids
			err = get_osiaf(day1985+1,h(1),2)
			err = get_osiaf(day1985+1,h(2),2)
		else
			! Replace both grids
			err = get_osiaf(day1985,h(1),1) .or. get_osiaf(day1985+1,h(1),2)
			err = get_osiaf(day1985,h(2),1) .or. get_osiaf(day1985+1,h(2),2)
		endif
		if (err) then
			call log_string ('Warning: No OSIAF field for current time')
		endif
		day1985old = day1985
	endif

! Compute the corresponding grid coordinates

	if (lat(i) > 0d0) then
		j = 1
	else
		j = 2
	endif

	call polar_stereo (q(j), lon(i), lat(i), wx, wy)
	wx = (wx - h(j)%x0) / h(j)%dx + 1d0
	wy = (wy - h(j)%y0) / h(j)%dy + 1d0

	ix = floor(wx)
	iy = floor(wy)

! If outside of the region covered by the grids, set sea ice concentation to 0

	if (ix < 1 .or. iy < 1 .or. ix >= h(j)%nx .or. iy >= h(j)%ny) then
		ice(i) = 0d0
		cycle
	endif

! Get the corresponding grid values

	g(:,:,:) = h(j)%grid(ix:ix+1,iy:iy+1,:)

! Set weights for bi-linear interpolation in space

	wx = wx - ix
	wy = wy - iy

	f(1,:,:) = (1d0-wx)
	f(2,:,:) = wx
	f(:,1,:) = f(:,1,:) * (1d0-wy)
	f(:,2,:) = f(:,2,:) * wy

! Set weights for linear interpolation in time

	wt = (time(i)/86400d0 - t1)/(t2 - t1)
	f(:,:,1) = f(:,:,1) * (1d0-wt)
	f(:,:,2) = f(:,:,2) * wt

! Where grid values are default (-999), set the weight to zero

	where (g == -999) f = 0d0

! If the sum of the weights if less than 0.5 assign NaN,
! otherwise the sea ice concentration becomes the weighted average.

	fsum = sum(f)
	if (fsum < 0.5d0) then
		ice(i) = nan
	else
		ice(i) = sum(f*g) / fsum * dz
	endif
enddo

! If requested, check for changes in sea ice first

if (update) then
	i = rads_verbose; rads_verbose = -1 ! Temporarily suspend warning
	call rads_get_var (S, P, 'seaice_conc_osiaf', tmp, .true.)
	rads_verbose = i
	do i = 1,n
!		write (*,'(2i4,f8.3,2f5.0)') i,nint(surf(i)),lat(i),tmp(i),ice(i)
		if (isnan_(tmp(i)) .and. isnan_(ice(i))) cycle
		if (isnan_(tmp(i))) exit
		if (nint(tmp(i)) /= nint(ice(i))) exit
	enddo
	if (i > n) then	! No changes
		call log_records (0)
		return
	endif
endif

! Store all data fields

call rads_put_history (S, P)

call rads_def_var (S, P, 'seaice_conc_osiaf')

call rads_put_var (S, P, 'seaice_conc_osiaf', ice)

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Polar sterographic projection
!-----------------------------------------------------------------------

subroutine polar_stereo_init (q, a, b, lat_ts, lon_0)
type(proj), intent(out) :: q
real(eightbytereal), intent(in) :: a, b, lat_ts, lon_0
real(eightbytereal) :: e, sp, mf, tf, k0

e = sqrt(1d0-(b/a)**2d0)
sp = abs(lat_ts)*rad
mf = cos(sp)/sqrt(1d0-(e*cos(sp))**2d0)
tf = tan(pi/4d0-sp/2d0)*(((1d0+e*sin(sp))/(1d0-e*sin(sp)))**(e/2d0))
k0 = mf/(2d0*tf)*sqrt((1d0+e)**(1d0+e)*(1d0-e)**(1d0-e))
k0 = 0.96985819d0 ! Added to force the right answers, aligned with PROJ9

q%akm1 = (2d0*a*k0)/sqrt((1d0+e)**(1d0+e)*(1d0-e)**(1d0-e))

q%e = e
q%lon_0 = lon_0
q%lat_ts = lat_ts
end subroutine polar_stereo_init

subroutine polar_stereo (q, lon, lat, x, y)
type(proj), intent(in) :: q
real(eightbytereal), intent(in) :: lon, lat
real(eightbytereal), intent(out) :: x, y
real(eightbytereal) :: phi, lam, coslam, sinphi, t

phi = lat*rad
lam = (lon-q%lon_0)*rad
coslam = cos(lam)
sinphi = sin(phi)

if (q%lat_ts < 0d0) then
	phi = -phi
	sinphi = -sinphi
	coslam = -coslam
endif

t = q%akm1 * tan(pi/4d0-phi/2d0)*((1d0+q%e*sinphi)/(1d0-q%e*sinphi))**(q%e/2d0)
x = t * sin(lam)
y = -t * coslam
end subroutine polar_stereo

!-----------------------------------------------------------------------
! Get the OSIAF sea ice concentration grid for "day1985"
!-----------------------------------------------------------------------

logical function get_osiaf (day1985, h, it)
integer(fourbyteint), intent(in) :: day1985
type(hemi), intent(inout) :: h
integer(fourbyteint), intent(in) :: it
character(len=rads_cmdl) :: filenm
integer(fourbyteint) ::	ncid,v_id,t_id,j,l,strf1985
real(eightbytereal) :: time
real(eightbytereal) :: z0, dz

600 format ('(',a,')')
1300 format (a,': ',a)

get_osiaf = .true.

! Determine file name

l = strf1985(filenm, trim(h%path), day1985*86400)

! Open input file

write (*,600,advance='no') trim(filenm)
if (nft(nf90_open(filenm,nf90_nowrite,ncid))) then
	write (*,1300) 'Error opening file',filenm(:l)
	h%grid(:,:,it) = -999
	return
endif

! Check if NetCDF file contains variable name time

if (nft(nf90_inq_varid(ncid,'time',t_id))) call fin('Error finding variable')

! Get time in 1980 seconds and check against input

if (nft(nf90_get_var(ncid,t_id,time))) call fin('Error reading time')
if (nint((time-220968000)/86400) /= day1985) call fin('Day does not match grid')

! Check if NetCDF file contains variable name sea_ice_fraction

if (nft(nf90_inq_varid(ncid,'ice_conc',v_id))) call fin('Error finding variable')

! Get scale factor, offset and missing value

if (nft(nf90_get_att(ncid,v_id,'add_offset',z0))) z0=0
if (nft(nf90_get_att(ncid,v_id,'scale_factor',dz))) dz=1

if (nft(nf90_get_var(ncid,v_id,h%grid(:,:,it),count = (/ h%nx, h%ny, 1 /)))) call fin('Error reading data grid')

j = nf90_close(ncid)
get_osiaf = .false.
end function get_osiaf

end program rads_add_seaice
