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

!*rads_add_sst -- Add SST temperature and ice concentration to RADS data
!+
! This program adds the Reynolds SST and sea ice concentration
! to the RADS data. These fields are only provided as reference
! or can be used for data screening or validation purposed.
!
! The Reynolds SST grids and sea ice concentration grids are produced
! weekly and are converted first to netCDF and stored in
! $ALTIM/data/sst. Note that the grids before 1990 were centered
! on Sundays, currently on Mondays.
!
! The SST grids have 1x1 degree resolution and contain no data gaps
! (invalid values or NaNs) and are interpolated linearly.
! Output units and resolution are the same as the input (0.01 degC).
!
! The sea ice concentation grids (1x1 deg) contain NaNs (over land);
! interpolation is linear with a minimum weight of 0.25.
! Output units and resolution are the same as the input (%).
!
! usage: rads_add_sst [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_sst

use rads
use rads_misc
use rads_devel
use rads_netcdf
use netcdf

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cyc, pass
logical :: lice=.false., lsst=.false., lmean=.false.

! Data variables

character(len=rads_cmdl) :: path
integer(fourbyteint) :: t1old=0, t2old=0, j
integer(fourbyteint), parameter :: nx=360, ny=180, secweek=7*86400, sec2000=473299200
integer(fourbyteint), parameter :: sun = 157982400	! Sun 1990-01-03 12:00
integer(fourbyteint), parameter :: wed = 157723200	! Wed 1989-12-31 12:00
real(eightbytereal), parameter :: x0=0.5d0, y0=-89.5d0, dx=1d0, dy=1d0
integer(twobyteint) :: grids(0:nx+1,ny,2,2)
real(eightbytereal) :: meangrid(0:nx+1,ny)

! Initialise

call synopsis ('--head')
call rads_set_options ('ism ice sst mean all')
call rads_init (S)

! Get template for path name

call parseenv ('${ALTIM}/data/sst/oisst.%Y.nc', path)

! Check all options
do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('i', 'ice')
		lice = .true.
	case ('s', 'sst')
		lsst = .true.
	case ('m', 'mean')
		lmean = .true.
	case ('all')
		lsst = .true.
		lmean = .true.
		lice = .true.
	end select
enddo

! Load the mean SST model if requested

if (lmean) call get_mean (meangrid)

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
if (rads_version ('$Revision$', 'Add SST temperature and ice concentration to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -i, --ice                 Add sea ice concentration' / &
'  -s, --sst                 Add sea surface temperature' / &
'  -m, --mean                Add local mean sea surface temperature' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, ix, iy
real(eightbytereal) :: time(n), lat(n), lon(n), ice(n), sst(n), meansst(n)
integer(fourbyteint) :: t1, t2
real(eightbytereal) :: wx, wy, wt, f(2,2,2)

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

! What was the last Wednesday or Sunday?

	if (time(i) < sun) then  ! Before Sun 1990-01-03 12:00
		t1 = floor((time(i) - wed)/secweek) * secweek + wed
	else
		t1 = floor((time(i) - sun)/secweek) * secweek + sun
	endif

! What is the next Wednesday of Sunday?

	if (time(i) < wed) then  ! Before Wed 1989-12-31 12:00
		t2 = floor((time(i) - wed)/secweek + 1) * secweek + wed
	else
		t2 = floor((time(i) - sun)/secweek + 1) * secweek + sun
	endif

! Load new grids when needed

	if (t1 /= t1old) then
		if (get_sst(t1,grids(:,:,1,:)) /= 0) then
			write (*,551,advance='no') 'Warning: No SST for current time.'
			write (*,552) 0
			return
		endif
		t1old = t1
	endif
	if (t2 /= t2old) then
		if (get_sst(t2,grids(:,:,2,:)) /= 0) then
			write (*,551,advance='no') 'Warning: No SST for current time.'
			write (*,552) 0
			return
		endif
		t2old = t2
	endif

! Set weights for bi-linear interpolation in space

	if (lon(i) < -360d0 .or. lat(i) < -360d0) then
		meansst(i) = nan
		sst(i) = nan
		ice(i) = nan
		cycle
	endif
	if (lon(i) < 0d0) lon(i) = lon(i) + 360d0
	wx = (lon(i)-x0)/dx + 1
	wy = (lat(i)-y0)/dy + 1
	ix = floor(wx)
	iy = floor(wy)
	wx = wx - ix
	wy = wy - iy

	f(1,:,:) = (1d0-wx)
	f(2,:,:) = wx
	f(:,1,:) = f(:,1,:) * (1d0-wy)
	f(:,2,:) = f(:,2,:) * wy

! Interpolate mean SST

	meansst(i) = sum(meangrid(ix:ix+1,iy:iy+1)*f(:,:,1))

! Set weights for linear interpolation in time

	wt = (time(i) - t1)/(t2 - t1)
	f(:,:,1) = f(:,:,1) * (1d0-wt)
	f(:,:,2) = f(:,:,2) * wt

! Interpolate SST

	sst(i) = sum(grids(ix:ix+1,iy:iy+1,:,1)*f)

! Interpolate ice concentration (has invalid values)

	ice(i) = mat_product(grids(ix:ix+1,iy:iy+1,:,2),f)
enddo

! Store all data fields

call rads_put_history (S, P)

if (lice)  call rads_def_var (S, P, 'seaice_conc')
if (lsst)  call rads_def_var (S, P, 'sst')
if (lmean) call rads_def_var (S, P, 'sst_mean')

if (lice)  call rads_put_var (S, P, 'seaice_conc', ice)
if (lsst)  call rads_put_var (S, P, 'sst', sst*1d-2)
if (lmean) call rads_put_var (S, P, 'sst_mean', meansst)

write (*,552) n
end subroutine process_pass

!-----------------------------------------------------------------------
! Get the mean SST
!-----------------------------------------------------------------------

subroutine get_mean (grid)
integer(fourbyteint) :: ncid, v_id, j
real(eightbytereal) :: grid(0:nx+1,ny)

! Open grid file and load mean SST

write (*,'(a)') '(Reading mean SST)'
if (nft(nf90_open(path(:len_trim(path)-4)//'mean.nc',nf90_nowrite,ncid))) call rads_exit ('No mean SST grid file found')
if (nft(nf90_inq_varid(ncid,'sst',v_id))) call rads_exit ('Could not find mean SST grid')
if (nft(nf90_get_var(ncid,v_id,grid(1:nx,:)))) call rads_exit ('Error reading mean SST grid')

! Copy left and right boundaries

grid(0,:) = grid(nx,:)
grid(nx+1,:) = grid(1,:)

! Close grid file

j = nf90_close(ncid)
end subroutine get_mean

!-----------------------------------------------------------------------
! Get the SST and ice concentration grids for time t.
!-----------------------------------------------------------------------

function get_sst (t, grid)
integer(fourbyteint) :: t, get_sst
integer(twobyteint) :: grid(0:nx+1,ny,2)
character(len=rads_cmdl) :: filename
integer(fourbyteint) ::	ncid, v_id, j, start(3)=1, nt
real(eightbytereal) :: time(60)

! Determine file name

call strf1985(filename,path,t)

! Open input file

if (nft(nf90_open(filename,nf90_nowrite,ncid))) then
	get_sst = 1
	return
endif

! Get hours since 2000 and convert to seconds since 1985

if (nft(nf90_inquire_dimension(ncid,3,len=nt))) call rads_exit ('Error reading time')
if (nt > 60) call rads_exit ('Time dimension too large')
if (nft(nf90_inq_varid(ncid,'time',v_id))) call rads_exit ('Error reading time')
if (nft(nf90_get_var(ncid,v_id,time(1:nt)))) call rads_exit ('Error reading time')
time = time * 3600d0 + sec2000

! What is the right grid index?

start(3) = nint((t - time(1)) / secweek + 1)
if (start(3) > nt) then
	get_sst = 3
	return
endif
call strf1985(filename,'(Reading %Y-%m-%d %H:%M)',nint(time(start(3))))
write (*,'(a)',advance='no') trim(filename)

! Read the required SST grid

if (nft(nf90_inq_varid(ncid,'sst',v_id))) call rads_exit ('Error reading SST grid')
if (nft(nf90_get_var(ncid,v_id,grid(1:nx,:,1),start))) call rads_exit ('Error reading SST grid')

! Read the required ice concentration grid

if (nft(nf90_inq_varid(ncid,'ice',v_id))) call rads_exit ('Error reading ice grid')
if (nft(nf90_get_var(ncid,v_id,grid(1:nx,:,2),start))) call rads_exit ('Error reading ice grid')

! Copy left and right boundaries

grid(0,:,:) = grid(nx,:,:)
grid(nx+1,:,:) = grid(1,:,:)

! Here normal termination

get_sst = 0
j = nf90_close(ncid)
end function get_sst

!-----------------------------------------------------------------------
! Weighted product with value checking on array a
!-----------------------------------------------------------------------

function mat_product (a, b)
real(eightbytereal) :: mat_product
integer(twobyteint), intent(in) :: a(8)
real(eightbytereal), intent(in) :: b(8)
real(eightbytereal) :: w
integer(fourbyteint) :: i
mat_product = 0d0
w = 0d0
do i = 1,8
	if (a(i) <= 100) then
	w = w + b(i)
	mat_product = mat_product + a(i)*b(i)
	endif
enddo
if (w < 0.5d0) then
	mat_product = 0d0
	w = 0d0
endif
mat_product = mat_product / w
end function mat_product

end program rads_add_sst
