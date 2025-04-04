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

!*rads2grd -- Make grid from RADS data
!+
program rads2grd

! A Quick'n'dirty gridding program to grid a single
! RADS data variable against two others.
!-----------------------------------------------------------------------
use rads
use rads_misc
use rads_time

character(len=rads_cmdl) :: grid_name=''
character(len=rads_strl) :: format_string=''
integer(fourbyteint) :: minnr=2, n(2), k(2), i, j, l, cycle, pass, ios, nvar, nsat
real(eightbytereal) :: lo(2), hi(2), res(2), x, y, limits(3,2)
logical :: c(2)=.false.
type :: stat
	integer(fourbyteint) :: nr
	real(eightbytereal) :: mean, sum2
endtype
type(stat), allocatable :: box(:,:,:)
real(eightbytereal), allocatable :: data(:,:)
integer(fourbyteint), parameter :: msat = 20
type(rads_sat) :: S(msat)
type(rads_pass) :: P

! Initialize RADS or issue help
call synopsis
call rads_set_options ('x:y:c::o: x: y: res: min: output: grd: line-format:')
call rads_init (S)
if (any(S%error /= rads_noerr)) call rads_exit ('Fatal error')

! If no sel= was given, use default (lon,lat,sla)
! If not three items, give error
nsat = 0
do j = 1,msat
	if (S(j)%sat == '') exit
	if (S(j)%nsel == 0) then
		call rads_parse_varlist (S(j), 'lon,lat,sla')
	else if (S(j)%nsel < 3) then
		call rads_exit ('-V|--var needs at least three elements')
	endif
	nsat = nsat + 1
enddo
nvar = S(1)%nsel - 2
if (nsat == 0) call rads_exit ('Specify at least one satellite')

! Set default cell sizes
limits(1:2,1) = S(1)%sel(1)%info%limits
limits(1:2,2) = S(1)%sel(2)%info%limits
limits(3:3,:) = 1d0

! Scan command line arguments
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('x', 'x:')
		call read_val (rads_opt(i)%arg, limits(:,1), '/', iostat=ios)
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('y', 'y:')
		call read_val (rads_opt(i)%arg, limits(:,2), '/', iostat=ios)
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('res')
		limits(3,2) = nan
		call read_val (rads_opt(i)%arg, limits(3,:), '/', iostat=ios)
		if (ios > 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
		if (isnan_(limits(3,2))) limits(3,2) = limits(3,1)
	case ('c')
		if (rads_opt(i)%arg == '' .or. rads_opt(i)%arg == 'a') then
			c = .true.
		else
			c(1) = (index (rads_opt(i)%arg, 'x') > 0)
			c(2) = (index (rads_opt(i)%arg, 'y') > 0)
		endif
	case ('min')
		read (rads_opt(i)%arg, *, iostat=ios) minnr
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('o', 'output', 'grd')
		grid_name = rads_opt(i)%arg
	case ('line-format')
		format_string = rads_opt(i)%arg
	end select
enddo

! Set up the grid cells
if (any(isnan_(limits(:,1)))) call rads_exit ('You have not set the range or resolution of the x-coordinate')
if (any(isnan_(limits(:,2)))) call rads_exit ('You have not set the range or resolution of the y-coordinate')
lo  = limits(1,:)
hi  = limits(2,:)
res = limits(3,:)
do j = 1,nsat
	do i = 1,2
		S(j)%sel(i)%info%limits = (/lo(i),hi(i)/)
	enddo
enddo
where (c)
	lo = lo + 0.5d0 * res
	hi = hi - 0.5d0 * res
endwhere
n = nint((hi-lo)/res+1d0)

allocate (box(n(1),n(2),nvar+2),stat=ios)
if (ios /= 0) then
	write (format_string, '(i0,"x",i0," points")') n, n
	call rads_exit ('Unable to allocate memory for '//trim(format_string))
endif

! Initialize statistics
box = stat(0, 0d0, 0d0)
if (grid_name == '') call write_header

! Fill the grid pass by pass and keep stats, using West (1979)
do j = 1,msat
	if (S(j)%sat == '') exit
	do cycle = S(j)%cycles(1), S(j)%cycles(2), S(j)%cycles(3)
		do pass = S(j)%passes(1), S(j)%passes(2), S(j)%passes(3)
			call rads_open_pass (S(j), P, cycle, pass)
			if (P%ndata > 0) then
				allocate (data(P%ndata,nvar+2))
				do l = 1,nvar+2
					call rads_get_var (S(j), P, S(j)%sel(l), data(:,l))
				enddo
				do i = 1,P%ndata
					if (any(isnan_(data(i,:)))) cycle
					k = nint((data(i,1:2) - lo)/res)+1
					k = max(1,min(k,n))
					do l = 1,nvar
						call update_stat (box(k(1),k(2),l), data(i,l+2))
					enddo
				enddo
				deallocate (data)
			endif
			call rads_close_pass (S(j), P)
		enddo
	enddo
enddo

! Post-process grid values

where (box%nr < minnr)
	box%mean = nan
	box%sum2 = nan
elsewhere (box%nr == 1)
	box%sum2 = nan
elsewhere
	box%sum2 = sqrt(box%sum2/(box%nr - 1))
endwhere

! Write out xyz or NetCDF grid
if (grid_name == '') then
	call write_xyz_grid
else
	call write_nc_grid (nvar)
endif

! Print statistics if requested
if (rads_verbose >= 1) call rads_stat (S)

! Deallocate memory
deallocate (box)
call rads_end (S)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('Quickly grid RADS data to xyz or NetCDF grid')) return
call rads_synopsis
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -V, --var X,Y,Z[,...]     Variables for x, y and z (default: lon,lat,sla); multiple z''s may be specified' / &
'  --x X0,X1[,DX]            Set x-range and interval (default: as set by default limits and --res)'/ &
'  --y Y0,Y1[,DY]            Set y-range and interval (default: as set by default limits and --res)'/ &
'  --res DX[,DY]             Set resolution in x and y (default: 1)'/ &
'  --min MINNR               Minimum number of points per grid cell (default: 2)'/ &
'  -o, --output, --grd GRIDNAME'/ &
'                            Create NetCDF grid (suppresses ASCII)'/ &
'  --line-format FORMAT      Format to be used for ASCII output (default is determined by variables)'/ &
'  -c                        Boundaries are cell oriented'/ &
'  -c[x|y]                   Only [x|y]-boundaries are cell oriented')
stop
end subroutine synopsis

!***********************************************************************
! Update box statistics

elemental subroutine update_stat (s, x)
type(stat), intent(inout) :: s
real(eightbytereal), intent(in) :: x
real(eightbytereal) :: q, r
s%nr = s%nr + 1
q = x - s%mean
r = q / s%nr
s%mean = s%mean + r
s%sum2 = s%sum2 + r * q * (s%nr - 1)
end subroutine update_stat

!***********************************************************************
! Write out the header

subroutine write_header
integer :: j

600 format ('# Grid of RADS variables'/'# Created: ',a,' UTC: ',a)
610 format ('#'/'# Satellite : ',a,'/',a/'# Cycles    :',i5,' -',i5/'# Passes    :',i5,' -',i5/'#')
620 format ('# Output columns per grid cell:'/ &
'#    ( 1) ',a,' [',a,']'/ &
'#    ( 2) ',a,' [',a,']')
621 format ('# (',i2,'-',i2,') mean and stddev of ',a,' [',a,']')
622 format ('#    (',i2,') nr of measurements')

write (*,600) timestamp(), trim(S(1)%command)
do j = 1,msat
	if (S(j)%sat /= '') write (*,610) trim(S(j)%sat),trim(S(j)%phase%name),S(j)%cycles(1:2),S(j)%passes(1:2)
enddo
write (*,620) (trim(S(1)%sel(j)%long_name),trim(S(1)%sel(j)%info%units),j=1,2)
do j = 1,nvar
	write (*,621) j*2+1,j*2+2,trim(S(1)%sel(j+2)%long_name),trim(S(1)%sel(j+2)%info%units)
enddo
write (*,622) nvar*2+3
end subroutine write_header

!***********************************************************************
! Write out xyz grid

subroutine write_xyz_grid
integer :: kx, ky
if (format_string == '') then
	format_string = '(' // trim(S(1)%sel(1)%info%format) // ',1x,' // trim(S(1)%sel(2)%info%format)
	do j = 1,nvar
		format_string = trim(format_string) // ',2(1x,' // trim(S(1)%sel(3)%info%format) // ')'
	enddo
	format_string = trim(format_string) // ',1x,i0)'
endif
do ky = 1,n(2)
	y = lo(2) + (ky-1)*res(2)
	do kx = 1,n(1)
		x = lo(1) + (kx-1)*res(1)
		if (box(kx,ky,1)%nr >= minnr) write (*,format_string) x,y,(box(kx,ky,j)%mean,box(kx,ky,j)%sum2,j=1,nvar),box(kx,ky,1)%nr
	enddo
enddo
end subroutine write_xyz_grid

!***********************************************************************
! Write out NetCDF grid

subroutine write_nc_grid (nvar)
use netcdf
use rads_netcdf
use rads_time
integer(fourbyteint), intent(in) :: nvar
integer(fourbyteint) :: ncid,varid(nvar*2+3),dimid(2),j,k
real(fourbytereal), parameter :: nan = transfer (not(0_fourbyteint),0e0)
character(len=rads_varl) :: ext
character(len=1) :: axis(2) = (/ 'X', 'Y' /)

call nfs (nf90_create(grid_name,nf90_write+nf90_nofill,ncid))
call nfs (nf90_put_att (ncid, nf90_global, 'Conventions', 'CF-1.7'))
call nfs (nf90_put_att(ncid,nf90_global,'title',grid_name))
call nfs (nf90_put_att(ncid,nf90_global,'history',timestamp()//' UTC: '//S(1)%command))
if (all(c)) call nfs (nf90_put_att(ncid,nf90_global,'node_offset',1))

do k = 1,2
	call nf90_def_axis(ncid,S(1)%sel(k)%name,S(1)%sel(k)%long_name,S(1)%sel(k)%info%units,n(k), &
		lo(k),hi(k),dimid(k),varid(k),axis(k),S(1)%sel(k)%info%standard_name)
enddo

ext = ''
do j = 1,nvar
	ext = trim(S(1)%sel(j+2)%name) // '_'
	call nfs (nf90_def_var(ncid,trim(ext) // 'mean'  ,nf90_real,dimid,varid(j*2+1)))
	call nfs (nf90_def_var(ncid,trim(ext) // 'stddev',nf90_real,dimid,varid(j*2+2)))
	call nfs (nf90_put_att(ncid,varid(j*2+1),'long_name','mean of '//trim(S(1)%sel(j+2)%long_name)))
	call nfs (nf90_put_att(ncid,varid(j*2+2),'long_name','std dev of '//trim(S(1)%sel(j+2)%long_name)))
	do k = 1,2
		call nfs (nf90_put_att(ncid,varid(j*2+k),'units',trim(S(1)%sel(j+2)%info%units)))
		call nfs (nf90_put_att(ncid,varid(j*2+k),'_FillValue',nan))
	enddo
enddo
call nfs (nf90_def_var(ncid,'nr',nf90_int,dimid,varid(nvar*2+3)))
call nfs (nf90_put_att(ncid,varid(nvar*2+3),'long_name','number of points per cell'))
call nfs (nf90_enddef(ncid))

call nf90_put_axis(ncid,varid(1))
call nf90_put_axis(ncid,varid(2))

do j = 1,nvar
	call nfs (nf90_put_var(ncid,varid(j*2+1),box(:,:,j)%mean))
	call nfs (nf90_put_var(ncid,varid(j*2+2),box(:,:,j)%sum2))
enddo
call nfs (nf90_put_var(ncid,varid(nvar*2+3),box%nr))
call nfs (nf90_close(ncid))
end subroutine write_nc_grid

!***********************************************************************

end program rads2grd
