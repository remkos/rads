!-----------------------------------------------------------------------
! $Id$
!
! Copyright (C) 2011  Remko Scharroo (Altimetrics LLC)
! See LICENSE.TXT file for copying and redistribution conditions.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!-----------------------------------------------------------------------

!*rads2grd -- Make grid from RADS data
!+
program rads2grd

! A Quick'n'dirty gridding program to grid a single
! RADS data variable against two others.
!
! usage: rads2grd sat=<sat> sel=x,y,z [RADS_options] [options]
!-----------------------------------------------------------------------
use rads
use rads_misc

character(len=rads_naml) :: arg, grid_name='', format_string='', option, value
integer(fourbyteint) :: minnr=2, n(2), k(2), i, j, l, cycle, pass, ios
real(eightbytereal) :: lo(2), hi(2), res(2)=1d0, x, y
logical :: c(2)=.false.
type :: stat
	integer(fourbyteint) :: nr
	real(eightbytereal) :: mean, sum2
endtype
type(stat), allocatable :: box(:,:)
real(eightbytereal), allocatable :: data(:,:)
integer(fourbyteint), parameter :: msat = 5
type(rads_sat) :: S(msat)
type(rads_pass) :: P

! Initialize RADS or issue help
call synopsis
call rads_init (S)
if (any(S%error /= rads_noerr)) call rads_exit ('Fatal error')

! If no sel= was given, use default (lon,lat,sla)
! If not three items, give error
do j = 1,msat
	if (S(j)%sat == '') exit
	if (S(j)%nsel == 0) then
		call rads_parse_varlist (S(j), 'lon,lat,sla')
	else if (S(j)%nsel /= 3) then
		call rads_exit ('-V|--var= needs exactly three elements')
	endif
enddo

! Scan command line arguments
do i = 1,iargc()
	call getarg (i,arg)
	call splitarg (arg, option, value)
	select case (option)
	case ('--x','-x')
		call chartrans (value, '/', ',')
		read (value,*,iostat=ios) S(1)%sel(1)%info%limits,res(1)
	case ('--y', '-y')
		call chartrans (value, '/', ',')
		read (value,*,iostat=ios) S(1)%sel(2)%info%limits,res(2)
	case ('--res')
		call chartrans (value, '/', ',')
		read (value,*,iostat=ios) res
	case ('-c')
		c(1) = (value /= 'y')
		c(2) = (value /= 'x') 
	case ('--min')
		read (value,*) minnr
	case ('--grd')
		grid_name = value
	case ('--fmt')
		format_string = value
	end select
enddo

! Set up the grid cells
forall (j=1:2)
	lo(j) = S(1)%sel(j)%info%limits(1)
	hi(j) = S(1)%sel(j)%info%limits(2)
end forall
where (c)
	lo = lo + 0.5d0 * res
	hi = hi - 0.5d0 * res
endwhere
n = nint((hi-lo)/res+0.999d0)

allocate (box(n(1),n(2)),stat=ios)
if (ios /= 0) then
	write (*,'("rads2grd: unable to allocate memory for ",i0,"x",i0," points")') n
	stop
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
				allocate (data(P%ndata,3))
				do l = 1,3
					call rads_get_var (S(j), P, S(j)%sel(l), data(:,l))
				enddo
				do i = 1,P%ndata
					if (any(isnan(data(i,:)))) cycle
					k = nint((data(i,1:2) - lo)/res)+1
					k = max(1,min(k,n))
					call update_stat (box(k(1),k(2)), data(i,3))
				enddo
				deallocate (data)
			endif
			call rads_close_pass (S(j), P)
		enddo
	enddo
enddo

! Post-process grid values

where (box%nr < minnr)
	box%mean = S(1)%nan
	box%sum2 = S(1)%nan
elsewhere (box%nr == 1)
	box%sum2 = S(1)%nan
elsewhere
	box%sum2 = sqrt(box%sum2/(box%nr - 1))
endwhere

! Write out xyz or netCDF grid
if (grid_name == '') then
	call write_xyz_grid
else
	call write_nc_grid
endif

! Print statistics if requested
if (any(S%debug >= 1)) call rads_stat (S)

! Deallocate memory
deallocate (box)
call rads_end (S)

contains

!***********************************************************************

subroutine synopsis ()
if (rads_version ('$Revision$','Quickly grid RADS data to xyz or netCDF grid')) return
call rads_synopsis ()
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  --x=x0,x1[,dx]    : set x-range and interval (def: as set by default limits and --res=)'/ &
'  --y=y0,y1[,dy]    : set y-range and interval (def: as set by default limits and --res=)'/ &
'  --res=dx,dy       : set resolution in x and y (def:1,1)'/ &
'  --var=x,y,z       : variables for x, y and z (def:lon,lat,sla)'/ &
'  --min=minnr       : minimum number of points per grid cell (def:2)'/ &
'  --grd=gridname    : create netCDF grid (suppresses ASCII)'/ &
'  --fmt=format      : format to be used for ASCII output (default is determined by variables)'/ &
'  -c                : boundaries are cell oriented'/ &
'  -c[x|y]           : only [x|y]-boundaries are cell oriented')
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

subroutine write_header ()
integer :: j

600 format ('# Grid of RADS variables'/'#')
610 format ('# Satellite : ',a,'/',a/'# Cycles    :',i5,' -',i5/'# Passes    :',i5,' -',i5/'#')
620 format ('# Output columns per grid cell:'/ &
'#   (1) ',a,' [',a,']'/ &
'#   (2) ',a,' [',a,']'/ &
'# (3-4) mean and stddev of ',a,' [',a,']'/ &
'#   (5) nr of measurements')

write (*,600)
do j = 1,msat
	if (S(j)%sat /= '') write (*,610) trim(S(j)%sat),trim(S(j)%phase%name),S(j)%cycles(1:2),S(j)%passes(1:2)
enddo
write (*,620) (trim(S(1)%sel(j)%info%long_name),trim(S(1)%sel(j)%info%units),j=1,3)
end subroutine write_header

!***********************************************************************
! Write out xyz grid

subroutine write_xyz_grid ()
integer :: kx, ky
if (format_string == '') format_string = '(' // trim(S(1)%sel(1)%info%format) // ',1x,' // &
	trim(S(1)%sel(2)%info%format) // ',2(1x,' // trim(S(1)%sel(3)%info%format) // '),1x,i0)'
do ky = 1,n(2)
	y = lo(2) + (ky-1)*res(2)
	do kx = 1,n(1)
		x = lo(1) + (kx-1)*res(1)
		if (box(kx,ky)%nr >= minnr) write (*,format_string) x,y,box(kx,ky)%mean,box(kx,ky)%sum2,box(kx,ky)%nr
	enddo
enddo
end subroutine write_xyz_grid

!***********************************************************************
! Write out netCDF grid

subroutine write_nc_grid ()
use netcdf
use rads_netcdf
use rads_time
integer(fourbyteint) :: ncid,varid(5),dimid(2),k

call nfs (nf90_create(grid_name,nf90_write+nf90_nofill,ncid))
call nfs (nf90_put_att (ncid, nf90_global, 'Conventions', 'CF-1.5'))
call nfs (nf90_put_att(ncid,nf90_global,'title',grid_name))
call nfs (nf90_put_att(ncid,nf90_global,'history',timestamp()//': '//S(1)%command))
if (all(c)) call nfs (nf90_put_att(ncid,nf90_global,'node_offset',1))

do k = 1,2
	call nf90_def_axis(ncid,S(1)%sel(k)%name,S(1)%sel(k)%info%long_name,S(1)%sel(k)%info%units,n(k), &
		S(1)%sel(k)%info%limits(1),S(1)%sel(k)%info%limits(2),dimid(k),varid(k))
enddo

call nfs (nf90_def_var(ncid,'mean'  ,nf90_real,dimid,varid(3)))
call nfs (nf90_def_var(ncid,'stddev',nf90_real,dimid,varid(4)))
call nfs (nf90_def_var(ncid,'nr'    ,nf90_int ,dimid,varid(5)))
call nfs (nf90_put_att(ncid,varid(3),'long_name','mean of '//trim(S(1)%sel(3)%info%long_name)))
call nfs (nf90_put_att(ncid,varid(4),'long_name','std dev of '//trim(S(1)%sel(3)%info%long_name)))
call nfs (nf90_put_att(ncid,varid(5),'long_name','number of points per cell'))
do k = 3,4
   call nfs (nf90_put_att(ncid,varid(k),'units',trim(S(1)%sel(3)%info%units)))
   call nfs (nf90_put_att(ncid,varid(k),'_FillValue',real(S(1)%nan)))
enddo
call nfs (nf90_enddef(ncid))

call nf90_put_axis(ncid,varid(1))
call nf90_put_axis(ncid,varid(2))

call nfs (nf90_put_var(ncid,varid(3),box%mean))
call nfs (nf90_put_var(ncid,varid(4),box%sum2))
call nfs (nf90_put_var(ncid,varid(5),box%nr))
call nfs (nf90_close(ncid))
end subroutine write_nc_grid

!***********************************************************************

end program rads2grd
