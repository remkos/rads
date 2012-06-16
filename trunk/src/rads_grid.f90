!-----------------------------------------------------------------------
! $Id$
!
! Copyright (C) 2012  Remko Scharroo (Altimetrics LLC)
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

!***********************************************************************
!*rads_grid -- A set of grid reading and interpolation routines
!+
module rads_grid
use typesizes
type grid
	character(128) :: filenm
	integer(fourbyteint) :: nx, ny, ntype, nxwrap
	real(eightbytereal) :: xmin, xmax, dx, ymin, ymax, dy, zmin, zmax, dz, z0, znan
	integer(onebyteint), allocatable :: grid_int1(:,:)
	integer(twobyteint), allocatable :: grid_int2(:,:)
	integer(fourbyteint), allocatable :: grid_int4(:,:)
	real(fourbytereal), allocatable :: grid_real(:,:)
	real(eightbytereal), allocatable :: grid_dble(:,:)
end type
!
! Use the following routines to load grids into memory and interpolate
! or query them:
! grid_load -- Load a grid into memory (NetCDF format)
! grid_free -- Free grid buffer
! grid_query -- Look-up value in buffered grid
! grid_lininter -- Bi-linear interpolation of buffered grid
! grid_splinter -- Bi-cubic spline interpolation of buffered grid
!
! Other practical functions are:
! grid_x -- Get x-coordinate of grid pixel
! grid_y -- Get y-coordinate of grid pixel
!
! Information about the grid is stored in a structure of type 'grid'.
! Use the module 'rads_grid' to define the structure.
!-----------------------------------------------------------------------
integer, parameter, private :: stderr = 0
private :: make_nan

contains

!***********************************************************************
!*grid_load -- Load a grid into memory (NetCDF format)
!+
function grid_load (filenm, info)
use netcdf
integer(fourbyteint) :: grid_load
character(len=*), intent(in) :: filenm
type(grid), intent(inout) :: info
!
! This routine allocates memory and loads the contents of a netCDF
! grid file into the allocated memory. To save memory, the grid is stored
! in its original representation, i.e., 2- or 4-byte integer or
! 4- or 8-byte floats.
!
! The argument filenm is the path name of the grid file. Optionally, one
! can append '?varname' to load a specific variable ('varname') from a netCDF
! grid. If not given, the first 2-dimensional grid is loaded.
!
! All information about the grid is stored in the structure info.
!
! After using <grid_load>, use <grid_query>, <grid_lininter> or
! <grid_splinter> to query the grid or do a bi-linear or bi-cubic spline
! interpolation.
!
! Use <grid_free> to free the allocated memory. Note that this does
! not free the grid info structure.
!
! When the allocation of memory or the loading of the grid was
! unsuccessful, this will be reflected in the returned function value.
!
! Input argument:
!  filenm   : Name of the file containing the grid.
!             Optionally append ?varname to indicate that the variable
!             varname must be read.
!
! Output argument:
!  info     : Structure containing information about the grid. Needs
!             to be declared or allocated in calling program.
!
! Return value:
!  grid_load: 0 = No error
!             1 = Could not find or open grid
!             2 = Variable not found or not 2-dimensional
!             3 = Illegal grid format
!             4 = Error loading grid
!-----------------------------------------------------------------------
integer(fourbyteint) :: i,ncid,z_id
character(80) :: units

call grid_load_nc
if (grid_load /= 0) then
	if (grid_load == -1) call grid_error (3, 'Error loading grid attributes from '//filenm)
	i = nf90_close (ncid)
	call grid_free (info)
	return
endif

! Load entire data block.
select case (info%ntype)
case (nf90_int1)
	allocate (info%grid_int1(info%nx,info%ny))
	i = nf90_get_var (ncid,z_id,info%grid_int1)
case (nf90_int2)
	allocate (info%grid_int2(info%nx,info%ny))
	i = nf90_get_var (ncid,z_id,info%grid_int2)
case (nf90_int4)
	allocate (info%grid_int4(info%nx,info%ny))
	i = nf90_get_var (ncid,z_id,info%grid_int4)
case (nf90_real)
	allocate (info%grid_real(info%nx,info%ny))
	i = nf90_get_var (ncid,z_id,info%grid_real)
case default
	allocate (info%grid_dble(info%nx,info%ny))
	i = nf90_get_var (ncid,z_id,info%grid_dble)
end select

if (i /= 0) then
	call grid_error (4, 'Error loading data from '//trim(filenm))
	i = nf90_close (ncid)
	call grid_free (info)
	return
endif

! Successful ending
grid_load = 0
i = nf90_close (ncid)
return

contains

subroutine grid_load_nc
integer(fourbyteint) :: i,l,noff,x_id,y_id,nvars,dims(2),ndims,dims2(1),start(1)=1,stride(1)=1
real(eightbytereal) :: dummy(2),x

! Check if a suffix '?varname' is used
l = index(filenm,'?')
grid_load = -1
info%filenm = filenm

if (l > 1) then  ! Varname is given

	if (nf90_open(filenm(:l-1),nf90_nowrite,ncid) /= 0) then
		call grid_error (1, 'No such netCDF file '//filenm(:l-1))
		return
	endif
	if (nf90_inq_varid(ncid,filenm(l+1:),z_id) /= 0) then
		call grid_error (2, 'No such variable '//trim(filenm(l+1:))//' in '//filenm(:l-1))
		return
	endif
	if (nf90_inquire_variable(ncid,z_id,ndims=ndims,xtype=info%ntype) /= 0) return
	if (ndims /= 2) then
		call grid_error (3, 'Variable is not 2-D in '//trim(filenm))
		return
	endif
	if (nf90_inquire(ncid,nvariables=nvars) /= 0) return

else ! Open grid file.

	if (nf90_open(filenm,nf90_nowrite,ncid) /= 0) then
		call grid_error (1, 'No such netCDF file '//trim(filenm))
		return
	endif

	! Look for first 2-dimensional (z) variable and determine
	if (nf90_inquire(ncid,nvariables=nvars) /= 0) return
	z_id = -1
	do i = 1,nvars
		if (nf90_inquire_variable(ncid,i,ndims=ndims,xtype=info%ntype) /= 0) return
		if (ndims == 2) then
			z_id = i
			exit
		endif
	enddo
	if (z_id < 0) then
		call grid_error (3, 'Could not find 2-D variable in '//trim(filenm))
		return
	endif
endif

! Get the data type
select case (info%ntype)
case (nf90_int1,nf90_int2,nf90_int4,nf90_real,nf90_double)
case default
	call grid_error (3, 'Unknown data type in '//trim(filenm))
	return
end select

! Get the ids of the x and y variables
if (nf90_inquire_variable(ncid,z_id,dimids=dims) /= 0) return
x_id = -1
y_id = -1
do i = 1,nvars
	if (nf90_inquire_variable(ncid,i,ndims=ndims) /= 0) return
	if (ndims == 1) then
		if (nf90_inquire_variable(ncid,i,dimids=dims2) /= 0) return
		if (dims2(1) == dims(1)) x_id=i
		if (dims2(1) == dims(2)) y_id=i
	endif
enddo
if (x_id < 0 .or. y_id < 0) then
	call grid_error (3, 'Could not find the x or y variables in '//trim(filenm))
	return
endif

! Fill grid structure
if (nf90_inquire_dimension(ncid,dims(1),len=info%nx) /= 0) return
if (nf90_inquire_dimension(ncid,dims(2),len=info%ny) /= 0) return

! Get z-range, -scale and -offset and missing value
if (nf90_get_att(ncid,z_id,'scale_factor',info%dz) /= 0) info%dz = 1d0
if (nf90_get_att(ncid,z_id,'add_offset'  ,info%z0) /= 0) info%z0 = 0d0

if (nf90_get_att(ncid,z_id,'actual_range',dummy) /= 0 .and. nf90_get_att(ncid,z_id,'valid_range' ,dummy) /= 0) then
	dummy(1) = -1d40
	dummy(2) = 1d40
endif
info%zmin = dummy(1)
info%zmax = dummy(2)

if (nf90_get_att(ncid,z_id,'_FillValue',info%znan) /= 0 .and. nf90_get_att(ncid,z_id,'missing_value',info%znan) /= 0) &
	info%znan = make_nan()

! Get x- and y-ranges
if (nf90_get_att(ncid,nf90_global,'node_offset',noff) /= 0) noff = 0
if (nf90_get_att(ncid,x_id,'actual_range',dummy) /= 0 .and. nf90_get_att(ncid,x_id,'valid_range',dummy) /= 0) then
	stride(1) = info%nx-1
	if (nf90_get_var(ncid,x_id,dummy,start=start,stride=stride) /= 0) return
	noff = 0
endif
info%dx = (dummy(2)-dummy(1))/(info%nx+noff-1)
info%xmin = dummy(1)+noff*info%dx/2
info%xmax = dummy(2)-noff*info%dx/2

if (nf90_get_att(ncid,nf90_global,'node_offset',noff) /= 0) noff = 0
if (nf90_get_att(ncid,y_id,'actual_range',dummy) /= 0 .and. nf90_get_att(ncid,y_id,'valid_range' ,dummy) /= 0) then
	stride(1)=info%ny-1
	if (nf90_get_var(ncid,y_id,dummy,start=start,stride=stride) /= 0) return
	noff = 0
endif
info%dy = (dummy(2)-dummy(1))/(info%ny+noff-1)
info%ymin = dummy(1)+noff*info%dy/2
info%ymax = dummy(2)-noff*info%dy/2

! Check if
! - grid is geographical
! - cellwidth is an integer fraction of 360
! - number of cells times cellwidth spans at least 360 degrees
! Determine how many cells needed to wrap
if (nf90_get_att(ncid,x_id,'units',units) == 0 .and. units == 'degrees_east') then
	x = 360d0/info%dx
	info%nxwrap = nint(x)
	if (info%nxwrap > info%nx .or. abs(x - info%nxwrap) > 0.01d0) info%nxwrap = 0
else
	info%nxwrap = 0
endif

grid_load = 0
end subroutine grid_load_nc

subroutine grid_error (err, string)
integer, intent(in) :: err
character(len=*), intent(in) :: string
character(len=80) :: prog
grid_load = err
call getarg (0,prog)
write (stderr, '(a,": ",a)') trim(prog),trim(string)
end subroutine grid_error

end function grid_load

!***********************************************************************
!*grid_free -- Free grid buffer
!+
elemental subroutine grid_free (info)
type (grid), intent(inout) :: info
!
! This routine frees up the memory allocated by <grid_load> to store a grid.
! Note that this routine does not deallocate the grid structure info
! itself, only the memory allocated to store the grid values.
!
! Input argument:
!  info : Grid info structure as returned by grid_load
!-----------------------------------------------------------------------
if (allocated(info%grid_int1)) deallocate(info%grid_int1)
if (allocated(info%grid_int2)) deallocate(info%grid_int2)
if (allocated(info%grid_int4)) deallocate(info%grid_int4)
if (allocated(info%grid_real)) deallocate(info%grid_real)
if (allocated(info%grid_dble)) deallocate(info%grid_dble)
info%ntype = 0
end subroutine grid_free

!***********************************************************************
!*grid_query -- Look-up value in buffered grid
!+
pure function grid_query (info, x, y)
use netcdf
type(grid), intent(in) :: info
real(eightbytereal), intent(in) :: x, y
real(eightbytereal) :: grid_query
!
! This function looks up a single value in a buffered grid that
! was previously loaded using <grid_load>. No interpolation is
! performed, the value at the closest grid point is returned.
!
! The location at which the grid is to be interpolated is given by
! <x> and <y>. When the grid is geographical, the routine will attemp
! to wrap <x> by a multiple of 360 degrees so that it is within the grid.
!
! Upon exit, the function value grid_query will be the value
! of the grid at a grid node closest to (<x>, <y>). When that point is
! undetermined or when (<x>, <y>) is outside the grid, even after
! wrapping, grid_query returns a NaN value.
!
! Input arguments:
!  info : Grid info structure as returned by grid_load
!  x, y : x- and y-coordinate of the point to be queried
!
! Output argument:
!  grid_query : Value at the location (x, y)
!-----------------------------------------------------------------------
real(eightbytereal) :: z
integer(fourbyteint) :: jx,jy

! If x or y are NaN or beyond allowed range, return NaN
if (.not.grid_inside (info, x, y)) then
	grid_query = make_nan()
	return
endif

! Determine jx (1->nx) of the closest point
jx = nint((x - info%xmin) / info%dx)
if (info%nxwrap == 0) then
	jx = jx + 1
else
	jx = modulo (jx, info%nxwrap) + 1
endif

! Determine jy (1->ny) of the closest point
jy = nint((y - info%ymin) / info%dy) + 1

! Lookup the value
select case (info%ntype)
case (nf90_int1)  ! BYTE
	z = info%grid_int1(jx,jy)
case (nf90_int2)  ! SHORT
	z = info%grid_int2(jx,jy)
case (nf90_int4)  ! INT
	z = info%grid_int4(jx,jy)
case (nf90_real)  ! FLOAT
	z = info%grid_real(jx,jy)
case default      ! DOUBLE
	z = info%grid_dble(jx,jy)
end select

! Check against missing value and scale
if (z == info%znan) then
	z = make_nan()
else if (isnan(z)) then
else
	z = z * info%dz + info%z0
endif
grid_query = z
end function grid_query

!***********************************************************************
!*grid_lininter -- Bi-linear interpolation of buffered grid
!+
pure function grid_lininter (info, x, y)
use netcdf
type(grid), intent(in) :: info
real(eightbytereal), intent(in) :: x, y
real(eightbytereal) :: grid_lininter
!
! This function interpolates a buffered grid that was previously loaded
! using <grid_load>. Bi-linear interpolation is used whenever possible.
!
! The location at which the grid is to be interpolated is given by
! <x> and <y>. When the grid is geographical, the routine will attemp
! to wrap <x> by a multiple of 360 degrees so that it is within the grid.
!
! Upon exit, the function value grid_lininter will be the interpolated
! value of the grid at the location (<x>, <y>). When (<x>, <y>) points
! directly to a grid point, the value at this grid point is returned,
! otherwise bi-linear interpolation is performed between the four grid points
! surrounding (<x>, <y>). When one of the grid points is not determined, a
! triangular interpolation is conducted.
!
! When (<x>, <y>) is outside the grid, even after wrapping, or when
! (<x>, <y>) is close to an undetermined grid point, <grid_lininter>
! returns a NaN value.
!
! Input arguments:
!  info : Grid info structure as returned by grid_load
!  x, y : x- and y-coordinate of the point to be interpolated
!
! Output argument:
!  grid_lininter : Interpolated value at the location (x, y)
!-----------------------------------------------------------------------
real(eightbytereal) :: xj,yj,z(2,2),weight(2,2),wtot,vtot
integer(fourbyteint) :: jx,jy,jx1,jy1

! If x or y are NaN or beyond allowed range, return NaN
if (.not.grid_inside (info, x, y)) then
	grid_lininter = make_nan()
	return
endif

! Determine jx,jx1 (1->nx) of the corners
xj = (x - info%xmin) / info%dx
if (info%nxwrap == 0) then
	jx = min(int(xj),info%nx-2)	! Use int() so we don't drop below zero
	xj = xj - jx
	jx = jx + 1
	jx1 = jx + 1
else
	jx = floor(xj)
	xj = xj - jx
	jx = modulo (jx, info%nxwrap) + 1
	jx1 = modulo (jx, info%nxwrap) + 1
endif

! Determine jy,jy1 (1->ny) of the corners
yj = (y - info%ymin) / info%dy
jy = min(int(yj),info%ny-2)
yj = yj - jy
jy = jy + 1
jy1 = jy + 1

! Lookup 4 corners
select case (info%ntype)
case (nf90_int1)  ! BYTE
	z(1,:) = info%grid_int1(jx ,jy:jy1)
	z(2,:) = info%grid_int1(jx1,jy:jy1)
case (nf90_int2)  ! SHORT
	z(1,:) = info%grid_int2(jx ,jy:jy1)
	z(2,:) = info%grid_int2(jx1,jy:jy1)
case (nf90_int4)  ! INT
	z(1,:) = info%grid_int4(jx ,jy:jy1)
	z(2,:) = info%grid_int4(jx1,jy:jy1)
case (nf90_real)  ! FLOAT
	z(1,:) = info%grid_real(jx ,jy:jy1)
	z(2,:) = info%grid_real(jx1,jy:jy1)
case default      ! DOUBLE
	z(1,:) = info%grid_dble(jx ,jy:jy1)
	z(2,:) = info%grid_dble(jx1,jy:jy1)
end select

! Set corner weights
weight(1,1) = (1-xj)*(1-yj)
weight(1,2) = (1-xj)*   yj
weight(2,1) =    xj *(1-yj)
weight(2,2) =    xj *   yj

! Add up weights
wtot = 0d0
vtot = 0d0
do jy = 1,2
	do jx = 1,2
		if (z(jx,jy) == info%znan .or. isnan(z(jx,jy))) then
		else
			wtot = wtot+weight(jx,jy)
			vtot = vtot+weight(jx,jy)*z(jx,jy)
		endif
	enddo
enddo
grid_lininter = vtot / wtot * info%dz + info%z0
end function grid_lininter

!***********************************************************************
!*grid_splinter -- Bi-cubic spline interpolation of buffered grid
!+
pure function grid_splinter (info, x, y)
use netcdf
type(grid), intent(in) :: info
real(eightbytereal), intent(in) :: x, y
real(eightbytereal) :: grid_splinter
!
! This function interpolates a buffered grid that was previously loaded
! using <grid_load>. Bi-cubic spline interpolation is used whenever possible.
!
! The location at which the grid is to be interpolated is given by
! <x> and <y>. When the grid is geographical, the routine will attemp
! to wrap <x> by a multiple of 360 degrees so that it is within the grid.
!
! Upon exit, the function value grid_splinter will be the interpolated
! value of the grid at the location (<x>, <y>). When (<x>, <y>) points
! directly to a grid point, the value at this grid point is returned,
! otherwise bi-cubic spline interpolation is performed between the
! 6-by-6 grid points surrounding (<x>, <y>).
!
! This routine DOES NOT take into account invalid values. Hence, it
! assumes that all values in the 6-by-6 subgrid are valid. This will
! work for most geoid and mean sea surface grids.
!
! When (<x>, <y>) is outside the grid, or too close to the boundary in order
! to perform bi-cubic spline interpolation, even after wrapping, or when
! (<x>, <y>) is close to an undetermined grid point, grid_splinter returns
! a NaN value.
!
! Input arguments:
!  info : Grid info structure as returned by grid_load
!  x, y : x- and y-coordinate of the point to be interpolated
!
! Output argument:
!  grid_splinter : Interpolated value at the location (x, y)
!-----------------------------------------------------------------------
integer(fourbyteint) :: jx,jy,j
real(eightbytereal) :: xj,yj,z(6,6),zy(6),w(6),u(6)

! If x or y are NaN or beyond allowed range, return NaN
if (.not.grid_inside (info, x, y)) then
	grid_splinter = make_nan()
	return
endif

! Determine coordinates of the lower left corner of the spline window (0->nx-1, 0->ny-1)
xj = (x - info%xmin) / info%dx
jx = int(xj)
xj = xj - jx
jx = jx - 2
yj = (y - info%ymin) / info%dy
jy = int(yj)
yj = yj - jy
jy = jy - 2

! Check if spline window is completely within y-range
if (jy < 0 .or. jy+5 >= info%ny) then
	grid_splinter = make_nan()
	return
endif

if (info%nxwrap == 0) then
	! Check if spline window is completely within x-range
	if (jx < 0 .or. jx+5 >= info%nx) then
		grid_splinter = make_nan()
		return
	endif

	! Load 6-by-6 subgrid
	select case (info%ntype)
	case (nf90_int1)  ! BYTE
		z = info%grid_int1(jx+1:jx+6,jy+1:jy+6)
	case (nf90_int2)  ! SHORT
		z = info%grid_int2(jx+1:jx+6,jy+1:jy+6)
	case (nf90_int4)  ! INT
		z = info%grid_int4(jx+1:jx+6,jy+1:jy+6)
	case (nf90_real)  ! FLOAT
		z = info%grid_real(jx+1:jx+6,jy+1:jy+6)
	case default      ! DOUBLE
		z = info%grid_dble(jx+1:jx+6,jy+1:jy+6)
	end select
else
	! Load 6-by-6 subgrid
	select case (info%ntype)
	case (nf90_int1)  ! BYTE
		forall (j=1:6) z(j,:) = info%grid_int1(modulo(jx+j-1,info%nxwrap)+1,jy+1:jy+6)
	case (nf90_int2)  ! SHORT
		forall (j=1:6) z(j,:) = info%grid_int2(modulo(jx+j-1,info%nxwrap)+1,jy+1:jy+6)
	case (nf90_int4)  ! INT
		forall (j=1:6) z(j,:) = info%grid_int4(modulo(jx+j-1,info%nxwrap)+1,jy+1:jy+6)
	case (nf90_real)  ! FLOAT
		forall (j=1:6) z(j,:) = info%grid_real(modulo(jx+j-1,info%nxwrap)+1,jy+1:jy+6)
	case default      ! DOUBLE
		forall (j=1:6) z(j,:) = info%grid_dble(modulo(jx+j-1,info%nxwrap)+1,jy+1:jy+6)
	end select
endif

! Conduct spline interpolation in horizontal direction
do jy = 1,6
	call grid_splintert(z(1,jy),6,w,u)
	zy(jy) = grid_splinteru(z(3,jy),z(4,jy),w(3),w(4),xj)
enddo

! Conduct spline interpolation in vertical direction
call grid_splintert(zy,6,w,u)
grid_splinter = grid_splinteru(zy(3),zy(4),w(3),w(4),yj) * info%dz + info%z0

contains

pure subroutine grid_splintert(y,n,w,u)
! Based on 'spline' from Numerical Recipes
integer(fourbyteint), intent(in) :: n
real(eightbytereal), intent(in) :: y(n)
real(eightbytereal), intent(out) :: w(n),u(n)
real(eightbytereal) :: p
integer(fourbyteint) :: k
w(1) = 0
u(1) = 0
do k = 2,n-1
   p = w(k-1) / 2 + 2
   w(k) = -0.5d0 / p
   u(k) = (3*(y(k+1)-2*y(k)+y(k-1))-u(k-1)/2) / p
enddo
w(n-1) = u(n-1)
do k = n-2,n/2,-1
   w(k) = w(k) * w(k+1) + u(k)
enddo
end subroutine grid_splintert

pure function grid_splinteru(y0,y1,w0,w1,b)
! Based on 'splint' from Numerical Recipes
real(eightbytereal), intent(in) :: y0,y1,w0,w1,b
real(eightbytereal) :: grid_splinteru
grid_splinteru = y0 + b * (y1-y0-w0/3-w1/6+b*(w0/2+b*(w1-w0)/6))
end function grid_splinteru

end function grid_splinter

!***********************************************************************
!*grid_x -- Get x-coordinate of grid pixel
!+
pure function grid_x (info, i)
type(grid), intent(in) :: info
integer(fourbyteint), intent(in) :: i
real(eightbytereal) :: grid_x
!
! This function returns the x-coordinate of a grid pixel.
! The input i is the index along the horizontal axis, running from 1 to
! <info%nx>. If <i> is out of range, NaN is returned.
!
! Input arguments:
!  info : Grid info structure as returned by grid_load
!  i    : Horizontal index of the grid pixel
!
! Output argument:
!  grid_x : x-coordinate of the grid pixel
!-----------------------------------------------------------------------
if (i < 1 .or. i > info%nx) then
	grid_x = make_nan()
else if (i == info%nx) then
	grid_x = info%xmax
else
	grid_x = info%xmin + (i-1) * info%dx
endif
end function grid_x

!***********************************************************************
!*grid_y -- Get y-coordinate of grid pixel
!+
pure function grid_y (info, i)
type(grid), intent(in) :: info
integer(fourbyteint), intent(in) :: i
real(eightbytereal) :: grid_y
!
! This function returns the y-coordinate of a grid pixel.
! The input i is the index along the vertical axis, running from 1 to
! <info%ny>. If <i> is out of range, NaN is returned.
!
! Input arguments:
!  info : Grid info structure as returned by grid_load
!  i    : Vertical index of the grid pixel
!
! Output argument:
!  grid_y : y-coordinate of the grid pixel
!-----------------------------------------------------------------------
if (i < 1 .or. i > info%ny) then
	grid_y = make_nan()
else if (i == info%ny) then
	grid_y = info%ymax
else
	grid_y = info%ymin + (i-1) * info%dy
endif
end function grid_y

!***********************************************************************
!*grid_inside -- Check if coordinates are within grid boundaries
!+
pure function grid_inside (info, x, y)
type(grid), intent(in) :: info
real(eightbytereal), intent(in) :: x,y
logical :: grid_inside
!
! This function determines if point (<x>, <y>) is within or on the
! boundaries of a grid designated by its struct info.
! If <x> or <y> are NaN, .false. is returned.
! If <x> is periodical, no boundary check is performed on <x>.
!-----------------------------------------------------------------------
if (y >= info%ymin .and. y <= info%ymax) then
	if (info%nxwrap == 0) then
		grid_inside = (x >= info%xmin .and. x <= info%xmax)
	else
		grid_inside = .not.isnan(x)
	endif
else
	grid_inside = .false.
endif
end function grid_inside

!***********************************************************************
!*make_nan -- Create a NaN value
!+
pure function make_nan ()
real(eightbytereal) :: make_nan
!
! This function returns a double float NaN scalar
!-----------------------------------------------------------------------
make_nan = 0d0
make_nan = make_nan / make_nan
end function make_nan

end module rads_grid
