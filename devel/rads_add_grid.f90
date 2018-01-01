!-----------------------------------------------------------------------
! Copyright (c) 2011-2018  Remko Scharroo
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

!*rads_add_grid -- Add interpolation of gridded data to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! with values computed by interpolating one or more grids.
!
! usage: rads_add_grid [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_grid

use rads
use rads_misc
use rads_grid
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P
type model
	character(rads_varl) :: xvar, yvar
	type(grid) :: info
	integer :: mode
	real(eightbytereal) :: rate
end type
type(model), allocatable :: grd(:)

! Command line arguments

character(rads_cmdl) :: path, filename
character(rads_naml), pointer :: src
integer(fourbyteint) :: i, j, k, ios, cyc, pass

! Parameters

real(eightbytereal), parameter :: t_2000 = 473299200d0, t_year = 365.25d0 * 86400d0

! Initialise

call synopsis ('--head')
call rads_init (S)

allocate (grd(S%nsel))

call parseenv ('${ALTIM}/data/', path)

! Load the selected grids

do k = 1,S%nsel
	src => S%sel(k)%info%parameters
	i = index(src, ' ')
	filename = src(:i-1)

	grd(k)%mode = rads_src_grid_lininter
	if (index(src, ' -s') > 0) grd(k)%mode = rads_src_grid_splinter
	if (index(src, ' -q') > 0) grd(k)%mode = rads_src_grid_query

	grd(k)%rate = 0d0
	i = index(src, ' -t')
	j = index(src(i+1:), ' ')
	if (i > 0) read (src(i+3:i+j), *, iostat=ios) grd(k)%rate

	grd(k)%xvar = 'lon'
	i = index(src, ' -x')
	j = index(src(i+1:), ' ')
	if (i > 0) grd(k)%xvar = src(i+3:i+j)

	grd(k)%yvar = 'lat'
	i = index(src, ' -y')
	j = index(src(i+1:), ' ')
	if (i > 0) grd(k)%yvar = src(i+3:i+j)

	call log_string ('Loading grid '//filename)
	if (filename(:1) == '/' .or. filename(:2) == './' .or. filename(:3) == '../') then
		if (grid_load(filename,grd(k)%info) /= 0) call rads_exit ('Error loading grid')
	else
		if (grid_load(trim(path)//filename,grd(k)%info) /= 0) call rads_exit ('Error loading grid')
	endif
	call log_string ('done', .true.)
enddo

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata,S%nsel)
		call rads_close_pass (S, P)
	enddo
enddo

! Free the allocated grids

do k = 1,S%nsel
	call grid_free(grd(k)%info)
enddo
deallocate (grd)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Add interpolation of gridded data to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -V, --var NAME[,...]      Select variable name(s) for interpolation (required)'// &
'All information about the grids and interpolation options are given by'/ &
'the <parameters> tags in the RADS configuration file.')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n, nmod)
integer(fourbyteint), intent(in) :: n, nmod
integer(fourbyteint) :: i, k
real(eightbytereal) :: x(n), y(n), z(n, nmod), years

call log_pass (P)

! Determine time in years since 2000

years = (P%equator_time - t_2000) / t_year

! Spline or linear interpolation of grid

do k = 1,nmod
	call rads_get_var (S, P, grd(k)%xvar, x, .true.)
	call rads_get_var (S, P, grd(k)%yvar, y, .true.)
	select case (grd(k)%mode)
	case (rads_src_grid_lininter)
		do i = 1,n
			z(i,k) = grid_lininter(grd(k)%info,x(i),y(i))
		enddo
	case (rads_src_grid_splinter)
		do i = 1,n
			z(i,k) = grid_splinter(grd(k)%info,x(i),y(i))
		enddo
	case default
		do i = 1,n
			z(i,k) = grid_query(grd(k)%info,x(i),y(i))
		enddo
	end select
	if (grd(k)%rate /= 0d0) z(:,k) = z(:,k) * grd(k)%rate * years
enddo

! Store all data fields.

call rads_put_history (S, P)
do k = 1,nmod
	call rads_def_var (S, P, S%sel(k))
enddo
do k = 1,nmod
	call rads_put_var (S, P, S%sel(k), z(:,k))
enddo

call log_records (n)
end subroutine process_pass

end program rads_add_grid
