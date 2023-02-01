program gridinter
use typesizes
use gridlib

! This program interpolates values in a grid

character(80) :: gridnm
integer(fourbyteint) :: ios
real(eightbytereal) :: lat,lon,z
logical :: linear,spline
type(grid) :: info

! Initialize

call getarg(0,gridnm)
linear=(index(gridnm,'gridinter') > 0)
spline=(index(gridnm,'gridspline') > 0)
call getarg(1,gridnm)

if (grid_load(gridnm,info) /= 0) call fin('Error reading input grid')

do
    read (*,*,iostat=ios) lat,lon
    if (ios /= 0) exit
    if (linear) then
	z = grid_lininter(info,lon,lat)
    else if (spline) then
	z = grid_splinter(info,lon,lat)
    else
	z = grid_query(info,lon,lat)
    endif
    write (*,'(f8.3,f9.3,f9.5)') lat,lon,z
enddo

call grid_free(info)
end program gridinter
