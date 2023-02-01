program basinstat

! Do a little statistics depending on the basin
!
! Syntax: basinstat gridfile
!
!-
! $Id: basinstat.f90,v 1.1 2009/12/15 05:53:52 rads Exp $
!-----------------------------------------------------------------------
use typesizes
use gridlib

type(grid) :: basin,info
integer(fourbyteint), parameter :: nbasin = 83
integer(fourbyteint) :: kx,ky,nx,i,l,ios,class(-3:nbasin)=0
real(eightbytereal), parameter :: rad = atan(1d0) / 45d0
real(eightbytereal) :: wgt(-3:nbasin)=0d0,mean(-3:nbasin)=0d0,rms(-3:nbasin)=0d0,x,y,w,val
character(80) :: gridnm,basinm(-3:nbasin)

! Set the basin names

basinm(-3) = '= Global'
basinm(-2) = '= Ocean'
basinm(-1) = '= Lakes'
basinm( 0) = 'Land'
basinm( 1) = 'Pacific Ocean'
basinm( 2) = 'Atlantic Ocean'
basinm( 3) = 'Indian Ocean'
basinm( 4) = 'Arctic Ocean'
basinm(10) = 'Bering Sea'
basinm(11) = 'Sea of Okhotsk'
basinm(12) = 'Sea of Japan'
basinm(13) = 'Yellow Sea'
basinm(14) = 'South China Sea'
basinm(15) = 'Indonesian Sea'
basinm(20) = 'Hudson Bay'
basinm(21) = 'Gulf of Mexico'
basinm(22) = 'Caribbean Sea'
basinm(23) = 'North Sea'
basinm(24) = 'Baltic Sea'
basinm(31) = 'Arabian Sea'
basinm(32) = 'Bay of Bengal'
basinm(33) = 'Andaman Sea'
basinm(34) = 'Persian Gulf'
basinm(35) = 'Red Sea'
basinm(41) = 'Great Slave Lake'
basinm(42) = 'Lake Winipeg'
basinm(43) = 'Lake Superior'
basinm(44) = 'Lake Michigan'
basinm(45) = 'Lake Huron'
basinm(46) = 'Lake Ontario'
basinm(47) = 'Lake Erie'
basinm(50) = 'Lake Titicaca'
basinm(60) = 'Mediterranan Sea'
basinm(61) = 'Adriatic Sea'
basinm(70) = 'Black Sea'
basinm(71) = 'Caspian Sea'
basinm(72) = 'Aral Sea'
basinm(73) = 'Lake Baikal'
basinm(74) = 'Lake Balkhash'
basinm(78) = 'Lake Ladoga'
basinm(80) = 'Lake Chad'
basinm(81) = 'Lake Malawi'
basinm(82) = 'Lake Tanganyika'
basinm(83) = 'Lake Victoria'
class(1:35) = -2    ! Ocean
class(41:50) = -1   ! Lake
class(60:70) = -2   ! Ocean
class(71:83) = -1   ! Lake

! Get the grid file name

call getarg (1,gridnm)
if (gridnm == ' ') then
    write (*,600)
    stop
endif
600 format ('basinstat - Compute grid statistics per basin'//'basinstat gridname')

! Load the data grid

ios = grid_load (gridnm,info)
if (ios /= 0) call fin('Error reading data grid')

! Load the basin grid

call checkenv ('ALTIM',gridnm,l)
gridnm(l+1:)='/data/basin_codes.nc'
ios = grid_load (gridnm,basin)
if (ios /= 0) call fin('Error reading basin grid')

! If the data grid is more than global, use only a subsection

if (abs(info%xmax - info%xmin) >= 360d0) then
    nx = nint(360d0 / info%dx)
else
    nx = info%nx
endif

! Now do the stats for each basin code

do ky = 1,info%ny
    y = grid_y(info,ky)
    w = cos(y * rad)
    do kx = 1,nx
	x = grid_x(info,kx)
	val = grid_query(info,x,y)
	if (isnan(val)) cycle
	if (x < basin%xmin) x = x + 360d0
	if (x > basin%xmax) x = x - 360d0
	i = nint(grid_query(basin,x,y))
	mean(i) = mean(i) + w * val
	rms(i) = rms(i) + w * val * val
	wgt(i) = wgt(i) + w
    enddo
enddo

! Sort in classes

do i = 0,nbasin
    l = class(i)
    if (l /= 0) then
	mean(l) = mean(l) + mean(i)
	rms(l) = rms(l) + rms(i)
	wgt(l) = wgt(l) + wgt(i)
    endif
    l = -3
    mean(l) = mean(l) + mean(i)
    rms(l) = rms(l) + rms(i)
    wgt(l) = wgt(l) + wgt(i)
enddo

! Output the stats

write (*,610)
mean = mean / wgt
rms = sqrt(max(0d0,rms / wgt - mean * mean))
do i = -3,nbasin
    if (wgt(i) > 0d0) write (*,620) basinm(i),mean(i),rms(i),wgt(i)/wgt(-3)*1d2
enddo
610 format ('# basin                     mean         rms  % area')
620 format (a20,2f12.6,f8.2)

end program basinstat
