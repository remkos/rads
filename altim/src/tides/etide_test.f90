program etide_test

use typesizes
use tides
integer :: i
real(eightbytereal) :: mjd,lat,lon,z

do
    read (*,*,iostat=i) mjd,lat,lon
    if (i/=0) exit
    do i=1,1000
        z=etide_ce(mjd,lat,lon)
    enddo
    write (*,'(4f15.6)') mjd,lat,lon,z
enddo
end program
