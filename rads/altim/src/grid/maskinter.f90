program maskinter
use typesizes
use masklib

! This program interpolates values in a mask

character(80) :: masknm
integer(fourbyteint) :: ios,z
real(eightbytereal) :: lat,lon
type(mask) :: info

! Initialize

call getarg(1,masknm)

if (mask_load(masknm,info) /= 0) call fin('Error reading mask file')

do
    read (*,*,iostat=ios) lat,lon
    if (ios /= 0) exit
    z = mask_query(info,lon,lat)
    write (*,'(f8.3,f9.3,i2)') lat,lon,z
enddo

call mask_free(info)
end program maskinter
