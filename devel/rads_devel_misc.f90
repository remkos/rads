module rads_devel_misc
use typesizes

! Struct to store the records from ORF file

type :: orfinfo
	integer(fourbyteint) :: cycle, pass
	real(eightbytereal) :: starttime, eqtime, eqlon
end type

contains

!*read_orf -- Read Orbital Revolution File (ORF)
!+
subroutine read_orf (sat, orf)
use rads_misc
character(len=3), intent(in) :: sat
type(orfinfo), intent(inout) :: orf(:)
character(len=320) :: line
integer :: hash, mjd, yy, mm, dd, hh, mn, ios, npass, orbitnr
real(eightbytereal) :: ss, lon, lat
!
! This routine reads an ORF file for the given 3-letter satellite
! abbreviation. Upon return the ORF structure will be filled.
!
! Argument:
!   sat  : 3-letter satellite abbreviation
!   orf  : structure containing the information from the ORF file
!-----------------------------------------------------------------------

! Open the equator crossing table

select case (sat)
case ('S3A')
	call parseenv ('${ALTIM}/data/ODR.SNTNL-3A/orf.txt', line)
case ('S3B')
	call parseenv ('${ALTIM}/data/ODR.SNTNL-3B/orf.txt', line)
case default
	stop 'Wrong satellite code'
end select
open (10, file=line, status='old')

! Initialise with dummy values

orf = orfinfo (-1, -1, nan, nan, nan)

! Skip until after the lines starting with #

hash = 0
do while (hash < 5)
	read (10,550,iostat=ios) line
	if (ios /= 0) stop 'Premature end of file'
	if (line(:1) == '#') hash = hash + 1
enddo

! Read the equator crossing table

npass = 1
do
	read (10,550,iostat=ios) line
	if (ios /= 0) exit
	read (line(:23),601,iostat=ios) yy,mm,dd,hh,mn,ss
	if (ios /= 0) exit
	read (line(24:),*,iostat=ios) orf(npass)%cycle,orf(npass)%pass,orbitnr,lon,lat
	if (ios /= 0) exit
	! Convert date and time to seconds since 1-1-2000
	call ymd2mjd(yy,mm,dd,mjd)
	ss = (mjd-51544)*86400d0 + hh*3600d0 + mn*60d0 + ss
	! Distinquish between rollover points (get starttime) and equator crossings (get eqtime and eqlon)
	if (abs(lat) > 1) then
		orf(npass)%starttime = ss
	else
		orf(npass)%eqtime = ss
		orf(npass)%eqlon = lon
		npass=npass + 1
	endif
enddo
close (10)
550 format (a)
601 format (i4,4(1x,i2),1x,f6.3)

end subroutine read_orf

end module rads_devel_misc
