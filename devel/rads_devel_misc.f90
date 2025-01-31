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

module rads_devel_misc
use typesizes

! Struct to store the records from ORF file

type :: orfinfo
	integer(fourbyteint) :: cycle, pass, abs_pass, abs_rev
	real(eightbytereal) :: starttime, eqtime, eqlon
end type

contains

!*read_orf -- Read Orbital Revolution File (ORF)
!+
subroutine read_orf (sat, orf)
use rads_misc
character(len=*), intent(in) :: sat
type(orfinfo), intent(inout) :: orf(:)
character(len=320) :: line
integer :: mjd, yy, mm, dd, hh, mn, ios, npass, unit, nr_passes, abs_pass_offset
real(eightbytereal) :: ss, lat, lon
!
! This routine reads an ORF file for the given 3-letter satellite
! abbreviation. Upon return the ORF structure will be filled.
!
! Argument:
!   sat  : 2- or 3-letter satellite abbreviation
!   orf  : structure containing the information from the ORF file
!-----------------------------------------------------------------------

! Open the equator crossing table

nr_passes = 254
abs_pass_offset = 0
select case (sat)
case ('ER1', 'er1', 'e1')
	call parseenv ('${ALTIM}/data/ODR.ERS-1/orf.txt', line)
case ('ER2', 'er2', 'e2')
	call parseenv ('${ALTIM}/data/ODR.ERS-2/orf.txt', line)
case ('EN1', 'en1', 'n1')
	call parseenv ('${ALTIM}/data/ODR.ENVISAT1/orf.txt', line)
case ('JA1', 'ja1', 'j1')
	call parseenv ('${RADSROOT}/ext/j1/JA1_ORF.txt', line)
case ('JA2', 'ja2', 'j2')
	call parseenv ('${RADSROOT}/ext/j2/JA2_ORF.txt', line)
case ('JA3', 'ja3', 'j3')
	call parseenv ('${RADSROOT}/ext/j3/JA3_ORF.txt', line)
case ('SWT', 'swt', 'sw')
	call parseenv ('${RADSROOT}/ext/sw/SWT_ORF.txt', line)
case ('CS_', 'CS2', 'cs2', 'c2')
	call parseenv ('${ALTIM}/data/ODR.CRYOSAT2/orf.txt', line)
case ('SRL', 'srl', 'sa')
	call parseenv ('${RADSROOT}/ext/sa/SRL_ORF.txt', line)
	nr_passes = 1024
case ('S3A', 's3a', '3a')
	call parseenv ('${ALTIM}/data/ODR.SNTNL-3A/orf.txt', line)
	nr_passes = 770
	abs_pass_offset = -54
case ('S3B', 's3b', '3b')
	call parseenv ('${ALTIM}/data/ODR.SNTNL-3B/orf.txt', line)
	nr_passes = 770
case ('S3C', 's3c', '3c')
	call parseenv ('${ALTIM}/data/ODR.SNTNL-3C/orf.txt', line)
	nr_passes = 770
case ('S6A', 's6a', '6a')
	call parseenv ('${ALTIM}/data/ODR.SNTNL-6A/orf.txt', line)
case ('S6B', 's6b', '6b')
	call parseenv ('${ALTIM}/data/ODR.SNTNL-6B/orf.txt', line)
case ('SWO')
	call parseenv ('${ALTIM}/data/ODR.SWOT/orf.txt', line)
case default
	stop 'Wrong satellite code: '//sat
end select
unit = getlun()
open (unit, file=line, status='old')

! Initialise with dummy values

orf = orfinfo (-1, -1, -1, -1, nan, nan, nan)

! Skip all lines before the first line starting with a hash

do
	read (unit, 550, iostat=ios) line
	if (ios /= 0) stop 'Premature end of file'
	if (line(:1) == '#') exit
enddo

! Read the equator crossing table
! Skip any lines starting with #

npass = 1
do
	read (unit, 550, iostat=ios) line
	if (ios /= 0) exit
	if (line(:1) == '#') cycle
	read (line(:23),601,iostat=ios) yy,mm,dd,hh,mn,ss
	if (ios /= 0) exit
	read (line(24:),*,iostat=ios) orf(npass)%cycle,orf(npass)%pass,orf(npass)%abs_rev,lon,lat
	orf(npass)%abs_pass = (orf(npass)%cycle - 1) * nr_passes + orf(npass)%pass + abs_pass_offset
	if (ios /= 0) exit
	! Convert date and time to seconds since 1-1-2000
	call ymd2mjd(yy,mm,dd,mjd)
	ss = (mjd-51544)*86400d0 + hh*3600d0 + mn*60d0 + ss
	! Distinquish between rollover points (get starttime) and equator crossings (get ready for next pass)
	if (abs(lat) > 1d0) then
		orf(npass)%starttime = ss
	else
		orf(npass)%eqtime = ss
		orf(npass)%eqlon = lon
		npass=npass + 1
	endif
enddo
close (unit)
550 format (a)
601 format (i4,4(1x,i2),1x,f6.3)

end subroutine read_orf

end module rads_devel_misc
