! Driver program for WebTide subroutines
! usage: wetide_f webtide_dir < input_file
!-
program webtide_f

real(8) :: sec85, ss, tide1, tide2, lat, lon, utc
integer :: yy, dd, hh, mm, id, ios
character(80) :: arg

call getarg (1, arg)
call webtideinit (arg, "s2c", id)

do
	read (*,*,iostat=ios) lon, lat, yy, dd, hh, mm, ss
	if (ios /= 0) exit
	utc = sec85 (4, yy*1d10+01*1d8+dd*1d6+hh*1d4+mm*1d2+ss)
	call webtide (id, utc, lat, lon, tide1, tide2)

	write (*,'(f7.4,2f14.8,i5,5i3,f6.2)') tide1, lon, lat, yy, 01, dd, dd, hh, mm, ss
enddo

call webtidefree (0)
end program webtide_f
