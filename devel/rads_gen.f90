!-----------------------------------------------------------------------
! Copyright (c) 2011-2026  Remko Scharroo
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

module rads_gen
use rads, only: fourbyteint, eightbytereal, rads_cmdl, rads_varl

! Command line arguments

integer(fourbyteint) :: cycles(2), nhz=0, nwvf=0, min_rec=0
real(eightbytereal) :: times(2)
character(len=rads_varl) :: sat

contains

subroutine rads_gen_getopt (satname, options)
use rads
use rads_misc
use rads_time
character(len=*), intent(in) :: satname
character(len=*), intent(in), optional :: options
character(len=rads_varl) :: optopt, optarg
integer(fourbyteint) :: ios

! Initialise

times = nan
cycles = (/0,999/)
sat = satname

! Scan command line for options

do
	if (present(options)) then
		call getopt ('vC:S:'//options//' verbose debug: sat: cycle: t: mjd: sec: ymd: doy:', optopt, optarg)
	else
		call getopt ('vC:S: verbose debug: sat: cycle: t: mjd: sec: ymd: doy:', optopt, optarg)
	endif
	select case (optopt)
	case ('!')
		exit
	case (':', '::')
		call rads_opt_error (optopt, optarg)
	case ('q', 'quiet')
		rads_verbose = -1
	case ('v', 'verbose')
		rads_verbose = rads_verbose + 1
	case ('debug')
		read (optarg, *, iostat=ios) rads_verbose
		if (ios /= 0) call rads_opt_error (optopt, optarg)
	case ('log')
		if (optarg == '-' .or. optarg == 'stdout') then
			rads_log_unit = stdout
		else if (optarg == '+' .or. optarg == 'stderr') then
			rads_log_unit = stderr
		else
			rads_log_unit = getlun()
			open (rads_log_unit, file=optarg, status='replace', iostat=ios)
			if (ios /= 0) rads_log_unit = stdout
		endif
	case ('C', 'cycle')
		cycles(2) = -1
		call read_val (optarg, cycles, '/-', iostat=ios)
		if (ios > 0) call rads_opt_error (optopt, optarg)
		if (cycles(2) < cycles(1)) cycles(2) = cycles(1)
	case ('S', 'sat')
		sat = optarg
	case ('time', 't:', 'sec', 'mjd', 'doy', 'ymd')
		call dateopt (optopt, optarg, times(1), times(2), iostat=ios)
		if (ios > 0) call rads_opt_error (optopt, optarg)
	case ('m', 'with-20hz')
		nhz = 20
	case ('w', 'with-wvf')
		nwvf = 256
		nhz = 20
	case ('min-rec')
		read (optarg, *, iostat=ios) min_rec
	end select
enddo

end subroutine rads_gen_getopt

end module rads_gen
