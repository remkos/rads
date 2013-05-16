!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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

!*rads_add_webtide -- Add WebTide ocean tide models to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! with values computed from one of the WebTide models. These models
! only provide the ocean tide.
!
! usage: rads_add_webtide [data-selectors] [options]
!
! [options] indicate the WebTide models to be used.
! See $ALTIM/data/WebTide for the various tide models that are
! available.
!-----------------------------------------------------------------------
program rads_add_webtide

use rads
use rads_misc
use rads_devel
use tides

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Command line arguments

integer(fourbyteint) :: cyc, pass
type(rads_var), pointer :: var
character(len=rads_naml) :: models, path

! Other variables

integer(fourbyteint) :: j, i0, i1, nmod = 0

! Initialise

call synopsis
call rads_set_options ('m: models:')
call rads_init (S)

! Get default coefficients from rads.xml file

var => rads_varptr (S, 'tide_ocean_webtide')
models = var%info%parameters

! Check for options

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('m', 'models')
		models = rads_opt(j)%arg
	end select
enddo

! Which models are to be used?

call getenv ('ALTIM',path)
path = trim(path) // '/data/WebTide/data/'
i1 = 0
do
	if (.not.next_word (models, i0, i1)) exit
	j = webtideinit (trim(path) // models(i0:i1-1), 's2c', nmod)
enddo

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Add WebTide ocean tide models to RADS data', flag=flag)) return
call synopsis_devel ('')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -m, --models=MODEL[,...]  Select WebTide models' / &
'                            (default models are defined in rads.xml)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n), lon(n), lat(n), z(n), tide2, lptide_eq, lptide_mf
integer(fourbyteint) :: i, j, nval

! Formats

551 format (a,' ...',i5,' records changed')

! Get time and location

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)

! Process data records

nval = 0
do i = 1,n

! Evaluate tide models.
! Take the first one that provides a non-NaN value.

	do j = 1,nmod
		if (webtide (j,time(i),lat(i),lon(i),z(i),tide2) < 0) cycle
		call lpetide (time(i),lat(i),1,lptide_eq,lptide_mf)
		z(i) = z(i) + lptide_eq
		nval = nval + 1
		exit
	enddo
enddo

! Store the values only when new ones are generated

if (nval == 0) return

! Store all data fields.

call rads_put_history (S, P)
call rads_def_var (S, P, var)
call rads_put_var (S, P, var, z)

write (*,551) trim(P%filename), nval
end subroutine process_pass

end program rads_add_webtide
