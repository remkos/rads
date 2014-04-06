!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2014  Remko Scharroo (Altimetrics LLC)
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

!*rads_add_ssb -- Add (hybrid) SSB model to RADS data
!
! This program adds (or overwrites) the hybrid SSB model to RADS data.
! This is done by interpolating a grid in sigma0-SWH space, of which the
! filename is given in the <parameter> field in the rads.xml file.
!
! High and low SWH and sigma0 values are clipped to the extremes of the
! grid first.
!
! usage: rads_add_ssb [data-selectors]
!-----------------------------------------------------------------------
program rads_add_ssb

use rads
use rads_misc
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

character(len=rads_cmdl) :: path
integer(fourbyteint) :: cyc, pass
type(grid) :: issb_hyb
type(rads_var), pointer :: var

! Initialise

call synopsis ('--head')
var => rads_varptr (S, 'ssb_hyb')

! Load SSB model

call parseenv ('${ALTIM}/data/models/' // var%info%parameters, path)
if (grid_load(path,issb_hyb) /= 0) call rads_exit ('Error loading '//trim(path))

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

call grid_free (issb_hyb)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Add SSB model to RADS data', flag=flag)) return
call synopsis_devel ('')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: x, y, sig0(n), swh(n), ssb(n)
integer :: i

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename(len_trim(S%dataroot)+2:))

! Compute SSB

call rads_get_var (S, P, 'sig0', sig0, .true.)
call rads_get_var (S, P, 'swh', swh, .true.)
do i = 1,n
	x = sig0(i)
	if (x < issb_hyb%xmin) x = issb_hyb%xmin
	if (x > issb_hyb%xmax) x = issb_hyb%xmax
	y = swh(i)
	if (y < issb_hyb%ymin) y = issb_hyb%ymin
	if (y > issb_hyb%ymax) y = issb_hyb%ymax
	ssb(i) = grid_lininter (issb_hyb, x, y)
enddo

! Write out all the data

call rads_put_history (S, P)
call rads_def_var (S, P, 'ssb_hyb')
call rads_put_var (S, P, 'ssb_hyb', ssb)

write (*,552) n

end subroutine process_pass

end program rads_add_ssb
