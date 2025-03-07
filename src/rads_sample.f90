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

program rads_sample
use rads                    ! Use rads module
type(rads_sat) :: S         ! Structure that contains all mission information
type(rads_pass) :: P        ! Structure for pass information
real(eightbytereal), allocatable :: lon(:), lat(:), sla(:)
integer :: cycle, pass, i
call rads_init (S, 'e1')    ! Initialize RADS. Select ERS-1. No mission phase needed!
do cycle = 1, 30            ! Loop through 30 cyles
    do pass = 1, 10, 2      ! List ascending passes 1-9
        call rads_open_pass (S, P, cycle, pass)              ! Open pass file
        allocate (lon(P%ndata), lat(P%ndata), sla(P%ndata))  ! Allocate memory
        call rads_get_var (S, P, 'lon', lon)                 ! Get longitude
        call rads_get_var (S, P, 'lat', lat)                 ! Get latitude
        call rads_get_var (S, P, 'sla', sla)                 ! Get sea level anomaly
        do i = 1,P%ndata
            write (*,*) lon(i), lat(i), sla(i)               ! Print the data
        enddo
        deallocate (lon, lat, sla)                           ! Deallocate memory
        call rads_close_pass (S, P)                          ! Close pass
    enddo
enddo
call rads_end (S)           ! End RADS; free data structure
end program rads_sample
