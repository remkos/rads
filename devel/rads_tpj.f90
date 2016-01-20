!-----------------------------------------------------------------------
! Copyright (c) 2011-2016  Remko Scharroo
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

!***********************************************************************
!*rads_get_var_tpj -- Create special TPJ data on the fly
!+
subroutine rads_get_var_tpj
!
! This routine produces data on the fly specific to TOPEX/Poseidon and
! Jason-1/2
!-----------------------------------------------------------------------
use tpj_subs
real(eightbytereal), save :: old_equator_time = 0d0, pbias
real(eightbytereal), save, allocatable :: alpha(:), betap(:), sapitch(:), scpitch(:), yaw(:)
logical, save, allocatable :: sunlit(:)
integer(fourbyteint), save, allocatable :: ymode(:)
integer(fourbyteint) :: i

if (old_equator_time /= P%equator_time) then
	deallocate (sunlit, alpha, betap, sapitch, scpitch, ymode, yaw, stat = i)
	allocate (sunlit(P%ndata), alpha(P%ndata), betap(P%ndata), sapitch(P%ndata), &
		scpitch(P%ndata), ymode(P%ndata), yaw(P%ndata))
	do i = 1,P%ndata
		call tpj_cgcorr (S%sat, P%equator_time, P%equator_lon, modulo(P%pass,2)==1, P%tll(i,1), &
			sunlit(i), alpha(i), betap(i), sapitch(i), scpitch(i), ymode(i), yaw(i), pbias)
	enddo
	pbias = cos(pbias*rad)**0.576d0 ! Correction for flexing due to reduced influx at non-zero SA pitch bias
	old_equator_time = P%equator_time
endif

S%error = rads_noerr

select case (var%name)
case ('tpj_sunlit')
	data = 0d0
	where (sunlit) data = 1d0
case ('tpj_alpha')
	data = alpha
case ('tpj_beta_prime')
	data = betap
case ('tpj_sa_pitch')
	data = sapitch
case ('tpj_sc_pitch')
	data = scpitch
case ('tpj_sin_pitch_sunlit')
	data = 0d0
	where (sunlit) data = pbias * sin(sapitch*rad)
case ('tpj_sin_pitch_shadow')
	data = 0d0
	where (.not.sunlit) data = sin(sapitch*rad)
case ('tpj_cos_pitch')
	data = cos(sapitch*rad)
case ('tpj_yaw_mode')
	data = ymode
case ('tpj_yaw_angle')
	data = yaw
case default
	S%error = rads_err_var
end select
end subroutine rads_get_var_tpj
