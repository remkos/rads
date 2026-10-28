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

!***********************************************************************
!*rads_get_var_tpj -- Dummy replacement for the real rads_get_var_tpj
!+
subroutine rads_get_var_tpj
!
! This routine does nothing but return error code rads_err_var
!-----------------------------------------------------------------------
S%error = rads_err_var
end subroutine rads_get_var_tpj
