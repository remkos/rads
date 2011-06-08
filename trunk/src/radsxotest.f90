!-----------------------------------------------------------------------
! $Id$
!
! Copyright (C) 2011  Remko Scharroo (Altimetrics LLC)
! See LICENSE.TXT file for copying and redistribution conditions.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!-----------------------------------------------------------------------

!*radsxotest -- Test with RADS xover file
!+
program radsxotest
use typesizes
use netcdf
use rads_netcdf

character(len=80) :: filename
integer(fourbyteint) :: ncid, varid

call getarg(1,filename)

call nfs(nf90_open(filename, nf90_write, ncid))

call nfs(nf90_redef(ncid))
call nfs(nf90_def_var(ncid, 'satellite', nf90_int, (/ 1 /), varid))
call nfs(nf90_enddef(ncid))
call nfs(nf90_close(ncid))

end program radsxotest
