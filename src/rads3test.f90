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

use rads3

integer maxdata
parameter (maxdata=3600)
real*8 eqtime,eqlon
real*8 data(maxdata,2),time(maxdata),dlon(maxdata),dlat(maxdata)
integer*4 cycle,pass,select(2),ndata,verbose,i
character*256 mission,metafile 

mission = 'e1/d'
verbose = 2
cycle = 111
select(1) = 0
select(2) = 16

call getraw_init(mission,verbose)
call getraw_limits(2,-90d0,0d0)
call getraw_options(16,11)

do pass = 1,10
	call getraw (cycle,pass,2,maxdata,select, &
	time,dlat,dlon,data,ndata,eqtime,eqlon,metafile)
	
	do i = 1,ndata
		write(*,*) i,time(i),dlat(i),dlon(i),data(i,1),data(i,2)
	enddo
enddo

call getraw_stat(0)
end
