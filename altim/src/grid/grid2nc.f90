!*GRID2NC -- Converts DEOS grid to COARDS-compatible NetCDF grid
!+
program grid2nc

! This program reads a DEOS grid and dumps it again as a NetCDF
! grid. No other transformations are made, i.e., the data remain
! stored as integers and the offset and scale remain unaltered.
!-
! 05-Dec-2009 - Converted to Fortran90
! 26-Jul-2006 - Avoiding use of %val
! 17-Aug-2005 - Created by Remko Scharroo, Altimetrics LLC
!-----------------------------------------------------------------------
use typesizes
use netcdf
use gridlib
use nf_subs
character(80) :: gridin=' ',gridout=' ',arg
character(256) :: history
integer(fourbyteint) :: i,ncid,x_id,y_id,z_id,dims(2),l,noff=0,iargc
logical :: geographic=.false.
real(eightbytereal), allocatable :: tmp(:)
real(eightbytereal) :: dummy(2),z0
type(grid) :: info

! Read input and output file name

call getarg(0,history)
l = len_trim(history)
do i = 1,iargc()
	call getarg(i,arg)
	history(:l+2)=arg
	l=len_trim(history)
	if (arg(:2) == '-F') then
		noff=1
	else if (arg(:2) == '-G') then
		geographic=.true.
	else if (gridin == ' ') then
		gridin=arg
	else
		gridout=arg
	endif
enddo

! If output grid is not specified, print usage message

if (gridout == ' ') then
	write (*,1300)
	goto 9999
endif
1300  format ('grid2nc -- Convert DEOS grid to NetCDF grid'// &
'syntax: grid2nc [-F] [-G] deos-grid netcdf-grid'// &
'deos-grid  : Name of DEOS grid file'/ &
'netcdf-grid: Name of NetCDF grid file'/ &
'-F         : Force pixel registration (Default: grid reg.)'/ &
'-G         : Make grid geographical')

! Load input grid and copy header

if (grid_load(gridin,info) /= 0) call fin ('Error reading input grid')

! Change grid boundaries when noff=1

info%xmin = info%xmin - noff * info%dx/2d0
info%xmax = info%xmax + noff * info%dx/2d0
info%ymin = info%ymin - noff * info%dy/2d0
info%ymax = info%ymax + noff * info%dy/2d0

! Open output file

call nfs(nf90_create(gridout,nf90_clobber,ncid))

! Set up grid variables and dimensions

if (geographic) then
	call nfs(nf90_def_dim(ncid,"lon",info%nx,dims(1)))
	call nfs(nf90_def_dim(ncid,"lat",info%ny,dims(2)))
	call nfs(nf90_def_var(ncid,"lon",nf90_double,dims(1),x_id))
	call nfs(nf90_def_var(ncid,"lat",nf90_double,dims(2),y_id))
	call nfs(nf90_put_att(ncid,x_id,"long_name","longitude"))
	call nfs(nf90_put_att(ncid,y_id,"long_name","latitude"))
	call nfs(nf90_put_att(ncid,x_id,"units","degrees_east"))
	call nfs(nf90_put_att(ncid,y_id,"units","degrees_north"))
else
	call nfs(nf90_def_dim(ncid,"x",info%nx,dims(1)))
	call nfs(nf90_def_dim(ncid,"y",info%ny,dims(2)))
	call nfs(nf90_def_var(ncid,"x",nf90_real,dims(1),x_id))
	call nfs(nf90_def_var(ncid,"y",nf90_real,dims(2),y_id))
endif
call nfs(nf90_def_var(ncid,"z",info%ntype,dims,z_id))

! Add attributes

z0 = (info%zmin + info%zmax) / 2d0
call nfs(nf90_put_att (ncid,z_id,"scale_factor",info%dz))
call nfs(nf90_put_att (ncid,z_id,"add_offset",z0))
if (info%ntype == 1) then
	call nfs(nf90_put_att(ncid,z_id,"_FillValue",nint(info%znan,onebyteint)))
else if (info%ntype == 3) then
	call nfs(nf90_put_att(ncid,z_id,"_FillValue",nint(info%znan,twobyteint)))
else if (info%ntype == 4) then
	call nfs(nf90_put_att(ncid,z_id,"_FillValue",nint(info%znan,fourbyteint)))
endif
call nfs(nf90_put_att (ncid,nf90_global,"node_offset",noff))
call nfs(nf90_put_att (ncid,nf90_global,"title",gridout))
call nfs(nf90_put_att (ncid,nf90_global,"history",history))
dummy(1) = info%xmin
dummy(2) = info%xmax
call nfs(nf90_put_att (ncid,x_id,"actual_range",dummy))
dummy(1) = info%ymin
dummy(2) = info%ymax
call nfs(nf90_put_att (ncid,y_id,"actual_range",dummy))
dummy(1) = info%zmin
dummy(2) = info%zmax
call nfs(nf90_put_att (ncid,z_id,"valid_range",dummy))

! Set up x- and y-columns

call nfs(nf90_enddef(ncid))
allocate (tmp(info%nx))
do i=1,info%nx
	tmp(i)=info%xmin+(i-1+noff/2d0)*info%dx
enddo
call nfs(nf90_put_var(ncid,x_id,tmp))
deallocate (tmp)
allocate (tmp(info%ny))
do i=1,info%ny
	tmp(i)=info%ymin+(i-1+noff/2d0)*info%dy
enddo
call nfs(nf90_put_var(ncid,y_id,tmp))
deallocate (tmp)

! Finally, dump the whole grid

if (info%ntype == 1) then
	call nfs(nf90_put_var(ncid,z_id,info%grid_int1))
else if (info%ntype == 3) then
	call nfs(nf90_put_var(ncid,z_id,info%grid_int2))
else if (info%ntype == 4) then
	call nfs(nf90_put_var(ncid,z_id,info%grid_int4))
else if (info%ntype == 5) then
	call nfs(nf90_put_var(ncid,z_id,info%grid_real))
else if (info%ntype == 6) then
	call nfs(nf90_put_var(ncid,z_id,info%grid_dble))
else
	call fin('Unknown data type')
endif

! Close grid

call nfs(nf90_close(ncid))
9999 continue

end program grid2nc
