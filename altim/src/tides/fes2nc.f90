program fes2nc

! This program convert the ASCII distributed files for the
! FES Tide Model to netCDF grid files in order to save
! disk space and time in reading the files.
! The conversion is 1:1 and should be completely reversible.
!-
! $Log: fes2nc.f90,v $
! Revision 1.7  2014/02/08 16:39:51  rads
! - Replaced cmplx by dcmplx to avoid loss of precision
! - Deflate netCDF files
!
! Revision 1.6  2011/11/04 00:24:18  rads
! - Let fes2nc read from standard input
! - Allow use of DTU10 grids
! - When concatenating two grids, subtract second from first
!
! Revision 1.5  2009/05/20 15:49:39  rads
! - Add one component at the time to a single netCDF file
!
! 22-Sep-2006 - Based on ascii2grd_fes
!-----------------------------------------------------------------------
use typesizes
use netcdf
use nf_subs
integer :: i,j,k,nx,ny,dn,mx,ncid,varid(4),dimid(2)
real(eightbytereal), parameter :: rad = atan(1d0) / 45d0
real(eightbytereal) :: xmin,ymin,xmax,ymax,dx,dy,mask,hmin,hmax,offset,fact,scale
real(eightbytereal), allocatable :: wra(:,:), wrg(:,:), ura(:), urg(:), c(:), s(:)
integer(twobyteint), allocatable :: amp(:,:), pha(:,:)
integer(twobyteint) :: i2min=-32768
character(80) :: file_bin
character(8) :: compo
logical :: new

500 format('fes2nc -- Convert FES Ascii to netCDF grid file'//'syntax: fes2nc netcdf_file compo < input')
call getarg(1,file_bin)
call getarg(2,compo)
if (compo == ' ') then
	write (*,500)
	stop
endif

! Read input parameters

hmin=1d30; hmax=-1d30; offset=0d0

write (*,*) 'Reading ASCII file ... component '//trim(compo)

! Read ascii header

read (*,*) xmin,xmax
read (*,*) ymin,ymax
read (*,*) dx,dy
read (*,*) nx,ny
read (*,*) mask,mask

! Allocate memory

mx=nx/2*2+1
allocate (wra(nx,ny),wrg(nx,ny),ura(nx),urg(nx),c(nx),s(nx),amp(mx,ny),pha(mx,ny))

! Read rest of ascii file
   
do j = 1,ny
	do k = 1,nx,30
		read (*,*) (wra(i,j),i=k,MIN(nx,k+29))
		read (*,*) (wrg(i,j),i=k,MIN(nx,k+29))
	enddo
enddo

! If there is more, then subtract the second part from the first

read (*,*,iostat=i) xmin,xmax
if (i == 0) then
	read (*,*) ymin,ymax
	read (*,*) dx,dy
	read (*,*) nx,ny
	read (*,*) mask,mask
	do j = 1,ny
		do k = 1,nx,30
			read (*,*) (ura(i),i=k,MIN(nx,k+29))
			read (*,*) (urg(i),i=k,MIN(nx,k+29))
		enddo
		c = wra(:,j) * cos(wrg(:,j)*rad) - ura * cos(urg*rad)
		s = wra(:,j) * sin(wrg(:,j)*rad) - ura * sin(urg*rad)
!		write (*,*) 1,j,wra(1,j),wrg(1,j),ura(1),urg(1),c(1),s(1)
		wra(:,j) = sqrt(c*c + s*s)
		wrg(:,j) = atan2(s, c) / rad
!		write (*,*) 1,j,wra(1,j),wrg(1,j)
	enddo
endif
do j = 1,ny
	do i = 1,nx
		if (wra(i,j) /= mask) then
			hmin = min(hmin,wra(i,j))
			hmax = max(hmax,wra(i,j))
		endif
		if (wrg(i,j) /= mask .and. wrg(i,j) >= 180d0) wrg(i,j) = wrg(i,j) - 360d0
	enddo
enddo
fact=1d-2
scale=1d-4
write (*,*) 'Conversion ready. Min/Max value = ',hmin/fact,hmax/fact
if (hmax/fact > 32767d0) offset=3d4

! Shift longitude to 0-360 degrees

if (xmin==-180) then
	dn=nx/2
	xmin=xmin+180d0
else
	dn=0
endif
xmax=xmin+360d0

! Create or open netCDF grid file

new = (nf90_open(file_bin,nf90_write,ncid) /= nf90_noerr)
if (new) then
	write (*,*) 'Creating netCDF file: '//trim(file_bin)
	call nfs(nf90_create(file_bin,nf90_write+nf90_hdf5,ncid))
	i = index(file_bin, 'FES')
	if (i == 0) i = index(file_bin, 'DTU')
	j = index(file_bin(i:), '/') + i - 2
	if (index(file_bin,'load') > 0) then
		call nfs(nf90_put_att(ncid,nf90_global,'model',file_bin(i:j)//' load tide'))
	else
		call nfs(nf90_put_att(ncid,nf90_global,'model',file_bin(i:j)//' ocean tide'))
	endif
	call nf90_def_axis(ncid,'lon','longitude','degrees_east',mx,xmin,xmax,dimid(1),varid(1))
	call nf90_def_axis(ncid,'lat','latitude','degrees_north',ny,ymin,ymax,dimid(2),varid(2))
else
	call nfs(nf90_inq_dimid(ncid,'lon',dimid(1)))
	call nfs(nf90_inq_dimid(ncid,'lat',dimid(2)))
	call nfs(nf90_redef(ncid))
endif

! Define new component

write (*,*) 'Writing netCDF file ... component '//trim(compo)
call nfs(nf90_def_var(ncid,'amp_'//trim(compo),nf90_int2,dimid,varid(3)))
call nfs(nf90_def_var_deflate(ncid,varid(3),1,1,9))
call nfs(nf90_put_att(ncid,varid(3),'long_name',trim(compo)//' amplitude'))
call nfs(nf90_put_att(ncid,varid(3),'units','m'))
if (scale /= 1d0) call nfs(nf90_put_att(ncid,varid(3),'scale_factor',scale))
if (offset /= 0d0) call nfs(nf90_put_att(ncid,varid(3),'add_offset',offset*scale))
call nfs(nf90_put_att(ncid,varid(3),'_FillValue',i2min))
call nfs(nf90_def_var(ncid,'pha_'//trim(compo),nf90_int2,dimid,varid(4)))
call nfs(nf90_def_var_deflate(ncid,varid(4),1,1,9))
call nfs(nf90_put_att(ncid,varid(4),'long_name',trim(compo)//' phase'))
call nfs(nf90_put_att(ncid,varid(4),'units','degrees'))
call nfs(nf90_put_att(ncid,varid(4),'scale_factor',1d-2))
call nfs(nf90_put_att(ncid,varid(4),'_FillValue',i2min))
call nfs(nf90_enddef(ncid))

! Write coordinate axes

if (new) then
	call nf90_put_axis(ncid,varid(1))
	call nf90_put_axis(ncid,varid(2))
endif

! Convert to integers

do j = 1,ny
	do i = 1,nx
		if (i <= dn) then
			k = i + dn
		else
			k = i - dn
		endif
		if (mask == 0d0 .and. wra(i,j) == mask .and. wrg(i,j) == mask) then
			amp(k,j) = i2min
			pha(k,j) = i2min
		else if (mask /= 0d0 .and. (wra(i,j) == mask .or. wrg(i,j) == mask)) then
			amp(k,j) = i2min
			pha(k,j) = i2min
		else
			amp(k,j) = nint(wra(i,j)/fact-offset, twobyteint)
			pha(k,j) = nint(wrg(i,j)/fact, twobyteint)
		endif
	enddo
	amp(mx,j) = amp(1,j)
	pha(mx,j) = pha(1,j)
enddo

! Write out grids

call nfs(nf90_put_var(ncid,varid(3),amp))
call nfs(nf90_put_var(ncid,varid(4),pha))

! Close grid file

call nfs(nf90_close(ncid))

deallocate (wra,wrg,amp,pha)

end program fes2nc
