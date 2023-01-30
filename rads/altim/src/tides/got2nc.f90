program got2nc

! This program convert the ASCII distributed files for the
! GOT Tide Model to NetCDF grid files in order to save
! disk space and time in reading the files.
! The conversion is 1:1 and should be completely reversible.
!-
! $Log: got2nc.f90,v $
! Revision 1.8  2014/02/08 16:39:51  rads
! - Replaced cmplx by dcmplx to avoid loss of precision
! - Deflate netCDF files
!
! Revision 1.7  2011/11/04 00:22:28  rads
! - Avoid warning on nint
!
! Revision 1.6  2011/03/02 21:32:29  rads
! - Get model name for first and S2 component (supporting GOT4.8 and GOT4.9)
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
real(eightbytereal) :: xmin,ymin,xmax,ymax,dx,dy,mask,hmin,hmax,offset,fact,scale
real(eightbytereal), parameter :: rad=atan(1d0)/45d0
real(eightbytereal), allocatable, dimension(:,:) :: wra, wrg
integer(twobyteint), allocatable, dimension(:,:) :: amp, pha
integer(twobyteint) :: i2min=-32768
character(80) :: file_bin,unit,format
character(160) :: head
character(8) :: compo
complex(eightbytereal) :: val
logical :: new

500 format('got2nc -- Convert GOT Ascii to netCDF grid file'//'syntax: got2nc netcdf_file < input')
call getarg(1,file_bin)
if (file_bin == ' ') then
	write (*,500)
	stop
endif

! Read input parameters

hmin=1d30; hmax=-1d30; offset=0d0
   
! Read ascii header

550 format (a)
read (*,550,iostat=i) head(1:80)
if (i /= 0) stop
i = index(head,' ')-1
compo = head(:i)
write (*,*) 'Reading ASCII file ... component '//trim(compo)
read (*,550) head(81:160)
read (*,*) ny,nx
read (*,*) ymin,ymax
read (*,*) xmin,xmax
read (*,*) mask
read (*,550) format
dx=(xmax-xmin)/(nx-1)
dy=(ymax-ymin)/(ny-1)

! Allocate memory

mx=nx/2*2+1
allocate (wra(nx,ny),wrg(nx,ny),amp(mx,ny),pha(mx,ny))

! Read rest of ascii file

do j=1,ny
	read(*,format) (wra(i,j),i=1,nx)
enddo
do j=1,7
	read(*,550)
enddo
do j=1,ny
	read(*,format) (wrg(i,j),i=1,nx)
	do i=1,nx
		if (wra(i,j)/=mask) then
			hmin=min(hmin,wra(i,j))
			hmax=max(hmax,wra(i,j))
		endif
		if (wrg(i,j)/=mask .and. wrg(i,j)>=180d0) wrg(i,j)=wrg(i,j)-360d0
	enddo
enddo

fact=1d-2
if (index(format,'.3')>0) fact=1d-3
if (index(head,'(cm)') > 0) then
	scale=fact*1d-2
	unit='m'
else if (index(head,'(mm)') > 0) then
	scale=fact*1d-3
	unit='m'
else
	fact=1d-1
	scale=1d-4
	unit='millibar'
endif
write (*,*) 'Conversion ready. Min/Max value = ',hmin/fact,hmax/fact
if (hmax/fact > 32767d0) offset=3d4

! Some grids have the North Pole undefined. Fix that here

if (wra(1,ny)==mask) then
	val = dot_product (wra(1:nx,ny-1),expj(wrg(1:nx,ny-1)*rad)) / nx
	wra(:,ny)=abs(val)
	wrg(:,ny)=atan2(dimag(val),real(val))/rad
	write (*,*) 'North pole values fixed = ',wra(1,ny),wrg(1,ny)
endif

! Shift longitude to 0-360 degrees

if (xmin==-180) then
	dn=nx/2
	xmin=xmin+180d0
else
	dn=0
endif
xmax=xmin+360d0

! Create or open netCDF grid

i=index(head,'GOT')
new = (nf90_open(file_bin,nf90_write,ncid) /= nf90_noerr)
if (new) then
	write (*,*) 'Creating netCDF file: '//trim(file_bin)
	call nfs(nf90_create(file_bin,nf90_write+nf90_hdf5,ncid))
	call nf90_def_axis(ncid,'lon','longitude','degrees_east',mx,xmin,xmax,dimid(1),varid(1))
	call nf90_def_axis(ncid,'lat','latitude','degrees_north',ny,ymin,ymax,dimid(2),varid(2))
else
	call nfs(nf90_inq_dimid(ncid,'lon',dimid(1)))
	call nfs(nf90_inq_dimid(ncid,'lat',dimid(2)))
	call nfs(nf90_redef(ncid))
endif
if (new .or. compo == 'S2') then
	if (index(head,'Load') > 0) then
		call nfs(nf90_put_att(ncid,nf90_global,'model',trim(head(i:i+6))//' load tide'))
	else
		call nfs(nf90_put_att(ncid,nf90_global,'model',trim(head(i:i+6))//' ocean tide'))
	endif
	call nfs(nf90_put_att(ncid,nf90_global,'description',trim(head(81:160))))
endif

! Define new component

write (*,*) 'Writing netCDF file ... component '//trim(compo)
call nfs(nf90_def_var(ncid,'amp_'//trim(compo),nf90_int2,dimid,varid(3)))
call nfs(nf90_def_var_deflate(ncid,varid(3),1,1,9))
call nfs(nf90_put_att(ncid,varid(3),'long_name',trim(compo)//' amplitude'))
call nfs(nf90_put_att(ncid,varid(3),'units',trim(unit)))
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

do j=1,ny
	do i=1,nx
		if(i<=dn) then
			k=i+dn
		else
			k=i-dn
		endif
		if (mask==0d0 .and. wra(i,j)==mask .and. wrg(i,j)==mask) then
			amp(k,j)=i2min
			pha(k,j)=i2min
		else if (mask/=0d0 .and. (wra(i,j)==mask .or. wrg(i,j)==mask)) then
			amp(k,j)=i2min
			pha(k,j)=i2min
		else
			amp(k,j)=nint(wra(i,j)/fact-offset,twobyteint)
			pha(k,j)=nint(wrg(i,j)/1d-2,twobyteint)
		endif
	enddo
	amp(mx,j)=amp(1,j)
	pha(mx,j)=pha(1,j)
enddo

! Write out grids

call nfs(nf90_put_var(ncid,varid(3),amp))
call nfs(nf90_put_var(ncid,varid(4),pha))

! Close grid file

call nfs(nf90_close(ncid))

deallocate (wra,wrg,amp,pha)
   
contains

elemental function expj (x)
real(eightbytereal), intent(in) :: x
complex(eightbytereal) :: expj
expj = dcmplx(cos(x),sin(x))
end function expj

end program got2nc
