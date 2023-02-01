program nc2fes

! This program converts a NetCDF file of amplitudes and phases to
! the the ASCII format used for the FES tide model distribution.
!-
! 26-Sep-2006 - Based on fes2nc
!-----------------------------------------------------------------------
use typesizes
use netcdf
use nf_subs
integer :: i,j,k,nx,ny,ncid
real(eightbytereal) :: xmin,ymin,xmax,ymax,dx,dy,mask=999.9d0,offset,fact,dummy(2),fill
real(eightbytereal), allocatable, dimension(:,:) :: wra, wrg
character(80) :: file_bin

500 format('fes2nc -- Convert NetCDF to FES Ascii grid file'//'syntax: nc2fes <NetCDF_files>')
if (iargc()==0) then
    write (*,500)
    stop
endif

! Read input file name

call getarg(1,file_bin)

write (0,*) 'Reading NetCDF file ... will be quick ...'
write (0,*) 'File : ', trim(file_bin)

! Open NetCDF file and read header info

call nfs(nf90_open(file_bin,nf90_nowrite,ncid))
call nfs(nf90_inquire_dimension(ncid,1,len=nx))
call nfs(nf90_inquire_dimension(ncid,2,len=ny))
call nfs(nf90_get_att(ncid,1,"actual_range",dummy))
xmin=dummy(1) ; xmax=dummy(2)
dx = (xmax - xmin) / (nx-1)
call nfs(nf90_get_att(ncid,2,"actual_range",dummy))
ymin=dummy(1) ; ymax=dummy(2)
dy = (ymax - ymin) / (ny-1)

! Allocate memory

allocate (wra(nx,ny),wrg(nx,ny))

! Write ascii header

write(*,'(2f10.3)') xmin,xmax
write(*,'(2f10.3)') ymin,ymax
write(*,'(f8.6,f13.6)') dx,dy
write(*,'(2i8)') nx,ny
write(*,'(2f9.3," (mask value)")') mask,mask

! Read NetCDF grids

call nfs(nf90_get_var(ncid,3,wra))
if (nf90_get_att(ncid,3,"scale_factor",fact)/=nf90_noerr) fact=1d0
if (nf90_get_att(ncid,3,"add_offset",offset)/=nf90_noerr) offset=0d0
fact=fact*1d2 ; offset=offset*1d2	! m to cm
wra=wra*fact+offset
if (nf90_get_att(ncid,3,"_FillValue",fill)==nf90_noerr) then
    fill=fill*fact+offset
    where (wra==fill) wra=mask
endif
where (wra/=wra) wra=mask

call getarg(2,file_bin)
if (file_bin==' ') then
    i=4
else
    i=3
    write (0,*) 'File : ', trim(file_bin)
    call nfs(nf90_close(ncid))
    call nfs(nf90_open(file_bin,nf90_nowrite,ncid))
endif

call nfs(nf90_get_var(ncid,i,wrg))
if (nf90_get_att(ncid,i,"scale_factor",fact)/=nf90_noerr) fact=1d0
if (nf90_get_att(ncid,i,"add_offset",offset)/=nf90_noerr) offset=0d0
wrg=wrg*fact+offset
if (nf90_get_att(ncid,i,"_FillValue",fill)==nf90_noerr) then
    fill=fill*fact+offset
    where (wrg==fill) wrg=mask
endif
where (wrg<0d0) wrg=wrg+360d0
where (wrg/=wrg) wrg=mask
call nfs(nf90_close(ncid))

! Write rest of ascii file
   
write (0,*) 'Writing Ascii file ... be patient ...'
do j=1,ny
    do k=1,nx,30
	write(*,'(30f7.2)') (wra(i,j),i=k,MIN(nx,k+29))
	write(*,'(30f7.2)') (wrg(i,j),i=k,MIN(nx,k+29))
    enddo
enddo

deallocate (wra,wrg)
   
end
