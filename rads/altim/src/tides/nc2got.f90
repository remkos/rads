program nc2got

! This program converts a NetCDF file of amplitudes and phases to
! the the ASCII format used for the GOT tide model distribution.
!-
! 26-Sep-2006 - Based on got2nc
!-----------------------------------------------------------------------
use typesizes
use netcdf
use nf_subs
integer :: i,j,nx,ny,ncid
real(eightbytereal) :: xmin,ymin,xmax,ymax,dx,dy,mask=999d0,offset,fact,dummy(2),fill
real(eightbytereal), allocatable, dimension(:,:) :: wra, wrg
character(80) :: file_bin

500 format('got2nc -- Convert NetCDF to GOT Ascii grid file'//'syntax: nc2got <NetCDF_files>')
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

write(*,'("#"/"#")')
write(*,'(2i21)') ny,nx-1
write(*,'(2f25.4)') ymin,ymax
write(*,'(2f25.4)') xmin,xmax-dx
write(*,'(2f25.2)') mask,mask
write(*,'("(11F7.2)")')

! Read NetCDF grids

call nfs(nf90_get_var(ncid,3,wra))
if (nf90_get_att(ncid,3,"scale_factor",fact)/=nf90_noerr) fact=1d0
if (nf90_get_att(ncid,3,"add_offset",offset)/=nf90_noerr) offset=0d0
fact=fact*1d3 ; offset=offset*1d3	! m to mm
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
    write (*,'(11f7.2)') wra(1:nx-1,j)
enddo
do j=1,7
    write (*,'("#")')
enddo
do j=1,ny
    write (*,'(11f7.2)') wrg(1:nx-1,j)
enddo

deallocate (wra,wrg)
   
end
