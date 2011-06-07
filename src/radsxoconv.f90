!*radsxoconv -- RADS crossover converter
!+
program radsxoconv

! This program converts RADS netCDF crossover files into the historical
! DEOS XAF, XTF and XXF/XXO format.
!-
! Created by Remko Scharroo, Altimetrics LLC, based in part on previous
! programs, max and max2, developed at DEOS.
!-----------------------------------------------------------------------
use rads
use netcdf
use rads_netcdf
integer(fourbyteint) :: ncid, i, nvars
real(eightbytereal), allocatable :: lat(:), lon(:), time(:,:), sla(:,:)
integer(fourbyteint), allocatable :: track(:,:)
integer(fourbyteint) :: id_sla = 0, id_alt_rate = 0, id_alt = 0, id_sig0 = 0, id_swh = 0, npar = 3, ios
type :: trk_
	real(eightbytereal) :: equator_lon, equator_time, start_time, end_time
	integer(twobyteint) :: nr_alt, nr_xover, satid, cycle, pass
endtype
type(trk_), allocatable :: trk(:)
type :: prod_
	character(len=1) :: fmt
	integer(fourbyteint) :: nr
endtype prod_
type(prod_) :: a, x, t
character(len=80) :: arg, filename = '', shortname
integer(fourbyteint), parameter :: msat = 11
real(eightbytereal), parameter :: period(msat) = (/6037.582d0,6037.582d0,6037.582d0,6035.928d0,6745.759d0,6745.759d0,6035.928d0, &
	6037.582d0,6745.759d0,6035.928d0,6745.759d0/)
real(eightbytereal), parameter :: altsig(msat) = (/0.2d0,0.07d0,0.07d0,0.05d0,0.02d0,0.04d0,0.05d0,0.05d0,0.04d0,0.05d0,0.04d0/)
real(eightbytereal), parameter :: orberr(msat) = (/0.8d0,0.8d0,0.15d0,0.05d0,0.035d0,0.035d0,0.05d0,0.05d0,0.035d0,0.05d0,0.035d0/)
real(eightbytereal), parameter :: inclination(msat) = (/60d0,108.05d0,108.05d0,98.54d0,66.04d0,66.04d0,98.54d0,108.05d0,66.04d0, &
	98.54d0,66.04d0/)
!character(len=2) :: sat(11) = (/'g3','ss','gs','e1','tx','pn','e2','g1','j1','n1','j2'/)
integer(twobyteint), parameter :: bound(4) = (/-180,180,-90,90/)
integer(fourbyteint), parameter :: maxtrk = 32768
logical :: l

! Initialize RADS or issue help
call synopsis
a = prod_ ('', 0)
t = prod_ ('', 0)
x = prod_ ('', 0)

! Start with this-is message
l = rads_version ('$Rev: 4$')

! Scan command line arguments
do i = 1,iargc()
	call getarg(i,arg)
	if (arg == '-xaf') then
		t%fmt = 'a'
	else if (arg == '-xtf') then
		t%fmt = 't'
	else if (arg == '-xxo') then
		x%fmt = 'o'
	else if (arg == '-xxf') then
		x%fmt = 'x'
	else if (arg(:2) == '-p') then
		read (arg(3:), *, iostat=ios) npar
	endif
enddo

! Now process each file
do i = 1,iargc()
	call getarg(i,filename)
	if (filename(:1) == '-') cycle
	call process
enddo

contains

!***********************************************************************

subroutine process
integer(fourbyteint) :: i

i = index(filename, '.nc')
if (i > 0) then
	shortname = filename(:i-1)
else
	shortname = filename
endif

! Open netCDF file
call nfs (nf90_open (filename, nf90_write, ncid))

! Check out which data variables are provided
call nfs (nf90_inquire (ncid, nvariables=nvars))
do i = 1,nvars
	call nfs (nf90_inquire_variable (ncid, i, name=arg))
	if (arg == 'sla') then
		id_sla = i
	else if (arg == 'alt_rate') then
		id_alt_rate = i
	else if (arg == 'sig0') then
		id_sig0 = i
	else if (arg == 'swh') then
		id_swh = i
	else if (arg == 'alt' .or. arg(:4) == 'alt_') then
		id_alt = i
	endif
enddo

! Check if we have required elements for requested format
if (id_sla == 0) then
	write (*,'(a)') 'Output format requires sla field. Program terminating.'
	stop
endif

! Get the number of xovers and number of tracks
call nfs (nf90_inquire_dimension (ncid, 2, len=x%nr))
call nfs (nf90_inquire_dimension (ncid, 3, len=t%nr))
if (t%nr >= maxtrk) then
	write (*,'(a)') 'Output format does not allow track numbers exceeding 32767. We will modulo them.'
endif

! Read all the "base variables" into memory
allocate (lat(x%nr), lon(x%nr), time(2,x%nr), track(2,x%nr), sla(2,x%nr), trk(t%nr))
call get_var_1d (get_varid('lat'), lat)
call get_var_1d (get_varid('lon'), lon)
call get_var_2d (get_varid('time'), time)
call nfs (nf90_get_var (ncid, get_varid('track'), track))
call nfs (nf90_get_var (ncid, get_varid('satid'), trk%satid))
call nfs (nf90_get_var (ncid, get_varid('cycle'), trk%cycle))
call nfs (nf90_get_var (ncid, get_varid('pass'), trk%pass))
call nfs (nf90_get_var (ncid, get_varid('equator_lon'), trk%equator_lon))
call nfs (nf90_get_var (ncid, get_varid('equator_time'), trk%equator_time))
call nfs (nf90_get_var (ncid, get_varid('start_time'), trk%start_time))
call nfs (nf90_get_var (ncid, get_varid('end_time'), trk%end_time))
call nfs (nf90_get_var (ncid, get_varid('nr_xover'), trk%nr_xover))
call nfs (nf90_get_var (ncid, get_varid('nr_alt'), trk%nr_alt))
call get_var_2d (id_sla, sla)

600 format (/ &
'File name              : ',a/ &
'Number of xovers read  : ',i9/ &
'Number of tracks found : ',i9)
write (*, 600) trim(filename), x%nr, t%nr

! Write out various file formats
if (a%fmt == 'a') call write_xaf
if (x%fmt == 'x') then
	call write_xxf
else if (x%fmt == 'o') then
	call write_xxo
endif
if (t%fmt == 't') call write_xtf

call nfs (nf90_close (ncid))
deallocate (lat, lon, time, track, sla, trk)
end subroutine process

!***********************************************************************

subroutine synopsis
if (rads_version ('$Revision$','RADS crossover file converter')) return
call rads_synopsis ()
write (0,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  filename          : Input file name (extension .nc)'/ &
'  -xaf              : Create XAF file (extension .xaf)'/ &
'  -xtf              : Create XTF file (extension .xtf)'/ &
'  -xxf              : Create XXF file (extension .xxf)'/ &
'  -xxo              : Create XXO file (extension .xxo, .xxo.aux)'// &
'The extensions are added to the input file after removing .nc')
stop
end subroutine synopsis

!***********************************************************************

function get_varid (name)
character(len=*), intent(in) :: name
integer(fourbyteint) :: get_varid
call nfs (nf90_inq_varid(ncid,name,get_varid))
end function get_varid

!***********************************************************************

subroutine get_var_1d (varid, array)
integer(fourbyteint), intent(in) :: varid
real(eightbytereal), intent(out) :: array(:)
real(eightbytereal) :: value
call nfs (nf90_get_var (ncid, varid, array))
if (nf90_get_att (ncid, varid, 'scale_factor', value) == nf90_noerr) array = array * value
if (nf90_get_att (ncid, varid, 'add_offset', value) == nf90_noerr) array = array + value
end subroutine get_var_1d

!***********************************************************************

subroutine get_var_2d (varid, array)
integer(fourbyteint), intent(in) :: varid
real(eightbytereal), intent(out) :: array(:,:)
real(eightbytereal) :: value
call nfs (nf90_get_var (ncid, varid, array))
if (nf90_get_att (ncid, varid, 'scale_factor', value) == nf90_noerr) array = array * value
if (nf90_get_att (ncid, varid, 'add_offset', value) == nf90_noerr) array = array + value
end subroutine get_var_2d

!***********************************************************************

subroutine write_xxf
type :: xxf_ ! 48 bytes
	integer(fourbyteint) :: lat, lon, time(2)
	integer(twobyteint) :: track(2)
	integer(fourbyteint) :: sla(2), sla_cor(2), omega(2)
	integer(twobyteint) :: sigma(2)
endtype
type(xxf_) :: xxf
integer(fourbyteint) :: i
real(eightbytereal) :: omega(2)
open (10, file=trim(shortname)//'.xxf', status='replace', access='direct', form='unformatted', recl=48)
write (10, rec=1) '@XXF',x%nr,bound,0
do i = 1,x%nr
	omega = (time(:,i) - trk(track(:,i))%equator_time) / period(trk(track(:,i))%satid)
	where (modulo(trk(track(:,i))%pass,2) == 0) omega = omega + 0.5d0
	xxf%lat = nint(lat(i)*1d6)
	xxf%lon = nint(lon(i)*1d6)
	xxf%time = nint(time(:,i))
	xxf%track = modulo(track(:,i),maxtrk)
	xxf%sla = nint(sla(:,i)*1d6)
	xxf%sla_cor = xxf%sla
	xxf%omega = nint(omega*360d6)
	xxf%sigma = nint(altsig(trk(track(:,i))%satid)*1d3)
	write (10,rec = i+1) xxf
enddo
close (10)
end subroutine write_xxf

!***********************************************************************

subroutine write_xxo
type :: xxo_ ! 44 bytes
	integer(fourbyteint) :: lat, lon, time(2,2)
	integer(twobyteint) :: track(2)
	integer(fourbyteint) :: sla(2), alt(2)
endtype
type(xxo_) :: xxo
integer(fourbyteint) :: i
real(eightbytereal) :: alt(2,x%nr)

! Determine specific for output file and load data
if (id_alt == 0) then
	write (*,'(a)') 'XXO format requires alt variable. Skipped.'
	return
endif

call get_var_2d (id_alt, alt(:,:))

! Open output file and write data
open (10, file=trim(shortname)//'.xxo', status='replace', access='direct', form='unformatted', recl=44)
write (10, rec=1) '@XXO',x%nr,bound,0
do i = 1,x%nr
	xxo%lat = nint(lat(i)*1d6)
	xxo%lon = nint(lon(i)*1d6)
	xxo%time(1,:) = floor(time(:,i))
	xxo%time(2,:) = nint((time(:,i)-xxo%time(1,:))*1d6)
	xxo%track = modulo(track(:,i),maxtrk)
	xxo%sla = nint(sla(:,i)*1d6)
	xxo%alt = nint(alt(:,i)*1d3)
	write (10,rec = i+1) xxo
enddo
close (10)

! Create aux file if needed
if (id_alt_rate > 0 .or. id_sig0 > 0 .or. id_swh > 0) then
	if (id_alt_rate == 0 .or. id_sig0 == 0 .or. id_swh == 0) then
		write (*,'(a)') 'Some of alt_rate, sig0 and swh are selected. Should select all or none. AUX file skipped'
		return
	endif
	call write_aux
endif
end subroutine write_xxo

!***********************************************************************

subroutine write_aux
integer(fourbyteint), parameter :: naux = 3
integer(fourbyteint) :: i
integer(twobyteint) :: aux(2,naux)
real(eightbytereal) :: data(2,x%nr,naux)

! Read auxiliary data
call get_var_2d (id_sig0, data(:,:,1))
call get_var_2d (id_swh, data(:,:,2))
call get_var_2d (id_alt_rate, data(:,:,3))

! Open output file and write data
open (10, file=trim(shortname)//'.xxo.aux', status='replace', access='direct', form='unformatted', recl=4*naux)
write (10, rec=1) '@AUX',x%nr,2*naux
do i = 1,x%nr
	aux = nint(data(:,i,:)*1d3)
	write (10,rec = i+1) aux
enddo
close (10)
end subroutine write_aux

!***********************************************************************

subroutine write_xtf
type :: xtf_
	integer(twobyteint) :: track, satid, nr_xover, nr_alt
	integer(fourbyteint) :: inclination, arglat, equator_time, equator_lon, start_time, end_time
endtype
type(xtf_) :: xtf
integer(fourbyteint) :: par(npar), sig(npar), i
integer(twobyteint) :: flags
par = 0
open (10, file=trim(shortname)//'.xtf', status='replace', access='direct', form='unformatted', recl=34+npar*8)
write (10, rec=1) '@XTB',t%nr,npar,bound,0
do i = 1,t%nr
	if (modulo(trk(i)%pass,2) == 0) then ! Move descending node to ascending node
		trk(i)%equator_lon = trk(i)%equator_lon + 180d0 * (period(trk(i)%satid) / 86400d0 - 1d0)
		trk(i)%equator_time = trk(i)%equator_time - 0.5d0 * period(trk(i)%satid)
	endif
	xtf%track = modulo(i,maxtrk)
	xtf%satid = trk(i)%satid
	xtf%nr_xover = trk(i)%nr_xover
	xtf%nr_alt = trk(i)%nr_alt
	xtf%inclination = nint(inclination(trk(i)%satid) * 1d6)
	xtf%arglat = nint((trk(i)%start_time - trk(i)%equator_time) / period(trk(i)%satid) * 360d6)
	xtf%equator_time = nint(trk(i)%equator_time)
	xtf%equator_lon = nint(trk(i)%equator_lon * 1d6)
	xtf%start_time = nint(trk(i)%start_time)
	xtf%end_time = nint(trk(i)%end_time)
	sig = nint(orberr(trk(i)%satid)*1d6)
	flags = 256 + modulo(trk(i)%pass,2)
	write (10,rec = i+1) xtf,par,sig,flags
enddo
close (10)
end subroutine write_xtf

!***********************************************************************

subroutine write_xaf
end subroutine write_xaf

!***********************************************************************

end program radsxoconv
