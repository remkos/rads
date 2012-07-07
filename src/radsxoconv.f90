!-----------------------------------------------------------------------
! $Id$
!
! Copyright (C) 2012  Remko Scharroo (Altimetrics LLC)
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

!*radsxoconv -- RADS crossover converter
!+
program radsxoconv

! This program converts RADS netCDF crossover files into the historical
! DEOS XAF, XTF and XXF/XXO format.
!-----------------------------------------------------------------------
use rads
use netcdf
use rads_netcdf
use rads_misc
integer(fourbyteint) :: ncid
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
character(len=rads_varl) :: optopt
character(len=rads_naml) :: optarg, shortname
integer(fourbyteint), parameter :: msat = 20
type :: sat_
	character(len=2) :: name
	real(eightbytereal) :: period, altsig, orberr, inclination
endtype sat_
type(sat_) :: sat(msat)
character(len=3*msat) :: satlist
character(len=*), parameter :: optlist = 'x:p:'
type(rads_sat) :: S
integer(twobyteint), parameter :: bound(4) = int((/-180,180,-90,90/), twobyteint)
integer(fourbyteint), parameter :: maxtrk = 32768
logical :: l

! Initialize RADS or issue help
call synopsis
a = prod_ ('', 0)
t = prod_ ('', 0)
x = prod_ ('', 0)

! Start with this-is message
l = rads_version ('$Revision$')

! Scan command line arguments
do
	call getopt (optlist, optopt, optarg)
	select case (optopt)
	case ('!') ! End of options
		exit
	case ('x')
		select case (optarg)
		case ('af')
			t%fmt = 'a'
		case ('tf')
			t%fmt = 't'
		case ('xo')
			x%fmt = 'o'
		case ('xf')
			x%fmt = 'x'
		end select
	case ('p')
		read (optarg, *, iostat=ios) npar
	case (' ')
		call process (optarg) ! Process each file
	end select
enddo

contains

!***********************************************************************

subroutine process (filename)
character(len=*), intent(in) :: filename
character(len=rads_varl) :: varname
integer(fourbyteint) :: i, nvars

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
	call nfs (nf90_inquire_variable (ncid, i, name=varname))
	if (varname == 'sla') then
		id_sla = i
	else if (varname == 'alt_rate') then
		id_alt_rate = i
	else if (varname == 'sig0') then
		id_sig0 = i
	else if (varname == 'swh') then
		id_swh = i
	else if (varname == 'alt' .or. varname(:4) == 'alt_') then
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

! Get essential satellite information
! Older files only have IDs, newer have the satellite abbreviations
if (nf90_get_att (ncid, get_varid('satid'), 'flag_meanings', satlist) == nf90_noerr) then
	do i = 1,len_trim(satlist)/3+1
		call rads_init (S, satlist(i*3-2:i*3-1))
		sat(S%satid) = sat_ (S%sat, 2*S%phase%pass_seconds, &
			S%xover_params(1), S%xover_params(2), S%inclination)
	enddo
else
	satlist = 'g3 ss gs e1 tx pn e2 g1 j1 n1 j2 c2'
	do i = 1,12
		if (.not.any(trk%satid == i)) cycle
		call rads_init (S, satlist(i*3-2:i*3-1))
		sat(S%satid) = sat_ (S%sat, 2*S%phase%pass_seconds, &
			S%xover_params(1), S%xover_params(2), S%inclination)
	enddo
endif

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

subroutine synopsis ()
if (rads_version ('$Revision$','RADS crossover file converter')) return
write (stderr,1300)
1300 format (/ &
'usage: radsxoconv [options] FILENAME ...' // &
'Required argument:' / &
'  FILENAME                  Name of input netCDF xover file'// &
'Optional arguments [options] are:'/ &
'  -xaf                      Create XAF file (extension .xaf)'/ &
'  -xtf                      Create XTF file (extension .xtf)'/ &
'  -xxf                      Create XXF file (extension .xxf)'/ &
'  -xxo                      Create XXO file (extension .xxo, .xxo.aux)'// &
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
	omega = (time(:,i) - trk(track(:,i))%equator_time) / sat(trk(track(:,i))%satid)%period
	where (modulo(trk(track(:,i))%pass,2) == 0) omega = omega + 0.5d0
	xxf%lat = nint(lat(i)*1d6)
	xxf%lon = nint(lon(i)*1d6)
	xxf%time = nint(time(:,i))
	xxf%track = int(modulo(track(:,i),maxtrk), twobyteint)
	xxf%sla = nint(sla(:,i)*1d6)
	xxf%sla_cor = xxf%sla
	xxf%omega = nint(omega*360d6)
	xxf%sigma = nint(sat(trk(track(:,i))%satid)%altsig*1d3, twobyteint)
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
	xxo%track = int(modulo(track(:,i),maxtrk), twobyteint)
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
	aux = nint(data(:,i,:)*1d3, twobyteint)
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
		trk(i)%equator_lon = trk(i)%equator_lon + 180d0 * (sat(trk(i)%satid)%period / 86400d0 - 1d0)
		trk(i)%equator_time = trk(i)%equator_time - 0.5d0 * sat(trk(i)%satid)%period
	endif
	xtf%track = int(modulo(i,maxtrk), twobyteint)
	xtf%satid = trk(i)%satid
	xtf%nr_xover = trk(i)%nr_xover
	xtf%nr_alt = trk(i)%nr_alt
	xtf%inclination = nint(sat(trk(i)%satid)%inclination * 1d6)
	xtf%arglat = nint((trk(i)%start_time - trk(i)%equator_time) / sat(trk(i)%satid)%period * 360d6)
	xtf%equator_time = nint(trk(i)%equator_time)
	xtf%equator_lon = nint(trk(i)%equator_lon * 1d6)
	xtf%start_time = nint(trk(i)%start_time)
	xtf%end_time = nint(trk(i)%end_time)
	sig = nint(sat(trk(i)%satid)%orberr*1d6)
	flags = int(256 + modulo(trk(i)%pass,2), twobyteint)
	write (10,rec = i+1) xtf,par,sig,flags
enddo
close (10)
end subroutine write_xtf

!***********************************************************************

subroutine write_xaf
write (*, 600)
600 format ('radsxoconv: writing of XAF files not implemented')
end subroutine write_xaf

end program radsxoconv
