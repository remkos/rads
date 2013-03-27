!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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

!*rads_gen_saral -- Converts SARAL data to RADS
!+
program rads_gen_saral

! This program reads SARAL files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/sa/a/sapPPPPcCCC.nc.
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: rads_gen_saral [options] < list_of_SARAL_file_names
!
! This program handles only OGDR, IGDR and GDR files in netCDF format.
!-----------------------------------------------------------------------
!
! Variables array fields to be filled are:
! time - Time since 1 Jan 85
! lat - Latitude
! lon - Longitude
! alt_gdrd - Orbit altitude
! alt_rate - Orbit altitude rate
! range_ka - Ocean range (retracked)
! range_rms_ka - Std dev of range
! range_numval_ka - Nr of averaged range measurements
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_rad - Radiometer wet tropo correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! iono_gim - GIM ionosphetic correction
! iono_nic09 - NIC09 ionospheric correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! tide_solid - Solid earth tide
! tide_ocean_fes04 - FES2008 ocean tide
! tide_ocean_got48 - GOT4.8 ocean tide
! tide_load_fes04 - FES2008 load tide
! tide_load_got48 - GOT4.8 load tide
! tide_pole - Pole tide
! ssb_bm3 - SSB
! mss_cnescls11 - CLS01 MSS
! geoid_egm96 - EGM96 geoid
! swh_ka - Significant wave height
! swh_rms_ka - Std dev of SWH
! sig0_ka - Sigma0
! sig0_rms_ka - Std dev of sigma0
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! wind_speed_alt - Altimeter wind speed
! tb_k - Brightness temperature (K-band)
! tb_ka - Brightness temperature (Ka-band)
! peakiness_ka - Peakiness
! flags - Engineering flags
! off_nadir_angle2_wf_ka - Mispointing from waveform squared
! liquid_water - Liquid water content
! water_vapor_content - Water vapor content
!-----------------------------------------------------------------------
use typesizes
use netcdf
use rads
use rads_misc
use rads_time
use rads_netcdf
use rads_devel

! Command line arguments

integer(fourbyteint) :: verbose=0, c0=0, c1=999, ios, i
real(eightbytereal) :: t0, t1
character(160) :: infile, arg
character(20) :: optopt, optarg
character(80), parameter :: optlist='vC: debug: sat: cycle: t: mjd: sec: ymd: doy:'

! Header variables

integer(fourbyteint) :: cyclenr, passnr, varid
real(eightbytereal) :: equator_time

! Data variables

integer(fourbyteint), parameter :: mrec=3500, mvar=50
integer(fourbyteint) :: nvar, nrec=0, ncid
real(eightbytereal), allocatable :: a(:), b(:), d(:,:)
logical, allocatable :: valid(:,:)
integer(twobyteint), allocatable :: flags(:)
type(rads_sat) :: S
type(rads_pass) :: P
type :: var_
	type(rads_var), pointer :: v ! Pointer to rads_var struct
	real(eightbytereal) :: d(mrec) ! Data array
	logical :: empty ! .true. if all NaN
endtype
type(var_) :: var(mvar)

! Other local variables

real(eightbytereal), parameter :: sec2000=473299200d0	! UTC seconds from 1 Jan 1985 to 1 Jan 2000

! Initialise

t0 = nan
t1 = nan
550 format (a)

! Scan command line for options

call synopsis
do
	call getopt (optlist, optopt, optarg)
	select case (optopt)
	case ('!')
		exit
	case ('v')
		verbose = 1
	case ('debug')
		read (optarg,*) verbose
	case ('C', 'cycle')
		c1 = -1
		read (optarg,*,iostat=ios) c0,c1
		if (c1 < c0) c1 = c0
	case default
		if (.not.dateopt (optopt, optarg, t0, t1)) then
			call synopsis ('--help')
			stop
		endif
	end select
enddo

! Initialise

call synopsis ('--head')
call rads_init (S, 'sa', verbose)

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

files: do
	read (*,550,iostat=ios) infile
	if (ios /= 0) exit files
	write (*,550,advance='no') trim(infile) // ' ...'

	if (nf90_open(infile,nf90_nowrite,ncid) /= nf90_noerr) then
	    write (*,550) 'error opening file'
		cycle files
	endif

! Read global attributes

	call nfs(nf90_inq_dimid(ncid,'time',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
	if (nrec == 0) then
		cycle files
	else if (nrec > mrec) then
		write (*,'("Error: Too many measurements:",i5)') nrec
		cycle files
	endif
	call nfs(nf90_get_att(ncid,nf90_global,'mission_name',arg))
	if (arg /= 'SARAL') then
		write (*,550) 'Error: Wrong misson-name found in header'
		cycle files
	endif

	call nfs(nf90_get_att(ncid,nf90_global,'history',arg))
	call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cyclenr))
	call nfs(nf90_get_att(ncid,nf90_global,'pass_number',passnr))
	call nfs(nf90_get_att(ncid,nf90_global,'equator_time',arg))
	equator_time = strp1985f (arg)

! Skip passes of which the cycle number or equator crossing time is outside the specified interval

	if (equator_time < t0 .or. equator_time > t1 .or. cyclenr < c0 .or. cyclenr > c1) then
		write (*,550) 'Skipped'
		cycle files
	endif

! Store relevant info

	call rads_init_pass_struct (S, P)
	P%cycle = cyclenr
	P%pass = passnr
	P%equator_time = equator_time
	call nfs(nf90_get_att(ncid,nf90_global,'equator_longitude',P%equator_lon))
	call nfs(nf90_get_att(ncid,nf90_global,'first_meas_time',arg))
	P%start_time = strp1985f(arg)
	call nfs(nf90_get_att(ncid,nf90_global,'last_meas_time',arg))
	P%end_time = strp1985f(arg)
	i = index(infile, '/', .true.)
	P%original = infile(i+1:)

! Allocate variables

	allocate (a(nrec),b(nrec),d(40,nrec),valid(40,nrec),flags(nrec))
	nvar = 0

! Compile flag bits

	flags = 0
	call nc2f ('qual_alt_1hz_off_nadir_angle_wf',1)	! bit  1: Quality off-nadir pointing
	call nc2f ('surface_type',2,val=2)			! bit  2: Continental ice
	call nc2f ('surface_type',4,lim=2)			! bit  4: Water/land
	call nc2f ('surface_type',5,lim=1)			! bit  5: Ocean/other
	call nc2f ('rad_surf_type',6,lim=2)			! bit  6: Radiometer land flag
	call nc2f ('ice_flag',7)					! bit  7: Ice flag
	call nc2f ('qual_rad_1hz_tb_k',9)			! bit  9: Quality 23.8 GHz channel
	call nc2f ('qual_rad_1hz_tb_ka',10)			! bit 10: Quality 37.0 GHz channel
	call nc2f ('qual_alt_1hz_range',11)			! bit 11: Quality of range
	call nc2f ('qual_alt_1hz_swh',12)			! bit 12: Quality of SWH
	call nc2f ('qual_alt_1hz_sig0',13)			! bit 13: Quality of sigma0
	call nc2f ('orb_state_flag_diode',15,lim=2)	! bit 15: Quality of DIODE tracking

! Convert all the necessary fields to RADS

	call get_var (ncid, 'time', a)
	call new_var ('time', a + sec2000)
	call cpy_var ('lat', 'lat')
	call cpy_var ('lon', 'lon')
	call cpy_var ('alt', 'alt_gdrd')
	call cpy_var ('orb_alt_rate', 'alt_rate')
	call cpy_var ('range', 'range_ka')
	call cpy_var ('range_numval', 'range_numval_ka')
	call cpy_var ('range_rms', 'range_rms_ka')
	call cpy_var ('model_dry_tropo_corr', 'dry_tropo_ecmwf')
	call cpy_var ('model_wet_tropo_corr', 'wet_tropo_ecmwf')
	call cpy_var ('rad_wet_tropo_corr', 'wet_tropo_rad')
	call cpy_var ('iono_corr_gim', 'iono_gim')
	call cpy_var ('sea_state_bias', 'ssb_bm3')
	call cpy_var ('swh', 'swh_ka')
	call cpy_var ('swh_rms', 'swh_rms_ka')
	call cpy_var ('sig0', 'sig0_ka')
	call cpy_var ('atmos_corr_sig0', 'dsig0_atmos_ka')
	call cpy_var ('off_nadir_angle_wf', 'off_nadir_angle2_wf_ka')
	call cpy_var ('tb_k', 'tb_238')
	call cpy_var ('tb_ka', 'tb_370')
	call cpy_var ('mean_sea_surface', 'mss_cnescls11')
	call cpy_var ('geoid', 'geoid_egm96')
	call cpy_var ('bathymetry', 'topo_dtm2000')
	call cpy_var ('inv_bar_corr', 'inv_bar_static')
	call cpy_var ('inv_bar_corr+hf_fluctuations_corr', 'inv_bar_mog2d')
	call cpy_var ('ocean_tide_sol1-load_tide_sol1', 'tide_ocean_got48')
	call cpy_var ('ocean_tide_sol2-load_tide_sol2', 'tide_ocean_fes04')
	call cpy_var ('load_tide_sol1', 'tide_load_got48')
	call cpy_var ('load_tide_sol2', 'tide_load_fes04')
	call cpy_var ('solid_earth_tide', 'tide_solid')
	call cpy_var ('pole_tide', 'tide_pole')
	call cpy_var ('wind_speed_model_u', 'wind_speed_ecmwf_u')
	call cpy_var ('wind_speed_model_v', 'wind_speed_ecmwf_v')
	call cpy_var ('wind_speed_alt', 'wind_speed_alt')
	call cpy_var ('rad_water_vapor', 'water_vapor_content')
	call cpy_var ('rad_liquid_water', 'liquid_water')
	call get_var (ncid, 'range_used_40hz', d)
	valid = (d == 0d0)
	call get_var (ncid, 'peakiness_40hz', d)
	call mean_1hz (d, valid, a, b)
	call new_var ('peakiness_ka', a)
	a = flags
	call new_var ('flags', a)

! Dump the data

	call put_rads
	deallocate (a, b, d, valid, flags)

enddo files

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Write SARAL data to RADS', flag=flag)) return
write (*,1310)
1310 format (/ &
'syntax: rads_gen_saral [options] < list_of_SARAL_file_names'// &
'This program converts SARAL OGDR/IGDR/GDR files to RADS data'/ &
'files with the name $RADSDATAROOT/data/sa/a/pPPPP/sapPPPPcCCC.nc.'/ &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Copy variable to RADS
!-----------------------------------------------------------------------

subroutine cpy_var (varin, varout)
character(len=*), intent(in) :: varin, varout
call get_var (ncid, varin, a)
call new_var (varout, a)
end subroutine cpy_var

!-----------------------------------------------------------------------
! Create new RADS variable
!-----------------------------------------------------------------------

subroutine new_var (varnm, data)
! Write variables one after the other to the output file
character(len=*), intent(in) :: varnm
real(eightbytereal), intent(in) :: data(:)
nvar = nvar + 1
if (nvar > mvar) stop 'Too many variables'
var(nvar)%v => rads_varptr (S, varnm)
var(nvar)%d = data
end subroutine new_var

!-----------------------------------------------------------------------
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads
integer(fourbyteint) :: i

! Check which variables are empty
do i = 1,nvar
	var(i)%empty = all(isnan_(var(i)%d(1:nrec)))
enddo
if (any(var(1:nvar)%empty)) then
	write (*,550,advance='no') '... No'
	do i = 1,nvar
		if (var(i)%empty) write (*,550,advance='no') trim(var(i)%v%name)
	enddo
endif

! Open output file
call rads_create_pass (S, P, nrec)

! Define all variables
do i = 1,nvar
	call rads_def_var (S, P, var(i)%v)
enddo

! Fill all the data fields
do i = 1,nvar
	call rads_put_var (S, P, var(i)%v, var(i)%d(1:nrec))
enddo

! Close the data file
write (*,552) nrec,trim(P%filename(len_trim(S%dataroot)+2:))
call rads_close_pass (S, P)

! Formats
550 format (a,1x)
552 format ('...',i5,' records written to ',a)

end subroutine put_rads

!-----------------------------------------------------------------------
! nc2f: Load flag field, then set corresponding bit in RADS
! varnm : source variable name
! bit   : RADS bit to be set when value is val
! lim   : set bit when value >= lim (optional)
! val   : set bit when value == val (optional, default = 1)
! neq   : set bit when value /= val (optional)

subroutine nc2f (varnm, bit, lim, val, neq)
character(*), intent(in) :: varnm
integer(fourbyteint), intent(in) :: bit
integer(fourbyteint), optional, intent(in) :: lim,val,neq
integer(twobyteint) :: flag(mrec),flag2d(1:1,1:nrec)
integer(fourbyteint) :: i,ival,ndims

if (nf90_inq_varid(ncid,varnm,varid) /= nf90_noerr) then
	write (*,'("No such variable: ",a)') trim(varnm)
	return
endif
call nfs(nf90_inquire_variable(ncid,varid,ndims=ndims))
if (ndims == 2) then
	call nfs(nf90_get_var(ncid,varid,flag2d(1:1,1:nrec)))
	flag(1:nrec) = flag2d(1,1:nrec)
else
	call nfs(nf90_get_var(ncid,varid,flag(1:nrec)))
endif
if (present(lim)) then
	do i = 1,nrec
		if (flag(i) >= lim) flags(i) = ibset(flags(i),bit)
	enddo
else if (present(neq)) then
	do i = 1,nrec
		if (flag(i) /= neq) flags(i) = ibset(flags(i),bit)
	enddo
else
	ival = 1
	if (present(val)) ival = val
	do i = 1,nrec
		if (flag(i) == ival) flags(i) = ibset(flags(i),bit)
	enddo
endif
end subroutine nc2f

end program rads_gen_saral
