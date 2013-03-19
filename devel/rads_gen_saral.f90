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

! This program reads REAPER files and converts them to the RADS format,
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
! alt_reaper - Orbit altitude
! alt_rate - Orbit altitude rate
! range_ku - Ocean range (retracked)
! dry_tropo_ecmwf - ECMWF dry tropospheric correction
! wet_tropo_rad - Radiometer wet tropo correction
! wet_tropo_ecmwf - ECMWF wet tropo correction
! iono_gim - GIM ionosphetic correction
! iono_nic09 - NIC09 ionospheric correction
! inv_bar_static - Inverse barometer
! inv_bar_mog2d - MOG2D
! tide_solid - Solid earth tide
! tide_ocean_fes04 - FES2008 ocean tide
! tide_ocean_got47 - GOT4.7 ocean tide
! tide_load_fes04 - FES2008 load tide
! tide_load_got47 - GOT4.7 load tide
! tide_pole - Pole tide
! ssb_bm3 - SSB
! mss_cls01 - CLS01 MSS
! geoid_egm2008 - EGM2008 geoid
! mss_ucl04 - UCL04 MSS
! swh_ku - Significant wave height
! sig0_ku - Sigma0
! wind_speed_ecmwf_u - ECMWF wind speed (U)
! wind_speed_ecmwf_v - ECMWF wind speed (V)
! range_rms_ku - Std dev of range
! range_numval_ku - Nr of averaged range measurements
! topo_macess - MACESS topography
! tb_238 - Brightness temperature (23.8 GHz)
! tb_365 - Brightness temperature (36.5 GHz)
! peakiness_ku - Peakiness
! flags - Engineering flags
! swh_rms_ku - Std dev of SWH
! sig0_rms_ku - Std dev of sigma0
! off_nadir_angle2_wf_ku - Mispointing from waveform squared
! liquid_water - Liquid water content
! water_vapor_content - Water vapor content
! tide_equil - Long-period equilibrium tide
! tide_non_equil - Long-period non-equilibrium tide
!-----------------------------------------------------------------------
use typesizes
use netcdf
use rads
use rads_misc
use rads_time
use rads_netcdf
use rads_devel

! Command line arguments

integer(fourbyteint) :: verbose=0, c0=0, c1=999, ios
real(eightbytereal) :: t0, t1
character(160) :: infile, old_infile
character(20) :: optopt, optarg
character(80), parameter :: optlist='vC: debug: sat: cycle: t: mjd: sec: ymd: doy:'

! Header variables

character(1) :: phasenm(2)
character(80) :: l2_proc_time, l2_version
logical :: meteo
real(eightbytereal) :: tnode(2), lnode(2)
integer(fourbyteint) :: orbitnr(2), cyclenr(2), passnr(2), varid

! Data variables

integer(fourbyteint), parameter :: mrec=15000, mvar=50
integer(fourbyteint) :: nvar, ndata=0, nrec=0, nout=0, ncid, ers=0
real(eightbytereal) :: start_time, end_time
real(eightbytereal), allocatable :: a(:), b(:), c(:), d(:,:), dh(:), sum_c_applied(:), sum_d_applied(:)
integer(twobyteint), allocatable :: flags(:)
integer(fourbyteint), allocatable :: f_error(:), f_applied(:)
logical, allocatable :: valid(:,:)
type(rads_sat) :: S
type(rads_pass) :: P
type :: var_
	type(rads_var), pointer :: v ! Pointer to rads_var struct
	real(eightbytereal) :: d(mrec) ! Data array
	logical :: empty ! .true. if all NaN
endtype
type(var_) :: var(mvar)

! Other local variables

real(eightbytereal), parameter :: sec1990=157766400d0	! UTC seconds from 1 Jan 1985 to 1 Jan 1990
real(eightbytereal), parameter :: picosec_to_m=0.5d-12*299792458d0	! picoseconds of 2-way range to mm 1-way
integer :: i
logical :: new

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

!----------------------------------------------------------------------
! Read all file names from standard input
!----------------------------------------------------------------------

! Start reading with at least first file

read (*,550,iostat=ios) infile
if (ios /= 0) then
	call synopsis ('--help')
else
	call synopsis ('--head')
endif
call get_reaper

files: do
	read (*,550,iostat=ios) infile
	if (ios /= 0) exit files
	write (*,551) trim(arg)

	if (nf90_open(arg,nf90_nowrite,ncid) /= nf90_noerr) then
	    write (*,550) 'error opening file'
		cycle files
	endif

! Read global attributes

	call nfs(nf90_inq_dimid(ncid,'time',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=datanr))
	if (datanr > maxmeas) then
		write (*,'("Error: Too many measurements:",i5)') datanr
		cycle files
	endif
	product=arg
	call nfs(nf90_get_att(ncid,nf90_global,'title',arg))
	ogdr = (arg(:4) == 'OGDR')
	call nfs(nf90_get_att(ncid,nf90_global,'mission_name',arg))
	if (arg /= 'SARAL') then
		write (*,550) 'Error: Wrong misson-name found in header'
		cycle files
	endif

	call nfs(nf90_get_att(ncid,nf90_global,'history',arg))
	createdate = arg(:19)
	call nfs(nf90_get_att(ncid,nf90_global,'cycle_number',cycnr))
	call nfs(nf90_get_att(ncid,nf90_global,'pass_number',passnr))
	call nfs(nf90_get_att(ncid,nf90_global,'equator_longitude',eq_lon))
	call nfs(nf90_get_att(ncid,nf90_global,'equator_time',arg))
	eq_time = strp1985f (arg)

! Skip passes of which the cycle number or equator crossing time is outside the specified interval

	if (eq_time < t0 .or. eq_time > t1 .or. cycnr < c0 .or. cycnr > c1) then
		write (*,550) 'Skipped'
		cycle files
	endif

! Store relevant info

	call rads_init_pass_struct (S, P)
	P%cycle = cyclenr
	P%pass = passnr
	P%equator_time = eq_time
	P%equator_lon = eq_lon
	call nfs(nf90_get_att(ncid,nf90_global,'first_meas_time',arg))
	P%start_time = strp1985f (arg)
	call nfs(nf90_get_att(ncid,nf90_global,'last_meas_time',arg))
	P%end_time = strp1985f (arg)

! Compile flag bits

	f2601 = 0
	call nc2f('alt_state_flag_oper',0)			! bit  0: Altimeter Side A/B
	if (.not.gdr_d) call nc2f('qual_alt_1hz_off_nadir_angle_wf_ku',1)	! bit  1: Quality off-nadir pointing
	call nc2f('surface_type',2,val=2)			! bit  2: Continental ice
	call nc2f('qual_alt_1hz_range_c',3)			! bit  3: Quality dual-frequency iono
	call nc2f('surface_type',4,lim=2)			! bit  4: Water/land
	call nc2f('surface_type',5,lim=1)			! bit  5: Ocean/other
	if (gdr_d) then
		call nc2f('rad_surf_type',6,lim=2)		! bit  6: Radiometer land flag
	else
		call nc2f('rad_surf_type',6)			! bit  6: Radiometer land flag
	endif
	call nc2f('ice_flag',7)						! bit  7: Ice flag
	call nc2f('rain_flag',8)					! bit  8: Rain flag
	call nc2f('qual_rad_1hz_tb187',9)			! bit  9: Quality 18.7 and 23.8 GHz channel
	call nc2f('qual_rad_1hz_tb238',9)
	call nc2f('qual_rad_1hz_tb340',10)			! bit 10: Quality 34.0 GHz channel

! Convert all the necessary fields to RADS

	call nc2r( 101,'time','time [seconds since 1985-01-01]',sec2000)
	t_start = array(1) + sec2000
	t_end = array(datanr) + sec2000
	call nc2r( 201,'lat','latitude [degrees_north]')
	call nc2r( 301,'lon','longitude [degrees_east]')
	call nc2r( 418,'alt','CNES GDR-D orbital altitude [m]')
	call nc2r( 501,'orb_alt_rate','orbital altitude rate [m/s]')
	call nc2r( 601,'range_ku','altimeter range corrected for instr. effects (Ku) [m]')
	if (gdr_d) call nc2r( 611,'range_ku_mle3','altimeter range corrected for instr. effects (Ku MLE3) [m]')
	call nc2r( 602,'range_c','altimeter range corrected for instr. effects (C) [m]')
	call nc2r( 701,'model_dry_tropo_corr','ECMWF dry tropospheric correction [m]')
	call nc2r( 801,'rad_wet_tropo_corr','AMR wet tropospheric correction [m]')
	call nc2r( 802,'model_wet_tropo_corr','ECMWF model wet tropospheric correction [m]')
	call nc2r( 901,'iono_corr_alt_ku','dual-frequency ionospheric correction (Ku) [m]')
	if (gdr_d) call nc2r( 911,'iono_corr_alt_ku_mle3','dual-frequency ionospheric correction (Ku MLE3) [m]')
	call nc2r( 906,'iono_corr_gim_ku','JPL GIM ionospheric correction [m]')
	call nc2r(1002,'inv_bar_corr','local+global inverse barometer correction [m]')
	if (.not.ogdr) call nc2r(1004,'inv_bar_corr','MOG2D total inverse barometer correction [m]', &
	varnm2='+hf_fluctuations_corr')
!	call nc2r(1101,'solid_earth_tide','solid earth tide [m]')
!	call nc2r(1207,'ocean_tide_sol1','GOT00.2 ocean tide [m]',varnm2='-load_tide_sol1')
!	call nc2r(1263,'ocean_tide_sol2','FES2004 ocean tide (GDR) [m]',varnm2='-load_tide_sol2')
!	call nc2r(1307,'load_tide_sol1','GOT00.2 load tide [m]')
!	call nc2r(1363,'load_tide_sol2','FES2004 load tide (GDR) [m]')
!	call nc2r(1401,'pole_tide','pole tide [m]')
	call nc2r(1502,'sea_state_bias_ku','CLS sea state bias (Ku) [m]')
	if (gdr_d) call nc2r(1512,'sea_state_bias_ku_mle3','CLS sea state bias (Ku MLE3) [m]')
	call nc2r(1505,'sea_state_bias_c','CLS sea state bias (C) [m]')
!	call nc2r(1601,'geoid','EGM96 geoid [m]')
!	call nc2r(1605,'mean_sea_surface','CLS01 mean sea surface [m]')
	call nc2r(1701,'swh_ku','significant wave height (Ku) [m]')
	if (gdr_d) call nc2r(1711,'swh_ku_mle3','significant wave height (Ku MLE3) [m]')
	call nc2r(1702,'swh_c','significant wave height (C) [m]')
	call nc2r(1801,'sig0_ku','backscatter coefficient (Ku) [dB]')
	if (gdr_d) call nc2r(1811,'sig0_ku_mle3','backscatter coefficient (Ku MLE3) [dB]')
	call nc2r(1802,'sig0_c','backscatter coefficient (C) [dB]')
	call nc2r(1901,'wind_speed_alt','altimeter wind speed [m/s]')
	if (gdr_d) call nc2r(1911,'wind_speed_alt_mle3','altimeter wind speed (MLE3) [m/s]')
	call nc2r(1902,'wind_speed_rad','radiometer wind speed [m/s]')
	call nc2r(1903,'wind_speed_model_u','U-component of model wind speed [m/s]')
	call nc2r(1904,'wind_speed_model_v','V-component of model wind speed [m/s]')
	call nc2r(2002,'range_rms_ku','std dev of range (20-Hz, Ku) [m]')
	if (gdr_d) call nc2r(2012,'range_rms_ku_mle3','std dev of range (20-Hz, Ku MLE3) [m]')
	call nc2r(2004,'range_rms_c','std dev of range (20-Hz, C) [m]')
	call nc2r(2101,'range_numval_ku','number of averaged 20-Hz range measurements (Ku) [count]')
	if (gdr_d) call nc2r(2111,'range_numval_ku_mle3','number of averaged 20-Hz range measurements (Ku MLE3) [count]')
	call nc2r(2102,'range_numval_c','number of averaged 20-Hz range measurements (C) [count]')
	call nc2r(2202,'bathymetry','DTM2000.1 topography [m]')
	call nc2r(2301,'tb_187','brightness temperature at 18.7 GHz [K]')
	call nc2r(2302,'tb_238','brightness temperature at 23.8 GHz [K]')
	call nc2r(2303,'tb_340','brightness temperature at 34.0 GHz [K]')
	call putraw('U',2601,datanr,f2601,2,1d0,0d0,2,1d0,0d0,'engineering flags')
	if (gdr_d) call putraw('U',2611,datanr,f2611,2,1d0,0d0,2,1d0,0d0,'engineering flags (MLE3)')
	call nc2r(2802,'swh_rms_ku','std dev of significant wave height (20-Hz, Ku) [m]')
	if (gdr_d) call nc2r(2812,'swh_rms_ku_mle3','std dev of significant wave height (20-Hz, Ku MLE3) [m]')
	call nc2r(2804,'swh_rms_c','std dev of significant wave height (20-Hz, C) [m]')
	call nc2r(2902,'sig0_rms_ku','std dev of backscatter coefficient (20-Hz, Ku) [dB]')
	if (gdr_d) call nc2r(2912,'sig0_rms_ku_mle3','std dev of backscatter coefficient (20-Hz, Ku MLE3) [dB]')
	call nc2r(2904,'sig0_rms_c','std dev of backscatter coefficient (20-Hz, C) [dB]')
	call nc2r(3002,'off_nadir_angle_wf_ku','off-nadir angle squared from waveforms (Ku) [degrees^2]')
	if (.not.gdr_d) call nc2r(3004,'off_nadir_angle_pf','off-nadir angle squared from platform [degrees^2]')
	call nc2r(3203,'atmos_corr_sig0_ku','backscatter coefficient atm. attenuation (added, Ku) [dB]')
	call nc2r(3204,'atmos_corr_sig0_c','backscatter coefficient atm. attenuation (added, C) [dB]')
	call nc2r(3301,'rad_liquid_water','liquid water content [kg/m^2]')
!   call nc2r(3901,'ocean_tide_equil','long period equilibrium tide [m]')
!   call nc2r(3902,'ocean_tide_non_equil','long period non-equilibrium tide [m]')
	! Write out the data
	nout = i - 1 ! Number of measurements to be written out
	call put_rads

	! Number of measurements remaining
	ndata = ndata - nout
	if (ios /= 0 .and. ndata == 0) exit ! We are out of data

	! Move the data to be beginning
	do i = 1,nvar
		var(i)%d(1:ndata) = var(i)%d(nout+1:nout+ndata)
	enddo
enddo

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
! Write content of memory to a single pass of RADS data
!-----------------------------------------------------------------------

subroutine put_rads
integer(fourbyteint) :: i
character(160) :: original

if (nout == 0) return	! Skip empty data sets
if (cyclenr(1) < c0 .or. cyclenr(1) > c1) return	! Skip chunks that are not of the selected cycle
if (tnode(1) < t0 .or. tnode(1) > t1) return	! Skip equator times that are not of selected range

! Update phase name if required
phasenm(1) = strtolower(phasenm(1))
if (S%phase%name /= phasenm(1)) S%phase => rads_get_phase (S, phasenm(1)//'.r')

! Store relevant info
call rads_init_pass_struct (S, P)
P%cycle = cyclenr
P%pass = passnr
P%start_time = t0
P%end_time = t1
P%equator_time = tnode(1)
P%equator_lon = lnode(1)

! Check which input files pertain
if (P%start_time >= start_time) then
	original = infile
else if (P%end_time < start_time) then
	original = old_infile
else
	original = trim(old_infile) //' '// infile
endif
P%original = trim(l2_version)//' data of '//l2_proc_time(:11)//': '//trim(original)

! Check which variables are empty
do i = 1,nvar
	var(i)%empty = all(isnan_(var(i)%d(1:nout)))
enddo
if (any(var(1:nvar)%empty)) then
	write (*,550,advance='no') '... No'
	do i = 1,nvar
		if (var(i)%empty) write (*,550,advance='no') trim(var(i)%v%name)
	enddo
endif

! Open output file
call rads_create_pass (S, P, nout)

! Define all variables
do i = 1,nvar
	call rads_def_var (S, P, var(i)%v)
enddo

! Fill all the data fields
do i = 1,nvar
	call rads_put_var (S, P, var(i)%v, var(i)%d(1:nout))
enddo

! Close the data file
write (*,552) nout,trim(P%filename(len_trim(S%dataroot)+2:))
call rads_close_pass (S, P)

! Formats
550 format (a,1x)
552 format ('...',i5,' records written to ',a)

end subroutine put_rads

!-----------------------------------------------------------------------
! nc2r: Load variable from source, then dump into RADS

subroutine nc2r (in, out)
character(len=*), intent(in) :: in
character(len=*), intent(out) :: out
integer(fourbyteint) :: varid,xtype,ntype,constant,i0,i1,l,ndims
real(eightbytereal) :: scale_factor,add_offset,array2(maxmeas),array2d(1,maxmeas),invalid

if (nf90_inq_varid(ncid,in,varid) /= nf90_noerr) then
	write (*,'("No such variable: ",a)') trim(in)
	return
endif

! If offset is specified, overrule add_offset from the file.

if (nf90_get_att(ncid,varid,'add_offset',add_offset) /= nf90_noerr) add_offset = 0d0
if (nf90_get_att(ncid,varid,'scale_factor',scale_factor) /= nf90_noerr) scale_factor = 1d0
if (present(offset)) add_offset = offset

! Convert netCDF type number to RADS type number

call nfs(nf90_inquire_variable(ncid,varid,xtype=xtype,ndims=ndims))
select case (xtype)
case (nf90_int1)
	ntype = 1
	invalid = huge(0_onebyteint)
case (nf90_int2)
	ntype = 2
	invalid = huge(0_twobyteint)
case (nf90_int4)
	ntype = 4
	invalid = huge(0_fourbyteint)
case (nf90_double)
	ntype = 8
	invalid = huge(0_eightbytereal)
case default
	write (*,'("Unknown type in ",a)') trim(varnm)
	return
end select

! Get the data array

if (ndims == 2) then
	call nfs(nf90_get_var(ncid,varid,array2d(1:1,1:datanr)))
	array(1:datanr) = array2d(1,1:datanr)
else
	call nfs(nf90_get_var(ncid,varid,array(1:datanr)))
endif

! If varnm2 is specified, add it to or subtract it from the current array

if (present(varnm2)) then
	i1 = 1
	l = len_trim(varnm2)
	do
		if (i1 > l) exit
		i0 = i1
		i1 = scan(varnm2(i0+1:),'+-') + i0
		if (i1 == i0) i1 = l + 1
		if (nf90_inq_varid(ncid,varnm2(i0+1:i1-1),varid) /= nf90_noerr) then
			write (*,'("No such variable: ",a)') trim(varnm2(2:))
			return
		endif
		if (ndims == 2) then
			call nfs(nf90_get_var(ncid,varid,array2d(1:1,1:datanr)))
			array2(1:datanr) = array2d(1,1:datanr)
		else
			call nfs(nf90_get_var(ncid,varid,array2(1:datanr)))
		endif
		constant = 0
		if (varnm2(i0:i0) == '-') constant = -1
		if (varnm2(i0:i0) == '+') constant = 1
		where (array == invalid)
		else where (array2 == huge(0_twobyteint))
		else where
			array = array + constant * array2
		end where
	enddo
endif

! For field 3004, we check if all may be just zero. If so, skip it.

if (sel /= 3004) then
else if (all(array(1:datanr) == 0d0)) then
	return
else
	write (*,551) 'warning: non-zero off-nadir angles found'
endif
551 format (a,' ...',$)

! Store the array

call putraw('U',sel,datanr,array,8,scale_factor,add_offset,ntype,scale_factor,add_offset,longname)
end subroutine nc2r

!-----------------------------------------------------------------------
! nc2f: Load flag field, then set corresponding bit in RADS
! varnm : source variable name
! bit   : RADS bit to be set when value is val
! lim   : set bit when value >= lim (optional)
! val   : set bit when value == val (optional, default = 1)
! neq   : set bit when value /= val (optional)

subroutine nc2f(varnm,bit,lim,val,neq)
character(*), intent(in) :: varnm
integer(fourbyteint), intent(in) :: bit
integer(fourbyteint), optional, intent(in) :: lim,val,neq
integer(twobyteint) :: flag(maxmeas),flag2d(1:1,1:datanr)
integer(fourbyteint) :: i,ival,ndims

if (nf90_inq_varid(ncid,varnm,varid) /= nf90_noerr) then
	write (*,'("No such variable :",a)') trim(varnm)
	return
endif
call nfs(nf90_inquire_variable(ncid,varid,ndims=ndims))
if (ndims == 2) then
	call nfs(nf90_get_var(ncid,varid,flag2d(1:1,1:datanr)))
	flag(1:datanr) = flag2d(1,1:datanr)
else
	call nfs(nf90_get_var(ncid,varid,flag(1:datanr)))
endif
if (present(lim)) then
	do i = 1,datanr
		if (flag(i) >= lim) f2601(i) = ibset(f2601(i),bit)
	enddo
else if (present(neq)) then
	do i = 1,datanr
		if (flag(i) /= neq) f2601(i) = ibset(f2601(i),bit)
	enddo
else
	ival = 1
	if (present(val)) ival = val
	do i = 1,datanr
		if (flag(i) == ival) f2601(i) = ibset(f2601(i),bit)
	enddo
endif
end subroutine nc2f

end program rads_gen_saral
