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

!*rads_gen_reaper -- Converts REAPER ERS-1/2 data to RADS
!+
program rads_gen_reaper

! This program reads REAPER files and converts them to the RADS format,
! written into files $RADSDATAROOT/data/eE/F.r/eEpPPPPcCCC.nc.
!     E = 1 or 2
!     F = mission phase
!  PPPP = relative pass number
!   CCC = cycle number
!
! syntax: reaperraw [options] < list_of_REAPER_file_names
!
! This program handles only the REAPER ALT_2M files in netCDF format.
!-----------------------------------------------------------------------
!
!  Variables array fields to be filled are:
!  101 - Time since 1 Jan 85
!  201 - Latitude
!  301 - Longitude
!  404 - Orbit altitude
!  501 - Orbit altitude rate
!  601 - Ocean range (retracked)
!  701 - ECMWF dry tropospheric correction
!  801 - Radiometer wet tropo correction
!  802 - ECMWF wet tropo correction
!  906 - GIM ionosphetic correction
!  908 - NIC09 ionospheric correction
! 1002 - Inverse barometer
! 1004 - MOG2D
! 1101 - Solid earth tide
! 1213 - FES2008 ocean tide
! 1217 - GOT4.7 ocean tide
! 1313 - FES2008 load tide
! 1317 - GOT4.7 load tide
! 1401 - Pole tide
! 1501 - SSB
! 1605 - CLS01 MSS
! 1610 - EGM2008 geoid
! 1613 - UCL04 MSS
! 1701 - Significant wave height
! 1801 - Sigma0
! 1901 - Altimeter wind speed
! 1903 - ECMWF wind speed (U)
! 1904 - ECMWF wind speed (V)
! 2002 - Std dev of range
! 2101 - Nr of averaged range measurements
! 2206 - MACESS topography
! 2302 - Brightness temperature (23.8 GHz)
! 2303 - Brightness temperature (36.5 GHz)
! 2401 - Peakiness (ALT_2M only)
! 2601 - Engineering flags
! 2702 - Internal calibration correction to range (appied) ((S)GDR only)
! 2704 - Doppler correction (applied) (ALT_2_/ALT_2S only)
! 2802 - Std dev of SWH
! 2902 - Std dev of sigma0
! 3002 - Mispointing from waveform squared
! 3203 - Sigma0 attenuation
! 3301 - Liquid water content
! 3302 - Water vapor content
! 3409 - Mean quadratic error of waveform fit (ALT_2_/ALT_2S only)
! 3901 - Long-period equilibrium tide
! 3902 - Long-period non-equilibrium tide
!-----------------------------------------------------------------------
use typesizes
use netcdf
use rads
use rads_misc
use rads_time
use rads_netcdf

! Command line arguments

integer(fourbyteint) :: verbose=0,c0=0,c1=999,ios
real(eightbytereal) :: t0, t1
character(160) :: arg
character(20) :: sat='',newsat='',oldsat='',optopt,optarg
character(80), parameter :: optlist='v debug: sat: cycle: t: mjd: sec: ymd: doy:'

! Header variables

character(2) :: mission
character(1) :: phasenm(2)
character(80) :: l2_proc_time, l2_version
logical :: alt_2m
real(eightbytereal) :: tnode(2),lnode(2)
integer(fourbyteint) :: datanr=0,nrec,ncid,varid,orbitnr(2),cyclenr(2),passnr(2),ers

! Data variables

integer(fourbyteint), parameter :: mrec = 9000
real(eightbytereal) :: f0101(mrec), f0201(mrec), f0301(mrec), f0404(mrec), f0501(mrec), f0601(mrec), &
	f0701(mrec), f0801(mrec), f0802(mrec), f0906(mrec), f0908(mrec), f1002(mrec), f1004(mrec), &
	f1101(mrec), f1213(mrec), f1217(mrec), f1313(mrec), f1317(mrec), f1401(mrec), f1501(mrec), &
	f1605(mrec), f1610(mrec), f1613(mrec), f1701(mrec), f1801(mrec), f1901(mrec), f1903(mrec), &
	f1904(mrec), f2002(mrec), f2101(mrec), f2206(mrec), f2302(mrec), f2303(mrec), f2401(mrec), &
	f2702(mrec), f2704(mrec), f2802(mrec), f2902(mrec), f3002(mrec), f3203(mrec), &
	f3301(mrec), f3302(mrec), f3409(mrec), f3901(mrec), f3902(mrec)
integer(twobyteint) :: f2601(mrec)

! Other local variables

real(eightbytereal), parameter :: sec1990=157766400d0
real(eightbytereal) :: nan, sum_d_applied
real(eightbytereal), allocatable :: a(:),b(:),c(:),d(:,:),sum_c_applied(:)
logical, allocatable :: valid(:,:)
integer(fourbyteint) :: i,k,i0,i1,flag
logical :: new, erspass

! Initialise

nan = 0d0
nan = nan / nan
t0 = nan
t1 = nan
550 format (a)
551 format (a,' ...')

! Scan command line for options

call synopsis ()
do
	call getopt (optlist, optopt, optarg)
	select case (optopt)
	case ('!')
		exit
	case ('v')
		verbose = 1
	case ('debug')
		read (optarg,*) verbose
	case ('sat')
		sat = optarg
	case ('cycle')
		c1 = -1
		read (optarg,*,iostat=ios) c0,c1
		if (c1 < c0) c1 = c0
	case default
		if (.not.dateopt (optopt, optarg, t0, t1)) then
			call synopsis ()
			stop
		endif
	end select
enddo

!----------------------------------------------------------------------
! Read all file names for standard input
!----------------------------------------------------------------------

files: do
	read (*,'(a)',iostat=ios) arg
	if (ios /= 0) exit files

! Check input file name

	i = index(arg,'ERS_ALT_2')
	if (i <= 0) then
		write (*,550) 'Error: Wrong input file'
		cycle files
	endif
	alt_2m = index(arg,'ERS_ALT_2M') > 0

! Open input file

	write (*,551) trim(arg)
	if (nf90_open(arg,nf90_nowrite,ncid) /= nf90_noerr) then
		write (*,550) 'Error opening file'
		cycle files
	endif

! Check for ERS-1 or -2
! Do not trust 'mission' attribute. It is always 'E1'.

	if (index(arg,'E1_') > 0) then
		ers = 1
		mission = 'e1'
	else
		ers = 2
		mission = 'e2'
	endif

! Read header records

	call nfs(nf90_inq_dimid(ncid,'Record',varid))
	call nfs(nf90_inquire_dimension(ncid,varid,len=nrec))
	if (nrec > mrec) then
		write (*,'("Error: Too many measurements:",i5)') nrec
		cycle files
	endif
	call nfs(nf90_get_att(ncid,nf90_global,'l2_proc_time',l2_proc_time))
	call nfs(nf90_get_att(ncid,nf90_global,'l2_software_ver',l2_version))

	i0 = datanr+1
	i1 = datanr+nrec
	if (i1 > mrec) then
		write (*,'("Error: Too many accumulated measurements:",i5)') i1
		cycle files
	endif

! Allocate arrays

	allocate (a(nrec),b(nrec),c(nrec),d(20,nrec),valid(20,nrec),sum_c_applied(nrec))

! Time and orbit: Low rate

	call get_var_1d ('time_day_1hz',a)
	call get_var_1d ('time_milsec_1hz',b)
	call get_var_1d ('time_micsec_1hz',c)
	f0101(i0:i1) = a * 86400d0 + b * 1d-3 + c * 1d-6 + sec1990
	call get_var_1d ('latitude_1hz', f0201(i0:i1))
	call get_var_1d ('longitude_1hz', f0301(i0:i1))
	call get_var_1d ('altitude_1hz', f0404(i0:i1))
	call get_var_1d ('altitude_rate_1hz', f0501(i0:i1))
	call get_var_1d ('wf_attitude_1hz', a)
	f3002(i0:i1) = a * a

! Range data: Low rate

	call get_var_1d ('ocean_range_1hz', f0601(i0:i1))
	call get_var_1d ('ocean_stdev_1hz', f2002(i0:i1))
	call get_var_1d ('ocean_valid_num_1hz', f2101(i0:i1))
	call invalidate (f2101(i0:i1) == 0, f0601(i0:i1))
	call invalidate (f2101(i0:i1) <= 1, f2002(i0:i1))
	call get_var_1d ('f_ocean_valid_bitmap_1hz', a)
	do i = 1,nrec
		flag = nint(a(i))
		do k = 1,20
			valid (k,i) = btest(flag,k-1)
		enddo
	enddo
	! For some reason ALT_2M has 1-Hz peakiness (not documented)
	! where ALT_2S and ALT_2_ have ocean_wind_1hz
	if (alt_2m) then
		call get_var_1d ('wf_pk_1hz', f2401(i0:i1))
	else
		call get_var_1d ('ocean_wind_1hz', f1901(i0:i1))
	endif

! Range data: High rate

	if (.not.alt_2m) then
		call get_var_2d ('ocean_mean_quadratic_error', d)
		call mean_1hz (d,f3409(i0:i1),a)
		call get_var_2d ('wf_pk', d)
		call mean_1hz (d,f2401(i0:i1),a)
		call get_var_2d ('inst_range_c', d)
		call mean_1hz (d,f2702(i0:i1),a)
		call get_var_2d ('dop_c+delta_dop_c', d)
		call mean_1hz (d,f2704(i0:i1),a)
		call get_var_2d ('sum_c_applied', d)
		sum_c_applied = d(1,:)
		do i = 1,nrec
			do k = 2,19
				if (d(k,i) /= d(1,i)) write (*,*) 'Error: sum_c_applied changed:',i,k,d(1,i),d(k,i)
			enddo
		enddo
	endif

! Sigma zero: Low rate

	call get_var_1d ('ocean_sig0_1hz', f1801(i0:i1))
	call get_var_1d ('ocean_sig0_stdev_1hz', f2902(i0:i1))
	call get_var_1d ('ocean_sig0_valid_num_1hz', a)
	call invalidate (a(1:nrec) == 0, f1801(i0:i1))
	call invalidate (a(1:nrec) == 0, f1901(i0:i1))
	call invalidate (a(1:nrec) == 0, f2401(i0:i1))
	call invalidate (a(1:nrec) <= 1, f2902(i0:i1))

! SWH: Low rate

	call get_var_1d ('swh_signed_1hz', f1701(i0:i1))
	call get_var_1d ('swh_stdev_1hz', f2802(i0:i1))
	call get_var_1d ('swh_valid_num_1hz', a)
	call invalidate (a(1:nrec) == 0, f1701(i0:i1))
	call invalidate (a(1:nrec) <= 1, f2802(i0:i1))
	call get_var_1d ('em_bias_1hz', f1501(i0:i1))

! MWR Flags: Low rate

	call get_var_1d ('f_sea_ice_flag_1hz', a)
	f2601(i0:i1) = 0
	call flag_set (a(1:nrec) == 1, f2601(i0:i1), 8)

! MWR: Low rate

	call get_var_1d ('tb_23_8_1hz', f2302(i0:i1))
	call get_var_1d ('tb_36_5_1hz', f2303(i0:i1))
	call get_var_1d ('f_mwr_srf_typ_1hz', a)
	call flag_set (a(1:nrec) == 1, f2601(i0:i1), 5)
	call get_var_1d ('f_mwr_interp_qual_1hz', a)
	call invalidate (a(1:nrec) == 3, f2302(i0:i1))
	call invalidate (a(1:nrec) == 3, f2303(i0:i1))
	call get_var_1d ('f_MWR_valid_1hz', a)
	call invalidate (a(1:nrec) == 1, f2302(i0:i1))
	call invalidate (a(1:nrec) == 1, f2303(i0:i1))

! Atmospheric and geophysical: Low rate

	call get_var_1d ('dry_c_1hz', f0701(i0:i1))
	call get_var_1d ('ib_c_1hz', f1002(i0:i1))
	call get_var_1d ('mog2d_c_1hz', f1004(i0:i1))
	call get_var_1d ('wet_c_mod_1hz', f0802(i0:i1))
	call get_var_1d ('wet_c_mwr_1hz', f0801(i0:i1))
	call get_var_1d ('water_vapor_content_1hz', f3302(i0:i1))
	call get_var_1d ('liquid_water_content_1hz', f3301(i0:i1))
	call get_var_1d ('u_wind_1hz', f1903(i0:i1))
	call get_var_1d ('v_wind_1hz', f1904(i0:i1))
	call get_var_1d ('iono_c_mod_1hz', f0908(i0:i1))
	call get_var_1d ('iono_c_gps_1hz', f0906(i0:i1))
	call get_var_1d ('h_mss_cls01_1hz', f1605(i0:i1))
	call get_var_1d ('h_mss_ucl04_1hz', f1613(i0:i1))
	call get_var_1d ('h_geo_1hz', f1610(i0:i1))
	call get_var_1d ('h_ot_1hz', f1217(i0:i1))
	call get_var_1d ('h_ot2_1hz', f1213(i0:i1))
	call get_var_1d ('h_olt_1hz', f1317(i0:i1))
	call get_var_1d ('h_olt2_1hz', f1313(i0:i1))
	call get_var_1d ('h_lpt_1hz', f3901(i0:i1))
	call get_var_1d ('h_lptne_1hz', f3902(i0:i1))
	call get_var_1d ('h_set_1hz', f1101(i0:i1))
	call get_var_1d ('h_pol_1hz', f1401(i0:i1))

	call get_var_1d ('f_srf_typ_1hz', a)
	call flag_set (a(1:nrec) >= 3, f2601(i0:i1), 2)
	call flag_set (a(1:nrec) >= 2, f2601(i0:i1), 4)
	call flag_set (a(1:nrec) >= 1, f2601(i0:i1), 5)

	call get_var_1d ('h_odle_1hz', f2206(i0:i1))
	if (.not.alt_2m) call get_var_1d ('sig0_attn_c_1hz', f3203(i0:i1))

	call get_var_1d ('f_corr_error_1hz', a)
	do i = 1,nrec
		k = datanr + i
		flag = nint(a(i))
		if (btest(flag, 1)) f0701(k) = nan
		if (btest(flag, 2)) f1002(k) = nan
		if (btest(flag, 3)) f1004(k) = nan
		if (btest(flag, 5)) f0802(k) = nan
		if (btest(flag, 6)) then
			f0801(k) = nan
			f3301(k) = nan
			f3302(k) = nan
		endif
		if (btest(flag, 7)) f1903(k) = nan
		if (btest(flag, 8)) f1904(k) = nan
		if (btest(flag, 9)) f0908(k) = nan
		if (btest(flag,10)) f0906(k) = nan
		if (btest(flag,11)) f1605(k) = nan
		if (btest(flag,12)) f1610(k) = nan
		if (btest(flag,13)) f1217(k) = nan
		if (btest(flag,14)) f1213(k) = nan
		if (btest(flag,15)) f1317(k) = nan
		if (btest(flag,16)) f1313(k) = nan
		if (btest(flag,17)) f3901(k) = nan
		if (btest(flag,18)) f3902(k) = nan
		if (btest(flag,19)) f1101(k) = nan
		if (btest(flag,20)) f1401(k) = nan
		if (btest(flag,22)) f2206(k) = nan
		if (btest(flag,23)) f1501(k) = nan
		if (btest(flag,27)) f1613(k) = nan
		if (btest(flag,28)) f3203(k) = nan
	enddo

! Remove applied corrections from range

	call get_var_1d ('f_corr_applied_1hz', a)
	do i = 1,nrec
		k = datanr + i
		flag = nint(a(i))
		sum_d_applied = 0d0
		if (btest(flag, 1)) sum_d_applied = sum_d_applied + f0701(k)
		if (btest(flag, 2)) sum_d_applied = sum_d_applied + f1002(k)
		if (btest(flag, 3)) sum_d_applied = sum_d_applied + f1004(k)
		if (btest(flag, 5)) sum_d_applied = sum_d_applied + f0802(k)
		if (btest(flag, 6)) sum_d_applied = sum_d_applied + f0801(k)
		if (btest(flag, 7)) sum_d_applied = sum_d_applied + f1903(k)
		if (btest(flag, 8)) sum_d_applied = sum_d_applied + f1904(k)
		if (btest(flag, 9)) sum_d_applied = sum_d_applied + f0908(k)
		if (btest(flag,10)) sum_d_applied = sum_d_applied + f0906(k)
		if (btest(flag,11)) sum_d_applied = sum_d_applied + f1605(k)
		if (btest(flag,12)) sum_d_applied = sum_d_applied + f1610(k)
		if (btest(flag,13)) sum_d_applied = sum_d_applied + f1217(k)
		if (btest(flag,14)) sum_d_applied = sum_d_applied + f1213(k)
		if (btest(flag,15)) sum_d_applied = sum_d_applied + f1317(k)
		if (btest(flag,16)) sum_d_applied = sum_d_applied + f1313(k)
		if (btest(flag,17)) sum_d_applied = sum_d_applied + f3901(k)
		if (btest(flag,18)) sum_d_applied = sum_d_applied + f3902(k)
		if (btest(flag,19)) sum_d_applied = sum_d_applied + f1101(k)
		if (btest(flag,20)) sum_d_applied = sum_d_applied + f1401(k)
		if (btest(flag,22)) sum_d_applied = sum_d_applied + f2206(k)
		if (btest(flag,23)) sum_d_applied = sum_d_applied + f1501(k)
		if (btest(flag,27)) sum_d_applied = sum_d_applied + f1613(k)
		f0601(k) = f0601(k) - sum_d_applied
		if (.not.alt_2m .and. .not.(sum_d_applied == sum_c_applied(i))) &
			write (*,*) "Error: sum_c_applied wrong: ",i,sum_d_applied,sum_c_applied(i)
	enddo

! Close this input file

	deallocate (a,b,c,d,valid)

	call nfs(nf90_close(ncid))

	datanr = i1

! Dump whatever data we can

	do
		new = erspass (ers, f0101(1), orbitnr(1), phasenm(1), cyclenr(1), passnr(1), tnode(1), lnode(1))
		do i = 2,datanr
			if (erspass (ers, f0101(i), orbitnr(2), phasenm(2), cyclenr(2), passnr(2), tnode(2), lnode(2))) exit
		enddo
		if (i > datanr) exit
		nrec = i - 1
		call write_data (nrec)
		datanr = datanr - nrec
		call move_data (nrec, datanr)
	enddo

enddo files ! Each file

! Dump whatever remains

new = erspass (ers, f0101(1), orbitnr(1), phasenm(1), cyclenr(1), passnr(1), tnode(1), lnode(1))
call write_data (datanr)

contains

subroutine synopsis ()
if (rads_version ('$Revision$', 'Write REAPER data to RADS')) return
write (*,1310)
1310 format (/ &
'syntax: reaperraw [options] < list_of_REAPER_file_names'// &
'This program converts REAPER ALT_2M files to RADS data'/ &
'files with the name $RADSDATAROOT/data/eE/F.r/pPPPP/eEpPPPPcCCC.nc.'/ &
'The directory is created automatically and old files are overwritten.')
stop
end subroutine synopsis

subroutine move_data (bottom, n)
integer(fourbyteint), intent(in) :: bottom, n
integer(fourbyteint) :: i0, i1
i0 = bottom + 1
i1 = bottom + n
f0101(1:n) = f0101(i0:i1)
f0201(1:n) = f0201(i0:i1)
f0301(1:n) = f0301(i0:i1)
f0404(1:n) = f0404(i0:i1)
f0501(1:n) = f0501(i0:i1)
f0601(1:n) = f0601(i0:i1)
f0701(1:n) = f0701(i0:i1)
f0801(1:n) = f0801(i0:i1)
f0802(1:n) = f0802(i0:i1)
f0906(1:n) = f0906(i0:i1)
f0908(1:n) = f0908(i0:i1)
f1002(1:n) = f1002(i0:i1)
f1004(1:n) = f1004(i0:i1)
f1101(1:n) = f1101(i0:i1)
f1213(1:n) = f1213(i0:i1)
f1217(1:n) = f1217(i0:i1)
f1313(1:n) = f1313(i0:i1)
f1317(1:n) = f1317(i0:i1)
f1401(1:n) = f1401(i0:i1)
f1501(1:n) = f1501(i0:i1)
f1605(1:n) = f1605(i0:i1)
f1610(1:n) = f1610(i0:i1)
f1613(1:n) = f1613(i0:i1)
f1701(1:n) = f1701(i0:i1)
f1801(1:n) = f1801(i0:i1)
f1901(1:n) = f1901(i0:i1)
f1903(1:n) = f1903(i0:i1)
f1904(1:n) = f1904(i0:i1)
f2002(1:n) = f2002(i0:i1)
f2101(1:n) = f2101(i0:i1)
f2206(1:n) = f2206(i0:i1)
f2302(1:n) = f2302(i0:i1)
f2303(1:n) = f2303(i0:i1)
f2401(1:n) = f2401(i0:i1)
f2601(1:n) = f2601(i0:i1)
f2702(1:n) = f2702(i0:i1)
f2704(1:n) = f2704(i0:i1)
f2802(1:n) = f2802(i0:i1)
f2902(1:n) = f2902(i0:i1)
f3002(1:n) = f3002(i0:i1)
f3203(1:n) = f3203(i0:i1)
f3301(1:n) = f3301(i0:i1)
f3302(1:n) = f3302(i0:i1)
f3409(1:n) = f3409(i0:i1)
f3901(1:n) = f3901(i0:i1)
f3902(1:n) = f3902(i0:i1)
end subroutine move_data

subroutine write_data (datanr)
integer(fourbyteint), intent(in) :: datanr
integer(fourbyteint) :: i
real(eightbytereal) :: dhellips
type(rads_sat) :: S
type(rads_pass) :: P

if (datanr == 0) return	! Skip empty data sets
if (cyclenr(1) < c0 .or. cyclenr(1) > c1) return	! Skip chunks that are not of the selected cycle
if (tnode(1) < t0 .or. tnode(1) > t1) return	! Skip equator times that are not of selected range

! Determine mission and phase
! Use sat from command line if given, otherwise autogenerate

if (sat /= '') then
	newsat = sat
else
	newsat = mission(1:2)//'/'//strtolower(phasenm(1))//'.r'
endif
if (newsat /= oldsat) then
	oldsat = newsat
	call rads_init (S, newsat,verbose)
endif

! Store relevant info

P%cycle=cyclenr(1)
P%pass=passnr(1)
P%start_time = f0101(1)
P%end_time = f0101(datanr)
P%equator_time = tnode(1)
P%equator_lon = lnode(1)
P%original = 'REAPER '//trim(l2_version)//' data of '//trim(l2_proc_time)

! Open output file

call rads_create_pass (S, P)

! Correct altitude WGS84 to TOPEX ellipsoid

do i = 1,datanr
	f0404(i) = f0404(i) + dhellips(1,f0201(i)*1d-6)*1d3
enddo

! Fill all the data fields

call put_var (S, P,  101, f0101 * 1d0 )
call put_var (S, P,  201, f0201 * 1d-6)
call put_var (S, P,  301, f0301 * 1d-6)
call put_var (S, P,  404, f0404 * 1d-3)
call put_var (S, P,  501, f0501 * 1d-3)
call put_var (S, P,  601, f0601 * 1d-3)
call put_var (S, P,  701, f0701 * 1d-3)
call put_var (S, P,  801, f0801 * 1d-3)
call put_var (S, P,  802, f0802 * 1d-3)
call put_var (S, P,  906, f0906 * 1d-3)
call put_var (S, P,  908, f0908 * 1d-3)
call put_var (S, P, 1002, f1002 * 1d-3)
call put_var (S, P, 1004, f1004 * 1d-3)
call put_var (S, P, 1101, f1101 * 1d-3)
call put_var (S, P, 1213, f1213 * 1d-3)
call put_var (S, P, 1217, f1217 * 1d-3)
call put_var (S, P, 1313, f1313 * 1d-3)
call put_var (S, P, 1317, f1317 * 1d-3)
call put_var (S, P, 1401, f1401 * 1d-3)
call put_var (S, P, 1501, f1501 * 1d-3)
call put_var (S, P, 1605, f1605 * 1d-3)
call put_var (S, P, 1610, f1610 * 1d-3)
call put_var (S, P, 1613, f1613 * 1d-3)
call put_var (S, P, 1701, f1701 * 1d-3)
call put_var (S, P, 1801, f1801 * 1d-2)
call put_var (S, P, 1903, f1903 * 1d-3)
call put_var (S, P, 1904, f1904 * 1d-3)
call put_var (S, P, 2002, f2002 * 1d-3)
call put_var (S, P, 2101, f2101 * 1d0 )
call put_var (S, P, 2206, f2206 * 1d-3)
call put_var (S, P, 2302, f2302 * 1d-2)
call put_var (S, P, 2303, f2303 * 1d-2)
call put_var (S, P, 2601, f2601 * 1d0 )
call put_var (S, P, 2702, f2702 * 1d-3)
call put_var (S, P, 2802, f2802 * 1d-3)
call put_var (S, P, 2902, f2902 * 1d-2)
call put_var (S, P, 3002, f3002 * 1d-4)
call put_var (S, P, 3301, f3301 * 1d-2)
call put_var (S, P, 3302, f3302 * 1d-2)
call put_var (S, P, 3901, f3901 * 1d-3)
call put_var (S, P, 3902, f3902 * 1d-3)
if (alt_2m) then
	call put_var (S, P, 2401, f2401 * 1d-3)
else
	call put_var (S, P, 1901, f1901 * 1d-3)
	call put_var (S, P, 3203, f3203 * 1d-2)
	call put_var (S, P, 2702, f2702 * 1d-3)
	call put_var (S, P, 2704, f2704 * 1d-3)
	call put_var (S, P, 3409, f3409 * 1d-4)
endif

! Write to the data file

552 format ('...',i5,' records written to ',a)
write (*,552) datanr,trim(P%filename)
call rads_close_pass (S, P)

end subroutine write_data

subroutine put_var (S, P, field, data)
type(rads_sat), intent(inout) :: S
type(rads_pass), intent(inout) :: P
integer(fourbyteint), intent(in) :: field
real(eightbytereal), intent(in) :: data(:)
character(len=4) :: name
integer :: i, start(1) = (/ 1 /)
! Look for field number in variable list
do i = 1,S%nvar
	if (any(S%var(i)%field == field)) then
		call rads_put_var (S, P, S%var(i), data(:P%ndata), start)
		return
	endif
enddo
write (name,'(i0)') field
call rads_error (S, rads_err_var, 'Variable with field number '//name//' not found')
end subroutine put_var

subroutine get_var_1d (varnm, array)
character(*), intent(in) :: varnm
real(eightbytereal), intent(out) :: array(:)
real(eightbytereal) :: array2(mrec)
integer(fourbyteint) :: i0,i1,l,varid,constant
i1 = 0
l = len_trim(varnm)
do
	if (i1 > l) exit
	i0 = i1
	i1 = scan(varnm(i0+1:), '+-') + i0
	if (i1 == i0) i1 = l + 1
	if (nf90_inq_varid(ncid,varnm(i0+1:i1-1),varid) /= nf90_noerr) then
		write (*,'("No such variable: ",a)') varnm(i0+1:i1-1)
		return
	endif
	if (i0 == 0) then
		call nfs(nf90_get_var(ncid,varid,array(1:nrec)))
	else
		call nfs(nf90_get_var(ncid,varid,array2(1:nrec)))
		constant = 0
		if (varnm(i0:i0) == '-') constant = -1
		if (varnm(i0:i0) == '+') constant = 1
		array(1:nrec) = array(1:nrec) + constant * array2(1:nrec)
	endif
enddo
end subroutine get_var_1d

subroutine get_var_2d (varnm, array)
character(*), intent(in) :: varnm
real(eightbytereal), intent(out) :: array(:,:)
real(eightbytereal) :: array2(20,mrec)
integer(fourbyteint) :: i0,i1,l,varid,constant
i1 = 0
l = len_trim(varnm)
do
	if (i1 > l) exit
	i0 = i1
	i1 = scan(varnm(i0+1:), '+-') + i0
	if (i1 == i0) i1 = l + 1
	if (nf90_inq_varid(ncid,varnm(i0+1:i1-1),varid) /= nf90_noerr) then
		write (*,'("No such variable: ",a)') varnm(i0+1:i1-1)
		return
	endif
	if (i0 == 0) then
		call nfs(nf90_get_var(ncid,varid,array(:,1:nrec)))
	else
		call nfs(nf90_get_var(ncid,varid,array2(:,1:nrec)))
		constant = 0
		if (varnm(i0:i0) == '-') constant = -1
		if (varnm(i0:i0) == '+') constant = 1
		array(:,1:nrec) = array(:,1:nrec) + constant * array2(:,1:nrec)
	endif
enddo
end subroutine get_var_2d

subroutine mean_1hz (y, mean, rms)
real(eightbytereal), intent(in) :: y(:,:)
real(eightbytereal), intent(out) :: mean(:), rms(:)
integer(fourbyteint) :: i, j, n
do j = 1,nrec
	mean(j) = 0d0
	rms(j) = 0d0
	n = 0
	do i = 1,20
		if (valid(i,j)) then
			n = n + 1
			mean(j) = mean(j) + y(i,j)
			rms(j) = rms(j) + y(i,j)**2
		endif
	enddo
	if (n < 1) then
		mean(j) = nan
	else
		mean(j) = mean(j) / n
	endif
	if (n < 2) then
		rms(j) = nan
	else
		rms(j) = sqrt ((rms(j) - n * mean(j)**2) / (n - 1))
	endif
enddo
end subroutine mean_1hz

subroutine flag_set (a, flags, bit)
logical, intent(in) :: a(:)
integer(twobyteint), intent(inout) :: flags(:)
integer(fourbyteint), intent(in) :: bit
integer(fourbyteint) :: i
integer(twobyteint) :: j
if (size(a) /= size(flags)) stop "Error in flag_set"
j = int(bit,twobyteint)
do i = 1,size(a)
	if (a(i)) flags(i) = ibset(flags(i),j)
enddo
end subroutine flag_set

subroutine invalidate (a, b)
logical, intent(in) :: a(:)
real(eightbytereal), intent(inout) :: b(:)
if (size(a) /= size(b)) stop "Error in invalidate"
where (a) b = nan
end subroutine invalidate

end program rads_gen_reaper
