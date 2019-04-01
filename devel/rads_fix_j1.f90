!-----------------------------------------------------------------------
! Copyright (c) 2011-2019  Remko Scharroo
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

!*rads_fix_j1 -- Patch RADS altimeter files of Jason-1 for various anomalies
!
! This program makes numerous patches to the Jason-1 RADS data processed
! by rads_gen_jason. These patches include:
!
! dry:
! - Take a jump out of the time series
!
! jmre:
! - Add enhanced JMR data
!
! range:
! - Adjust C-band range for bias and adjust iono accordingly
!
! sig0:
! - Adjust backscatter for bias
!
! wind:
! - Recompute wind speed from adjusted sigma0 based on Collard model
!
! usage: rads_fix_j1 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_j1

use rads
use rads_devel
use rads_misc
use meteo_subs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

real(eightbytereal), parameter :: dsig0_ku = -2.40d0, dsig0_c = -0.73d0	! Ku- and C-band Sigma0 bias of Jason-1
real(eightbytereal), parameter :: range_c_bias = -0.0023d0
integer(fourbyteint) :: i, cyc, pass
logical :: ldry = .false., ljmre = .false., lrange = .false., lsig0 = .false., lwind = .false.

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' dry jmre range sig0 all wind')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('dry')
		ldry = .true.
	case ('jmre')
		ljmre = .true.
	case ('range')
		lrange = .true.
	case ('sig0')
		lsig0 = .true.
	case ('all')
		ldry = .true.
		ljmre = .true.
		lrange = .true.
		lsig0 = .true.
	case ('wind')
		lwind = .true.
	end select
enddo

! If nothing selected, stop here

if (.not.(ldry .or. ljmre .or. lrange .or. lsig0 .or. lwind)) stop

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Patch Jason-1 data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --dry                     Take a jump out of the time series' / &
'  --jmre                    Add enhanced JMR data' / &
'  --range                   Adjust C-band range for bias and adjust iono accordingly' / &
'  --sig0                    Adjust backscatter coefficient for bias' / &
'  --all                     All of the above' / &
'  --wind                    Recompute wind speed from adjusted sigma0 based on Collard model')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: sig0_ku(n), sig0_c(n), swh_ku(n), range_c(n), dry_tropo_ecmwf(n), &
	wet_tropo_rad(n), flags(n), u(n), iono_alt(n)
integer(twobyteint) :: jmre(n,4), flag
integer(fourbyteint) :: i
logical :: cjmre

call log_pass (P)

! Adjust backscatter for bias

if (lsig0) then
	call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	call rads_get_var (S, P, 'sig0_c', sig0_c, .true.)
	sig0_ku = sig0_ku + dsig0_ku
	sig0_c  = sig0_c  + dsig0_c
endif

! Compute wind speed from Collard model (using unadjusted sig0_ku)

if (lwind) then
	call rads_get_var (S, P, 'swh_ku', swh_ku, .true.)
	if (.not.lsig0) call rads_get_var (S, P, 'sig0_ku', sig0_ku, .true.)
	do i = 1,n
		u(i) = wind_j1 (1, sig0_ku(i) - dsig0_ku, swh_ku(i))
	enddo
endif

! Jason manual suggests to add 131 mm to range.
! Also add bias to C-band to fix offset in iono delay.

if (lrange) then
	call rads_get_var (S, P, 'range_c', range_c, .true.)
	call rads_get_var (S, P, 'iono_alt', iono_alt, .true.)
	range_c = range_c + range_c_bias
	iono_alt = iono_alt - 0.17984d0 * range_c_bias
endif

! Remove a 0.3 mm bias from the ECMWF dry tropo correction after 31 Jan 2006

if (ldry .and. P%equator_time > 665280000d0) then
	call rads_get_var (S, P, 'dry_tropo_ecmwf', dry_tropo_ecmwf)
	dry_tropo_ecmwf = dry_tropo_ecmwf + 3d-4
endif

! Add the enhanced JMR products (if available)

cjmre = ljmre .and. get_jmr(n,jmre)
if (cjmre) then
	call rads_get_var (S, P, 'flags', flags, .true.)
	do i = 1,n
		flag = nint(flags(i),twobyteint)

		! Bit 8: Radiometer rain/ice flag
		if (jmre(i,2) == 1 .or. jmre(i,3) == 1) then
			flag = ibset(flag,8)
		else
			flag = ibclr(flag,8)
		endif

		! Bit 6: Radiometer land flag
		if (jmre(i,4) == 2) then
			flag = ibset(flag,6)
		else
			flag = ibclr(flag,6)
		endif

		flags(i) = flag
		wet_tropo_rad(i) = jmre(i,1) * 1d-4
	enddo
endif

! Update history

call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Write out all the data

if (ldry .and. P%equator_time > 665280000d0) &
	call rads_put_var (S, P, 'dry_tropo_ecmwf', dry_tropo_ecmwf)
if (lsig0) then
	call rads_put_var (S, P, 'sig0_ku', sig0_ku)
	call rads_put_var (S, P, 'sig0_c' , sig0_c)
endif
if (lwind) call rads_put_var (S, P, 'wind_speed_alt', u)
if (lrange) then
	call rads_put_var (S, P, 'range_c', range_c)
	call rads_put_var (S, P, 'iono_alt' , iono_alt)
endif
if (cjmre) then
	call rads_put_var (S, P, 'flags', flags)
	call rads_put_var (S, P, 'wet_tropo_rad', wet_tropo_rad)
endif

call log_records (n)
end subroutine process_pass

function get_jmr (n, jmre)
use netcdf
use rads_netcdf
integer(fourbyteint), intent(in) :: n
integer(twobyteint), intent(out) :: jmre(:,:)
logical :: get_jmr
character(128) :: filenm
integer(fourbyteint) :: ncid, i

get_jmr = .false.

! Get filename like JA1_GPN_JMR_EXP_2PcP300_093_20100225_150810_20100225_160422.nc
! And convert to $RADSROOT/ext/j1/jmr_enh/c300/JA1_GPN_JMR_EXP_2PcP300_093_20100225_150810_20100225_160422.nc

i = index(P%original, ' ') - 1
call parseenv ('${RADSROOT}/ext/j1/jmr_enh/c'//P%original(13:15)// &
    '/JA1_GPN_JMR_EXP_2PcP'//P%original(13:i), filenm)

if (nf90_open(filenm,nf90_nowrite,ncid) /= nf90_noerr) return
i = index(filenm, '/j1/jmr_enh')
call log_string (filenm(i+1:))

! Check file size

call nfs(nf90_inquire_dimension(ncid,1,len=i))
if (i /= n) then
	write (rads_log_unit, 555, advance='no') i
    return
endif
555 format ('wrong file size:',i5,' ... ')

! Now get variables

do i = 1,4
    call nfs(nf90_get_var(ncid,i+3,jmre(:,i)))
enddo

! Succesful completion

i = nf90_close(ncid)
get_jmr = .true.
end function get_jmr

end program rads_fix_j1
