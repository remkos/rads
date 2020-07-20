!-----------------------------------------------------------------------
! Copyright (c) 2011-2020  Remko Scharroo
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

!*rads_fix_c2 -- Patch RADS altimeter files of CryoSat for various anomalies
!
! This program makes numerous patches to the CryoSat RADS data processed
! by rads_gen_c2_l1r. These patches include:
!
! meteo:
! - Set dry tropo and wet tropo to NaN when both are zero
! - Set iono to NaN when it is zero as well
!
! sig0 (FDM L2 and LRM L2 only):
! - Reduce backscatter by 1.4 dB during period 2010-08-13 00:00 to 2010-09-03 00:00
! - Reduce backscatter by 1.4 dB during period 2010-10-04 00:00 to 2010-10-21 00:00
! - Reduce backscatter by 40. dB during period 2011-01-28 00:00 to 2011-02-27 23:38
! - Increase FDM L2 backscatter by 4.8 dB
! - ... Then flip the backscatter around (sig0 := 21 - sig0)
!
! sig0 (LRM and PLRM L1):
! - Increase backscatter by 1.5 dB during period 2010-08-13 00:00 to 2010-09-03 00:00
! - Increase backscatter by 1.5 dB during period 2010-10-04 00:00 to 2010-10-21 00:00
! - Reduce backscatter by 3.02 dB (LRM) or 3.04 dB (PLRM)
!
! drift (sig0):
! - Increase backscatter by 0.22 dB/year (LRM) or 0.27 dB/year (PLRM), based on
!   Thales documentation and own analysis 18-19 Dec 2013
! - Reference time (time with no correction) is 2011-05-01, the middle of the period
!   over which original sigma0 bias and SSB was determined
!
! usage: rads_fix_c2 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_c2

use rads
use rads_devel
use rads_misc

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

! Sigma0 drift and bias with reference time 1-MAY-2011. See notes of 18-DEC-2013
real(eightbytereal) :: sig0_drift_lrm = 0.22d0 / 365.25d0 / 86400d0, sig0_bias_lrm = -3.02d0, &
	sig0_drift_sar = 0.27d0 / 365.25d0 / 86400d0, sig0_bias_sar = -3.04d0, time_drift = 830822400d0
integer(fourbyteint) :: i, cyc, pass
integer(twobyteint) :: flag
logical :: ldrift = .false., lmeteo = .false., lsig0 = .false.
character(len=1) :: baseline

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' drift meteo sig0 all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('drift')
		ldrift = .true.
	case ('meteo')
		lmeteo = .true.
	case ('sig0')
		lsig0 = .true.
	case ('all')
		ldrift = .true.
		lmeteo = .true.
		lsig0 = .true.
	end select
enddo

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
if (rads_version ('Patch CryoSat-2 data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --drift                   Correct sigma0 for apparent drift (Baseline B only)' / &
'  --meteo                   Set dry, wet, IB (and iono) to NaN when zero' / &
'  --sig0                    Correct sigma0 for biases' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n), dry(n), wet(n), iono(n), sig0(n), flagword(n), &
	ib(n), sig0_20hz(20,n), dsig0
integer(fourbyteint) :: i
logical :: lrm_l2, fdm_l2, sar, cmeteo, ciono, cdrift

call log_pass (P)

! Do we have FDM or LRM? Do we have an older version?
lrm_l2 = index(P%original,'IPF2LRM') > 0
fdm_l2 = index(P%original,'IPF2FDM') > 0

! Determine baseline version
i = index(P%original,'CS_')
baseline = P%original(i+51:i+51)

! These may be set to true later on
cmeteo = .false.
ciono = .false.
cdrift = ldrift .and. baseline < 'C'

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'sig0_ku', sig0, .true.)
call rads_get_var (S, P, 'flags', flagword, .true.)
if (P%n_hz > 0) call rads_get_var (S, P, 'sig0_20hz_ku', sig0_20hz, .true.)
if (lmeteo) then
	call rads_get_var (S, P, 'dry_tropo_ecmwf', dry, .true.)
	call rads_get_var (S, P, 'wet_tropo_ecmwf', wet, .true.)
	call rads_get_var (S, P, 'inv_bar_static', ib, .true.)
	call rads_get_var (S, P, 'iono_gim', iono, .true.)
endif

! Process data records

do i = 1,n
	flag = nint(flagword(i),twobyteint)
	sar = btest(flag,0)

! Set dry and wet to NaN when both are zero (not available)

	if (lmeteo .and. dry(i) == 0d0 .and. wet(i) == 0d0) then
		dry(i) = nan
		wet(i) = nan
		ib(i) = nan
		cmeteo = .true.
		if (iono(i) == 0d0) then
			iono(i) = nan
			ciono = .true.
		endif
	endif

! For data older than Baseline C: Correct for sigma0 drift

	if (.not.cdrift) then
		dsig0 = 0
	else if (sar) then
		dsig0 = sig0_drift_sar * (time(i)-time_drift)
	else
		dsig0 = sig0_drift_lrm * (time(i)-time_drift)
	endif

! For all data: Correct sigma0 for biases

	if (.not.lsig0) then
		! Skip
	else if (fdm_l2) then
		dsig0 = dsig0 - 1.4d0
	else if (lrm_l2) then
		dsig0 = dsig0 - 7.0d0
	else	! All L1 data (LRM and PLRM)
		if (baseline < 'C' .and. ((time(i) > 808272000d0 .and. time(i) < 810086400d0) .or. &	! 2010-08-13 00:00 - 2010-09-03 00:00
			 (time(i) > 812764800d0 .and. time(i) < 814233600d0))) &	! 2010-10-04 00:00 - 2010-10-21 00:00
			dsig0 = dsig0 + 1.5d0 ! Baseline B data for these periods were biased low by 1.5 dB, so we add 1.5 dB
		if (sar) then
			dsig0 = dsig0 + sig0_bias_sar
		else
			dsig0 = dsig0 + sig0_bias_lrm
		endif
	endif

	sig0(i) = sig0(i) + dsig0
	if (P%n_hz > 0) sig0_20hz(:,i) = sig0_20hz(:,i) + dsig0
enddo

! Baseline B LRM L2 data: Prior to Cycle 5 pass 333 all iono is bogus

if (lmeteo .and. baseline < 'C' .and. lrm_l2 .and. P%equator_time < 808700000d0) then
	iono = nan
	ciono = .true.
endif

! If nothing changed, stop here

if (.not.(cmeteo .or. ciono .or. lsig0 .or. cdrift)) then
	call log_records (0)
	return
endif

! Update history and define new variables (if required)

P%start_time = time(1)
P%end_time = time(n)
call rads_put_passinfo (S, P)
call rads_put_history (S, P)

! Write out all the data

if (cmeteo) then
	call rads_put_var (S, P, 'dry_tropo', dry)
	call rads_put_var (S, P, 'wet_tropo', wet)
	call rads_put_var (S, P, 'inv_bar_static', ib)
endif
if (ciono) call rads_put_var (S, P, 'iono_gim', iono)
if (lsig0 .or. cdrift) then
	call rads_put_var (S, P, 'sig0_ku', sig0)
	if (P%n_hz > 0) call rads_put_var (S, P, 'sig0_20hz_ku', sig0_20hz)
endif

call log_records (n)
end subroutine process_pass

end program rads_fix_c2
