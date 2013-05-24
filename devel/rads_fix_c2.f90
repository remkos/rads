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

!*rads_fix_c2 -- Patch RADS altimeter files of CryoSat for various anomalies
!
! This program makes numerous patches to the CryoSat RADS data processed
! by c2lrmraw. These patches include:
!
! meteo:
! - Set dry tropo and wet tropo to NaN when both are zero
! - Set iono to NaN when it is zero as well
!
! range:
! - Add -3.33 m range bias for LRM L2 data prior to Feb 2011
! - Supposedly we need to add 7 mm (1/64 of a range gate) to FDM/LRM L1 and
!   FDM/LRM L2 range because of error in CAL1
!
! sig0 (FDM L2 and LRM L2 only):
! - Reduce backscatter by 1.4 dB during period 2010-08-13 00:00 to 2010-09-03 00:00
! - Reduce backscatter by 1.4 dB during period 2010-10-04 00:00 to 2010-10-21 00:00
! - Reduce backscatter by 40. dB during period 2011-01-28 00:00 to 2011-02-27 23:38
! - Increase FDM L2 backscatter by 4.8 dB
! - ... Then flip the backscatter around (sig0 := 21 - sig0)
!
! drift (sig0):
! - Increase backscatter by 0.32 dB per year, based on information provided by ESA (Marco)
!   on 15 April 2013
! - Reference time (time with no correction) is 2011-05-01, the middle of the period
!   over which original sigma0 bias and SSB was determined
!
! ssb:
! - Interpolate hybrid model
!
! swh (FDM L2 and LRM L2 only):
! - Set SWH to NaN when SWH quality flag is raised (otherwise reports 0)
! - Correct for error in L2 SWH algorithm
! - Adjust pseudo-LRM SWH to match LRM data (see notes of 16 Sep 2012)
!
! tbias:
! - Apply timing bias according to Marco's table (different for SAR and LRM data)
! - Account for the change in timing bias in SIR1FDM/2.4 (See Ruby's e-mail of 22 Apr 2013)
! - Adjust altitude from the product accordingly (using altitude rate)
! - Note that this does NOT change the equator time or longitude!
!
! wind:
! - Add wind speed according to ECMWF model
!
! usage: rads_fix_c2 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_c2

use rads
use rads_misc
use rads_grid
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

character(len=rads_cmdl) :: path
real(eightbytereal), parameter :: fai = 7.3d-3, &
	sig0_drift = 0.32d0 / 365.25d0 / 86400d0, time_drift = 830822400d0, & ! Sigma0 drift is 0.32 dB per year since 1-MAY-2011
	swh_adjustment = 1.8737**2 * (0.513**2 - 0.383**2) ! Adjustment to be added to SWH squared
integer(fourbyteint) :: i,cyc,pass
integer(twobyteint) :: flag
logical :: ldrift=.false.,lmeteo=.false.,lrange=.false.,lssb=.false.,lswh=.false.,ltbias=.false., &
	lwind=.false.,lsig0=.false.,cswh,cmeteo,ciono,fdm_l1_v24
type(grid) :: issb_hyb

! Scan command line for options

call synopsis
call rads_set_options (' drift meteo range sig0 ssb swh tbias wind all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('drift')
		ldrift = .true.
	case ('meteo')
		lmeteo = .true.
	case ('range')
		lrange = .true.
	case ('sig0')
		lsig0 = .true.
	case ('ssb')
		lssb = .true.
	case ('swh')
		lswh = .true.
	case ('tbias')
		ltbias = .true.
	case ('wind')
		lwind = .true.
	case ('all')
		ldrift = .true.
		lmeteo = .true.
		lrange = .true.
		lsig0 = .true.
		lssb = .true.
		lswh = .true.
		ltbias = .true.
		lwind = .true.
	end select
enddo

! Load SSB model if requested

if (.not.lssb) then
	! Do nothing
else if (index(S%phase%dataroot, '/a.l2') > 0) then
	call parseenv ('${ALTIM}/data/models/c2_l2i_hyb.nc?ssb_hyb', path)
	if (grid_load(path,issb_hyb) /= 0) call rads_exit ('Error loading '//trim(path))
else
	call parseenv ('${ALTIM}/data/models/c2_l1c_hyb.nc?ssb_hyb', path)
	if (grid_load(path,issb_hyb) /= 0) call rads_exit ('Error loading '//trim(path))
endif

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
if (rads_version ('$Revision$', 'Patch CryoSat-2 data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --drift                   Correct sigma0 for apparent drift' / &
'  --meteo                   Set dry, wet, IB (and iono) to NaN when zero' / &
'  --range                   Correct range for biases' / &
'  --sig0                    Correct sigma0 for biases and reversal (L2 only)' / &
'  --ssb                     Add hybrid SSB model' / &
'  --swh                     Correct SWH' / &
'  --tbias                   Correct time and orbital altitude for timing bias' / &
'  --wind                    Add wind speed' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n),alt(n),alt_rate(n),dry(n),wet(n), &
	iono(n),sig0(n),swh(n),ssb(n),wind(n),flagword(n),range(n),ib(n),ecmwf_ws,t,tbias
integer(fourbyteint) :: i
logical :: lrm_l2,fdm_l2,old_version_a,version_a,sar

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename)

! Do we have FDM or LRM? Do we have an older version?
lrm_l2 = index(P%original,'IPF2LRM') > 0
fdm_l2 = index(P%original,'IPF2FDM') > 0
version_a = index(P%original,'_A00') > 0
fdm_l1_v24 = index(P%original,'SIR1FDM/2.4') > 0

old_version_a = (index(P%original,'-2010') > 0 .or. &
	index(P%original,'-JAN-2011') > 0 .or. index(P%original,'-FEB-2011') > 0)

! These may be set to true later on
cmeteo = .false.
ciono = .false.

! SWH corrections that apply to LRM L2 and FDM L2 (version A)
cswh = lswh .and. (lrm_l2 .or. fdm_l2) .and. version_a

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'alt_cnes', alt, .true.)
call rads_get_var (S, P, 'alt_rate', alt_rate, .true.)
call rads_get_var (S, P, 'range_ku', range, .true.)
call rads_get_var (S, P, 'dry_tropo_ecmwf', dry, .true.)
call rads_get_var (S, P, 'wet_tropo_ecmwf', wet, .true.)
call rads_get_var (S, P, 'inv_bar_static', ib, .true.)
call rads_get_var (S, P, 'iono_gim', iono, .true.)
call rads_get_var (S, P, 'swh_ku', swh, .true.)
call rads_get_var (S, P, 'sig0_ku', sig0, .true.)
call rads_get_var (S, P, 'flags', flagword, .true.)

! Process data records

do i = 1,n
	flag = nint(flagword(i),twobyteint)
	sar = btest(flag,0)

! Apply range bias
! All LRM L2 data processed prior to MAR-2011 have a 3.33 m range bias
! According to Marco Fornari:
! All FDM/LRM L1/L2 data have a CAL1 which is off by 1 FAI (1/64 gate). Thus add 7.3 mm to range.

	if (lrange) then
		if (.not.sar) range(i) = range(i) + fai
		if (lrm_l2 .and. old_version_a) range(i) = range(i) - 3.33d0
	endif

! Apply timing bias.
! Different values apply for SAR and LRM separately.
! See IPF_datation_biases_v4.xlsx by Marco Fornari.

	if (ltbias) then
		if (sar) then
			tbias = +0.520795d-3
		else
			tbias = -4.699112d-3
			if (fdm_l1_v24) tbias = tbias + 4.4436d-3 ! Partial correction of timing bias (See Ruby's e-mail of 22 Apr 2013)
		endif
		time(i) = time(i) + tbias
		alt(i) = alt(i) + tbias * alt_rate(i)
	endif

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

! For LRM/FDM L2 version A data:
! Set SWH to NaN when SWH quality flag (bit 12) is on. Is always zero in product.
! Correct SWH for error in L2 computation. See notes 4 Oct 2011.
! Analysed for LRM L2 data only. Could not be done for FDM, where SWH = NaN always.

	if (cswh) then
		if (btest(flag,12)) then
			swh(i) = nan
		else
			t = swh(i) * abs(swh(i))
			t = (t + 2.7124d0) / 0.5777d0
			swh(i) = sign(sqrt(abs(t)), t)
			if (t < 0d0) swh(i) = -swh(i)
		endif
	endif

! For all data: Correct sigma0 for biases
! For L2 version A: Flip sigma around

	if (.not.lsig0) then
		! Skip
	else if (fdm_l2) then
		if (version_a) then
			sig0(i) = sig0(i) + 4.8d0
			sig0(i) = 21d0 - sig0(i)
		else
			sig0(i) = sig0(i) - 1.4d0
		endif
	else if (lrm_l2) then
		if (version_a) then
			if ((time(i) > 808272000d0 .and. time(i) < 810086400d0) .or. &	! 2010-08-13 00:00 - 2010-09-03 00:00
				(time(i) > 812764800d0 .and. time(i) < 814233600d0)) then	! 2010-10-04 00:00 - 2010-10-21 00:00
				sig0(i) = sig0(i) - 1.4d0
			else if (time(i) > 822787200d0 .and. time(i) < 825464280d0) then! 2011-01-28 00:00 - 2011-02-27 23:38
				sig0(i) = sig0(i) - 40d0
			endif
			sig0(i) = 21d0 - sig0(i)
		else
			sig0(i) = sig0(i) - 7.0d0
		endif
	else	! All L1 data (LRM and PLRM)
		sig0(i) = sig0(i) - 3.0d0
	endif

! Correct for (ad hoc) sigma0 drift

	if (ldrift) sig0(i) = sig0(i) + sig0_drift * (time(i)-time_drift)

! Compute wind speed with ECMWF model

	if (lwind) wind(i) = ecmwf_ws(sig0(i))

! Apply hybrid SSB

	if (lssb) then
		t = swh(i)
		if (t < 0d0) t = 0d0
		if (t > 12d0) t = 12d0
		ssb(i) = grid_lininter (issb_hyb,sig0(i),t)
	endif

enddo

! Prior to Cycle 5 pass 333 all iono is bogus

if (lmeteo .and. lrm_l2 .and. P%equator_time < 808700000d0) then
	iono = nan
	ciono = .true.
endif

! If nothing changed, stop here

if (.not.(ltbias .or. lrange .or. cmeteo .or. ciono .or. lssb .or. cswh .or. lsig0 .or. lwind)) then
	write (*,552) 0
	return
endif

! Update history and define new variables (if required)

P%start_time = time(1)
P%end_time = time(n)
call rads_put_passinfo (S, P)
call rads_put_history (S, P)

if (lssb) call rads_def_var (S, P, 'ssb_hyb')
if (lwind) call rads_def_var (S, P, 'wind_speed_alt')

! Write out all the data

if (ltbias) then
	call rads_put_var (S, P, 'time', time)
	call rads_put_var (S, P, 'alt_cnes', alt)
endif
if (lrange) call rads_put_var (S, P, 'range_ku', range)
if (cmeteo) then
	call rads_put_var (S, P, 'dry_tropo', dry)
	call rads_put_var (S, P, 'wet_tropo', wet)
	call rads_put_var (S, P, 'inv_bar_static', ib)
endif
if (ciono) call rads_put_var (S, P, 'iono_gim', iono)
if (lssb) call rads_put_var (S, P, 'ssb_hyb', ssb)
if (cswh) call rads_put_var (S, P, 'swh_ku', swh)
if (lsig0 .or. ldrift) call rads_put_var (S, P, 'sig0_ku', sig0)
if (lwind) call rads_put_var (S, P, 'wind_speed_alt', wind)

write (*,552) n
end subroutine process_pass

end program rads_fix_c2
