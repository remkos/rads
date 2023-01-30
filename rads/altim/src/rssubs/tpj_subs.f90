module tpj_subs
contains

subroutine tpj_yawmode (sat, time, ymode, pbias, pitch)
use typesizes
character(2) :: sat
real(eightbytereal), intent(in) :: time
integer(fourbyteint), intent(out) :: ymode
real(eightbytereal), intent(out) :: pbias, pitch
!
! This function returns the yaw steering mode for a given time and satellite
!
! Input:
!  sat   : Use 'tx', 'pn', 'j1', or 'j2' for TOPEX, Poseidon, Jason-1 or Jason-2
!  time  : Time in seconds since 1985
!
! Output:
!  ymode : yaw steering mode
!          0 = unknown
!          1 = fixed yaw at low beta prime (< 15 deg)
!          3 = sinusoidal yaw
!          5 = fixed yaw at high beta prime (> 80 deg)
!          (ymode is negative for negative beta', flying backward)
!  pbias : bias of the solar array pitch (to be added to nominal pitch angle) (deg)
!  pitch : space craft pitch angle (deg)
!-
! 15-Apr-2010 : Created by Remko Scharroo (Altimetrics LLC)
!-----------------------------------------------------------------------
character(2) :: oldsat = '  '
character(120) :: line
integer(fourbyteint), parameter :: maxmodes = 600
integer(fourbyteint) :: i, nmodes, yawmode(maxmodes), yy, dd, hh, mm, ss, unit, ios, &
	cyc, pass
real(eightbytereal) :: yawtime(maxmodes), sapitch(maxmodes), scpitch(maxmodes), t, t0
logical :: opened
save oldsat, i, nmodes, yawmode, yawtime, sapitch, scpitch

! Load satellite attitude file

if (sat /= oldsat) then
	! Start with a dummy record
	nmodes = 1
	yawmode(1) = 0
	yawtime(1) = -1d30
	sapitch(1) = 0d0
	scpitch(1) = 0d0
	ymode = 1
	pbias = 0d0
	pitch = 0d0
	t0 = 0d0

	! Open the appropriate attitude file
	call getenv ('ALTIM',line)
	i = len_trim (line)
	if (sat == 'tx' .or. sat == 'pn') then
		line(i+1:) = '/data/tables/topatt.txt'
	else if (sat == 'j1') then
		line(i+1:) = '/data/tables/ja1att.txt'
	else if (sat == 'j2') then
		line(i+1:) = '/data/tables/ja2att.txt'
	else
		stop 'Unknown satellite'
	endif
	do unit=99,7,-1
		inquire (unit=unit,opened=opened)
		if (.not.opened) exit
	enddo
	open (unit, file=line, status='old')

	! Look for mode changes in the attitude file
	do
		read (unit,'(a)',iostat=ios) line
		if (ios /= 0) exit

		i = index(line,':') ! Find delimiter between hh:mm
		if (i <= 0) cycle
		read (line(i-11:),'(i4,1x,i3,3(1x,i2),1x,2i6)',iostat=ios) yy,dd,hh,mm,ss,cyc,pass
		if (ios /= 0) cycle
		yy = yy - 1985
		dd = dd + yy * 365 + yy / 4 - 1
		t = dd*86400d0 + hh*3600d0 + mm*60d0 + ss

		if (index(line,'Yaw -90 deg') > 0) then
			ymode = -5
		else if (index(line,'Start Switch to 90 Degree Fixed') > 0 .or. &
			index(line,'Start 90 fixed yaw') > 0) then
			ymode = +5
		else if (index(line,'to Fixed') > 0 .and. index(line,'P=-') > 0) then
			ymode = -1
		else if (index(line,'to Fixed') > 0 .and. index(line,'P=') > 0) then
			ymode = +1
		else if ((index(line,'to Sinusoid') > 0 .and. index(line,'P=-') > 0) .or. &
			index(line,'Stop  Fixed Yaw') > 0) then
			ymode = -3
		else if (index(line,'to Sinusoid') > 0 .and. index(line,'P=') > 0 .or. &
			index(line,'Start Sinusoidal Yaw Steering') > 0) then
			ymode = +3
		else if (index(line,'Start Yaw Flip') > 0) then
			t0 = t  ! Just save the start time
			cycle
		else if (index(line,'Stop  Yaw Flip') > 0) then
			if (abs(ymode) /= 1) stop 'Yaw flip occurred when not in fixed yaw'
			t = (t0 + t) / 2d0	! Average start and stop time
			ymode = -ymode
		else if (index(line,'Stop  Fixed') > 0) then
			ymode = -3	! Case for TX not picked up by previous
		else if (index(line,'Normal') > 0) then
			! Return to normal operations, we assume same yaw mode and pitch bias as before
		else if (index(line,'Stop  Solar Array Ramp to Bias Position') > 0 .or. &
			index(line,'Stop  Set Solar Array') > 0) then
			! Extract the solar array pitch bias
			i = index(line,'deg')
			if (i > 0) line(i:i+2) = ' ' ! Remove the word 'deg'
			read (line(40:47),*) pbias
		else if (index(line,'Stop  Solar Array at 0 Deg Bias') > 0 .or. &
			index(line,'Stop  Solar Array Ramp to Zero deg Bias') > 0 .or. &
			index(line,'Stop  Solar Panel Ramp to Null Position') > 0) then
			pbias = 0d0
		else if (line(:10) == 'Pitch Bias') then
			i = index(line,'=')
			read (line(i+1:i+6),*) pitch
		else
			cycle
		endif
		nmodes = nmodes + 1
		if (nmodes > maxmodes) stop 'To many yaw mode switches'
		if (nmodes > 2 .and. abs(ymode-yawmode(nmodes-1)) > 2) then
			! Changed too much since last. Correct.
			if (abs(ymode) == 3) yawmode(nmodes-1) = ymode / abs(ymode)	! Change +/-3 to +/-1
		endif
		yawtime(nmodes) = t
		yawmode(nmodes) = ymode
		sapitch(nmodes) = pbias
		scpitch(nmodes) = pitch
		oldsat = sat
	enddo
	close (unit)

	! Add extra dummy record
	nmodes = nmodes + 1
	yawmode(nmodes) = 0
	yawtime(nmodes) = 1d30
	sapitch(nmodes) = 0d0
	scpitch(nmodes) = 0d0
	i = 1
endif

! Look for correct time interval
do
	if (time < yawtime(i)) then
		i = i - 1
		cycle
	else if (time > yawtime(i + 1)) then
		i = i + 1
		cycle
	else
		exit
	endif
enddo

! Return yaw mode, SA pitch bias and SC pitch angle
ymode = yawmode(i)
pbias = sapitch(i)
pitch = scpitch(i)
end subroutine tpj_yawmode

subroutine tpj_steering (alpha, betap, ymode, yaw, pitch)
use typesizes
real(eightbytereal), intent(in) :: alpha, betap
integer(fourbyteint), intent(in) :: ymode
real(eightbytereal), intent(out) :: yaw, pitch
!
! This function computes the yaw of the satellite and the pitch of the solar
! array of the TOPEX and Jason spacecraft, depending on the beta' angle of
! the sun and the in-orbit angle since "06:00" solar time along the orbit.
! The results do not include any optional yaw or pitch biases.
!
! Input:
!  alpha  : in-orbit angle since the point of 6:00 LST (radians)
!  betap  : angle of the sun above the orbital plane (radians)
!  ymode  : yaw steering mode
!           0 = determine automatically from beta prime
!           1 = fixed yaw at low beta prime (< 15 deg)
!           3 = sinusoidal yaw
!           5 = fixed yaw at high beta prime (> 80 deg)
!           (use negative for negative betap, flying backward)
!
! Output:
!  yaw    : yaw angle of the satellite (radians)
!  pitch  : pitch angle of the solar array (radians)
!           (not including a possible pitch bias)
!-
! 15-Apr-2010 : Created by Remko Scharroo (Altimetrics LLC)
!-----------------------------------------------------------------------
real(eightbytereal), parameter :: pi = 4d0*atan(1d0), tpi = 2d0*pi, hpi = pi/2d0, rad = pi/180d0
real(eightbytereal) :: x
integer(fourbyteint) :: is

! Determine yaw steering mode automatically when ymode == 0

if (ymode == 0) then
	x = abs(betap/rad)
	if (x < 15d0) then
		is = 1
	else if (x > 80d0) then
		is = 5
	else
		is = 3
	endif
	if (betap < 0d0) is = -is
else
	is = ymode
endif

! Determine yaw

select case (is)
case (1)
	yaw = 0d0
case (-1)
	yaw = pi
case (5)
	yaw = hpi
case (-5)
	yaw = -hpi
case (3)
	yaw = hpi + (hpi - betap) * cos(alpha)
case (-3)
	yaw = -hpi - (hpi + betap) * cos(alpha)
end select

! Determine pitch of the solar array

select case (is)
case (1)
	pitch = modulo (alpha, tpi)
case (-1)
	pitch = modulo (pi - alpha, tpi)
case (3, 5)	! Formulas suggested by J2 telemetry (without the + pi)
	pitch = - (hpi - betap) * sin(alpha) + pi
case (-3, -5)
	pitch = - (hpi + betap) * sin(alpha) + pi
case default
	! This is the "official" formula from Marshall, 1992 (without the + pi)
	! However, J2 telemetry suggests the simpler formulas above
	x = cos(betap)
	pitch = atan ((sin(alpha) * x) / (cos(alpha) * x * cos(yaw) - sin(yaw) * sin(betap))) + pi
end select
end subroutine tpj_steering

subroutine tpj_cgcorr (sat, eqtime, eqlon, asc, time, sunlit, alpha, betap, sapitch, scpitch, ymode, yaw, pbias)
use typesizes
use solar_subs
character(*), intent(in) :: sat
real(eightbytereal), intent(in) :: eqtime, eqlon, time
logical, intent(in) :: asc
logical, intent(out) :: sunlit
real(eightbytereal), intent(out) :: alpha, betap, sapitch, scpitch, yaw, pbias
integer(fourbyteint), intent(out) :: ymode
!
! This function computes the ingredients needed for the center-of-gravity
! correction for the movement of the solar array and spacecraft pitch
! of TOPEX/Poseidon or Jason-1 or Jason-2, based on
! their respective yaw mode, steering laws and solar array pitch angle.
!
! Input:
!  sat     : Satellite abbreviation ('tx', 'pn', 'j1', 'j2')
!  eqtime  : Time of equator crossing (sec85)
!  eqlon   : Longitude of equator crossing (deg)
!  asc     : TRUE if equator crossing is ascending
!  time    : Time (sec85)
!
! Output:
!  sunlit  : TRUE if satellite is in sunlight
!  alpha   : in-orbit angle since the point of 6:00 LST (deg)
!  betap   : angle of the sun above the orbital plane (deg)
!  sapitch : Solar array pitch angle (deg)
!  scpitch : Spacecraft pitch angle (deg)
!  ymode   : yaw steering mode
!            0 = unknown
!            1 = fixed yaw at low beta prime (< 15 deg)
!            3 = sinusoidal yaw
!            5 = fixed yaw at high beta prime (> 80 deg)
!            (ymode is negative for negative beta', flying backward)
!  yaw     : yaw angle (deg)
!  pbias   : Solar array pich bias (deg)
!-
! 26-Apr-2010 : Created by Remko Scharroo (Altimetrics LLC)
!-----------------------------------------------------------------------
real(eightbytereal), parameter :: pi = 4d0*atan(1d0), rad = pi/180d0, tpi = 2d0*pi, &
	incl = 66.04d0*rad, tpass = 3373d0, re = 6375d0, r = re + 1350d0, fact = sqrt(1d0 - (re/r)**2)
real(eightbytereal) :: mjd, alphap, raan

! Convert time to MJD
mjd = eqtime / 86400d0 + 46066d0

! Convert equator longitude to raan
raan = eqlon * rad + gmst (mjd)
if (.not.asc) raan = raan + pi

! Get alpha' and beta' angles
call sun_orbit_geometry (mjd, raan, incl, alphap, betap)

! Determine in-orbit angle since "6AM" position (alpha)
alpha = (time - eqtime) / tpass * pi - alphap
if (.not.asc) alpha = alpha + pi	! Adjust descending to ascending node
alpha = modulo(alpha, tpi)

! Determine of satellite is in sunlight
sunlit = (sin(alpha) > -fact / cos(betap))

! Get yaw and pitch angles
call tpj_yawmode (sat, time, ymode, pbias, scpitch)
call tpj_steering (alpha, betap, ymode, yaw, sapitch)
sapitch = modulo(sapitch/rad + pbias, 360d0)
alpha = alpha / rad
betap = betap / rad
yaw = yaw / rad

end subroutine tpj_cgcorr
end module tpj_subs

!----------------------------------------------------------------------
!*tpj_cgcorr_f77 -- Fortran 77 interface to tp_cgcorr
!+
subroutine tpj_cgcorr_f77 (sat, eqtime, eqlon, asc, time, sunlit, alpha, betap, sapitch, scpitch, ymode, yaw, pbias)
use typesizes
use tpj_subs
character(*), intent(in) :: sat
real(eightbytereal), intent(in) :: eqtime, eqlon, time
logical, intent(in) :: asc
logical, intent(out) :: sunlit
real(eightbytereal), intent(out) :: alpha, betap, sapitch, scpitch, yaw, pbias
integer(fourbyteint), intent(out) :: ymode
call tpj_cgcorr (sat, eqtime, eqlon, asc, time, sunlit, alpha, betap, sapitch, scpitch, ymode, yaw, pbias)
end subroutine tpj_cgcorr_f77
