module meteo_subs
contains

!*pp_isa - International Standard Atmosphere correction to pressure
!+
elemental function pp_isa (h)
use typesizes
real(eightbytereal) :: pp_isa
real(eightbytereal), intent(in) :: h
!
! This routine computes the scale factor to sea level pressure in order
! to convert it to pressure at altitude. In other words:
!
! pressure(h) = pressure(0) * pp_isa(h)
!
! References:
! G. Gyatt (2011), The Standard Atmosphere, http://www.atmosculator.com/The%20Standard%20Atmosphere.html
! U.S. Standard Atmosphere, 1976 (1976) NOAA-S/T 76-1562, NOAA/NASA/USAF, Washington, D.C.
!
! Input:
!  h      : altitude (m)
!
! Output:
!  pp_isa : scale factor to sea level pressure to account for altitude
!-
! 17-Nov-2011 - Remko Scharroo (c) Altimeterics LLC
!-----------------------------------------------------------------------
real(eightbytereal), parameter :: h0 = 0d0        ! Reference altitude (m)
real(eightbytereal), parameter :: l0 = -0.0065d0  ! Temperature lapse rate (K/m)
real(eightbytereal), parameter :: T0 = 288.15d0   ! Mean temperature at sea level (K)
real(eightbytereal), parameter :: R = 8.31432d0   ! Universal gas constant (J/mole/K)
real(eightbytereal), parameter :: M = 0.0289644d0 ! Mean molucular mass of air (kg/mole)
real(eightbytereal), parameter :: g0 = 9.80665d0  ! Gravitational acceleration at sea level (m/s^2)
real(eightbytereal), parameter :: l0T0 = l0 / T0, coeff = -g0 / (l0 * R / M)

pp_isa = (1d0 + (h - h0) * l0t0) ** coeff
end function pp_isa

!*pp_hop - Correction to pressure according to Hopfield (1969)
!+
elemental function pp_hop (h, lat, t2m)
use typesizes
real(eightbytereal) :: pp_hop
real(eightbytereal), intent(in) :: h, lat, t2m
!
! This routine computes the scale factor to sea level pressure in order
! to convert it to pressure at altitude. This takes into account variations
! in gravity as function of latitude as well as variations of sea level
! temperature (t2m)
!
! In other words:
!
! pressure(h) = pressure(0) * pp_hop (h, lat, t2m)
!
! References:
! Fernandes, M. J., N. Pires, C. Lázaro, and A L. Nunes (2012), Tropospheric delays fro
! GNSS for application in coastal altimetry, Adv. Space Res., 10.1016/j.asr.2012.04.025, in press.
!
! Input:
!  h      : altitude (m)
!  lat    : latitude (deg)
!  t2m    : temperature at 2 meters above the sea surface (K)
!
! Output:
!  pp_hop : scale factor to sea level pressure to account for altitude
!-
! 23-Jul-2012 - Remko Scharroo (c) Altimeterics LLC
!-----------------------------------------------------------------------
real(eightbytereal), parameter :: h0 = 0d0        ! Reference altitude (m)
real(eightbytereal), parameter :: l0 = -0.0065d0  ! Temperature lapse rate (K/m)
real(eightbytereal), parameter :: R = 8.31432d0   ! Universal gas constant (J/mole/K)
real(eightbytereal), parameter :: M = 0.0289644d0 ! Mean molucular mass of air (kg/mole)
real(eightbytereal), parameter :: g0 = 9.784d0    ! Gravitational acceleration at sea level (m/s^2)
real(eightbytereal), parameter :: rad = atan(1d0)/45d0 ! Conversion from degrees to radians
real(eightbytereal), parameter :: RM = R/M, rad2 = 2d0 * rad
real(eightbytereal) :: gm, tm

tm = t2m + 0.5d0 * l0 * h ! Mean temperature
gm = g0 * (1d0 - 266d-5 * cos (lat * rad2) - 28d-8 * h) ! Mean gravity

pp_hop = exp (-gm * (h-h0) / RM / tm)
end function pp_hop

!*wind_ecmwf -- ECMWF wind speed model
!+
elemental function wind_ecmwf (sigma0, ka_band)
use typesizes
real(eightbytereal), intent(in) :: sigma0
logical, intent(in), optional :: ka_band
real(eightbytereal) :: wind_ecmwf
!
! This is a 1-D wind speed model relating backscatter (sigma0) to wind speed at
! 10 meters altitude. It was developed at ECMWF by Saleh Abdalla to support the
! Envisat mission.
!
! To support also the SARAL processing, the optional parameter <ka_band> was added,
! which will select parameters appropriate for Ka-band instead of Ku-band.
!
! References:
! Abdalla, S., Ku-band radar altimeter surface wind speed algorithm,
!  in Proc. of the 2007 Envisat Symposium, Montreux, Switzerland, 23-27 April 2007,
!  Eur. Space Agency Spec. Publ., ESA SP-636, 2007.
! Abdalla, S., Ku-band radar altimeter surface wind speed algorithm,
!  Mar. Geod., 35(sup1), 276–298, doi: 10.1080/01490419.2012.718676, 2012.
! Lillibridge, J. L., R. Scharroo, S. Abdalla, and D. C. Vandemark, One- and
!  two-dimensional wind speed models for Ka-band altimetry, J. Atmos. Oceanic
!  Technol., 31(3), 630-638, doi:10.1175/JTECH-D-13- 00167.1, 2014.
!
! Input:
!  sigma0 : backscatter coefficient (dB)
!  ka_band: backscatter is from Ka-band instead of Ku-band
!
! Output:
!  wind_ecmwf : wind speed at 10 meters altitude (m/s)
!-
real(eightbytereal) :: a, b, c, d, sigmab
real(eightbytereal) :: um, ux
logical :: ka

if (present(ka_band)) then
	ka = ka_band
else
	ka = .false.
endif
if (ka) then ! Ka-band values (c and sigmab modified for continuity)
	a = 34.2d0
	b = 2.48d0
	c = 711.6d0
	d = 0.42d0
	sigmab = 11.409d0
else ! Ku-band values
	a = 46.5d0
	b = 3.6d0
	c = 1690d0
	d = 0.5d0
	sigmab = 10.917d0
endif

! Step 1: First approximation

if (sigma0 /= sigma0) then
	! If SIGMA0 is NaN, return NaN
	wind_ecmwf = sigma0
	return
else if (sigma0 <= sigmab) then
	! Linear portion of the model
	um = a - b * sigma0
else
	! Exponential portion of the model
	um = c * exp (-d * sigma0)
endif

! Step 2: Fine tuning

ux = um ** 0.096d0
wind_ecmwf = um + 1.4d0 * ux * exp (-0.32d0 * um * ux)

end function wind_ecmwf

!*nn_l2_mwr -- Neural Network Envisat MWR wet tropospheric model
!+
function nn_l2_mwr (tb23, tb36, sigma0, qty)
use typesizes
real(eightbytereal) :: nn_l2_mwr
real(eightbytereal), intent(in) :: tb23, tb36, sigma0
integer(fourbyteint), intent(in) :: qty

! This routine computes each of the 5 radiometric level 2 (L2)
! quantities in the same way as they are computed in Envisat products,
! i.e., using the coefficients obtained from a neural network study
! of the dependence of these 5 radiometric level 2 quantities on
! both brightness temperatures and Ku-band Sigma0.
!
! Input to the routine are 23.8 GHz and 36.5 GHz brightness temperature
! in Kelvin and (uncorrected) altimeter Ku-band backscatter coefficient
! in dB: TB23, TB36 and SIGMA0. The uncorrected backscatter coefficient
! is obtained by subtracting the backscatter atmospheric attenuation
! from the backscatter coefficient.
!
! Depending on the argument QTY, the routine returns atmospheric
! attenuation on the Ku- or S-band sigma0, wet tropospheric correction,
! water vapour conent or liquid water content.
!
! In case of an input error, the routine returns 1D30.
!
! Input arguments:
!  tb23   : 23.8 GHz brightness temperature (K)
!  tb36   : 36.5 GHz brightness temperature (K)
!  sigma0 : (Uncorrected) altimeter backscatter coefficient (dB)
!  qty    : Output quantity:
!           =1 : Atmospheric attenuation on Ku-band sigma0 (dB)
!           =2 : Atmospheric attenuation on S-band sigma0 (dB)
!           =3 : Radiometer-derived wet tropospheric correction (m)
!           =4 : Radiometer-derived water vapour content (g/cm^2)
!           =5 : Radiometer-derived liquid water content (mg/cm^2)
!
! Return value:
!  nn_l2_mwr : One of the quantities as specified by QTY, or 1D30
!              upon error.
!
! Notes:
! 1) The original routine required 1 dB to be added to the backscatter
!    coefficient before calling the routine. This 1 dB bias is now included
!    in routine by lowering dw_InputMean(3) by 1.
! 2) To use this routine with ERS-2, the values need to be biased
!    in order to match Envisat values. Simply use
!
!    NN_L2_MWR ( TB23' + 3.37d0, TB36 + 2.36d0, SIGMA0 + 0.21d0 )
!
!    where TB23' is the corrected brightness temperature according to NW6:
!
!    TB23' = 19.203 + 0.93895 TB23 (beyond 5 years after launch)
!
! Reference:
! S. Labroue and E. Obligis,
! Neural network retrieval algorithms for the Envisat/MWR,
! Tech. Rep. CLS/DOS/NT/03.848, CLS, Ramonville St.Agne, 2003
!-
! 10-Apr-2014 - Integrated into meteo_subs
! 14-Aug-2008 - Make sure NaN on input results in NaN on output
! 24-Jun-2008 - Converted from Fortran 77 to 90
! 02-Feb-2005 - Return NaN on out-of-bound, not 1d30
! 13-Sep-2004 - Extended manual, no chance in functionality
! 10-Feb-2004 - Changed manual only
! 29-Dec-2003 - Decreased CPU by replacing exp() calls with tanh()
! 11-Dec-2003 - Corrected, simplified and improved by Remko Scharroo (NOAA)
! 13-Mar-2003 - Original NN_L2_MWR by M. van den Bossche (CLS)
!-----------------------------------------------------------------------

! Neural Network Coefficients

real(eightbytereal) :: dw_InputMean(3) = (/ 177.90316158d0, 173.45177492d0, 11.97140805d0 /)

real(eightbytereal) :: dw_InputStd(3) = (/ 21.16395751d0,  13.33997377d0,  1.08309444d0 /)

real(eightbytereal) :: dw_Weights1(8,3,5) = reshape((/ &
      0.44696310d0,   0.72220056d0,  -0.80716305d0,  0.20325626d0, & ! 1
     -0.62522626d0,  -0.23312843d0,   0.33566595d0,  1.35234027d0, &
     -0.26862238d0,   2.05446692d0,   0.11175954d0,  0.05766287d0, &
      1.23925023d0,   0.34652949d0,  -0.24864995d0, -1.49118535d0, &
      0.14130567d0,   0.36519738d0,  -0.10754916d0,  0.05025406d0, &
     -0.00462589d0,   0.02103813d0,  -0.10220144d0,  1.85671196d0, &
      0.90454448d0,   0.10071430d0,   1.68904800d0,  2.81384595d0, & ! 2
      0.52859072d0,   3.94392507d0,   6.77282432d0, -0.62388637d0, &
     -2.79581500d0,  -0.14916371d0,  -0.43745905d0, -1.11521141d0, &
     -0.90535679d0,  -3.98263270d0,  -9.61069818d0, -0.58687942d0, &
     -0.01235963d0,  -0.01202587d0,   0.65426505d0, -0.78274514d0, &
      0.03443161d0,   5.10835392d0,   0.03524043d0, -0.07082674d0, &
     -0.89552820d0,  -0.03297866d0,   0.25545500d0,  0.52267161d0, & ! 3
      0.72511737d0,  -0.47927862d0,  -0.69027030d0,  0.78903778d0, &
      1.34534555d0,   0.28625740d0,  -0.70072999d0, -0.14315618d0, &
     -1.37795694d0,   0.31775696d0,   1.14886824d0, -0.34057041d0, &
     -0.42812062d0,  -0.00633370d0,   0.00449982d0, -0.01506284d0, &
      0.70438201d0,  -0.03733946d0,  -0.54311885d0, -0.21026933d0, &
     -0.54139486d0,   0.37818313d0,  -0.26754756d0,  0.44947685d0, & ! 4
      1.75577079d0,  -0.32956621d0,   0.02547023d0, -1.19510396d0, &
      0.55832606d0,   0.02217675d0,  -0.10928559d0, -0.05580770d0, &
     -2.37684429d0,   0.38835569d0,   0.27645473d0,  0.52061306d0, &
     -0.03588430d0,  -0.03220448d0,  -0.22025970d0,  0.06771010d0, &
      1.03585041d0,  -0.04628975d0,   0.00015090d0,  0.20231264d0, &
      0.72631347d0,   0.96958809d0,   1.14603101d0, -0.12408819d0, & ! 5
      0.52902469d0,   0.33981985d0,  -0.10472735d0, -0.40687924d0, &
     -0.29259677d0,  -0.46276556d0,  -0.59018200d0,  0.33034467d0, &
     -1.31313582d0,  -1.00623029d0,   0.27860956d0,  0.16612519d0, &
      0.01502955d0,   0.13737705d0,   0.11857027d0,  0.00012563d0, &
     -0.02835980d0,   0.04447786d0,  -0.32599584d0,  0.11527573d0 /),(/8,3,5/))

real(eightbytereal) :: dw_Weights2(8,5) = reshape((/ &
     -6.01521638d0,   0.05673385d0,   0.86818533d0,  4.31398179d0, &
      8.80848855d0,   6.22459180d0,   7.77021534d0, -6.46198696d0, &
    -19.10855900d0, -21.48264157d0,  -0.13681912d0,  0.14754153d0, &
    -25.75343810d0,  -7.13295436d0, -17.55866459d0,  6.54912774d0, &
     -5.54265727d0,   3.21172540d0,   0.85245131d0,  6.26893361d0, &
      0.80490100d0,  -2.92601319d0,   2.95899176d0, -1.53422135d0, &
     -3.88755082d0,   1.82713759d0,   1.83008551d0,  6.45017814d0, &
      3.11388202d0,  -1.41100836d0,   2.02299006d0,  0.22977112d0, &
     -6.32733950d0,   2.12955323d0,  -2.04922147d0,  7.19951870d0, &
     -0.35391656d0,  -9.57243059d0,  -1.04047189d0, -3.04045926d0 /),(/8,5/))

real(eightbytereal) :: dw_Bias1(8,5) = reshape((/ &
     -2.11115274d0,   1.17852961d0,  -0.16826516d0, -0.06023521d0, &
     -7.13894675d0,  -1.42349247d0,   1.93617157d0, 14.24175815d0, &
     -4.58268849d0,   1.08658091d0,   0.38275346d0, -1.61448894d0, &
      5.69683964d0,  36.99539002d0,  -8.77354719d0, -3.15340641d0, &
     -7.98464772d0,   0.53395655d0,  -0.23145209d0, -1.40254307d0, &
      7.80384869d0,  -1.25261261d0,  -6.85063315d0, -1.95532644d0, &
     -2.64524836d0,   0.49913280d0,   1.58748752d0, -1.56560119d0, &
     15.52012851d0,  -0.59742264d0,  -1.27245596d0, -0.69086003d0, &
     -2.27087343d0,   0.29229745d0,   0.22048456d0, -1.31376626d0, &
      0.27622941d0,   6.16186925d0,  -1.98966909d0, -1.94396380d0 /),(/8,5/))

real(eightbytereal) :: dv_Bias2(5) = (/ 7.84544738d0, 19.83645336d0,  -3.25179309d0,  -2.54502683d0,  5.64953122d0 /)

real(eightbytereal) :: dv_OutputMean(5) = (/ 0.23336214d0, 0.07527797d0, -18.14563614d-2,  2.99711125d0, 13.44185845d0 /)

real(eightbytereal) :: dv_OutputStd(5) = (/ 0.09925092d0, 0.00498339d0,  -9.47044771d-2,  1.59653331d0, 25.06408418d0 /)

! Local variables

integer(fourbyteint) :: j
real(eightbytereal) :: dv_X, dw_P(3), dw_A(8)
save

! Input data normalisation

dw_P(1) = (tb23   - dw_InputMean(1)) / dw_InputStd(1)
dw_P(2) = (tb36   - dw_InputMean(2)) / dw_InputStd(2)
dw_P(3) = (sigma0 - dw_InputMean(3)) / dw_InputStd(3)

! First neural network layer summation

do j=1,8
	dv_X = dw_Bias1(j,qty) + dot_product(dw_Weights1(j,:,qty), dw_P(:))
	if (dv_X >= -700d0 .and. dv_X <= 700d0) then
		dw_A(j) = tanh(dv_X)
	else
		! Also captures NaNs on input
		nn_l2_mwr = 0d0
		nn_l2_mwr = nn_l2_mwr / nn_l2_mwr   ! Create NaN
		return
	endif
enddo

! Second neural network layer summation

dv_X = dv_Bias2(qty) + dot_product(dw_A(:), dw_Weights2(:,qty))

! Denormalization of output quantity

nn_l2_mwr = dv_OutputMean(qty) + dv_OutputStd(qty) * dv_X
end function nn_l2_mwr

!*wind_j1 -- Convert sigma0 to wind and vv using Collard 2-parameter model
!+
function wind_j1 (i, val, swh)
use typesizes
real(eightbytereal) :: wind_j1
real(eightbytereal), intent(in) :: val, swh
integer(fourbyteint), intent(in) :: i
!
! This function computes the wind speed (U10, referenced to 10 m above
! sea level) from the JASON altimeter backscatter coefficient (SIGMA0) and
! significant wave height (SWH). The reverse is also possible:
! converting U10 and SWH to SIGMA0.
! The parameter i selects whether the forward (i=1) or backward
! (i=2) conversion is requested. NOTE: backward conversion is NOT
! CORRECTLY IMPLEMENTED.
!
! The unit for SIGMA0 is dB, for SWH is m and for U10 is m/s.
!
! The respective models F1(SIGMA0,SWH) and F2(U10,SWH) are described in
! Gourrion et al. (2002) and adapted by Collard.
!
! Input arguments:
!  i   : Select direction of conversion
!        1 = Convert SIGMA0 and SWH to U10 (f1)
!        2 = Convert U10 and SWH to SIGMA0 (f2)
!  val : Either SIGMA0 in dB (I=1) or U10 in m/s (I=2)
!  swh : Significant wave height in m
!
! Returned value:
!  wind_j1: Either U10 in m/s (I=1) or SIGMA0 in dB (I=2)
!
! References:
!
! Gourrion, J., D. Vandemark, S. Bailey, B. Chapron, G. P. Gommenginger,
! P. G. Challenor, and M. A. Srokosz,
! A two-parameter wind speed algorithm for Ku-band altimeters,
! J. Atmos. Ocean. Technol., 19(12), 2030-2048, 2002.
!
! Collard, F. Algo Etude BOOST. Rapport BO-021-CLS-0407-RF.
!-
! $Log: meteo_subs.f90,v $
! Revision 1.11  2017/09/27 18:46:36  rads
! Add proper reference
!
! Revision 1.10  2015/12/04 11:16:20  rads
! Documentation update only.
!
! Revision 1.9  2014/09/08 16:05:08  rads
! - Add wind_j1
!
! Revision 1.1  2008/11/18 17:26:41  rads
! - Initial revision based on f1f2.f
!
!-----------------------------------------------------------------------
real(eightbytereal) :: a(3) = (/-0.6322581d0,0.0869731d0,0.0991799d0/)
real(eightbytereal) :: b(3) = (/ 0.0806452d0,0.0321528d0,0.0683395d0/)
! Note: the order of the second and third element of a and b is intensionally switched
real(eightbytereal) :: wx(2,2,2) = reshape((/58.76821d0, 22.08226d0, 3.05118d0, 0.18608d0, &
	-43.39541d0, -6.92550d0, 2.78612d0, 1.22293d0/),(/2,2,2/))
real(eightbytereal) :: wy(2,2) = reshape((/ -0.20367d0, -22.35999d0, 1.18281d0,-3.30096d0/),(/2,2/))
real(eightbytereal) :: bx(2,2) = reshape((/ -33.61161d0, 1.06796d0, 7.83459d0,-1.46489d0/),(/2,2/))
real(eightbytereal) :: by(2) = (/ 20.10259d0, 1.13906d0/)
save a,b,wx,wy,bx,by
real(eightbytereal) :: x(2),y,p(2)
logical :: isnan

! If input values are NaN, or incomprehensible input, return NaN

if (isnan(val) .or. isnan(swh) .or. i < 1 .or. i > 2) then
    wind_j1 = 0d0
    wind_j1 = wind_j1 / wind_j1
    return
endif

! Normalise input values

p(1) = a(i) + b(i) * val
p(2) = a(3) + b(3) * swh

! Determine X vector

x = wx(1,:,i) * p(1) + wx(2,:,i) * p(2) + bx(:,i)
x = 1d0 / (1d0 + exp(-x))

! Determine normalized output Y

y = wy(1,i) * x(1) + wy(2,i) * x(2) + by(i)
y = 1d0 / (1d0 + exp(-y))

! Return unnormalized value

wind_j1 = (y - a(3-i)) / b(3-i)

end function wind_j1

end module meteo_subs
