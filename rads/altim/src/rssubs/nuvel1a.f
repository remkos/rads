**NUVEL1A -- Compute station velocities according to NUVEL1A
*+
      SUBROUTINE NUVEL1A (PLATE, XYZ, UVW)
      INTEGER*4 PLATE
      REAL*8    XYZ(3), UVW(3)

* This routine determines station velocities according to the NUVEL1A
* NNR plate motion model.
*
* Arguments:
*   PLATE (input) : Number of the reference plate
*   XYZ   (input) : XYZ coordinates of station (m)
*   UVW  (output) : XYZ velocities according to NUVEL1A (m/s)
*-
*  9-Sep-1998 - Created by Remko Scharroo
* 13-Mar-2001 - Version without external table
*-----------------------------------------------------------------------
      integer mplate
      real*8 omega(3,16),fact
      parameter (mplate=16,fact=1d9*365.25d0*86400d0)

* Following values are rotational velocities for each plate in
* 1D-9 radians per year

      data omega /
     |    -0.992,   -2.387,    3.153,	! EURA
     |     0.249,   -3.591,   -0.153,	! NOAM
     |     7.830,    5.132,    6.283,	! AUST
     |    -1.519,    4.847,   -9.970,	! PCFC
     |     0.882,   -3.092,    3.922,	! AFRC
     |    -1.047,   -1.508,   -0.870,	! SOAM
     |    -0.830,   -1.694,    3.707,	! ANTA
     |     6.676,   -0.514,    6.761,	! ARAB
     |    -0.187,   -3.377,    1.580,	! CARB
     |   -10.434,  -21.597,   10.925,	! COCO
     |    -1.542,   -8.570,    9.609,	! NAZC
     |     6.661,    0.047,    6.790,	! INDI
     |     5.442,    8.746,   -5.678,	! JUFU
     |    10.158,   -7.166,  -10.320,	! PHIL
     |    -9.392,  -30.951,   12.054,	! RIVR
     |    -0.391,   -2.609,   -1.263/	! SCOT
      save omega

* If plate number out of bounds, return all zeros

      if (plate.lt.1 .or. plate.gt.mplate) then
         uvw(1)=0d0
	 uvw(2)=0d0
	 uvw(3)=0d0
      else

* Compute omega * xyz and scale by 1D9 and seconds in a year to obtain m/s

         uvw(1)=(omega(2,plate)*xyz(3)-omega(3,plate)*xyz(2))/fact
         uvw(2)=(omega(3,plate)*xyz(1)-omega(1,plate)*xyz(3))/fact
         uvw(3)=(omega(1,plate)*xyz(2)-omega(2,plate)*xyz(1))/fact
      endif
      end
