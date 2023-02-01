**J2000 -- Transform ephemerides between J2000 and TOD
*+
      SUBROUTINE J2000 (MJD, A, B, NDIM, ICH)
      implicit none
      REAL*8 MJD, A(*), B(*)
      INTEGER*4 NDIM, ICH
*
* This routine transforms ephemerides from J2000 to True-Of-Date or vice
* versa. A and B are state vectors (position and/or velocity)
* in the inertial reference frame. NDIM indicates whether the vectors
* contain 3 or 6 coordinates. ICH regulates the transformation direction
* (to or from J2000). The epoch is given by MJD (UTC).
* 
* A and B may share the same array in the calling routine.
*
* Arguments:
*  MJD  (input): Mean Julian Date in UTC.
*  A    (input): Input state vector.
*  B   (output): Transformed state vector.
*  NDIM (input): Dimension of the state vector (3 or 6).
*  ICH  (input): Direction of the transformation:
*                ICH=1: J2000 -> TOD
*                ICH=2: TOD -> J2000
*-
*  2-May-1995 - Created from routines by Arthur Smith and Thomas
*               Schildknecht (AIUB)
*-----------------------------------------------------------------------
      real*8 epst,dpsi,epsm,za,ta,xa
*     real*8 tdb,dt
*     parameter (dt=32.186d0/86400d0)
*     integer i
*     tdb=mjd-dt

      call nutate (mjd,epst,dpsi,epsm)
      call precess(mjd,  za,  ta,  xa)

      if (ich.eq.2) then

* True-of-Date to J2000

         call rotate(1, epst,a(1),b(1))
         call rotate(3, dpsi,b(1),b(1))
         call rotate(1,-epsm,b(1),b(1))
         call rotate(3,   za,b(1),b(1))
         call rotate(2,  -ta,b(1),b(1))
         call rotate(3,   xa,b(1),b(1))
         if (ndim.eq.6) then
            call rotate(1, epst,a(4),b(4))
            call rotate(3, dpsi,b(4),b(4))
            call rotate(1,-epsm,b(4),b(4))
            call rotate(3,   za,b(4),b(4))
            call rotate(2,  -ta,b(4),b(4))
            call rotate(3,   xa,b(4),b(4))
         endif

      else

* J2000 to True-of-Date

         call rotate(3,  -xa,a(1),b(1))
         call rotate(2,   ta,b(1),b(1))
         call rotate(3,  -za,b(1),b(1))
         call rotate(1, epsm,b(1),b(1))
         call rotate(3,-dpsi,b(1),b(1))
         call rotate(1,-epst,b(1),b(1))
         if (ndim.eq.6) then
            call rotate(3,  -xa,a(4),b(4))
            call rotate(2,   ta,b(4),b(4))
            call rotate(3,  -za,b(4),b(4))
            call rotate(1, epsm,b(4),b(4))
            call rotate(3,-dpsi,b(4),b(4))
            call rotate(1,-epst,b(4),b(4))
         endif

      endif
      end


      SUBROUTINE PRECESS (MJD, ZA, TA, XA)
      implicit none
      REAL*8 MJD, ZA, TA, XA
*
* Computation of precession-angles for rotation from J2000
* to True-Of-Date. See USNO CIRCULAR NO 163, 1981.
*
* Arguments:
*  MJD (input): Modified Julian Date in UTC days
*  ZA (output): Z angle (radians)
*  TA (output): THETA angle (radians)
*  XA (output): ZETA angle (radians)
*-----------------------------------------------------------------------
      real*8 tl
C
C  TIME INTERVAL (IN JUL. CENTURIES) BETWEEN EPOCH AND J2000.0
      TL=(MJD-51544.5D0-32.186D0/86400D0)/36525.D0
C
C  ROTATION ANGLES XI,Z AND THETA
      XA=2306.2181  D0*TL
     1     +0.30188 D0*TL**2
     2     +0.017998D0*TL**3
      ZA=2306.2181  D0*TL
     1     +1.09468 D0*TL**2
     2     +0.018203D0*TL**3
      TA=2004.3109  D0*TL
     1     -0.42665 D0*TL**2
     2     -0.041833D0*TL**3
      XA=XA/206264.8D0
      ZA=ZA/206264.8D0
      TA=TA/206264.8D0
      END




      SUBROUTINE NUTATE (MJD, EPST, DPSI, EPSM)
      implicit none
      REAL*8 MJD, EPST, DPSI, EPSM
*
* Compute nutation angles at epoch MJD. Nutation is the difference
* between the Mean Equinox and Equator and the True Equinox and
* Equator at epoch.
* See USNO CIRCULAR NO 163, 1981.
*
* Arguments:
*  MJD (input): Modified Julian Date (in UTC)
*  EPSM (output): Mean obliquity of date (radians)
*  DPSI (output): Nutation in longitude (radians)
*  EPST (output): True obliquity of date (radians)
*-----------------------------------------------------------------------
      REAL*8 ARG1(9,22),ARG2(7,84),ALP(5)
      real*8 pi,r
      integer i,j
      real*8 roh,tu,tu4,eps,deps,phi
C
C DATA
C ----
      DATA PI/3.141592653589793D0/,R/1296000.D0/
C
C 22 COEFFICIENTS FOR THE STONGEST NUTATIONPERIODS (TIMEDEPENDENT)
C
C UNIT: ARCSEC/10000
C
      DATA ((ARG1(I,J),I=1,9),J=1,11)/
     1           0., 0., 0., 0., 1., -171996., -174.2, 92025., 8.9,
     2           0., 0., 0., 0., 2., 2062.   , .2    , -895. , .5 ,
     3          -2., 0., 2., 0., 1.,   46.   , .0    ,  -24. , .0 ,
     9           0., 0., 2.,-2., 2., -13187. ,-1.6   , 5736. ,-3.1,
     -           0., 1., 0., 0., 0.,   1426. ,-3.4   ,   54. ,-0.1,
     1           0., 1., 2.,-2., 2.,   -517. , 1.2   ,  224. ,-0.6,
     2           0.,-1., 2.,-2., 2.,    217. ,-0.5   ,  -95. , 0.3,
     3           0., 0., 2.,-2., 1.,    129. , 0.1   ,  -70. , 0.0,
     4           2., 0., 0.,-2., 0.,     48. , 0.0   ,    1. , 0.0,
     6           0., 2., 0., 0., 0.,     17. ,-0.1   ,    0. , 0.0,
     8           0., 2., 2.,-2., 2.,    -16. , 0.1   ,    7. , 0.0/
C
      DATA ((ARG1(I,J),I=1,9),J=12,22)/
     1           0., 0., 2., 0., 2.,  -2274. ,-0.2   ,  977. ,-0.5,
     2           1., 0., 0., 0., 0.,    712. , 0.1   ,   -7. , 0.0,
     3           0., 0., 2., 0., 1.,   -386. ,-0.4   ,  200. , 0.0,
     4           1., 0., 2., 0., 2.,   -301. , 0.0   ,  129. ,-0.1,
     5           1., 0., 0.,-2., 0.,   -158. , 0.0   ,   -1. , 0.0,
     6          -1., 0., 2., 0., 2.,    123. , 0.0   ,  -53. , 0.0,
     7           0., 0., 0., 2., 0.,     63. , 0.0   ,   -2. , 0.0,
     8           1., 0., 0., 0., 1.,     63. , 0.1   ,  -33. , 0.0,
     9          -1., 0., 0., 0., 1.,    -58. ,-0.1   ,   32. , 0.0,
     -          -1., 0., 2., 2., 2.,    -59. , 0.0   ,   26. , 0.0,
     1           1., 0., 2., 0., 1.,    -51. , 0.0   ,   27. , 0.0/
C
C 84 ADDITIONAL AMPLITUDES WITH NO TIMEDEPENDENT COEFFICIENTS
C
      DATA ((ARG2(I,J),I=1,7),J=1,19)/
     4           2., 0.,-2., 0., 0.,  11. ,   0. ,
     5          -2., 0., 2., 0., 2.,  -3. ,   1. ,
     6           1.,-1., 0.,-1., 0.,  -3. ,   0. ,
     7           0.,-2., 2.,-2., 1.,  -2. ,   1. ,
     8           2., 0.,-2., 0., 1.,   1. ,   0. ,
     5           0., 0., 2.,-2., 0., -22. ,   0. ,
     7           0., 1., 0., 0., 1., -15. ,   9. ,
     9           0.,-1., 0., 0., 1., -12. ,   6. ,
     -          -2., 0., 0., 2., 1.,  -6. ,   3. ,
     1           0.,-1., 2.,-2., 1.,  -5. ,   3. ,
     2           2., 0., 0.,-2., 1.,   4. ,  -2. ,
     3           0., 1., 2.,-2., 1.,   4. ,  -2. ,
     4           1., 0., 0.,-1., 0.,  -4. ,   0. ,
     5           2., 1., 0.,-2., 0.,   1. ,   0. ,
     6           0., 0.,-2., 2., 1.,   1. ,   0. ,
     7           0., 1.,-2., 2., 0.,  -1. ,   0. ,
     8           0., 1., 0., 0., 2.,   1. ,   0. ,
     9          -1., 0., 0., 1., 1.,   1. ,   0. ,
     -           0., 1., 2.,-2., 0.,  -1. ,   0. /
C
      DATA ((ARG2(I,J),I=1,7),J=20,37)/
     2           0., 0., 2., 2., 2., -38. ,  16. ,
     3           2., 0., 0., 0., 0.,  29. ,  -1. ,
     4           1., 0., 2.,-2., 2.,  29. , -12. ,
     5           2., 0., 2., 0., 2., -31. ,  13. ,
     6           0., 0., 2., 0., 0.,  26. ,  -1. ,
     7          -1., 0., 2., 0., 1.,  21. , -10. ,
     8          -1., 0., 0., 2., 1.,  16. ,  -8. ,
     9           1., 0., 0.,-2., 1., -13. ,   7. ,
     -          -1., 0., 2., 2., 1., -10. ,   5. ,
     1           1., 1., 0.,-2., 0.,  -7. ,   0. ,
     2           0., 1., 2., 0., 2.,   7. ,  -3. ,
     3           0.,-1., 2., 0., 2.,  -7. ,   3. ,
     4           1., 0., 2., 2., 2.,  -8. ,   3. ,
     5           1., 0., 0., 2., 0.,   6. ,   0. ,
     6           2., 0., 2.,-2., 2.,   6. ,  -3. ,
     7           0., 0., 0., 2., 1.,  -6. ,   3. ,
     8           0., 0., 2., 2., 1.,  -7. ,   3. ,
     9           1., 0., 2.,-2., 1.,   6. ,  -3. /
C
      DATA ((ARG2(I,J),I=1,7),J=38,47)/
     -           0., 0., 0.,-2., 1.,  -5. ,   3. ,
     1           1.,-1., 0., 0., 0.,   5. ,   0. ,
     2           2., 0., 2., 0., 1.,  -5. ,   3. ,
     3           0., 1., 0.,-2., 0.,  -4. ,   0. ,
     4           1., 0.,-2., 0., 0.,   4. ,   0. ,
     5           0., 0., 0., 1., 0.,  -4. ,   0. ,
     6           1., 1., 0., 0., 0.,  -3. ,   0. ,
     7           1., 0., 2., 0., 0.,   3. ,   0. ,
     8           1.,-1., 2., 0., 2.,  -3. ,   1. ,
     9          -1.,-1., 2., 2., 2.,  -3. ,   1. /
C
      DATA ((ARG2(I,J),I=1,7),J=48,57)/
     -          -2., 0., 0., 0., 1.,  -2. ,   1. ,
     1           3., 0., 2., 0., 2.,  -3. ,   1. ,
     2           0.,-1., 2., 2., 2.,  -3. ,   1. ,
     3           1., 1., 2., 0., 2.,   2. ,  -1. ,
     4          -1., 0., 2.,-2., 1.,  -2. ,   1. ,
     5           2., 0., 0., 0., 1.,   2. ,  -1. ,
     6           1., 0., 0., 0., 2.,  -2. ,   1. ,
     7           3., 0., 0., 0., 0.,   2. ,   0. ,
     8           0., 0., 2., 1., 2.,   2. ,  -1. ,
     9          -1., 0., 0., 0., 2.,   1. ,  -1. /
C
      DATA ((ARG2(I,J),I=1,7),J=58,67)/
     -           1., 0., 0.,-4., 0.,  -1. ,   0. ,
     1          -2., 0., 2., 2., 2.,   1. ,  -1. ,
     2          -1., 0., 2., 4., 2.,  -2. ,   1. ,
     3           2., 0., 0.,-4., 0.,  -1. ,   0. ,
     4           1., 1., 2.,-2., 2.,   1. ,  -1. ,
     5           1., 0., 2., 2., 1.,  -1. ,   1. ,
     6          -2., 0., 2., 4., 2.,  -1. ,   1. ,
     7          -1., 0., 4., 0., 2.,   1. ,   0. ,
     8           1.,-1., 0.,-2., 0.,   1. ,   0. ,
     9           2., 0., 2.,-2., 1.,   1. ,  -1. /
C
      DATA ((ARG2(I,J),I=1,7),J=68,77)/
     -           2., 0., 2., 2., 2.,  -1. ,   0. ,
     1           1., 0., 0., 2., 1.,  -1. ,   0. ,
     2           0., 0., 4.,-2., 2.,   1. ,   0. ,
     3           3., 0., 2.,-2., 2.,   1. ,   0. ,
     4           1., 0., 2.,-2., 0.,  -1. ,   0. ,
     5           0., 1., 2., 0., 1.,   1. ,   0. ,
     6          -1.,-1., 0., 2., 1.,   1. ,   0. ,
     7           0., 0.,-2., 0., 1.,  -1. ,   0. ,
     8           0., 0., 2.,-1., 2.,  -1. ,   0. ,
     9           0., 1., 0., 2., 0.,  -1. ,   0. /
C
      DATA ((ARG2(I,J),I=1,7),J=78,84)/
     -           1., 0.,-2.,-2., 0.,  -1. ,   0. ,
     1           0.,-1., 2., 0., 1.,  -1. ,   0. ,
     2           1., 1., 0.,-2., 1.,  -1. ,   0. ,
     3           1., 0.,-2., 2., 0.,  -1. ,   0. ,
     4           2., 0., 0., 2., 0.,   1. ,   0. ,
     5           0., 0., 2., 4., 2.,  -1. ,   0. ,
     6           0., 1., 0., 1., 0.,   1. ,   0. /
C
      ROH = PI/648000.D0
C
C  TIME INTERVAL (IN JUL. CENTURIES) BETWEEN XTDB AND J2000.0
      TU =(MJD-51544.5D0-32.186D0/86400D0)/36525.D0
C TU4 = MOD JUL. DATUM OF REFERENCE EPOCH ( XMOD = XTDB FOR EPSM)
      TU4=TU
C
C  FUNDAMENTAL ARGUMENTS (IN RAD)
      ALP(1)=(485866.733D0+(1325.D0*R+715922.633D0)*TU+31.31D0*TU*TU+
     1       .064D0*TU*TU*TU)*ROH
      ALP(2)=(1287099.804D0+(99.D0*R+1292581.224D0)*TU-0.577D0*TU*TU-
     1       .012D0*TU*TU*TU)*ROH
      ALP(3)=(335778.877D0+(1342.D0*R+295263.137D0)*TU-13.257D0*TU*TU+
     1       .011D0*TU*TU*TU)*ROH
      ALP(4)=(1072261.307D0+(1236.D0*R+1105601.328D0)*TU-6.891D0*TU*TU+
     1       .019D0*TU*TU*TU)*ROH
      ALP(5)=(450160.28D0-(5.D0*R+482890.539D0)*TU+7.455D0*TU*TU+
     1       .008D0*TU*TU*TU)*ROH
C
C EPSM-ZERO EPOCH J2000 (RAD)
      EPSM=(84381.448D0-46.815D0*TU4-0.00059D0*TU4*TU4+
     1       .001813D0*TU4*TU4*TU4)*ROH
C
      DPSI=0.D0
      DEPS=0.D0
      DO 10 J=1,22
        PHI = 0.D0
        DO 5 I=1,5
          PHI = PHI + ARG1(I,J)*ALP(I)
5       CONTINUE
        PHI = DMOD (PHI,2.D0*PI)
C       NUTATION IN ARCSEC*10000
        DPSI = DPSI+(ARG1(6,J)+ARG1(7,J)*TU)*DSIN(PHI)
        DEPS = DEPS+(ARG1(8,J)+ARG1(9,J)*TU)*DCOS(PHI)
10    CONTINUE
      DO 20 J=1,84
        PHI = 0.D0
        DO 15 I=1,5
          PHI = PHI + ARG2(I,J)*ALP(I)
15      CONTINUE
        PHI = DMOD (PHI,2.D0*PI)
C SERIES FOR NUTATION IN ARCSEC*10000 IN LONGITUDE AND OBLIQUITY
        DPSI = DPSI+ARG2(6,J)*DSIN(PHI)
        DEPS = DEPS+ARG2(7,J)*DCOS(PHI)
20    CONTINUE
C
C DPSI AND DEPS IN RADIANS, TRUE OBLIQUITY EPST
      DEPS=DEPS*ROH/10000.D0
      DPSI=DPSI*ROH/10000.D0
      EPST=EPSM+DEPS
      END
