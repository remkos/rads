**J2000 -- Transform ephemerides between J2000 and TOD
*+
      SUBROUTINE J2000 (ICH, MJD, DPSI, DEPS, NDIM, A, B)
      REAL*8 MJD, DPSI, DEPS, A(*), B(*)
      INTEGER*4 NDIM, ICH
*
* This routine transforms ephemerides from J2000 to True-Of-Date
* (TOD) reference frame or vice versa.
*
* A and B are state vectors (position and/or velocity) in the inertial
* reference frame.
* A and B may share the same array in the calling routine.
*
* NDIM indicates whether the vectors contain 3 or 6 coordinates.
* ICH determines the transformation direction (to or from J2000).
* The epoch of the TOD system is given by MJD (Terrestial Time).
*
* Note 1: J2000 is the system of reference date 1.5 Jan 2000.
* Note 2: Following recommendations by IERS [McCarthy, 1996] the
*         expressions are to be evaluated as a function of Terrestial
*         Time and not Barycentric Dynamical Time.
* Note 3: Routine is based on IAU 1980 nutation theory [Seidelmann, 1982;
*         Wahr, 1981]. The coefficient values from IERS 1996
*         [McCarthy, 1996] are used to compute the nutation angles in
*         longitude and obliquity.
*         These angles can be corrected for the observed celestial pole
*         offsets (DPSI and DEPS) reported in the IERS Bulletins.
*
* Arguments:
*  ICH  (input): Direction of the transformation:
*                 1 = J2000 -> TOD
*                -1 = TOD -> J2000
*  MJD  (input): Reference epoch of the TOD system (Mean Julian Date)
*  DPSI (input): Correction on nutation in longitude (milliarcseconds)
*  DEPS (input): Correction on nutation in obliquity (milliarcseconds)
*  NDIM (input): Dimension of the state vector:
*                 3 = position vector only
*                 6 = full state vector (position and velocity)
*  A    (input): Input state vector
*  B   (output): Transformed state vector
*
* References:
*
* D. D. McCarthy, "IERS Conventions (1996)", IERS Technical Note 21,
* Central Bureau of IERS, Observatoire de Paris, 1996.
*
* O. Montenbruck and E. Gill, "Satellite Orbits: Models, Methods,
* Applications", Springer Verlag, 2000.
*
* P. K. Seidelmann, "1980 IAU Nutation: The final report of the IAU
* Working Group on Nutation", Celest. Mech., 27, 79-106, 1982.
*
* J. M. Wahr, "The forced nutations of an elliptical, rotating, elastic,
* and oceanless Earth", Geophys. J. Roy. Astron. Soc., 64, 705-727, 1981.
*-
*  2-May-1995 - Created from routines by Arthur Smith (DEOS)
*               and Thomas Schildknecht (AIUB)
* 12-Jan-2001 - Add D0 to all real values (also small ones)
*-----------------------------------------------------------------------
      REAL*8 N(3,3),P(3,3),C(3)

* Compute precession and nutation matrices to transform ephemerides
* from J2000 to TOD

      CALL J2000P(MJD,P)
      CALL J2000N(MJD,DPSI,DEPS,N)

      IF (ICH.EQ.1) THEN

* J2000 -> TOD :    B = N * P * A

	 call matmmv(3,3,3,P,A(1),C)
	 call matmmv(3,3,3,N,C,B(1))
         if (ndim.eq.6) then
	    call matmmv(3,3,3,P,A(4),C)
	    call matmmv(3,3,3,N,C,B(4))
	 endif
      ELSE

* TOD -> J2000 :    B = P' * N' * A

	 call mattmv(3,3,3,N,A(1),C)
	 call mattmv(3,3,3,P,C,B(1))
         if (ndim.eq.6) then
	    call mattmv(3,3,3,N,A(4),C)
	    call mattmv(3,3,3,P,C,B(4))
	 endif
      ENDIF
      RETURN
      END

**J2000P -- Compute precession matrix
*+
      SUBROUTINE J2000P (XMJD,PREC)
      REAL*8 XMJD, PREC(3,3)
*
* Computation of precession matrix from J2000 to epoch XMJD
* (in modified julian date, terrestial time).
*
* Arguments:
*  XMJD  (input): Epoch in MJD, Terrestial Time
*  PREC (output): 3x3 Precession matrix
*
*-
*  2-May-1985 - Adapted from BERNESE software (T. Schildknecht, AIUB)
* 28-Feb-2001 - Changed reference to Barycentric time to Terrestial
*               Time according to IERS Conventions (1996)
*-----------------------------------------------------------------------

      REAL*8 XA,ZA,TA,T,T2,T3,PI,ARCSEC
      REAL*8 SINXA,SINZA,SINTA
      REAL*8 COSXA,COSZA,COSTA
      PARAMETER (PI=3.141592653589793D0,ARCSEC=PI/180/3600)

* Time interval (in Julian centuries) between epoch (XMJD) and
* J2000 (1.5 Jan 2000 = MJD 51544.5)

      T=(XMJD-51544.5D0)/36525.D0
      T2=T*T
      T3=T*T2

* Rotation angles zeta, z and theta (in arcseconds)

      XA=2306.2181  D0*T
     |     +0.30188 D0*T2
     |     +0.017998D0*T3
      ZA=2306.2181  D0*T
     |     +1.09468 D0*T2
     |     +0.018203D0*T3
      TA=2004.3109  D0*T
     |     -0.42665 D0*T2
     |     -0.041833D0*T3

* Conversion from arcseconds to radians

      XA=XA*ARCSEC
      ZA=ZA*ARCSEC
      TA=TA*ARCSEC

* Compute SIN and COS of angles

      COSXA=DCOS(XA)
      SINXA=DSIN(XA)
      COSZA=DCOS(ZA)
      SINZA=DSIN(ZA)
      COSTA=DCOS(TA)
      SINTA=DSIN(TA)

* Construct rotation matrix

      PREC(1,1)= COSXA*COSZA*COSTA-SINXA*SINZA
      PREC(2,1)= COSXA*SINZA*COSTA+SINXA*COSZA
      PREC(3,1)= COSXA*SINTA
      PREC(1,2)=-SINXA*COSZA*COSTA-COSXA*SINZA
      PREC(2,2)=-SINXA*SINZA*COSTA+COSXA*COSZA
      PREC(3,2)=-SINXA*SINTA
      PREC(1,3)=-COSZA*SINTA
      PREC(2,3)=-SINZA*SINTA
      PREC(3,3)= COSTA

      RETURN
      END

**J2000N -- Compute nutation matrix
*+
      SUBROUTINE J2000N (XMJD,dDPSI,dDEPS,NUTN)
      REAL*8 XMJD, dDPSI, dDEPS, NUTN(3,3)
*
* Computation of nutation matrix from J2000 to epoch XMJD
* (in modified julian date, terrestial time).
* Formulas according to IERS 1996 conventions (adaptation of IAU 1980).
* The nutation according to the IAU 1980 theory can be corrected with
* observed celestial pole offsets (dDPSI and dDEPS).
*
* Arguments:
*  XMJD   (input): Epoch in MJD, Terrestial Time
*  dDPSI  (input): Correction on nutation in longitude (milliarcseconds)
*  dDEPS  (input): Correction on nutation in obliquity (milliarcseconds)
*  NUTN  (output): 3x3 Nutation matrix
*
*-
*  2-May-1985 - Adapted from BERNESE software (E. Brockmann, AIUB)
*  2-Mar-2001 - Changed reference to Barycentric time to Terrestial
*               Time according and to IERS Conventions (1996).
*               Additional corrections dDPSI and dDEPS included.
*-----------------------------------------------------------------------
      REAL*8 ARGA(0:4,6),ARGT(2,15),FARG(6)
      REAL*8 PI,TWOPI
      INTEGER*4 I,J,ARGI(7,106)
      REAL*8 ARCSEC,T(4),DPSI,DEPS,PHI,EPSTRU
      REAL*8 SINDP,COSDP,SINEM,COSEM,SINET,COSET
      PARAMETER (PI=3.141592653589793D0,ARCSEC=PI/180/3600,
     |TWOPI=2*PI)

* Coefficients of the 5 fundamental arguments and for obliquity.
* Coefficients are in arcseconds: the first is a constant, the remaining 4
* have to be multiplied by T, T**2, T**3 and T**4 (T in Julian Centuries).

      DATA ARGA /
     | 485868.249036D0,1717915923.2178D0, 31.8792D0,51635D-6,-24470D-8,
     |1287104.793048D0, 129596581.0481D0, -0.5532D0,  136D-6, -1149D-8,
     | 335779.526232D0,1739527262.8478D0,-12.7512D0,-1037D-6,   417D-8,
     |1072260.703692D0,1602961601.2090D0, -6.3706D0, 6593D-6, -3169D-8,
     | 450160.398036D0,  -6962890.2665D0,  7.4722D0, 7702D-6, -5939D-8,
     |  84381.448   D0,       -46.8150D0,-0.00059D0, 1813D-6,     0D-8/

* 15 Coefficients for the time-dependent nutation angles
* in 0.1 milliarcseconds (1/36 000 000 degrees)

      DATA ARGT/
     |      -174.2D0, 8.9D0,
     |         0.2D0, 0.5D0,
     |        -1.6D0,-3.1D0,
     |        -3.4D0,-0.1D0,
     |         1.2D0,-0.6D0,
     |        -0.5D0, 0.3D0,
     |         0.1D0, 0.0D0,
     |        -0.1D0, 0.0D0,
     |         0.1D0, 0.0D0,
     |        -0.2D0,-0.5D0,
     |         0.1D0, 0.0D0,
     |        -0.4D0, 0.0D0,
     |         0.0D0,-0.1D0,
     |         0.1D0, 0.0D0,
     |        -0.1D0, 0.0D0/
      DATA ((ARGI(I,J),I=1,7),J=1,15)/
     |      0, 0, 0, 0, 1,-171996,92025,
     |      0, 0, 0, 0, 2,   2062, -895,
     |      0, 0, 2,-2, 2, -13187, 5736,
     |      0, 1, 0, 0, 0,   1426,   54,
     |      0, 1, 2,-2, 2,   -517,  224,
     |      0,-1, 2,-2, 2,    217,  -95,
     |      0, 0, 2,-2, 1,    129,  -70,
     |      0, 2, 0, 0, 0,     17,    0,
     |      0, 2, 2,-2, 2,    -16,    7,
     |      0, 0, 2, 0, 2,  -2274,  977,
     |      1, 0, 0, 0, 0,    712,   -7,
     |      0, 0, 2, 0, 1,   -386,  200,
     |      1, 0, 2, 0, 2,   -301,  129,
     |      1, 0, 0, 0, 1,     63,  -33,
     |     -1, 0, 0, 0, 1,    -58,   32/

* 91 additional amplitudes without timedependent coefficients

      DATA ((ARGI(I,J),I=1,7),J=16,30)/
     |     -2, 0, 2, 0, 1,  46, -24,
     |      2, 0,-2, 0, 0,  11,   0,
     |     -2, 0, 2, 0, 2,  -3,   1,
     |      1,-1, 0,-1, 0,  -3,   0,
     |      0,-2, 2,-2, 1,  -2,   1,
     |      2, 0,-2, 0, 1,   1,   0,
     |      2, 0, 0,-2, 0,  48,   1,
     |      0, 0, 2,-2, 0, -22,   0,
     |      0, 1, 0, 0, 1, -15,   9,
     |      0,-1, 0, 0, 1, -12,   6,
     |     -2, 0, 0, 2, 1,  -6,   3,
     |      0,-1, 2,-2, 1,  -5,   3,
     |      2, 0, 0,-2, 1,   4,  -2,
     |      0, 1, 2,-2, 1,   4,  -2,
     |      1, 0, 0,-1, 0,  -4,   0/
      DATA ((ARGI(I,J),I=1,7),J=31,45)/
     |      2, 1, 0,-2, 0,   1,   0,
     |      0, 0,-2, 2, 1,   1,   0,
     |      0, 1,-2, 2, 0,  -1,   0,
     |      0, 1, 0, 0, 2,   1,   0,
     |     -1, 0, 0, 1, 1,   1,   0,
     |      0, 1, 2,-2, 0,  -1,   0,
     |      1, 0, 0,-2, 0,-158,  -1,
     |     -1, 0, 2, 0, 2, 123, -53,
     |      0, 0, 0, 2, 0,  63,  -2,
     |     -1, 0, 2, 2, 2, -59,  26,
     |      1, 0, 2, 0, 1, -51,  27,
     |      0, 0, 2, 2, 2, -38,  16,
     |      2, 0, 0, 0, 0,  29,  -1,
     |      1, 0, 2,-2, 2,  29, -12,
     |      2, 0, 2, 0, 2, -31,  13/
      DATA ((ARGI(I,J),I=1,7),J=46,60)/
     |      0, 0, 2, 0, 0,  26,  -1,
     |     -1, 0, 2, 0, 1,  21, -10,
     |     -1, 0, 0, 2, 1,  16,  -8,
     |      1, 0, 0,-2, 1, -13,   7,
     |     -1, 0, 2, 2, 1, -10,   5,
     |      1, 1, 0,-2, 0,  -7,   0,
     |      0, 1, 2, 0, 2,   7,  -3,
     |      0,-1, 2, 0, 2,  -7,   3,
     |      1, 0, 2, 2, 2,  -8,   3,
     |      1, 0, 0, 2, 0,   6,   0,
     |      2, 0, 2,-2, 2,   6,  -3,
     |      0, 0, 0, 2, 1,  -6,   3,
     |      0, 0, 2, 2, 1,  -7,   3,
     |      1, 0, 2,-2, 1,   6,  -3,
     |      0, 0, 0,-2, 1,  -5,   3/
      DATA ((ARGI(I,J),I=1,7),J=61,75)/
     |      1,-1, 0, 0, 0,   5,   0,
     |      2, 0, 2, 0, 1,  -5,   3,
     |      0, 1, 0,-2, 0,  -4,   0,
     |      1, 0,-2, 0, 0,   4,   0,
     |      0, 0, 0, 1, 0,  -4,   0,
     |      1, 1, 0, 0, 0,  -3,   0,
     |      1, 0, 2, 0, 0,   3,   0,
     |      1,-1, 2, 0, 2,  -3,   1,
     |     -1,-1, 2, 2, 2,  -3,   1,
     |     -2, 0, 0, 0, 1,  -2,   1,
     |      3, 0, 2, 0, 2,  -3,   1,
     |      0,-1, 2, 2, 2,  -3,   1,
     |      1, 1, 2, 0, 2,   2,  -1,
     |     -1, 0, 2,-2, 1,  -2,   1,
     |      2, 0, 0, 0, 1,   2,  -1/
      DATA ((ARGI(I,J),I=1,7),J=76,90)/
     |      1, 0, 0, 0, 2,  -2,   1,
     |      3, 0, 0, 0, 0,   2,   0,
     |      0, 0, 2, 1, 2,   2,  -1,
     |     -1, 0, 0, 0, 2,   1,  -1,
     |      1, 0, 0,-4, 0,  -1,   0,
     |     -2, 0, 2, 2, 2,   1,  -1,
     |     -1, 0, 2, 4, 2,  -2,   1,
     |      2, 0, 0,-4, 0,  -1,   0,
     |      1, 1, 2,-2, 2,   1,  -1,
     |      1, 0, 2, 2, 1,  -1,   1,
     |     -2, 0, 2, 4, 2,  -1,   1,
     |     -1, 0, 4, 0, 2,   1,   0,
     |      1,-1, 0,-2, 0,   1,   0,
     |      2, 0, 2,-2, 1,   1,  -1,
     |      2, 0, 2, 2, 2,  -1,   0/
      DATA ((ARGI(I,J),I=1,7),J=91,106)/
     |      1, 0, 0, 2, 1,  -1,   0,
     |      0, 0, 4,-2, 2,   1,   0,
     |      3, 0, 2,-2, 2,   1,   0,
     |      1, 0, 2,-2, 0,  -1,   0,
     |      0, 1, 2, 0, 1,   1,   0,
     |     -1,-1, 0, 2, 1,   1,   0,
     |      0, 0,-2, 0, 1,  -1,   0,
     |      0, 0, 2,-1, 2,  -1,   0,
     |      0, 1, 0, 2, 0,  -1,   0,
     |      1, 0,-2,-2, 0,  -1,   0,
     |      0,-1, 2, 0, 1,  -1,   0,
     |      1, 1, 0,-2, 1,  -1,   0,
     |      1, 0,-2, 2, 0,  -1,   0,
     |      2, 0, 0, 2, 0,   1,   0,
     |      0, 0, 2, 4, 2,  -1,   0,
     |      0, 1, 0, 1, 0,   1,   0/

* Time interval (in Julian centuries) between epoch (XMJD) and
* J2000 (1.5 Jan 2000 = MJD 51544.5)

      T(1)=(XMJD-51544.5D0)/36525.D0	! t**1
      T(2)=T(1)*T(1)			! t**2
      T(3)=T(1)*T(2)			! t**3
      T(4)=T(2)*T(2)			! t**4

* 5 Luni-solar fundamental arguments (in arcseconds): l, l', F, D, Omega
* 6th parameter is obliquity: epsilon
* Last step converts to radians

      DO J=1,6
	FARG(J)=ARGA(0,J)
        DO I=1,4
          FARG(J)=FARG(J)+ARGA(I,J)*T(I)
        ENDDO
        FARG(J)=FARG(J)*ARCSEC
      ENDDO

* Compute DPsi and Deps (in 0.0001 arcsec) for the 15 time dependent terms

      DPSI=0D0
      DEPS=0D0
      DO J=1,15
        PHI = 0D0
        DO I=1,5
          PHI = PHI + ARGI(I,J)*FARG(I)
        ENDDO
        PHI = DMOD (PHI,TWOPI)
        DPSI = DPSI+(ARGI(6,J)+ARGT(1,J)*T(1))*DSIN(PHI)
        DEPS = DEPS+(ARGI(7,J)+ARGT(2,J)*T(1))*DCOS(PHI)
      ENDDO

* Add the remaining terms without time dependency

      DO J=16,106
        PHI = 0D0
        DO I=1,5
          PHI = PHI + ARGI(I,J)*FARG(I)
        ENDDO
        PHI = DMOD (PHI,TWOPI)
        DPSI = DPSI+ARGI(6,J)*DSIN(PHI)
        DEPS = DEPS+ARGI(7,J)*DCOS(PHI)
      ENDDO

* Add celestial pole offsets (in arcseconds) and convert DPsi and Deps to radians

      DPSI=(DPSI/1D4+dDPSI)*ARCSEC
      DEPS=(DEPS/1D4+dDEPS)*ARCSEC

* Determine true obliquity of the ecliptic (eps')

      EPSTRU=FARG(6)+DEPS

* Store SIN and COS of DPsi, eps and eps'

      SINDP=DSIN(DPSI)
      COSDP=DCOS(DPSI)
      SINEM=DSIN(FARG(6))
      COSEM=DCOS(FARG(6))
      SINET=DSIN(EPSTRU)
      COSET=DCOS(EPSTRU)

* Construct rotation matrix

      NUTN(1,1)= COSDP
      NUTN(2,1)= COSET*SINDP
      NUTN(3,1)= SINET*SINDP
      NUTN(1,2)=-COSEM*SINDP
      NUTN(2,2)= COSEM*COSET*COSDP + SINEM*SINET
      NUTN(3,2)= COSEM*SINET*COSDP - SINEM*COSET
      NUTN(1,3)=-SINEM*SINDP
      NUTN(2,3)= SINEM*COSET*COSDP - COSEM*SINET
      NUTN(3,3)= SINEM*SINET*COSDP + COSEM*COSET

      RETURN
      END

**MATMMV - Multiply Matrix with Vector
*+
      SUBROUTINE MATMMV(M,N,L,A,X,B)
*-
      INTEGER M,N,L,I,J
      REAL*8 A(L,N),X(N),B(M),T
      DO I=1,M
         T=0D0
         DO J=1,N
            T=T+A(I,J)*X(J)
         ENDDO
         B(I)=T
      ENDDO
      END

**MATTMV - Multiply Transposed Matrix with Vector
*+
      SUBROUTINE MATTMV(M,N,L,A,X,B)
*-
      INTEGER M,N,L,I,J
      REAL*8 A(L,M),X(N),B(M),T
      DO I=1,M
         T=0D0
         DO J=1,N
            T=T+A(J,I)*X(J)
         ENDDO
         B(I)=T
      ENDDO
      END
