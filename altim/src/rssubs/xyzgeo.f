**XYZGEO -- Convert ECF coordinates to geodetic longitude and latitude.
*+
      SUBROUTINE XYZGEO (XYZ, R, LAT, LON, HEIGHT)
      REAL*8 XYZ(3), R, LAT, LON, HEIGHT
*
* This subroutine converts the Earth Centered Fixed coordinates X, Y, and Z
* to the geodetic latitude (LAT) and longitude (LON) and the height (HEIGHT)
* above the reference ellipsoid GRS80.
* An iterative method is used, where a parameter DZ (to be added to the Z
* coordinate) has to converge up to 1.D-4 m. This results in an accuracy in
* LAT of better than 1.D-11 radian and HEIGHT of 1.D-4 metres.
* Generally this routine converges in 4-5 iterations, after which LAT
* and HEIGHT are computed and control is given back to the calling
* (sub)routine.
* In any event, the iteration is stopped after 10 iterations, after which
* a warning message is printed and results are returned.
*
* Reference:
* O. Montenbruck and E. Gill, "Satellite Orbits: Models, Methods,
* Applications", Springer Verlag, p. 188, 2000.
*
* Arguments:
* XYZ(3)  (input) : Earth Centered Fixed coordinates X, Y, and Z (m).
* R      (output) : Distance to geocenter (m).
* LAT    (output) : Geodetic latitude (rad).
* LON    (output) : Geodetic longitude (rad).
* HEIGHT (output) : Height above the reference ellipsoid (m).
*-
* 15-Mar-1991: Created - Remko Scharroo
* 18-Nov-1991: New Version, new input
* 17-Feb-1994: Set bounds to number of iterations
*  2-Aug-1998: Loop changed to avoid error with Silicon compiler optimasation
* 18-May-1999: Algorithm improved: much faster iteration and more accurate
*  6-Mar-2001: Yet another new algorithm. CPU time reduced by factor 3-4.
* 10-May-2001: Output of R was broken. Repaired.
* 16-Jan-2002: Version with adjustable ellipsoid coefficients
*-----------------------------------------------------------------------
      real*8  x2y2,dz,dz0,zdz,sinphi,rn,flat,ffact,ae,finv
      integer*4 i

      call getearth(ae,finv)
      flat=1d0/finv
      ffact=flat*(flat-2)

* Values for later usage:
*   x2y2 = x**2+y**2
*   e=sqrt(1-(1-f)**2) => -e**2 = ffact

      x2y2=xyz(1)*xyz(1)+xyz(2)*xyz(2)
      r=sqrt(x2y2+xyz(3)*xyz(3))

* Initialize Delta z = e**2 z

      dz = -ffact*xyz(3)

* Start iteration

      do i=1,10
         dz0=dz
	 zdz=xyz(3)+dz
         sinphi=zdz/sqrt(x2y2+zdz**2)
	 rn=ae/sqrt(1+ffact*sinphi**2)
	 dz=-rn*ffact*sinphi
	 if (abs(dz-dz0).lt.1d-4) goto 10
      enddo

* Print warning when convergence is not yet achieved

      write (*,'("WARNING: xyzgeo did not fully converge in 10 steps")')

* Compute relevant output

10    lon=atan2(xyz(2),xyz(1))
      lat=atan2(zdz,sqrt(x2y2))
      height=sqrt(x2y2+zdz**2)-rn
      return
      end
