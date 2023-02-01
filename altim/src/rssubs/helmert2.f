**HELMERT2 -- Perform Helmert transformation of a number of points
*+
      SUBROUTINE HELMERT2 (COEFF, N, XYZFROM, XYZTO)
      INTEGER N
      REAL*8 COEFF(7), XYZFROM(3,N), XYZTO(3,N)
*
* This routine performs the Helmert transformation of N points of which
* the Cartesian coordinates are given in the array XYZFROM. The result
* ends up in the array XYZTO. XYZFROM and XYZTO may share the same
* memory space.
*
* The seven Helmert transformation coefficients are stored in the array
* COEFF. They are:
*
* COEFF(1..3) : Translation vector (X,Y,Z) in meters
* COEFF(4)    : Scale difference
* COEFF(5..7) : Rotations around the X, Y and Z axes in radians
*
* Note that the position coordinates in XYZFROM and XYZTO have to be in
* the same units as COEFF(1..3).
*
* Arguments:
*  COEFF   (input) : Helmert transformation coefficients (see above)
*  N       (input) : Number of points to transform
*  XYZFROM (input) : XYZ-Coordinates of the points to transform
*  XYZTO  (output) : Transformed XYZ-Coordinates
* 
*-
*  4-Apr-2001 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer i
      real*8 x,y,z
      do i=1,n
         x=xyzfrom(1,i)+coeff(1)
     |+coeff(4)*xyzfrom(1,i)-coeff(7)*xyzfrom(2,i)+coeff(6)*xyzfrom(3,i)
         y=xyzfrom(2,i)+coeff(2)
     |+coeff(7)*xyzfrom(1,i)+coeff(4)*xyzfrom(2,i)-coeff(5)*xyzfrom(3,i)
         z=xyzfrom(3,i)+coeff(3)
     |-coeff(6)*xyzfrom(1,i)+coeff(5)*xyzfrom(2,i)+coeff(4)*xyzfrom(3,i)
         xyzto(1,i)=x
	 xyzto(2,i)=y
	 xyzto(3,i)=z
      enddo
      end
