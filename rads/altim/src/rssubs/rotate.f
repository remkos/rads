**ROTATE -- Rotate coordinate system
*+
      SUBROUTINE ROTATE (AXIS, ANGLE, A, B)
      INTEGER*4 AXIS
      REAL*8 ANGLE, A(3), B(3)
*
* This routine computes the object of a coordinate vector while
* rotating the coordinate system over an angle 'ANGLE' around the
* coordinate axis 'AXIS'. Vector A in the original system becomes
* vector B in the new system. In the call to this subroutine arrays
* A and B may be the same array.
*
* Arguments:
*  AXIS  (input): Number of the coordinate axis (1,2,3)
*  ANGLE (input): Rotation angle (right-handed) in radians
*  A     (input): Vector in the original system
*  B    (output): The same vector, but now in the rotated coordinate
*                 system
*-
*  1-Dec-1989: Created (Remko Scharroo)
* 12-Nov-1991: Revised
*-----------------------------------------------------------------------
      integer*4 i1,i2,i3
      real*8 bitwo

      i1=mod(axis+2,3)+1
      i2=mod(axis+3,3)+1
      i3=mod(axis+4,3)+1
      b(i1)=a(i1)
      bitwo=a(i2)*cos(angle)+a(i3)*sin(angle)
      b(i3)=a(i3)*cos(angle)-a(i2)*sin(angle)
      b(i2)=bitwo
      end
