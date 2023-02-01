**ANTTOM -- Convert True anomaly to Mean anomaly.
*+
      FUNCTION ANTTOM (THETA, E)
      REAL*8 THETA, E, ANTTOM
*
* Computes the mean anomaly out of the true anomaly (and eccentricity)
* using straightforward algorithms. No approximation. No iteration.
*
* Arguments:
*  THETA   (input): True anomaly (radians)
*  E       (input): Eccentricity (unity)
*  ANTTOM (output): Mean anomaly (radians)
*-
*  1-Aug-2004: Created (Remko Scharroo)
*-----------------------------------------------------------------------
      real*8 ecc	! eccentric anomaly
      ecc=atan2(sqrt(1d0-e*e)*sin(theta),(cos(theta)+e))
      anttom=ecc-e*sin(ecc)
      end
