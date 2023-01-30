**ANMTOT -- Convert Mean anomaly to True anomaly.
*+
      FUNCTION ANMTOT (M, E)
      REAL*8 M, E, ANMTOT
*
* Computes the true anomaly out of the mean anomaly (and eccentricity)
* with an accuracy of (eccentricity**7)
*
* Arguments:
*  M       (input): Mean anomaly (radians)
*  E       (input): Eccentricity (unity)
*  ANMTOT (output): True anomaly (radians)
*-
* 12-Nov-1991: Redesigned (Remko Scharroo)
*-----------------------------------------------------------------------
      real*8 e2,e4,e6
      e2=e*e
      e4=e2*e2
      e6=e2*e4
      anmtot=m+sin(m)*e*(2-e2/4+5*e4/96+107*e6/4608)
     |        +sin(2*m)*(5*e2/4-11*e4/24+17*e6/192)
     |        +sin(3*m)*e*(13*e2/12-43*e4/64+95*e6/512)
     |        +sin(4*m)*(103*e4/96-451*e6/480)
     |        +sin(5*m)*e*(1097*e4/960-5957*e6/4608)
     |        +sin(6*m)*1223*e6/960
     |        +sin(7*m)*e*47273*e6/32256
      end
