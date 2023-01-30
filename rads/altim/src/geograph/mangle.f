**MANGLE -- Compute sine and cosine of m*angle
*+
      SUBROUTINE MANGLE (L, ANGLE, COSMA, SINMA)
      INTEGER L
      REAL*8  ANGLE, COSMA(0:L), SINMA(0:L)

* This routine computes the sine and cosine of all multiples of
* ANGLE. The arrays COSMA and SINMA are filled with COS(M*ANGLE) and
* SIN(M*ANGLE), respectively, for M=0...L and L >= 0
*
* Arguments:
*  L      (input) : maximum multiple of ANGLE
*  ANGLE  (input) : angle in radians
*  COSMA (output) : array with COS(M*ANGLE) for M=0...L
*  SINMA (output) : array with SIN(M*ANGLE) for M=0...L
*-
*  6-Aug-2001 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer m

      cosma(0)=1d0
      sinma(0)=0d0
      if (l.eq.0) return
      cosma(1)=dcos(angle)
      sinma(1)=dsin(angle)
      do m=2,l
         cosma(m)=cosma(m-1)*cosma(1)-sinma(m-1)*sinma(1)
         sinma(m)=cosma(m-1)*sinma(1)+sinma(m-1)*cosma(1)
      enddo
      end
