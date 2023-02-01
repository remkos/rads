**GRXHEX -- Convert RGB values to index code
*+
      SUBROUTINE GRXHEX (R, G, B, HEX)
      REAL R, G, B
*     INTEGER HEX
*
* This routine converts the red, green, and blue values (in the range
* 0.0-1.0) to the appropriate negative `true color' index code.
* This code consists of four bytes:
*   Byte 1:   255 (fixed)
*   Byte 2: 0-255 (red)
*   Byte 3: 0-255 (green)
*   Byte 4: 0-255 (blue)
* Index codes produced by this routine can be used to set the color
* with the PGSCI routine or produce a true color pixel dump with the
* PGPIXL routine. Obviously, this does only work for `true color' devices
* like /PPM.
*
* Arguments:
*  R    (input): Red value   (0.0 to 1.0)
*  G    (input): Green value (0.0 to 1.0)
*  B    (input): Blue value  (0.0 to 1.0)
*  HEX (output): True color index (-256^3 to -1)
*-
*  9-Aug-1993 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      INTEGER*1 HEX(4)
      HEX(1)=255
      HEX(2)=NINT(R*255)
      HEX(3)=NINT(G*255)
      HEX(4)=NINT(B*255)
      END
