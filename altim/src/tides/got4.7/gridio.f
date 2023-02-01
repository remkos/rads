      SUBROUTINE gridin( G,NDIM,NX,NY,TITLE,LU,
     *                   LATMIN,LATMAX,LONMIN,LONMAX,UNDEF )
*
*  Reads an ascii gridfile in the "standard" format.
*  This format uses 80-byte card-image records in ascii.
*
*  Calling arguments:
*
*   G     - O - two-dimensional array to contain the gridded data.
*               The grid will be loaded in the following way:
*               G(i,j), with i going west to east and j south to north.
*               G(1,1) is thus the southwest corner of grid.
*               G(NX,NY) is the northeast corner.
*               User must make certain G is dimensioned to at least
*               the expected values of NX,NY.
*
*   NDIM  - I - first dimension of G in calling program.
*               NDIM must not be less than the expected value of NX.
*
*   NX,NY - O - size of the grid G.
*
*   TITLE - O - 160-byte character-string title that describes data.
*
*   LU    - I - fortran unit to use for reading input data.
*
*   LATMIN,LATMAX,LONMIN,LONMAX - O - area limits for grid,
*               in decimal degrees.
*               Note: these are REAL variables.
*             The grid intervals are therefore:
*             DX = (LONMAX-LONMIN)/(NX-1)  and
*             DY = (LATMAX-LATMIN)/(NY-1).
*
*   UNDEF - O - value to denote a null or missing value in G.
*
*
*  Written by R. Ray      Aug. 1990
*  Revised 8 Apr 1998 
*    - change name.
*    - handle 2 nulls on input file.
*

      DIMENSION   G(NDIM,*)
      REAL        LATMIN,LATMAX,LONMIN,LONMAX
      CHARACTER   TITLE*160, FORMAT*80

      READ(LU,1) TITLE(1:80)
      READ(LU,1) TITLE(81:160)
    1 FORMAT(A80)

      READ(LU,2) NY,NX
    2 FORMAT(16X,I5,16X,I5)
      IF (NX.GT.NDIM) THEN
         WRITE(6,*) 'Increase NDIM in call to UTCSRI'
         STOP
      ENDIF
      READ(LU,3) LATMIN,LATMAX
      READ(LU,3) LONMIN,LONMAX
    3 FORMAT(16X,F9.0,16X,F9.0)
      READ(LU,3) UNDEF1, UNDEF2
      UNDEF = UNDEF2
      READ(LU,1) FORMAT
      DO 4 J=1,NY
         READ(LU,FORMAT) (G(I,J),I=1,NX)
    4 CONTINUE
      if (undef1.ne.undef2) then
         do 5 j=1,ny
         do 5 i=1,nx
            if (g(i,j).eq.undef1) g(i,j) = undef
    5    continue
      endif
      RETURN
      END

      SUBROUTINE gridut( G,NDIM,NX,NY,TITLE,LU,
     *                   LATMIN,LATMAX,LONMIN,LONMAX,UNDEF,FORMAT )

*  Outputs an ascii gridfile in "standard" format.
*  This format uses 80-byte card-image records in ascii.
*
*  Calling arguments:
*
*   G     - I - input two-dimensional array containing gridded data.
*               The grid should be loaded in the following way:
*               G(i,j), with i going west to east and j south to north.
*               (same as for DIUTIL contouring package).
*               G(1,1) is southwest corner of grid.
*
*   NDIM  - I - first dimension of G in calling program.
*
*   NX,NY - I - size of G.  Thus G(NX,NY) is northeast corner of grid.
*
*   TITLE - I - 160-byte character-string title for image.
*
*   LU    - I - fortran unit to use for writing output data.
*
*   LATMIN,LATMAX,LONMIN,LONMAX - I - area limits for grid,
*               in decimal degrees.
*               Note: these are REAL variables.
*
*   UNDEF - I - value to denote a null or missing value in G.
*
*   FORMAT- I - character string containing fortran format to use
*               for writing G.
*               Example:   '(10F8.2)'
*               Warning: If G contains values too large for FORMAT,
*                        '*'s can be written to the output dataset.
*                        User should check that this hasn't happened.
*
*
*  Written by R. Ray      Aug. 1990
*
*  Revised - R. Ray  5/4/91
*    1. Included UNDEF in calling arguments.
*    2. Removed REWIND & error checking (for '*'s) because of this
*       routine's possible use in multiple calls for many outputs.
*    3. Allowed FORMAT variable length character.
*  Revised 3/8/98 - changed name.

      DIMENSION   G(NDIM,*)
      REAL        LATMIN,LATMAX,LONMIN,LONMAX
      CHARACTER   TITLE*160, FORMAT*(*)
      INTRINSIC   LEN,MIN

      WRITE(LU,1) TITLE(1:80)
      WRITE(LU,1) TITLE(81:160)
    1 FORMAT(A80)
   11 FORMAT(80A1)

      WRITE(LU,2) NY,NX
    2 FORMAT(16X,I5,16X,I5)
      WRITE(LU,3) LATMIN,LATMAX
      WRITE(LU,3) LONMIN,LONMAX
    3 FORMAT(16X,F9.4,16X,F9.4)
   13 FORMAT(16X,F9.2,16X,F9.2)
      WRITE(LU,13) UNDEF,UNDEF
      LFL = MIN(LEN(FORMAT),80)
      WRITE(LU,11) (FORMAT(JJ:JJ),JJ=1,LFL)
      DO 4 J=1,NY
         WRITE(LU,FORMAT) (G(I,J),I=1,NX)
    4 CONTINUE
      RETURN
      END
