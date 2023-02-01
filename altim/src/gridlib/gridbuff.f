**GRIDBUFF -- Load a grid into memory
*+
      FUNCTION GRIDBUFF (FILENM, POINTER)
      INTEGER*4 GRIDBUFF, POINTER
      CHARACTER FILENM*(*)
*
* This routine allocates memory and loads the contents of a grid file
* into the allocated memory. In contrast to the GRIDRD4 and GRIDRD8
* routines, the grid values are not directly available.
*
* Other information on the grid, like grid dimension, can be obtained
* with the GRIDBINF routine.
*
* After using GRIDBUFF, GRIDBINT or GRIDBSPL can be used to interpolate
* the grid, based on bi-linear or bi-cubic spline interpolation.
*
* The GRIDBUFF routine handles grids of either 2-, 3- or 4-byte integers,
* in little or big endian notation, on little or big endian machines.
* These grids can be recognised by their @GRID, @GR2L, @GR2B, @GR3L,
* @GR3B, @GR4L or @GR4B preamble.
*
* GRIDBUFF is also more efficient in the memory use than GRIDRD4 and
* GRIDRD8. In case of 2-byte grids, it uses almost half the memory
* GRIDRD4 uses, and about a quarter of the memory GRIDRD8 needs.
* This is achieved by loading the original 2-byte-per-cell grid into
* memory without further manipulation. The value returned by GRIDBUFF
* is a pointer to a structure that includes both the grid dimensions
* and the grid itself.
*
* When the allocation of memory or the loading of the grid was
* unsuccessful, this will be reflected in the returned pointer value.
*
* Input argument:
*  FILENM   : Name of the file containing the grid.
*
* Output argument:
*  POINTER  : Pointer to the grid structure
*
* Return value:
*  GRIDBUFF : 0 = No error
*             1 = Could not find or open grid
*             2 = Ilegal grid format
*             3 = Memory allocation was unsuccessful
*             4 = Error loading grid
*             5 = Wrong specifier
*-
* 22-Nov-2000 - Created by Remko Scharroo
* 21-Aug-2005 - Completely redeveloped
*-----------------------------------------------------------------------
      integer*4 gridbuff_deos,gridbuff_nc,l,lnblnk

* Try to open as DEOS grid. If wrong format, try NetCDF

      gridbuff=gridbuff_deos(filenm,pointer)
      if (gridbuff.eq.5) gridbuff=gridbuff_nc(filenm,pointer)
      if (gridbuff.eq.0)  return

      l=lnblnk(filenm)
1300  format ('GRIDBUFF: ',a,1x,a)
      if (gridbuff.eq.1) then
	 write(0,1300) 'file not found:',filenm(:l)
      else if (gridbuff.eq.2) then
	 write(0,1300) 'illegal grid format in',filenm(:l)
      else if (gridbuff.eq.3) then
	 write(0,1300) 'error allocating memory for',filenm(:l)
      else if (gridbuff.eq.4) then
	 write(0,1300) 'error loading grid',filenm(:l)
      else if (gridbuff.eq.5) then
	 write(0,1300) 'unknown file type in',filenm(:l)
      endif
      end
