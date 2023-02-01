**LAND -- Check whether point is over land
*+
      FUNCTION LAND (LAT, LON)
      LOGICAL LAND
      REAL*4  LAT, LON
*
* This function queries the direct-access binary data file "ut_land12.bin"
* to determine if a given latitude/longitude location is on land or over the
* ocean.
* It returns a .TRUE. value if the point is on land or a .FALSE. value if the
* point is over the ocean.  The landmask file is a 1/12th degree grid of the
* world, with single bits set to 1 for Land, 0 for Ocean.
*
* The file "ut_land12.bin" must be located in directory $ALTIM/data.
*
* Arguments:
*  LAT   (input): Latitude of the point (deg)
*  LON   (input): Longitude of the point (deg) (non restricted)
*  LAND (output): .TRUE. if point is over land, .FALSE. otherwise.
*-
* John L. Lillibridge @ NOAA/NOS/C&GS/NGSD: 9/24/91.
* Remko Scharroo @ DUT/DEOS:
* 12-Nov-1993 - All output to stderr through routine warning of "fastio.c"
* 22-Nov-1993 - Include byte swapping for Little Endian. Use cyclic lon.
* 11-Jan-1996 - Include getenv of LAND_DIR. files.h removed.
* 10-Aug-1998 - PERRORF removed, FIN introduced.
* 17-Jun-2003 - Remove sysdep.h
*-----------------------------------------------------------------------
* Note the use of base 0 addressing of the 2-D mask array to eliminate
* arithmetic steps below.

      integer*4 irec,iwrd,ibit,mask(0:134,0:2159)
      integer*4 ios,nbytes,rbytes,l
      parameter (nbytes=2160*540)
      logical*4	btest,ltlend

* Need to specify a character array, and equivalence it to the I*4 mask
* array in order to use the fast, low-level UNIX read(2) routine for ingestion
* of the landmask file on the first call.

      character*1 cbuf(nbytes)
      equivalence (cbuf,mask)

* Functions required for the fast low-level read.

      integer*4 openf,closef,readf,fd

* Get data directory name from command line

      character*80 land_file

* Save the mask data between calls of the function, and use the logical
* init flag so the file is read in on the first call only.

      logical init
      save init,mask
      data init/.True./

* Read in the landmask array on the first call, open the file read only.

      if(init)then
        init=.False.
	land_file='/user/altim'
	call checkenv('ALTIM',land_file,l)
        land_file(l+1:)='/data/ut_land12.bin'
        fd=openf(land_file,'r')
        if(fd.lt.0)call fin("land: Error opening landmask")

* Use a single read to ingest the entire file & check that we get the
* right number of bytes: 180*12*360*12/8 = 2160*135*4 = 1,166,400.
* Place the data in the character array cbuf, which is equivalenced to
* the integer*4 2-D mask array.
* Swap integers if machine has Little Endian notation.

        rbytes=readf(fd,nbytes,cbuf)
        if(rbytes.ne.nbytes) call fin("land: Error reading landmask")

        ios=closef(fd)

        if (ltlend()) call i4swap(nbytes/4,mask)

      endif
 
* Change from E/W longitudes to positive only E longitude.

      lon=mod(lon,360.0)
      if(lon.lt.0)lon=360.0+lon

* Calculate the record number from 1/12th degree latitude grid starting at 
* the S. Pole. Note that the first row in the mask array is number 0.

      irec=int((lat+90.0)*12.0)

* Calculate the total bit offset within the strip from the longitude,
* and convert that to the 32-bit word number (column) and the bit offset 
* within that word.  Again, the first column in the mask array is number 0.

      ibit=int(lon*12.0)

* Check that the lat,lon input values are in-bounds.

      if(irec.lt.0)then
        irec=0
        call warning('land: Warning, lat < -90')
      elseif(irec.gt.2159)then
        irec=2159
        call warning('land: Warning, lat >= 90')
      endif

* Convert total bit offset to column/word number and bit offset (0->31).

      iwrd=ibit/32
      ibit=mod(ibit,32)

* Finally, if this gridcell`s bit is set, the point is over Land and we
* return a true value.

      land=btest(mask(iwrd,irec),ibit)

      return
      end
