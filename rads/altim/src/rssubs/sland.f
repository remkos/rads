**SLAND -- Check whether point is over deep ocean
*+
      FUNCTION SLAND (LAT, LON)
      LOGICAL SLAND
      REAL*4  LAT, LON
*
* This function queries the direct access binary file "sland3.bin" to
* determine if a given latitude/longitude location is over land or "shallow"
* seas ( < 2251 m), or is over "deep" seas ( > 2251 m).
* The function returns a logical value of .TRUE. for land/shallow sea and a
* value of .FALSE. for deep sea points.
*
* The data file "sland3.bin" is derived from Russ Agreen`s sequential binary 
* version.  One-degree squares in the region 79 S -> 81 N; 0 -> 360 E are
* represented by bits within the mask file. Hence there are 360*160=57600 bits,
* which are contained within 7200 bytes, or 1800 integer*4 words.  A bit "set"
* to "1" indicates a land/shallow sea 1-degree square, while a bit "cleared" 
* to "0" indicates a deep sea square.  The cycling of the indices has the bit
* number cycling fastest from 0->31, the longitude cycles next fastest from
* 0->360, and the latitude cycles slowest from -79->81. Since 360 bits do not 
* divide evenly by the 32-bit I*4 word boundary, the array is 1-D vs. 2-D 
* (lat/long). The element or word number within the mask array is determined
* from a combination of latitude and longitude: ilat=int(lat+79); 
* ilon=int(lon); ibit=ilat*360+ilon; iwrd=ibit/32+1; ibit=mod(ibit,32).  
* The zero bit of the first word of the array corresponds to 79S,0E and bit 31
* of the 1800th word corresponds to 81N,360E. Note that the grid cells are 
* centered on 1/2 degree lat/long points... it is the EDGES that fall on the 
* whole degree boundaries.  
*
* Any locations south of 79S or north of 81N are considered land/shallow sea 
* and a .TRUE. value is returned.
*
* The file "sland3.bin" must be located in directory $ALTIM/data.
*
* Arguments:
*  LAT    (input): Latitude of the point (deg)
*  LON    (input): Longitude of the point (deg)
*  SLAND (output): Flag indicating whether point is over deep ocean (.TRUE.
*                  if depth > 2251 m)
*-
* John L. Lillibridge @ NOAA/NOS/CGS/NGSD: 9/24/91.
* Remko Scharroo @ DUT/DEOS:
* 22-Nov-1993 - Use byte swapping for Little Endian machines
*  7-Jan-1996 - Include getenv of ALTIM_DATA. files.h rmoved.
* 10-Aug-1998 - PERRORF removed, FIN introduced.
* 17-Jun-2003 - Remove sysdep.h
*-----------------------------------------------------------------------
      integer*4   nbytes
      parameter   (nbytes=7200)
      integer*4   mask(nbytes/4),irec,iwrd,ibit,ios,l
      logical*4   btest,ltlend

      character*1 cbuf(nbytes)
      equivalence (cbuf,mask)

* Low-Level fast-IO routines for ingesting data file

      integer*4 fd,openf,readf,closef

* Get data directory name

      character*80 sea_file

* Save the mask array and use logical flag init so we only have to read
* in the mask file on the first call.

      logical init
      save init,mask
      data init/.True./

* Read in the shallow sea mask on the first call using a single low-level
* "read(2)" C-routine, read only.
* Swap integers if machine has Big Endian notation.

      if(init)then
        init=.False.
	sea_file='/user/altim'
        call checkenv('ALTIM',sea_file,l)
        sea_file(l+1:)='/data/sland3.bin'
        fd=openf(sea_file,'r')
        if(fd.lt.0)call fin("sland: error opening sea_file")

        ios=readf(fd,nbytes,cbuf)
        if(ios.ne.nbytes)call fin("sland: error reading sea_file")

        ios=closef(fd)

        if (ltlend()) call i4swap(nbytes/4,mask)

      endif

* Convert the latitude to a "record" number starting at 79S. Points out of
* bounds to the N or S are considered land/shallow sea.

      irec=int(lat+79.0)
      if(irec.lt.0.or.irec.gt.159)then
        sland=.True.
        return
      endif

* Convert the longitude to a bit offset within the "record".
* Use a cyclic representation of longitude for either E/W (-180 -> 180)
* or positive E (0 -> 360) values of longitude.

      ibit=int(lon)
      ibit=mod(ibit,360)
      if(ibit.lt.0)ibit=360+ibit

* Now combine the latitude "record" with the longitude "offset" to get the
* total bit offset within the file, and compute the 32-bit word and bit
* offset within that word from the total # of bits.

      ibit=irec*360+ibit
      iwrd=ibit/32+1
      ibit=mod(ibit,32)

* Finally, test the bit for this gridcell: a "1" (set) corresponds to
* land or shallow sea and a TRUE value is returned; a "0" (clear) corresponds
* to blue water ocean and a FALSE value is returned.

      sland=btest(mask(iwrd),ibit)
      return

      end
