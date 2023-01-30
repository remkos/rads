**GIMREAD -- Read GPS Ionosphere Map (IONEX version)
*+
      FUNCTION GIMREAD (MJD, TYPE, TECMAP)
      INTEGER*4 MJD, GIMREAD
      CHARACTER TYPE*4
      INTEGER*2 TECMAP(73,71,*)

* This routine loads the TEC map stored in an IONEX file into
* memory. Some minor checks are made on the sanity of the data.
* It is assumed that all grids cover the same area and
* have the same grid spacing and the each daily file contains
* twelve or thirteen TEC maps.
*
* Longitude: -180 to +180 with intervals of 5 degrees (73 grid lines)
* Latitude:  -87.5 to +87.5 with intervals of 2.5 degrees (71 grid lines)
* TEC:       In units of 0.1 TEC Units, or  10**15 electrons/m**2
*
* IONEX files should be stored in $ALTIM/data/gim/YYYY, where ALTIM
* is an environment variable containing a directory name chosen by
* the user, with /user/altim as default. YYYY is the year (e.g. 2002).
* The IONEX files are named either:
* - JPLGddd0.yyI.gz (JPL analysis files, Gnu compressed)
* - JPLQddd0.yyI.gz (JPL quick-look files, Gnu compressed)
* - CODGddd0.yyI.Z  (CODE analysis files, Unix compressed)
* where yy is the year (e.g. 02) and ddd the day number (e.g. 011)
*
* Temporary uncompressed files are written to /tmp.
* Requires "gunzip" program.
* Files are opened on a self-chosen free unit number.
*
* Arguments:
*  MJD      (input) : Modified Julian Date of IONEX file
*  TYPE     (input) : Source of IONEX file (JPLG, JPLQ or CODG)
*  TECMAP  (output) : Integer array with the 12 or 13 TEC maps [in 0.1 TECU]
*  GIMREAD (output) : Status report: 1 = file does not exist,
*                     2 = error opening file, 3 = format error,
*                     0 = no error.
*-
* 23-May-2002 - Created by Remko Scharroo
* 15-Nov-2002 - Adjusted to allow new GIM files (13 maps/file)
*  4-Dec-2002 - Added JPLQ. Updated manual.
*-----------------------------------------------------------------------
      character*80 line
      character*256 filenm,dirnm,cmd,ext*3
      integer*4 unit,freeunit,l,lnblnk,i,j/0/,k,yyyy,ddd,mdate,ldir,
     |		next,nmaps
      logical   exist

* Construct the pathname of the IONEX file

      dirnm='/user/altim'
      call checkenv('ALTIM',dirnm,ldir)
      yyyy=mdate(3,mjd)/10000
      ddd=mjd-mdate(2,yyyy*10000+0101)+1
      if (type(:3).eq.'COD') then
         ext='.Z'
      else
         ext='.gz'
      endif
      next=lnblnk(ext)
      nmaps=0
      write (filenm,620) dirnm(:ldir),yyyy,type,ddd,mod(yyyy,100),ext

* Open the file. Uncompress to /tmp first.

      l=lnblnk(filenm)
      inquire (file=filenm,exist=exist)
      if (.not.exist) goto 1300
      write (cmd,600) filenm(:l),filenm(l-11-next:l-next)
      filenm='/tmp/'//filenm(l-11-next:l-next)
      call system(cmd)
      unit=freeunit()
      open (unit,file=filenm,status='old',err=1310)

* Read the file. Load the data into memory. Check the sanity of the file

100   read (unit,550,end=1320) line
      if (line(61:).eq.'# OF MAPS IN FILE') then
	 read (line(:6),*) nmaps
	 if (nmaps.lt.12 .or. nmaps.gt.13) goto 1320
      else if (line(61:).eq.'HGT1 / HGT2 / DHGT') then
         if (line(:14).ne.'   450.0 450.0') goto 1320
      else if (line(61:).eq.'LAT1 / LAT2 / DLAT') then
         if (line(:20).ne.'    87.5 -87.5  -2.5') goto 1320
      else if (line(61:).eq.'LON1 / LON2 / DLON') then
         if (line(:20).ne.'  -180.0 180.0   5.0') goto 1320
      else if (line(61:).eq.'EXPONENT') then
         if (line(:6).ne.'    -1') goto 1320
      else if (line(61:).eq.'START OF TEC MAP') then
         read (line(:6),*) k
	 j=71
      else if (line(61:).eq.'LAT/LON1/LON2/DLON/H') then
         read (unit,610) (tecmap(i,j,k),i=1,73)
         j=j-1
      else if (line(61:).eq.'END OF TEC MAP') then
         if (j.ne.0) goto 1320
	 if (k.ge.nmaps) goto 200
      endif
      goto 100

* Close up

200   close (unit)
      call system('rm -f '//filenm(:l))
      gimread=0
      return

* Formats

550   format (a)
600   format ('gunzip < ',a,' > /tmp/',a)
610   format ((16i5))
620   format (a,'/data/gim/',i4.4,'/',a4,i3.3,'0.',i2.2,'I',a)

* Error exits

1300  write (0,550) ' GIMREAD: file does not exist: '//filenm(:l)
      gimread=1
      return
1310  write (0,550) ' GIMREAD: error opening file: '//filenm(:l)
      gimread=2
      return
1320  write (0,550) ' GIMREAD: error in format of file: '//filenm(:l)
      write (0,550) line
      gimread=3
      return
      end
