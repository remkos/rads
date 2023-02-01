**SYSDEP -- Do a number of system dependent checks
*+
      program sysdep
*
* This program creates files related to system dependencies
* - sysdep.h   : System dependent stuff as C-header file
* - sysdep.mk  : System dependent stuff as Makefile inclusion
*-
* 13-Feb-1998 - Remko Scharroo
*-----------------------------------------------------------------------

      integer*4 ios,underscore
      integer*4 itest(2)
      character*4 ctest(2)
      equivalence (ctest,itest)
      character*1 c
      data ctest /'1234','4321'/

* Open the output files

      open (20,file='sysdep.h')
      open (22,file='sysdep.mk')

* Check byte representation

      if (itest(1).gt.itest(2)) then
	 write (20,600)
      else
	 write (20,601)
      endif

* Check whether bytes or words are used in recl= argument of open

      open (10,file='sysdep.f',form='unformatted',status='old',
     |	recl=1,access='direct')
      read (10,rec=2) c
      if (c.eq.'*') then
         write (22,620)
      else
         write (22,621)
      endif

* Check the use of underscores in Fortran objects

      ios=underscore()
      if (ios.eq.1) then
         write (20,631)
      else if (ios.eq.2) then
         write (20,632)
      else if (ios.eq.3) then
         write (20,633)
      else
         write (20,630)
      endif

600   format('#define SWAP')
601   format('#undef SWAP')
620   format(50('#')/
     |'# The record length is in bytes. You are not SGI #'/
     |50('#')/
     |'FFLAGS2 =')
621   format(51('#')/
     |'# The record length is in words. You are like SGI #'/
     |51('#')/
     |'FFLAGS2 = -bytereclen')
630   format('#undef CAPITALS'/
     |'#undef UNDERSCOREAFTER'/'#undef UNDERSCOREBEFORE')
631   format('#define CAPITALS'/
     |'#undef UNDERSCOREAFTER'/'#undef UNDERSCOREBEFORE')
632   format('#undef CAPITALS'/
     |'#define UNDERSCOREAFTER'/'#undef UNDERSCOREBEFORE')
633   format('#undef CAPITALS'/
     |'#undef UNDERSCOREAFTER'/'#define UNDERSCOREBEFORE')

      close (10)
      close (20)
      close (22)
      end
