      program max

* MAX main program
*
* 16-Oct-1994 - 9410.2 - Blocked reading of ADR file (both in SORTRK and XOFIND)
*   to minimize disk access. Old binary formats removed. Xover guessing
*   and determination combined.
* 20-Oct-1994 - 9410.3 - -nNSKIP flag introduced. Argument line scanning in
*   max.f. Bug fixed in DAREAD.
* 21-Oct-1994 - 9410.4 - Improved splitting of tracks
*               9410.5 - Minimum angle can be specified on command line
*               9410.6 - Attempt to speed up xovdet
* 22-Nov-1994 - 9411.0 - New version of getsat
*  9-Aug-1995 - 9508.0 - Limits on crossover time difference through dtxo
*   and dtxoff
* 14-Sep-1995 - 9509.0 - Set parameter 1 to -altbias. Include boundaries and
*   iteration number in file headers.
* 19-Sep-1995 - 9509.1 - Adopted new satcat.
*  3-Oct-1995 - 9510.0 - LINFIT completely removed from XOVDET
*  9-Oct-1995 - 9510.1 - Auxiliary fields
*  7-Dec-1995 - 9512.0 - File headers properly initialized. Bug removed that
*   led to crash when auxiliary data where processed.
* 16-Jul-1996 - 9607.0 - "nomerge" option to avoid crossing of data from
*   different ADR input files.
* 14-Aug-1996 - 9608.0 - Include option for cubic polynomial interpolation
*   and make it default.
* 13-Sep-1996 - 9609.0 - Small adjustment to xovdet to have rlonx be between
*   minlon and maxlon
* 15-Feb-1998 - 9802.0 - Obliterate getsat. Bug fixes in xofind.
*  2-Jul-1998 - 9807.0 - Increase number of tracks to 65536. Tracknumber is
*   automatically stored as negative for itrk > 32767. Maxblocks increased
*   to 400.
*  8-Jul-2000 - 2000.0 - Make sure MAX does not crash when nrxfnd=0
*  3-Jun-2002 - 0206.0 - Minor changes

      character*6 version
      parameter (version="0206.0")

      include "maxcom.inc"
      character*80 usestr,arg
      integer*4 i,iargc,j
      integer ids(maxsat),mjd0(maxsat),mjd1(maxsat)
      real*4  orbalt(maxsat),distmin(maxsat),maxerr(maxsat)

      namelist /nml/ npar,specif,timper,nskip,shortl,usestr,angmin,
     |		dtxoff,dtxo
      namelist /satcat_nml/ name,orberr,altsig,ids,abias,mjd0,mjd1,
     |          orbinc,orbalt,distmin,maxerr

* Initialisation and defaults.
* dtxoff(i,j) is mean time difference between satellites i and j.
* dtxo(i,j) is maximum time difference between satellites i and j around
* dtxoff(i,j).
* Example: allow combination of GEOSAT and ERS-1 data with time difference
* of 35 days around a mean of 3 years: dtxoff(3,4)=1095, dtxo(3,4)=35
* (to be set in max.nml). Automatically dtxoff(4,3) will be set to -1095
* and dtxo(4,3) to 35.
 
      dtrms=0
      dtmax=0
      xrms=0

      do i=1,maxsat
	 do j=1,maxsat
	    dtxoff(i,j)=0
	    dtxo(i,j)=1d35
	 enddo
      enddo

      arg='/user/altim'
      call checkenv('ALTIM',arg,i)
      open (9,file=arg(:i)//'/nml/max.nml',status='old',err=10)
      read (9,nml)
      close (9)
10    open (9,file='max.nml',status='old',err=11)
      read (9,nml)
      close (9)
11    continue

      open (9,file=arg(:i)//'/nml/satcat.nml',status='old',err=20)
      read (9,satcat_nml)
      close (9)
20    open (9,file='satcat.nml',status='old',err=21)
      read (9,satcat_nml)
      close (9)
21    continue

* More initialisations

      nrfil=0
      nrxins=0
      nrxdt=0
      nrxinc=0
      nrxiter=0
      nrxdat=0
      nrxfnd=0

      afile='-'
      xfile='-'
      tfile='-'
      nomerge=.false.

* Read Input Parameters, partly stored in common /inparm/

      do i=1,iargc()
	 call getarg(i,arg)
	 if (arg(1:4).eq.'nml=') then
	    open (9,file=arg(5:),status='old')
	    read (9,nml)
	    close (9)
	 else if (arg(1:2).eq.'-p') then
	    read (arg(3:),*) npar
	 else if (arg(1:2).eq.'-a') then
	    read (arg(3:),*) angmin
	 else if (arg(1:2).eq.'-t') then
	    read (arg(3:),*) timper
	 else if (arg(1:4).eq.'-nom') then
	    nomerge=.true.
	 else if (arg(1:2).eq.'-n') then
	    read (arg(3:),*) nskip
	 else if (arg(1:2).eq.'-m') then
	    read (arg(3:),*) shortl
	 else if (arg(1:2).eq.'-s') then
	    specif=arg(3:)
	 else if (arg(1:2).eq.'-u') then
	    usestr=arg(3:)
	 else if (arg.eq.'-xxa') then
	    specif='AC-A-'
	 else if (arg.eq.'-xxo') then
	    specif='AC-O-'
	 else if (arg.eq.'-xxs') then
	    specif='AC-S-'
	 else if (arg(1:2).eq.'-x') then
	    read (arg(3:),*) naux
	 else if (arg(1:1).eq.'-') then
	    goto 1300
	 else
	    nrfil=nrfil+1
	    if (nrfil.gt.maxfil) then
	       write (*,*) 'Sorry: no more files'
	       stop
	    endif
	    infile(nrfil)=arg
	 endif
      enddo
      angmin=angmin*rad

      if (nrfil.eq.0) goto 1300
      if (nrfil.gt.1) nrfil=nrfil-1
      i=index(arg,' ')-1
      if (specif(3:3).eq.'A') afile=arg(:i)//'.xaf'
      if (specif(4:4).eq.'X') xfile=arg(:i)//'.xxf'
      if (specif(4:4).eq.'A') xfile=arg(:i)//'.xxa'
      if (specif(4:4).eq.'O') xfile=arg(:i)//'.xxo'
      if (specif(4:4).eq.'S') xfile=arg(:i)//'.xxs'
      if (specif(5:5).eq.'T') tfile=arg(:i)//'.xtf'

      do i=1,maxpar
         use(i)=.not.(usestr(i:i).eq.'-')
      enddo

* Scan the dtxoff and dtxo.
* Choose for dtxoff the largest one, for dtxo the one that is smallest.
* Convert to seconds.

      do i=1,maxsat
	 do j=1,maxsat
	    if (abs(dtxoff(i,j)).gt.abs(dtxoff(j,i))) then
	       dtxoff(j,i)=-dtxoff(i,j)
	    endif
	    if (dtxo(i,j).lt.dtxo(j,i)) then
	       dtxo(j,i)=dtxo(i,j)
	    endif
	 enddo
      enddo
      do i=1,maxsat
	 do j=1,maxsat
	    dtxo(i,j)=dtxo(i,j)*86400
	    dtxoff(i,j)=dtxoff(i,j)*86400
	 enddo
      enddo

      write (*,100) version,timper,specif,usestr,nint(shortl),
     |              nskip,angmin/rad,afile,xfile,tfile

* Start Crossover generation by calling the subprogrammes
*
* SORTRK sorts the observations in tracks, determines inclination,
* and produces crossover location guesses

      call sortrk

* XOFIND is a routine to find locations of crossovers

      if (specif(4:4).ne.'-') call xofind

* WRIXTF writes the track catalogue file

      if (nrxfnd.eq.0) then
	 write (*,550) '   No xovers found'
      else
         if (specif(5:5).eq.'T') call wrixtf
         write (*,110)
      endif
      goto 9999

* Manual

1300  write (*,1301) version
1301  format(
     |'MAX version ',a//
     |'usage: max [ options ] adr(s) output_prefix'//
     |'with'/
     |'  adr(s)   : Input ADR file(s)'/
     |'  prefix   : Optional name for XAF, XXF and XTF files',
     |' (excl. extension)'/
     |'             If not specified: (last) ADR filename is',
     |' used for this purpose.'//
     |'and [options]'/
     |'-aANGMIN : Minimum angle between tracks at crossover (degrees)'/
     |'-pNPAR   : Number of orbit error parameters per track'/
     |'-nNSKIP  : Skip one crossover out of every NSKIP guessed'/
     |'-sSPECIF : specify 5 switches as character*5. Default = APAXT.'/
     |'       1 : Crossover type: A = all, S = singles only,',
     |' D = duals only'/
     |'       2 : Interpolation : C = cubic spline,'/
     |'                           Q = quadratic poly, P = cubic poly'/
     |'       3 : Generate XAF  : A = yes, - = no'/
     |'       4 : Generate XXF  : X = standard, S = short,',
     |' A = ascii,'/
     |'                                  O = incl orbit, - = no.'/
     |'       5 : Generate XTF         : T = yes, - = no'/
     |'-tTIMPER : Maximum time interval between first or last',
     |' observation used for'/
     |'           interpolation and crossover epoch in seconds.',
     |' Default = 12.'/
     |'-mMINLEN : Minimum length in kilometers for ''long'' track'/
     |'-uUSE    : Parameters to be used for short track. Example:'/
     |'           -u-XX-- <- only 2nd and 3rd parameters is',
     |' estimated if track is short')
* Formats

  100 format(1x,78('*')/' *',t38,'MAX ',a,t79,'*'/1x,78('*')//
     |'Options'/
     |'   Maximum time interval over 6 observations   -> ',f9.3,' sec.'/
     |'   Process specifier                           -> ',4x,a5/
     |'   Short track orbit model specifier           -> ',4x,a5/
     |'   Short track length                          -> ',i9,' km.'/
     |'   Interval for skipping crossovers (0=none)   -> ',i9/
     |'   Minimum angle between tracks                -> ',f9.3,' deg'//
     |'Output files'/
     |'   Altimeter file       -> ',a50/
     |'   Crossover file       -> ',a50/
     |'   Tracks catalogue     -> ',a50//
     |'Reading ADR files ...')
  110 format(/'Crossover processing succesfully ended')
  550 format(a)
 9999 end
