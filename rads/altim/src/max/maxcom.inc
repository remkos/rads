**MAXCOM -- Parameters and common blocks for MAX
*+
* Parameters:
* pi, twopi, rad : pi, 2*pi, and conversion factor degrees -> radians
* maxtrk : maximum number of tracks
* maxobs : maximum number of observations per track
* maxsat : maximum number of satellites
* maxpar : maximum number of orbit error parameters per track
* maxobsblk : maximum number of observations per block
* maxfil : maximum number of input (ADR) files
* maxint2, maxint4 : maximum integer*2 and integer*4
* rotate : earth rotation velocity (rad/s)
* rdunit : unit number for reading (ADR)
* wrunit : unit number for writing (XAF, XXF, XTF)
*-
      REAL*8 PI, RAD, TWOPI, MURAD
      PARAMETER (PI=3.14159265358979D0,RAD=PI/180D0,
     |   TWOPI=2*PI,MURAD=RAD/1D6)

      integer*4  maxtrk,maxobs,maxpar,maxfil,maxsat
      integer*4  maxint2,maxint4,rdunit,wrunit
      integer*4  maxobsblk,mblocks,maux,naux,ladr,nblocks,nskip
      real*8     rotate
      parameter (maxint2=32767,maxint4=214748364)
      parameter (maxtrk=65536,maxobs=5000,maxsat=12,
     |		 maxobsblk=300000,mblocks=400,
     |           maxpar=5,maxfil=200,maux=3)
      parameter (rotate=7.292115855d-5)
      parameter (rdunit=10,wrunit=20)
*+
* Arrays and variables:
* locat(4,maxtrk): lat and lon of start and end point of track (rad)
* node (maxtrk)  : longitude of nearest/ascending node (rad)
* inc(maxtrk)    : inclination (rad)
* theta(maxtrk)  : orbital velocity (rad/s)
* time(3,maxtrk) : time of start and end and node of track (sec85)
* track(5,maxobs): time(sec85), lat, lon(rad), orb(km) and height(m) of obs
* atrack(3,maxobs): auxiliary data per obs
* u0(maxtrk)     : argument of latitude of start of track (rad)
* istrec(2,maxtrk): record number and filenr for start of track
* iobs(maxtrk)   : number of observation per track
* ixpt(maxtrk)   : number of crossovers per track
* timper         : maximum time period over 6 observations around crossover
* shortl         : short track for tracklength < shortl
* distmn         : mean track length
* satid(maxtrk)  : satellite id per track
* filenr(maxtrk) : file number from which track originated
* npar           : number of orbit error parameters
* nrxins         : number of track combinations do not intersect
* nrxdt          : number of xover rejected because of time difference
* nrxinc         : number of xover rejected on small angle
* nrxiter        : number of xover rejected because of nr of iterations
* nrxdat         : number of xovers with too few data.
* nrxfnd         : number of processed crossovers
* nrget          : number of times that blocks are read from file
* ishrt          : number of long tracks
* nrfil          : number of input files
* itrknr         : current track nr or nr of tracks
* minlon,maxlon  : minimum and maximum longitude
* bound(4)	 : area boundaries
* specif         : selection word 'ACXAF' etc
* dtrms,dtmax    : RMS and maximum time correction on xover guess
* xrms           : xover difference RMS
* angmin         : minimum angle between tracks at crossover
* cosmax         : maximum cos of angle between tracks at xover
* nskip          : skip one per nskip crossovers
*-
      real*8 locat(4,maxtrk),node(maxtrk),inc(maxtrk),
     |       theta(maxtrk),time(3,maxtrk),track(5,maxobs),
     |       u0(maxtrk),angmin,hfac
      integer*4 mtrkblk(0:mblocks),blkfil(0:mblocks),
     |       blkrec(0:mblocks),blkobs(0:mblocks)
      integer*4 bdata ( 6,maxobsblk,2),bdatx ( 7,maxobsblk,2)
      integer*2 bdata2(12,maxobsblk,2),bdatx2(14,maxobsblk,2)
      equivalence (bdatx,bdata),(bdata,bdata2),(bdatx,bdatx2)
      real*8 timper,shortl,distmn,dtrms,dtmax,xrms,cosmax
      integer*2 atrack(maux,maxobs)
      integer*2 satid(maxtrk),iobs(maxtrk),ixpt(maxtrk),filenr(maxtrk)
      integer*4 istrec(maxtrk),npar,jiter
      integer*4 nrxiter,nrxinc,nrxins,nrxdt,nrxdat,nrxfnd,nrget,
     |	     ishrt,nrfil,itrknr
      real*8 minlon,maxlon
      integer*2 bound(4)
      character*5 specif
      character*80 infile(maxfil),afile,xfile,tfile
      logical asc(maxtrk),good(maxtrk),nomerge,
     |        short(maxtrk),progrd(maxtrk),use(maxpar),fast

      character*8 name(maxsat)
      real*4 orberr(maxsat),altsig(maxsat),orbinc(maxsat),abias(maxsat)

      common /trax/ u0,track,time,locat,node,inc,theta,		! real*8
     |		    istrec,npar,asc,good,short,progrd,		! int*4, log*4
     |              atrack,satid,filenr,iobs,ixpt		! int*2
      common /inparm/ angmin,cosmax,timper,shortl,minlon,maxlon,! real*8
     |		      nskip,fast,nomerge,use,			! int*4, log*4
     |		      bound,					! int*2
     |                infile,afile,xfile,tfile,specif		! char
      common /cnstnt/ nrfil,itrknr,nrxiter,nrxinc,nrxins,nrxdt,
     |		      nrxdat,nrxfnd,nrget
      common /cstat/  distmn,dtrms,dtmax,xrms,jiter,ishrt
      common /satcom/ name,orberr,altsig,orbinc,abias
      common /blocks/ hfac,bdatx2,mtrkblk,blkfil,blkrec,blkobs,
     |		      nblocks,naux,ladr

*+
* Delta-Time limits on xovers
*-
      real*8 dtxoff(maxsat,maxsat),dtxo(maxsat,maxsat)
      common /dtxo/ dtxoff,dtxo
