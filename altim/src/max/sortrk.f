**SORTRK -- MAX module to generate XAF, devide data in tracks
*+
* SORTRK reads ADR files and generates an XAF data set from these records.
* Meanwhile data are clipped into tracks running from pole to pole.
* Information on the tracks is later stored in an XTF file by module
* WRIXTF
*-
* 16-Oct-1994 - XOGUESS split off from SORTRK. Old XAF formats removed.
* 21-Oct-1994 - Better check of track splitting.
* 16-Feb-1998 - Initialize counters
*-----------------------------------------------------------------------
      subroutine sortrk
      include "maxcom.inc"
      real*8    dist,dta,dtb,dtim/3.5d3/
      real*8    ttop,toppar,lata,latb,lat
      integer*4 sid(maxsat),npoint,nobsblk
      parameter (npoint=4)
      logical xaf,prograde,contlon
      character*8 satnam
      character filspf*4
      integer*4 ifilnr,isatel,nameid,nrrecs,nrxaf
      integer*4 i,ik,irec,isub,l,null/0/
*
* Open XAF and write header
*
      xaf=(specif(3:3).eq.'A')
      l=index(afile,' ')-1
      if (xaf) then
         open (wrunit,file=afile,access='direct',form='unformatted',
     .      status='unknown',recl=28)
	 filspf='@XAB'
         write (wrunit,rec=1) filspf,null
	 if (naux.gt.0) then
	    open (wrunit+1,file=afile(:l)//'.aux',
     |		access='direct',form='unformatted',recl=2*naux)
	 endif
      endif
*
* Initialize counters
*
      distmn=0
      nrxaf=0
      nblocks=1
      ifilnr=1
      isatel=0
      ishrt=0
*
* Open the altimeter input files one by one. If last file is read jump
* to crossover guessing (1000).
*
    4 if (ifilnr.gt.nrfil) goto 1000
      blkrec(nblocks)=1
      nobsblk=0
      open (rdunit,file=infile(ifilnr),access='direct',status='old',
     .      form='unformatted',recl=24)
      read (rdunit,rec=1) filspf,satnam,bound,nrrecs

      close (rdunit)
*
* Determine whether satellite occurred before
*
      sid(isatel+1)=nameid(satnam)
      do ik=1,isatel
         if (sid(isatel+1).eq.sid(ik)) goto 2
      enddo
      isatel=isatel+1
    2 continue

      write (6,100) infile(ifilnr),bound
  100 format (/3x,'File currently read  -> ',a50/3x,
     .       'longitude boundaries -> ',i4,' , ',i4,' deg.'/3x,
     .       'latitude  boundaries -> ',i4,' , ',i4,' deg.')
      irec=0
      minlon=bound(1)*rad
      maxlon=bound(2)*rad
      call statbar(0,nrrecs,'Reading and sorting')
*
* And since empty ADRs do exist ...
*
      if (nrrecs.eq.0) then
         ifilnr=ifilnr+1
         goto 4
      endif
*
* Read altimeter input file IFILNR. Each record is read with the routine
* DAREAD.
* The data are split into tracks. The end of the track is determined by
* means of
* 1. the time interval between measurements isub-1 and isub: If it is larger
*    than DTIM, isub does not belong to the present track.
* 2. a polynomial fit through latitude and longitude coordinates:
*    If the top latitude is at a longitude between the measurements
*    isub-1 and isub, measurement isub does not belong to the track, but
*    to the next.
*
* First read 2 points. Shift longitude such that longitude always
* increases for prograde tracks and decreases for retrograde tracks.
*
5     isub=1
      call daread(ifilnr,irec+isub,nrrecs,track(1,isub),atrack(1,isub))
      if (irec+isub.ge.nrrecs) goto 50
      prograde=(orbinc(sid(ik)).lt.90)
      isub=2
      call daread(ifilnr,irec+isub,nrrecs,track(1,isub),atrack(1,isub))
      if (contlon(prograde,track(3,isub-1),track(3,isub))) goto 49

10    if (irec+isub.ge.nrrecs) goto 50
      isub=isub+1
      if (isub.gt.maxobs) goto 1145
      call daread(ifilnr,irec+isub,nrrecs,track(1,isub),atrack(1,isub))
      if (contlon(prograde,track(3,isub-1),track(3,isub))) goto 49
*
* Check whether track is too long (> dtim seconds)
*
      if (track(1,isub)-track(1,1).gt.dtim) goto 49
*
* Check if top is passed in (t,lat) space
*
      dta=track(1,isub-1)-track(1,isub-2)
      dtb=track(1,isub  )-track(1,isub-1)
      lata=track(2,isub-2)
      lat =track(2,isub-1)
      latb=track(2,isub  )
      ttop=toppar(track(1,isub-2)-track(1,isub-1),lata,
     .   0d0,lat,dtb,latb)
      if (ttop.ge.0d0 .and. ttop.le.dtb) then
*         write (6,'(i9,3f10.6)') itrknr+1,-dta,ttop,dtb
	 goto 49
      endif
*
* Check if track changes from descending to ascending or v.v.
* If the largest time gap was already before mesurement (isub-1),
* then (isub-2) is probably the last measurement on the track.
*
      if ((lata.lt.lat).neqv.(lat.lt.latb)) then
	 if (dta.gt.dtb) isub=isub-1
*        write (6,'(i9,3f10.3)') itrknr+1,lata/rad,lat/rad,latb/rad
	 goto 49
      endif

      goto 10
*
* Close current track, and reject last read point.
*
   49 isub=isub-1
   50 itrknr=itrknr+1
      if (itrknr.ge.maxtrk) go to 1150
      asc(itrknr)=(track(2,1).lt.track(2,isub))
      satid(itrknr)=sid(ik)
      if (isub.ge.npoint) then
         call inclin(ik,isub)
         good(itrknr)=.true.
      endif
      iobs(itrknr)=isub
      time(1,itrknr)=track(1,1)
      time(2,itrknr)=track(1,isub)
      filenr(itrknr)=ifilnr
*
* Determine whether track is short or long
*
      dist=(track(1,isub)-track(1,1))*6.5d0
      distmn=distmn+dist
      short(itrknr)=(dist.lt.shortl)
      if (short(itrknr)) ishrt=ishrt+1
*
* Dump array TRACK in XAF
*
      if (xaf) call wrixaf(wrunit,isub)
      nrxaf=nrxaf+isub
      istrec(itrknr)=nobsblk+1
      nobsblk=nobsblk+isub
      irec=irec+isub
      call statbar(1,irec,' ')
*     write (6,*) itrknr,isub,nint(track(1,1)),
*    .	nint(track(1,isub)-track(1,1)),track(3,1)/rad,
*    .  track(3,isub)/rad
      if (nobsblk.gt.maxobsblk-maxobs .or. irec.ge.nrrecs) then
         blkfil(nblocks)=ifilnr
	 blkobs(nblocks)=nobsblk
	 blkrec(nblocks+1)=irec+1
	 mtrkblk(nblocks)=itrknr	! last track number in block
	 nblocks=nblocks+1
	 if (nblocks.gt.mblocks) goto 1160
	 nobsblk=0
      endif
*
* If end of file, close it.
*
      if (irec.lt.nrrecs) goto 5
      ifilnr=ifilnr+1
      goto 4
*
* Start intersection search
*
 1000 continue
      nblocks=nblocks-1

* Give message on number of tracks and satellites

      write (*,610) nrxaf,itrknr,nblocks,
     |		nint(distmn/itrknr),nint(1d2*ishrt/itrknr),
     |		isatel,(name(sid(i)),i=1,isatel)
610   format(/
     |'   Number of measurements processed            -> ',i9/
     |'   Number of tracks encountered                -> ',i9/
     |'   Number of blocks created                    -> ',i9/
     |'   Mean track length                           -> ',i9,' km.'/
     |'   Percentage of short tracks                  -> ',i9,' %'/
     |'   Number of altimeters encountered            -> ',i9/
     |'   ... being',20(1x,a8))
*
* Write XAF header and close altimeter file
*
      if (xaf) then
	 filspf='@XAB'
         write (wrunit,rec=1) filspf,nrxaf,bound,0
         close (wrunit)
	 if (naux.gt.0) then
            close (wrunit+1)
      	    l=index(afile,' ')-1
	    open (wrunit+1,file=afile(:l)//'.aux',
     |		access='direct',form='unformatted',recl=12)
	    write (wrunit+1,rec=1) '@AUX',nrxaf,naux
            close (wrunit+1)
	 endif
      endif
      return

 1145 write(6,*) 'warning message : number of datapoints per pass is exc
     .eeding ',maxobs
      stop

 1150 write(6,*) 'warning message : number of passes is exceeding ',
     .maxtrk
      stop

 1160 write(6,*) 'warning message : number of blocks is exceeding ',
     .mblocks
      stop

      end
