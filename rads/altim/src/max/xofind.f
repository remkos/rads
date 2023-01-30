**XOFIND -- Find all crossover locations
*+
* If XXF has to be created, loop through all block-block and all track-track
* combinations. Guess and determine xover locations.
*-
* 16-Oct-1994 - Created, splitted off from SORTRK
* 20-Oct-1994 - Skip one out of NSKIP crossovers
* 16-Feb-1998 - Little big with dtxoff repaired
*-----------------------------------------------------------------------
      subroutine xofind

      include "maxcom.inc"
      integer i,l0,l,k0,k,ios,satk,satl,nrxtype,nrxskip
      real*8  dt
      character*4 filspf
*
* Initialization
*
      cosmax=(cos(angmin))**2
      nrget=0
*
* Open crossover file
*
      l=index(xfile,' ')-1
      if (specif(4:4).eq.'X') then
         open (wrunit,file=xfile,status='new',form='unformatted',
     |      access='direct',recl=48)
         filspf='@XXB'
      else if (specif(4:4).eq.'O') then
         open (wrunit,file=xfile,status='new',form='unformatted',
     |      access='direct',recl=44)
         filspf='@XXO'
      else if (specif(4:4).eq.'S') then
         open (wrunit,file=xfile,status='new',form='unformatted',
     |      access='direct',recl=36)
         filspf='@XXS'
      else if (specif(4:4).eq.'A') then
	 open (wrunit,file=xfile,status='new')
      endif
      if (naux.gt.0) then
	 open (wrunit+1,file=xfile(:l)//'.aux',
     |		access='direct',form='unformatted',recl=4*naux)
      endif
*
* Start intersection search
*
* Guess crossovers if both tracks are marked good.
* Check if time interval is within bounds.
* Check if only singles or only duals are requested.
* Opening of temporary file transported to insect.
*
      write (*,600)
      call statbar(0,itrknr*(itrknr-1)/2,'Track combinations')

      i=0
      nrxtype=0
      nrxskip=0
      if (nskip.eq.0) nskip=maxint4

      do l0=1,nblocks
	 do k0=l0,nblocks
	    do l=mtrkblk(l0-1)+1,mtrkblk(l0)
	       do k=max(l+1,mtrkblk(k0-1)+1),mtrkblk(k0)
		  i=i+1
		  satk=satid(k)
		  satl=satid(l)
		  dt=abs(time(3,k)-time(3,l)+dtxoff(satk,satl))
		  if (.not.good(k).or..not.good(l)) then
		  else if (dt.gt.dtxo(satk,satl)+3600d0) then
		     nrxdt=nrxdt+1
		  else if (nomerge .and. filenr(k).ne.filenr(l)) then
		     nrxdt=nrxdt+1
		  else if (specif(1:1).eq.'D' .and. satk.eq.satl) then
		     nrxtype=nrxtype+1
		  else if (specif(1:1).eq.'S' .and. satk.ne.satl) then
		     nrxtype=nrxtype+1
		  else if (mod(i,nskip).eq.0) then
		     nrxskip=nrxskip+1
		  else
		     call blload(1,l0)
		     call blload(2,k0)
		     call insect(l,k)
		  endif
		  call statbar(1,i,' ')
	       enddo
	    enddo
	 enddo
      enddo
*
* WRITE XXF header and close crossover file
*
      if (specif(4:4).ne.'A') write (wrunit,rec=1) filspf,nrxfnd,
     |		bound,0
      close (wrunit)
      if (naux.gt.0) then
         close (wrunit+1)
         l=index(xfile,' ')-1
	 open (wrunit+1,file=xfile(:l)//'.aux',
     |		access='direct',form='unformatted',recl=12)
	 write (wrunit+1,rec=1) '@AUX',nrxfnd,naux*2
         close (wrunit+1)
      endif

      xrms=dsqrt(xrms/nrxfnd)
      write (*,601,iostat=ios) i,nrxdt,nrxtype,nrxskip,
     |   nrxins,nrxinc,nrxdat,nrxiter,nrxfnd,
     |   dble(nrget)/nblocks,sqrt(dtrms/2/nrxfnd),dtmax,xrms*1d2
      if (ios.ne.0) write (*,*) 'Error writing crossover :',ios

600   format(//'Guessing and determining crossovers ...')
601   format(/
     |'   Track combinations                          -> ',i9//
     |'   - Time interval out of limit                -> ',i9/
     |'   - Wrong xover type (dual/single)            -> ',i9/
     |'   - Skipped                                   -> ',i9/
     |'   - No intersection                           -> ',i9/
     |'   - Shallow angle                             -> ',i9/
     |'   - Too few data around                       -> ',i9/
     |'   - Too many iterations required              -> ',i9//
     |'   Nr of xovers found and processed            -> ',i9/
     |'   Average number of track readings            -> ',f9.3/
     |'   RMS crossover time correction on guess [s]  -> ',f9.3/
     |'   Max. crossover time correction on guess [s] -> ',f9.3/
     |'   Crossover difference RMS [cm]               -> ',f9.3)
      return
      end
