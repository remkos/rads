*XOVDET -- Determine xover location
*
* XOVDET determines the real crossover location from the estimate
* (indices a, en b represent the two tracks involved)
* (indices 3, en 4 represent the two nearest points to the crossover)
*-
* 16-Oct-1994 - Old binary formats removed. Do XOVDET for each guess.
* 21-Oct-1994 - Faster (?)
* 13-Sep-1996 - Adjust rlonx to always lay between monlon and maxlon
* 16-Feb-1998 - Bug fixed with regard to dtxoff.
* 10-May-1999 - Changed order of tracks in case both descending
*  4-Jun-2002 - Linear interpolation of auxiliary fields, instead of
*               plain average of two epochs. Delayed time interval checks
*               to final computed xover.
*-----------------------------------------------------------------------
      subroutine xovdet(ja,jb,xdecl,xlam)
      implicit none
      include "maxcom.inc"
      real*8 small,xdecl,xlam,xlam0,xlam1,arglat
      integer*4 maxitr,sata,satb

      integer*4 ireca,ireca0,ireca1,ja
      integer*4 irecb,irecb0,irecb1,jb
      real*8 utcxa,orbhxa,seahxa,coefa1,coefa2,timxa,ta0,ta1
      real*8 utcxb,orbhxb,seahxb,coefb1,coefb2,timxb,tb0,tb1

      parameter  (small=1d-15,maxitr=10)
      real*8     rdata(5,6),rdatb(5,6)
      integer*2  adata(maux,6),adatb(maux,6)
      logical    nexta,nextb,firsta,intorb,blnew

      real*8 rlatx,rlonx,dt
      integer*4 i,itera

      integer*2 auxa(maux),auxb(maux)

      intorb=(specif(4:4).eq.'O' .or. specif(4:4).eq.'S')
*     write (*,*) "tracks,xdecl,xlam:",ja,jb,xdecl/rad,xlam/rad
*
* Determine guessed time of crossover for each track.
* If guessed time is far beyond track limits, forget about it.
*
      ta0=time(1,ja)
      ta1=time(2,ja)
      timxa=time(3,ja)+arglat(xdecl,inc(ja))/theta(ja)
      if (timxa.lt.ta0-timper .or. timxa.gt.ta1+timper) goto 140

      tb0=time(1,jb)
      tb1=time(2,jb)
      timxb=time(3,jb)+arglat(xdecl,inc(jb))/theta(jb)
      if (timxb.lt.tb0-timper .or. timxb.gt.tb1+timper) goto 140

* Reject xover if time interval exceeds dtxo by some seconds

      sata=satid(ja)
      satb=satid(jb)
      dt=abs(timxa-timxb+dtxoff(sata,satb))
      if (dt.gt.dtxo(sata,satb)+timper) goto 120

      xlam0=xlam-pi
      xlam1=xlam+pi

      ireca0=istrec(ja)
      ireca1=ireca0+iobs(ja)-1
      irecb0=istrec(jb)
      irecb1=irecb0+iobs(jb)-1
*
* Start search for two points nearest to crossover guess epoch
* for track JA
*
      call bltween(1,rdata(1,3),adata(1,3),
     |		ireca,ireca0,ireca1,ta0,ta1,timxa)
      call grouplon(xlam0,xlam1,rdata(3,3))
      call grouplon(xlam0,xlam1,rdata(3,4))
*
* If guessed crossover is far from any data, just forget about it.
*
      if (timxa-rdata(1,3).gt.timper .and. rdata(1,4)-timxa.gt.timper)
     .     goto 100
*
* Start search for two points nearest to crossover guess epoch
* for track JB
*
      call bltween(2,rdatb(1,3),adatb(1,3),
     |		irecb,irecb0,irecb1,tb0,tb1,timxb)
      call grouplon(xlam0,xlam1,rdatb(3,3))
      call grouplon(xlam0,xlam1,rdatb(3,4))
*
* If guessed crossover is far from any data, just forget about it.
*
      if (timxb-rdatb(1,3).gt.timper .and. rdatb(1,4)-timxb.gt.timper)
     .     goto 100
*
* Start computing the crossover location in geodetic coordinates
*
* Keep track of the number of iterations: problems may occur at
* very high latitudes
*
      itera=-1
    2 itera=itera+1
      if (itera.ge.maxitr) goto 130
      coefa1=(rdata(2,3)-rdata(2,4))/(rdata(3,3)-rdata(3,4))
      coefa2=(rdata(2,3)*rdata(3,4)-rdata(2,4)*rdata(3,3))/
     .       (rdata(3,4)-rdata(3,3))
      coefb1=(rdatb(2,3)-rdatb(2,4))/(rdatb(3,3)-rdatb(3,4))
      coefb2=(rdatb(2,3)*rdatb(3,4)-rdatb(2,4)*rdatb(3,3))/
     .       (rdatb(3,4)-rdatb(3,3))
      if (abs(coefa1-coefb1).lt.small) goto 110
      rlonx=(coefb2-coefa2)/(coefa1-coefb1)
      rlatx=coefa1*rlonx+coefa2
*
* Check whether computed crossover point is still between the two
* datapoints
*
      nexta=(rlatx.lt.min(rdata(2,3),rdata(2,4)) .or. rlatx.gt.
     .       max(rdata(2,3),rdata(2,4)))
      nextb=(rlatx.lt.min(rdatb(2,3),rdatb(2,4)) .or. rlatx.gt.
     .       max(rdatb(2,3),rdatb(2,4)))
      if (nexta) then
         if (blnew(1,asc(ja),rdata(1,3),adata(1,3),
     .		ireca,ireca0,ireca1,rlatx)) goto 140
         call grouplon(xlam0,xlam1,rdata(3,3))
         call grouplon(xlam0,xlam1,rdata(3,4))
      endif
      if (nextb) then
         if (blnew(2,asc(jb),rdatb(1,3),adatb(1,3),
     .		irecb,irecb0,irecb1,rlatx)) goto 140
         call grouplon(xlam0,xlam1,rdatb(3,3))
         call grouplon(xlam0,xlam1,rdatb(3,4))
      endif
      if (nexta.or.nextb) goto 2
*
* Check whether track boundaries are not exceeded and read 2 more
* points on either side of the crossover (a total of six points
* will be used for creating the quasi observation)
*
      if (ireca+3.gt.ireca1 .or. ireca-2.lt.ireca0 .or.
     .       irecb+3.gt.irecb1 .or. irecb-2.lt.irecb0) goto 100
*
* Read additional points
*
      call blread(1,ireca-2,rdata(1,1),adata(1,1))
      call blread(1,ireca-1,rdata(1,2),adata(1,2))
      call blread(1,ireca+2,rdata(1,5),adata(1,5))
      call blread(1,ireca+3,rdata(1,6),adata(1,6))
      call blread(2,irecb-2,rdatb(1,1),adatb(1,1))
      call blread(2,irecb-1,rdatb(1,2),adatb(1,2))
      call blread(2,irecb+2,rdatb(1,5),adatb(1,5))
      call blread(2,irecb+3,rdatb(1,6),adatb(1,6))
*
* Check time interval between the two points furthest from crossover
*
      if (rdata(1,6)-rdata(1,1).gt.timper) goto 100
      if (rdatb(1,6)-rdatb(1,1).gt.timper) goto 100
*
* Check time interval between the two closest points
*
      if (rdata(1,4)-rdata(1,3).gt.timper/4) goto 100
      if (rdatb(1,4)-rdatb(1,3).gt.timper/4) goto 100
*
* Convert longitude to fit within minimum and maximum longitude
*
      call grouplon(minlon,maxlon,rlonx)
*
* Determine pseudo observation at crossover location
*
      call quasob (ja,rdata(1,1),rlatx,utcxa,seahxa,intorb,orbhxa)
      call quasob (jb,rdatb(1,1),rlatx,utcxb,seahxb,intorb,orbhxb)

* Finally check the actual time interval and see if it exceeds dtxo

      dt=abs(utcxa-utcxb+dtxoff(sata,satb))
      if (dt.gt.dtxo(sata,satb)) goto 120

* Interpolate auxiliary data

      dt=(utcxa-rdata(1,3))/(rdata(1,4)-rdata(1,3))
      do i=1,naux
	 auxa(i)=nint((1d0-dt)*adata(i,3)+dt*adata(i,4))
      enddo
      dt=(utcxb-rdatb(1,3))/(rdatb(1,4)-rdatb(1,3))
      do i=1,naux
	 auxb(i)=nint((1d0-dt)*adatb(i,3)+dt*adatb(i,4))
      enddo
*
* Keep statistics on time difference with crossover guess
*
      dtmax=max(dtmax,abs(timxa-utcxa),abs(timxb-utcxb))
      dtrms=dtrms+(timxa-utcxa)**2+(timxb-utcxb)**2
*
* Write record to xover file
*
      nrxfnd=nrxfnd+1
*
* Put ascending track first.
* If both are ascending the lowest track number comes first.
* If both are descending the lowest track number comes last.
*
      if (asc(ja).and.asc(jb)) then
         firsta=(ja.lt.jb)
      else if (.not.asc(ja).and..not.asc(jb)) then
         firsta=(ja.gt.jb)
      else
         firsta=asc(ja)
      endif

      if (firsta) then
         call wrixxf(rlatx,rlonx,ja,jb,utcxa,utcxb,
     .		seahxa,seahxb,orbhxa,orbhxb,auxa,auxb)
      else
         call wrixxf(rlatx,rlonx,jb,ja,utcxb,utcxa,
     .		seahxb,seahxa,orbhxb,orbhxa,auxb,auxa)
      endif
*
* Determine sum of squared crossover residuals to determine the RMS
* of the crossover differences
*
      ixpt(ja)=ixpt(ja)+1
      ixpt(jb)=ixpt(jb)+1
      xrms=xrms+(seahxa-seahxb)**2

      return

100   continue
      nrxdat=nrxdat+1
      return
110   continue
      nrxinc=nrxinc+1
      return
120   continue
      nrxdt=nrxdt+1
      return
130   continue
      nrxiter=nrxiter+1
      return
140   continue
      nrxins=nrxins+1
      return
      end
