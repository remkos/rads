**XSTAT -- determine crossover statistics
*
      program xstat
*
* XSTAT computes the mean, rms, and std dev of crossovers stored in one or
* several files. Optionally, edit levels can be defined, or datasets can be
* sampled, selected on area, or deep/shallow waters.
*-
* $Log: xstat.f,v $
* Revision 1.14  2012/05/24 12:52:57  rads
* - Higher precision
*
* Revision 1.13  2011/03/08 01:20:52  rads
* - Start bins consistently
*
* Revision 1.12  2011/02/10 18:40:24  rads
* - Widened columns
* - Allow up to 5000000 xovers
*
* (c) Remko Scharroo - TU Delft, Altimetrics
*----------------------------------------------------------------------
      implicit none

* General

      integer maxfiles
      parameter (maxfiles=200)
      character filenm(maxfiles)*80,spec*4,arg*80,text*80,table*80/' '/,
     |	title*80/'-'/
      integer*4 i,j,j1,j2,nrec,n1,n2,npar,iarg,iargc,lnblnk
      integer*4 t_n1/0/,t_n2/0/,l,ftype/0/,recl
      integer*4 dt,dorb/0/,dtrk,if/0/,nf/0/,sort/0/,chbias(2)/2*0/
      real*8 avg1,t_avg1/0d0/,rms1,sgm1,t_rms1/0d0/,t_sgm1,
     |	tbias(2)/2*0d0/,ssb(2)/2*0d0/,dssh,cor0,cor1,cor2,
     |  avg2,t_avg2/0d0/,rms2,sgm2,t_rms2/0d0/,t_sgm2,ssvar,
     |	dtmin/0d0/,dtmax/1d35/,edit/3.5d0/,xmax/1d35/,bias(2)/2*0d0/
      real*8 rlon,rlat
      real*8 latmin/-1e35/,latmax/1e35/,lonmin/-1e35/,lonmax/1e35/
      logical ocean/.false./,sland,shallow/.false./,datearg,swap/.false./
      real*8 smin/0d0/,smax/1d35/,slo/1d35/,shi/-1d35/,dum

* For timing error model

      real*8	tppar(3)/15d3,0.048d0,-2d3/
      logical	tpalt/.false./

* File handling

      integer*4 openf,readf,closef,ios,fd/-1/,fdaux/-1/,fddlt/-1/
      integer*4 xxf4(16),txo(2)
      integer*4 maux,naux/0/,fldn/0/,fld/-1/
      parameter (maux=20)
      integer*2 xxf2(32),aux(maux)/maux*0/,dlt(2)/2*0/
      equivalence (xxf2,xxf4)

* Storage

      integer mxover
      parameter (mxover=4000000)
      real*4 xover(mxover)

* For variability

      character*80 varnm
      integer nx,ny
      real*4 x0,x1,y0,y1,z0,z1,fx/1d0/,fy/1d0/,vmax/1e10/
      integer maxgrd
      parameter (maxgrd=360*180)
      real*4 var(maxgrd),v(mxover)
      logical novar/.false./,dovar/.false./

* For plotting

      real t(mxover,2)
      real*8 tmin,tmax,hmin,hmax,day,bin/1d0/
      parameter (day=86400d0)
      integer yref,tref,mdate,itmin,itmax
      character*80 dev/'-'/
      logical pgopen/.false./,nodots/.false./,norms/.false./,
     |	dtplot/.false./,combi/.false./,noarcnr/.false./,
     |	norejected/.false./,years/.false./
      integer it,nt/0/,maxt,lw/1/,pgbeg
      parameter (maxt=10000)
      real*8 trms(maxt),tmean(maxt),trej(maxt),rmax/18/
      real xp0,xp1,yp0,yp1
      integer tnr(maxt)

* Scan argument line

      do iarg=1,iargc()
      call getarg(iarg,arg)
      if (arg(:1).ge.'A' .and. arg(:1).le.'Z' .and.
     |		datearg(arg,slo,shi,dum)) then
      else if (datearg(arg,smin,smax,dum)) then   
      else if (arg(:4).eq.'bin=') then
	 read(arg(5:),*) bin
      else if (arg(:3).eq.'dt=') then
	 dtmax=0d0
	 read (arg(4:),*,iostat=ios) dtmin,dtmax
	 if (dtmax.eq.0d0) then
	    dtmax=dtmin
	    dtmin=0d0
	 endif
	 dtmin=dtmin*86400
	 dtmax=dtmax*86400
      else if (arg(:3).eq.'lw=') then
	 read (arg(4:),*) lw
      else if (arg(:4).eq.'tab=') then
	 table=arg(5:)
	 open (10,file=table)
	 write (10,550) '# Date Mean(cm) Std(cm) Rej(%) Nr'
      else if (arg(:4).eq.'dev=') then
	 dev=arg(5:)
      else if (arg(:4).eq.'lon=') then
	 read (arg(5:),*) lonmin,lonmax
      else if (arg(:4).eq.'lat=') then
	 read (arg(5:),*) latmin,latmax
      else if (arg(:5).eq.'edit=') then
	 read (arg(6:),*) edit
      else if (arg(1:7).eq.'chbias=') then
	 chbias(2)=9999
	 read (arg(8:),*,iostat=ios) chbias
	 if (chbias(2).eq.9999) chbias(2)=chbias(1)
      else if (arg(:2).eq.'-S') then
      	 swap=.true.
      else if (arg(:2).eq.'-v') then
	 read (arg(3:),*,iostat=ios) vmax
	 vmax=vmax/100
	 novar=.true.
      else if (arg(:2).eq.'+v') then
	 read (arg(3:),*,iostat=ios) vmax
	 vmax=vmax/100
	 dovar=.true.
      else if (arg(:2).eq.'-c') then
	 combi=.true.
      else if (arg(:2).eq.'-n') then
	 nodots=.true.
      else if (arg(:2).eq.'-N') then
	 nodots=.true.
	 norms=.true.
      else if (arg(:2).eq.'-y') then
	 years=.true.
      else if (arg(:2).eq.'-r') then
	 norejected=.true.
      else if (arg(:2).eq.'-f') then
	 noarcnr=.true.
      else if (arg(:5).eq.'rmax=') then
	 read (arg(6:),*) rmax
      else if (arg(:5).eq.'xmax=') then
	 read (arg(6:),*) xmax
      else if (arg(:6).eq.'-tpalt') then
         tpalt=.true.
      else if (arg(:2).eq.'-o') then
	 write (*,550) '(Ocean only)'
	 ocean=.true.
	 shallow=.false.
      else if (arg(:4).eq.'fld=') then
	 read (arg(5:),*) fld
         write (*,550) '(auxiliary field statistics) '
	 fld=(fld-1)*2
	 fdaux=0
      else if (arg(:4).eq.'-swh') then
         write (*,550) '(SWH statistics)'
	 fld=0
	 fdaux=0
      else if (arg(:2).eq.'-s') then
	 write (*,550) '(Shallow waters)'
	 shallow=.true.
	 ocean=.false.
      else if (arg(:2).eq.'-a' .or. arg(:2).eq.'-e') then
	 write (*,550) '(asc - des, ers - t/p)'
	 sort=0
      else if (arg(:2).eq.'-t') then
	 write (*,550) '(first - last)'
	 sort=+1
      else if (arg(:2).eq.'+t') then
	 write (*,550) '(last - first)'
	 sort=-1
      else if (arg(:2).eq.'-T') then
	 write (*,550) '(sat A - sat B)'
	 sort=+2
      else if (arg(:2).eq.'+T') then
	 write (*,550) '(sat B - sat A)'
	 sort=-2
      else if (arg(:3).eq.'-dt') then
	 dtplot=.true.
      else if (arg(:5).eq.'bias=') then
	 bias(2)=1d30
	 read (arg(6:),*,iostat=ios) bias
	 if (bias(2).gt.1d20) bias(2)=bias(1)
      else if (arg(:6).eq.'tbias=') then
	 tbias(2)=1d30
	 read (arg(7:),*,iostat=ios) tbias
	 if (tbias(2).gt.1d20) tbias(2)=tbias(1)
* tbias in milliseconds. Positive when tags are late.
* Will be multiplied by orbit rate in mm/s to give micrometres
	 fdaux=0
	 nf=nf+1
	 filenm(nf)=arg
      else if (arg(:4).eq.'ssb=') then
	 ssb(2)=1d30
	 read (arg(5:),*,iostat=ios) ssb
	 if (ssb(2).gt.1d20) ssb(2)=ssb(1)
* Addition ssb in percentage. Positive when ssb should be higher.
	 ssb(1)=ssb(1)*1d1
	 ssb(2)=ssb(2)*1d1
* ssb in permille will be multiplied by SWH in mm to give micrometres
	 fdaux=0
      else if (arg(:6).eq.'title=') then
	 title=arg(7:)
      else
	 nf=nf+1
	 if (nf.gt.maxfiles) call fin("xstat: too many files")
         filenm(nf)=arg
      endif
      enddo
      if (nf.eq.0) goto 1300

* Intialise periods

      if (slo.lt.shi) then
         smin=max(smin,slo)
         smax=min(smax,shi)
      endif
      tmin=slo/bin/day
      tmax=shi/bin/day

* Initialise statistics

      do i=1,maxt
         trej(i)=0
	 trms(i)=0
	 tmean(i)=0
	 tnr(i)=0
      enddo

* If variability grid is requested: load grid

      if (novar .or. dovar) then
         nx=0
         ny=maxgrd
         varnm='/user/altim'
         call checkenv('ALTIM',varnm,l)
         varnm(l+1:)='/data/world/ssvar95a.grd'
         call gridrd4(varnm,nx,ny,var,x0,x1,y0,y1,z0,z1)
         fx=(x1-x0)/(nx-1)
         fy=(y1-y0)/(ny-1)
      endif

* Process all xover files

      write (*,602)

* Reset counters

10    continue
      n1=0

* Open xover file and  determine file type
* (XXO, XXS, XXE, XXF, XXB, and D_H supported)

11    if=if+1
      arg=filenm(if)
      if (filenm(if)(:6).eq.'tbias=') then
	 tbias(2)=1d30
	 read (arg(7:),*,iostat=ios) tbias
	 if (tbias(2).gt.1d20) tbias(2)=tbias(1)
	 goto 11
      else if (arg(:5).eq.'bias=') then
	 bias(2)=1d30
	 read (arg(6:),*,iostat=ios) bias
	 if (bias(2).gt.1d20) bias(2)=bias(1)
	 goto 11
      endif
      fddlt=-1
12    fd=openf(arg,'r')
      if (fd.lt.0) goto 90
      ios=readf(fd,4,spec)
      if (spec.eq.'@D_H') then
	 fddlt=fd
	 call seekf(fddlt,24,0)
	 l=index(arg,'.xxo')
	 if (l.le.0) call fin("Error in DLT filename")
	 arg=arg(:l+3)
	 goto 12
      endif
      if (ios.ne.4) goto 90
      ios=readf(fd,4,nrec)
      if (swap) call i4swap(1,nrec)
      if (ios.ne.4 .or. nrec.eq.0) goto 90
      ios=readf(fd,4,npar)
      if (swap) call i4swap(1,npar)
      if (ios.ne.4) goto 90
      if (spec.eq.'@XXO') then
	 recl=44
         ftype=1
      else if (spec.eq.'@XXS') then
	 recl=36
	 ftype=2
      else if (spec.eq.'@XXE' .or. spec.eq.'@XXF'
     |		.or. spec.eq.'@XXB') then
	 recl=48
	 if (npar.eq.5) recl=64
	 ftype=3
      else
	 call fin('illegal filetype')
      endif
      call seekf(fd,recl,0)

* If needed, open auxiliary file

      l=index(arg,' ')-1
      if (fdaux.ge.0) then
         fdaux=openf(arg(:l)//'.aux','r')
	 if (fdaux.lt.0) goto 90
	 ios=readf(fdaux,4,spec)
	 if (spec.ne.'@AUX') call fin ('illegal file type on aux file')
	 ios=readf(fdaux,4,i)
	 if (swap) call i4swap(1,i)
	 if (i.ne.nrec) call fin('aux file does not match main file')
	 ios=readf(fdaux,4,naux)
	 if (swap) call i4swap(1,naux)
	 if (naux.gt.maux) call fin('too many parameters in aux file')
	 call seekf(fdaux,naux*2,0)
	 fldn=naux-2
      endif

* Read all xovers and store height difference in xover()

      do i=1,nrec
	 ios=readf(fd,recl,xxf4)
	 if (fdaux.ge.0) then
	    ios=readf(fdaux,naux*2,aux)
	    if (swap) call i2swap(naux*2,aux)
	 endif

	 if (ftype.eq.1) then
	    if (swap) then
	    	call i4swap(6,xxf4)
		call i2swap(2,xxf2(13))
		call i4swap(4,xxf4(8))
	    endif
	    txo(1)=xxf4(3)
	    txo(2)=xxf4(5)
	    dtrk=xxf2(13)-xxf2(14)
	    dssh=xxf4(8)-xxf4(9)
	    dorb=xxf4(10)-xxf4(11)
	 else if (ftype.eq.2) then
	    if (swap) then
	    	call i4swap(6,xxf4)
		call i2swap(2,xxf2(13))
		call i4swap(2,xxf4(8))
	    endif
	    txo(1)=xxf4(3)
	    txo(2)=xxf4(5)
	    dtrk=xxf2(13)-xxf2(14)
	    dssh=xxf4(8)-xxf4(9)
	 else
	    if (swap) then
	    	call i4swap(4,xxf4)
		call i2swap(2,xxf2(9))
		call i4swap(4,xxf4(6))
	    endif
	    txo(1)=xxf4(3)
	    txo(2)=xxf4(4)
	    dtrk=xxf2(9)-xxf2(10)
	    dssh=xxf4(8)-xxf4(9)
         endif
	 if (fddlt.ge.0) then
	    ios=readf(fddlt,4,dlt)
	    if (swap) call i2swap(2,dlt)
	    if (dlt(1).eq.-32768 .or. dlt(2).eq.-32768) goto 20
	    dssh=dssh+(dlt(1)-dlt(2))*1d3
	 endif

	 dt=abs(txo(1)-txo(2))
	 rlat=xxf4(1)/1e6
	 rlon=xxf4(2)/1e6

* Model T/P altitude rate if requested

	 if (tpalt) then
	    aux(fldn+1)=
     |		tppar(1)*sin(tppar(2)*rlat)+tppar(3)*sin(2*tppar(2)*rlat)
	    aux(fldn+2)=-aux(fldn+1)
	 endif

* Determine whether we need (1)-(2) or (2)-(1)

	 if ((sort.eq.+1 .and. txo(1).gt.txo(2)) .or.
     |	     (sort.eq.-1 .and. txo(1).lt.txo(2)) .or.
     |	     (sort.eq.+2 .and. dtrk.gt.0) .or.
     |	     (sort.eq.-2 .and. dtrk.lt.0) .or.
     |	     (sort.eq.0  .and. dorb.gt.100 000 000)) then
	    j1=2
	    j2=1
	    dssh=-dssh
	 else
	    j1=1
	    j2=2
	 endif

	 if (fld.ge.0) then
* Replace SSH difference with auxiliary field if requested
	    dssh=(aux(j1+fld)-aux(j2+fld))*1000
	 else
* Apply timing bias and/or sea state bias correction
	    dssh=dssh+(ssb(1)*aux(j1)-ssb(2)*aux(j2))-
     |		 (tbias(1)*aux(j1+fldn)-tbias(2)*aux(j2+fldn))
	 endif

* Add bias to SSH when requested

	 dssh=dssh+(bias(1)-bias(2))*1d4

* Change SPTR bias if requested from new to old

	 if (chbias(1).ne.0) then
	    call altbias(+chbias(1),txo(j1),cor0,cor1,cor2)
	    dssh=dssh-(cor0-cor1-cor2)*1d6
	    call altbias(-chbias(1),txo(j1),cor0,cor1,cor2)
	    dssh=dssh+(cor0-cor1-cor2)*1d6
         endif
	 if (chbias(2).ne.0) then
	    call altbias(+chbias(2),txo(j2),cor0,cor1,cor2)
	    dssh=dssh+(cor0-cor1-cor2)*1d6
	    call altbias(-chbias(2),txo(j2),cor0,cor1,cor2)
	    dssh=dssh-(cor0-cor1-cor2)*1d6
         endif

* Do data screening

	 if (dt.lt.dtmin .or. dt.gt.dtmax) then
	 else if (rlat.lt.latmin .or. rlat.gt.latmax) then
	 else if (rlon.lt.lonmin .or. rlon.gt.lonmax) then
	 else if (ocean .and. sland(rlat,rlon)) then
	 else if (shallow .and.	.not.sland(rlat,rlon)) then
	 else if (txo(1).lt.smin .or. txo(1).gt.smax) then
	 else if (txo(2).lt.smin .or. txo(2).gt.smax) then
	 else
            if (n1.gt.mxover) goto 1310
	    if (novar) then
	       l=nint((rlat-y0)/fy)*nx+nint((rlon-x0)/fx)+1
	       if (var(l).gt.vmax) goto 20
	       n1=n1+1
	       v(n1)=var(l)*100
	    else if (dovar) then
	       l=nint((rlat-y0)/fy)*nx+nint((rlon-x0)/fx)+1
	       if (var(l).lt.vmax) goto 20
	       n1=n1+1
	       v(n1)=var(l)*100
	    else
	       n1=n1+1
	       v(n1)=0
	    endif
	    xover(n1)=dssh*1d-4
	    if (dtplot) then
	    ! If time difference is to be plotted, store difference
	       t(n1,1)=dt/bin/day
	       t(n1,2)=dt/bin/day
	    else
	       t(n1,1)=txo(j1)/bin/day
	       t(n1,2)=txo(j2)/bin/day
	    endif
	    tmin=min(tmin,t(n1,1),t(n1,2))
	    tmax=max(tmax,t(n1,1),t(n1,2))
	    ! Do not use second time tag when sorting or T/P
	    if (sort.ne.0 .or. abs(dorb).gt.100000000) t(n1,2)=1e30
	 endif
20       continue
      enddo
90    continue
      if (fd.gt.0) ios=closef(fd)
      if (fdaux.gt.0) ios=closef(fdaux)
      if (fddlt.gt.0) ios=closef(fddlt)
      if (combi .and. if.ne.nf) goto 11

* Compute statistics and then edit.

      itmin=int(tmin)
      itmax=int(tmax+1)
      call meanrms(n1,xover,avg1,rms1)
      avg2=avg1
      rms2=sqrt(rms1**2-avg1**2)
      n2=n1
      call editor(n2,xover,v,avg2,rms2,ssvar,edit,xmax)

* Update overall statistics

      t_n1=t_n1+n1
      t_avg1=t_avg1+avg1*n1
      t_rms1=t_rms1+rms1*rms1*n1
      t_n2=t_n2+n2
      t_avg2=t_avg2+avg2*n2
      t_rms2=t_rms2+rms2*rms2*n2
      sgm1=sqrt(rms1**2-avg1**2)
      sgm2=sqrt(rms2**2-avg2**2)
      hmin=avg2-min(xmax,sgm2*edit)
      hmax=avg2+min(xmax,sgm2*edit)
      if (combi) then
	 arg='* combi *'
      else
         arg=filenm(if)
      endif
      l=lnblnk(arg)

      if (dtplot) then

* Generate xover RMS and mean statistics as function of time difference

         nt=itmax
         if (nt.gt.maxt) call fin("xstat: too many days")
         do i=1,nt
	    tnr(i)=0
            trms(i)=0
	    tmean(i)=0
         enddo
         do i=1,n1
	    it=int(t(i,1)+1)
	    if (xover(i).lt.hmin .or. xover(i).gt.hmax) then
	       trej(it)=trej(it)+1
            else
	       trms(it)=trms(it)+xover(i)**2
	       tmean(it)=tmean(it)+xover(i)
	       tnr(it)=tnr(it)+1
	    endif
         enddo
	 tmin=0
	 tmax=nt
	 tref=0
	 itmin=0
	 itmax=nt

      else

* Generate daily RMS and mean statistics
* First determine reference day ( = beginning of the year )

         yref=mdate(1,int(tmin*bin+46066))/10000
         tref=mdate(2,yref*10000+0101)-1-46066
	 yref=yref+1900
	 if (yref.lt.1950) yref=yref+100
   
         nt=itmax-itmin
         if (nt.gt.maxt) call fin("xstat: too many days")
         do i=1,nt
	    tnr(i)=0
            trms(i)=0
	    tmean(i)=0
         enddo
         do i=1,n1
	    do j=1,2
               if (t(i,j).lt.1e20) then
	          it=int(t(i,j)-itmin+1)
	          if (xover(i).lt.hmin .or. xover(i).gt.hmax) then
	             trej(it)=trej(it)+1
                  else
	             trms(it)=trms(it)+xover(i)**2
	             tmean(it)=tmean(it)+xover(i)
	             tnr(it)=tnr(it)+1
	          endif
	       endif
	    enddo
         enddo

      endif

* Continue with daily statistics

      do i=1,nt
	 if (tnr(i).gt.1) then
	    tmean(i)=tmean(i)/tnr(i)
	    trms(i)=sqrt(trms(i)/tnr(i)-tmean(i)**2)
	    trej(i)=100.*trej(i)/(tnr(i)+trej(i))
	    if (table.eq.' ') then
	    else if (dtplot) then
	       write (10,'(f6.2,2f8.2,f7.2,i8)')
     |		(itmin+i-1)*bin,tmean(i),trms(i),trej(i),tnr(i)
	    else
	       write (10,'(i6.6,2f8.2,f7.2,i8)')
     |		mdate(1,int(46066+(itmin+i-1)*bin)),
     |		tmean(i),trms(i),trej(i),tnr(i)
	    endif
	 else
	    trms(i)=1e30
	    tmean(i)=1e30
	    trej(i)=1e30
	 endif
      enddo

* If plotting is requested: plot xover height differences

      if (dev.ne.'-') then
      if (norms) then
	 hmin=-5
	 hmax=5
      else if (nodots) then
	 hmin=-rmax/2
	 hmax=rmax
      endif
      if (.not.pgopen) then
         if (pgbeg(0,dev,1,1).ne.1) 
     |		call fin("Error opening plot device")
         call pgscr(0,0.8,0.8,0.8)
         call pgscr(1,0.0,0.0,0.0)
         call pgscr(2,1.0,0.0,0.0)
         call pgscr(3,0.0,0.5,0.0)
	 pgopen=.true.
      else
         call pgpage
      endif
      call pgslw(1)
      if (.not.dtplot .and. .not.noarcnr) then
         call pgsci(2)
         call pgsch(0.7)
         call pgswin(real((tmin*bin-2381)/3.5d0),real((tmax*bin-2381)/3.5d0),
     |		0.0,1.0)
         call pgbox('csti',1.0,0,' ',0.0,0)
         if (tmin*bin.gt.3761) then
            call pgswin(real((tmin*bin-2381)/3.5d0-394.5d0),
     |      real((tmax*bin-2381)/3.5d0-394.5d0),0.0,1.0)
         else
            call pgswin(real((tmin*bin-2381)/3.5d0-0.5d0),
     |		real((tmax*bin-2381)/3.5d0-0.5d0),0.0,1.0)
         endif
         call pgbox('m',0.0,0,' ',0.0,0)
      endif
      call pgsci(1)
      call pgsch(1.0)
      call pgswin(real(tmin*bin-tref),real(tmax*bin-tref),0.0,1.0)
      call pgqvp(0,xp0,xp1,yp0,yp1)
      if (years) then
	 call pgswin(real(tmin*bin/365.25d0),real(tmax*bin/365.25d0),0.0,1.0)
	 call pgbox('cst',1.0,12,' ',0.0,0)
	 call pgbox('bsti',1.0,12,' ',0.0,0)
	 call pgswin(real(tmin*bin/365.25d0+1984.5d0),
     |		real(tmax*bin/365.25d0+1984.5d0),0.0,1.0)
	 call pgbox('n',1.0,0,' ',0.0,0)
         text='year'
      else
         call pgbox('cst',0.0,0,' ',0.0,0)
         call pgbox('bnsti',0.0,0,' ',0.0,0)
         text='time difference (days)'
         if (.not.dtplot) write (text,700) yref
700      format ('time (days of the year ',i4,')')
      endif
      if (norejected) then
         call pgswin(real(tmin),real(tmax),real(hmin),real(hmax))
         call pgbox(' ',0.0,0,'cmsti',5.0,5)
      else
         call pgsvp(xp0,xp1,yp0,yp0+(yp1-yp0)/5)
         call pgswin(real(tmin*bin-tref),real(tmax*bin-tref),0.0,10.0)
         call pgmtxt('R',2.5,0.5,0.5,'rejected (%)')
         call pgbox(' ',0.0,0,'cmsti',2.0,4)
         call pgsvp(xp0,xp1,yp0+(yp1-yp0)/4,yp1)
         call pgswin(real(tmin*bin-tref),real(tmax*bin-tref),
     |		real(hmin+(hmax-hmin)/4),real(hmax))
         call pgbox(' ',0.0,0,'cmsti',5.0,5)
         call pgsvp(xp0,xp1,yp0,yp1)
      endif
      call pgswin(real(tmin),real(tmax),real(hmin),real(hmax))
      call pgbox('a',0.0,0,'bnsti',5.0,5)
      if (nodots) then
         call pgmtxt('R',2.5,0.33,0.5,'mean')
         call pgmtxt('R',2.5,0.70,0.5,'std dev')
      else
         call pgmtxt('R',2.5,0.50,0.5,'mean')
	 if (abs(sort).eq.2) then
            call pgmtxt('R',2.5,0.65,0.5,'std dev')
	 else
            call pgmtxt('R',2.5,0.65,0.5,'RMS')
	 endif
      endif
      if (title.eq.'-') then
      call pglab(text,'xover height difference (cm)',arg)
      else
      call pglab(text,'xover height difference (cm)',title)
      endif

      if (nodots) then
      else if (dtplot) then

* Plot xover height differences (asc-des or first-last)

      call pgsci(3)
      call pgslw(lw)
      do i=1,n1
	 t(i,1)=abs(t(i,1)-t(i,2))
      enddo
      call pgpt(n1,t(1,1),xover,-1)
      call pgslw(1)
      call blabel(sort,2.5)

      else

* Plot xover height differences (asc-des or first-last)

      call pgsci(2)
      call pgslw(lw)
      call pgpt(n1,t(1,1),xover,-1)
      call pgslw(1)
      call blabel(sort,2.5)

* Plot xover height differences (des-asc or last-first)

      do i=1,n1
	 xover(i)=-xover(i)
      enddo
      call pgsci(3)
      call pgslw(lw)
      call pgpt(n1,t(1,2),xover,-1)
      call pgslw(1)
      call blabel(-sort,3.5)

      endif

* Plot daily statistics

      call pgsci(1)
      call pgslw(lw)
      call pgdaily(nt,real(trms),real(itmin),1.0)
      call pgdaily(nt,real(tmean),real(itmin),1.0)
      if (.not.norejected) then
         call pgswin(real(tmin),real(tmax),0.0,20.0)
         call pgdaily(nt,real(trej),real(itmin),1.0)
      endif

      endif ! END plotting

* Write statistics per file

      write (*,600)
     |  n1,avg1,rms1,sgm1,n2,avg2,rms2,sgm2,ssvar,arg(:l)
      if (table.ne.' ') write (10,603)
     |  n1,avg1,rms1,sgm1,n2,avg2,rms2,sgm2,ssvar,arg(:l)

      if (if.ne.nf) goto 10

  550 format (a)
  600 format (i8,3f7.2,1x,i8,4f7.2,2x,a)
  601 format (i8,3f7.2,1x,i8,3f7.2,8x,'* Total *')
  602 format ('----------- Input ----------- ---------- Edited -----------'/
     |2('      nr   mean    rms    std '),' ssvar  filenm')
  603 format ('# ----------- Input ----------- ---------- Edited -----------'/
     |'# ',2('      nr   mean    rms    std '),' ssvar  filenm'/
     |'# ',i8,3f7.2,1x,i8,4f7.2,2x,a)

* Print overall statistics

      t_avg1=t_avg1/t_n1
      t_avg2=t_avg2/t_n2
      t_rms1=sqrt(t_rms1/t_n1)
      t_rms2=sqrt(t_rms2/t_n2)
      t_sgm1=sqrt(t_rms1**2-t_avg1**2)
      t_sgm2=sqrt(t_rms2**2-t_avg2**2)
      write (*,601) t_n1,t_avg1,t_rms1,t_sgm1,
     |   t_n2,t_avg2,t_rms2,t_sgm2
      goto 9999

 1300 write (*,1301)
 1301 format (
     |'xstat: compute crossover height difference statistics'//
     |'usage: xstat [options] file(s)'//
     |'required:'/
     |'  file(s)        : xover file names (XXF, XXO, XXS)'//
     |'options:'/
     |'  edit=mult      : use mult*sigma as edit level (def=3.5)'/
     |'  step=step      : step through file by "step"'/
     |'  xmax=xmax      : use xmax as edit level (cm)'/
     |'  t=tmin,tmax    :',
     |' specify time interval ([yy]yymmdd[hhmmss],mjd,doy)'/
     |'              .... or use mjd=, doy=, ymd=, sec='/
     |'  T=tmin,tmax    :',
     |' specify plot interval ([yy]yymmdd[hhmmss],mjd,doy)'/
     |'              .... or use MJD=, DOY=, YMD=, SEC='/
     |'  dt=dtmax       : only xovers with dt < dtmax (days)'/
     |'  dt=dtmin,dtmax : only xovers with dtmin < dt < dtmax (days)'/
     |'  bin=dt         : bin size for plotting (days)'/
     |'  lat=lat0,lat1  : only xovers with lat0<lat<lat1'/
     |'  lon=lon0,lon1  : only xovers with lon0<lon<lon1'/
     |'  -s             : only xovers over shallow waters'/
     |'  -o             : only xovers over deep waters'/
     |'  -vVAR          : only xovers with ssvar < VAR (cm) (def:100)'/
     |'  +vVAR          : only xovers with ssvar > VAR (cm) (def:100)'/
     |'  -a             : compute asc - des (default)'/
     |'  -e             : compute ERS - T/P (default)'/
     |'  -t             : compute first - last (time sort)'/
     |'  +t             : compute last - first (time sort)'/
     |'  -T             : compute sat A - sat B (tracknr sort)'/
     |'  +T             : compute sat B - sat A (tracknr sort)'/
     |'  -n             : do not plot xover height differences'/
     |'  -N             : as -n, also no xover RMS, only mean'/
     |'  -c             : combine statistics of all xover files'/
     |'  -f             : plot no arc numbers'/
     |'  -y             : label time axis with years'/
     |'  -tpalt         : model T/P altitude rate'/
     |'  rmax=rmax      : maximum xover rms in plot (cm, def=18)'/
     |'  -swh           : use SWH field i.s.o. SSH'/
     |'  fld=num        : use auxiliary field num i.s.o. SSH'/
     |'  tbias=tb1(,tb2): remove time tag bias',
     |' (ms, only with AUX file)'/
     |'  bias=b1(,b2)   : add bias to SSH (cm)'/
     |'  chbias=i1,i2   : exchange new for old SPTR correction'/
     |'  ssb=ssb1(,ssb2): apply extra ssb (%, only with AUX file)'/
     |'  tab=table      : generate table'/
     |'  title=title    : plot title above plot i.s.o. filename'/
     |'  dev=device     : plot xover differences on device'/
     |'  -dt            : when plot selected: xover as function of dt')
      goto 9999

 1310 write (*,1311) n1,mxover
 1311 format('xstat: too many xovers: ',i7,'>',i7)
 9999 if (pgopen) then
	 call pgask(.false.)
	 call pgend
      endif
      end

      subroutine meanrms (n,x,mean,rms)
      integer n,i
      real*4 x(n)
      real*8 mean,rms

      mean=0d0
      rms=0d0
      do i=1,n
	 mean=mean+x(i)
	 rms =rms +x(i)**2
      enddo
      if (n.gt.0) then
         mean=mean/n
         rms =sqrt(rms/n)
      else
	 mean=0d0
	 rms=0d0
      endif
      end

      subroutine editor (n,x,v,mean,rms,var,edit,xmax)
      integer n,i,j,jold
      real*4 x(n),v(n)
      real*8 mean,rms,var,edit,xmax,x0,x1,sgm,sgmold

      j=n
      sgm=rms
10    jold=j
      sgmold=sgm
      x0=mean-min(xmax,sgm*edit)
      x1=mean+min(xmax,sgm*edit)
      mean=0d0
      rms=0d0
      var=0d0
      j=0

      do i=1,n
	 if (x(i).ge.x0 .and. x(i).le.x1) then
	    mean=mean+x(i)
	    rms =rms +x(i)**2
	    var =var +v(i)**2
	    j=j+1
	 endif
      enddo

      if (j.eq.0) then
	 mean=0d0
	 rms=0d0
	 var=0d0
	 sgm=0d0
	 n=0
	 return
      else
         mean=mean/j
         rms =sqrt(rms/j)
         var =sqrt(var/j)
         sgm =sqrt(rms**2-mean**2)
      endif
      if (j.ne.jold .or. abs(sgmold/sgm-1).gt.1d-2) goto 10
      n=j
      end

      subroutine pgdaily(n,y,x0,dx)
      integer*4 n,i
      real*4 y(*),x0,dx,x
      do i=1,n
	 x=x0+(i-1)*dx
         if (i.eq.1) then
	    call pgmove(x,y(i))
	 else if (y(i-1).gt.1e20) then
	    call pgmove(x,y(i))
         else if (y(i).lt.1e20) then
	    call pgdraw(x,y(i))
         endif
         if (y(i).lt.1e20) then
	    call pgdraw(x+dx,y(i))
         endif
      enddo
      end

      subroutine blabel(sort,offset)
      integer sort
      real offset
      if (sort.eq.+1) then
         call pgmtxt('B',offset,1.0,1.0,'1 - 2')
      else if (sort.eq.-1) then
         call pgmtxt('B',offset,1.0,1.0,'2 - 1')
      else if (sort.eq.+2) then
         call pgmtxt('B',offset,1.0,1.0,'A - B')
      else if (sort.eq.-2) then
         call pgmtxt('B',offset,1.0,1.0,'B - A')
      else
         call pgmtxt('B',offset,1.0,1.0,'asc - des')
      endif
      end
