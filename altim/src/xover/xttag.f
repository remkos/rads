**XTTAG -- Determine time tag bias from crossovers
*+
      program xttag
*
* This program uses XXO files to determine time tag biases in altimeter
* data. This is mainly based on the 2-cpr signal that is the result of such
* a time tag bias.
*-
* $Log: xttag.f,v $
* Revision 1.12  2012/05/24 12:52:58  rads
* - Higher precision
*
* Revision 1.11  2012/04/23 19:35:02  rads
* - Update from 3.5 to 5 million xovers maximum
*
* Revision 1.10  2011/03/08 01:20:52  rads
* - Start bins consistently
*
* (c) Remko Scharroo - TU Delft, Altimetrics
*-----------------------------------------------------------------------
      implicit none

* General

      character filenm(50)*80,spec*4,arg*80,text*80,table*80/' '/
      integer*4 i,j,k,nrec,n1,npar,iarg,iargc,j1,j2
      integer*4 l,ftype/0/,recl,chbias(2)/2*0/
      integer*4 lat,lon,dt,dorb/0/,dtrk,if,nf/0/,sort/0/,tmod(2)/2*0/
      real*8 	dtmin/-1d40/,dtmax/1d40/,edit/3.5d0/,xmax/1d40/,
     |	mean,rms,ssb/0d0/,tbias(2)/2*0d0/,bmin,bmax,dssh,
     |	cor0,cor1,cor2
      real*8 rlon,rlat
      real*8 latmin/-1e35/,latmax/1e35/,lonmin/-1e35/,lonmax/1e35/
      logical ocean/.false./,sland,shallow/.false./,
     |	ersasc/.false./,ersdes/.false./,datearg,zero/.false./
      real*8 smin/-1d40/,smax/1d40/,slo/1d40/,shi/-1d40/,dum

* For timing error model

      real*8 pi,rad,u,u_e,u_k,sin_i,sin_e,mjd,epharg,
     |	b0(2)/-1.468,-1.354/,
     |	b1(2)/ 0.483, 0.739/,
     |	b2(2)/ 0.286, 0.284/,
     |	u1(2)/ 102.6, 104.5/,
     |	u2(2)/  15.3,  11.7/
      real*8	tppar(3)/15d3,0.048d0,-2d3/
      logical	tpalt/.false./

* File handling

      integer*4 openf,readf,closef,ios,fd/-1/,fdaux/-1/,fddlt/-1/
      integer*4 xxf4(16),txo(2)
      integer*4 maux,naux/0/,fldn/0/,fld/-1/
      parameter (maux=20)
      integer*2 xxf2(32),aux(maux)/maux*0/,dlt(2)/2*0/
      equivalence (xxf2,xxf4)
      integer   mode/0/

* Storage

      integer mxover
      parameter (mxover=5000000)
      real*8 xover(mxover)

* For variability

      character*80 varnm
      integer nx,ny
      real*8 x0,x1,y0,y1,z0,z1,fx/1e0/,fy/1e0/,vmax/1e10/
      integer maxgrd
      parameter (maxgrd=360*180)
      real*8 var(maxgrd),v(mxover)
      logical dovar/.false./

* For plotting

      real*8 t(mxover,2),day,bin/1.0/
      parameter (day=86400.)
      real*8 tmin,tmax,tref
      integer yref,mdate
      character*80 dev/'-'/,title/'-'/
      logical pgactive/.false./,combi/.false./,years/.false./
      integer it1,it2,nt,maxt,lw/1/,pgbeg
      parameter (maxt=4000)
      real*8 zmin,zmax
      integer tnr(maxt),mrk/17/

* For time tag estimation

      real*8 atr(maxt),ata(maxt*(maxt+1)/2),
     |       atr_save(maxt),ata_save(maxt*(maxt+1)/2),
     |       rtr0,rtr1(maxt),rtr2,sum,cor,s
      integer h(maxt),itmin,itmax
      real*8 dxdt(mxover,2)

* Determine mode

      call getarg(0,arg)
      if (index(arg,'ttag-h').gt.0) then
         mode=4
	 fdaux=0
      else if (index(arg,'ttag').gt.0) then
	 mode=1
	 fdaux=0
      else if (index(arg,'ssb').gt.0) then
	 mode=2
	 fdaux=0
      else if (index(arg,'bias').gt.0) then
	 mode=3
      else
	 call fin("xttag: unrecognized call")
      endif

* Scan argument line

      zmin=0
      zmax=0
      do iarg=1,iargc()
      call getarg(iarg,arg)
      if (arg(:1).ge.'A' .and. arg(:1).le.'Z'
     |	.and. datearg(arg,slo,shi,dum)) then
      else if (datearg(arg,smin,smax,dum)) then
      else if (arg(1:2).eq.'z=') then
         read (arg(3:),*) zmin,zmax
      else if (arg(1:4).eq.'bin=') then
	 read(arg(5:),*) bin
      else if (arg(1:3).eq.'dt=') then
	 dtmax=0d0
	 read (arg(4:),*,iostat=ios) dtmin,dtmax
	 if (dtmax.eq.0d0) then
	    dtmax=dtmin
	    dtmin=0d0
	 endif
	 dtmin=dtmin*86400
	 dtmax=dtmax*86400
      else if (arg(1:3).eq.'lw=') then
	 read (arg(4:),*) lw
      else if (arg(1:4).eq.'tab=') then
         table=arg(5:)
         open (10,file=table)
         if (mode.eq.1) then
	    write (10,550) '# Date Tbias(msec) Nr Var'
	 else if (mode.eq.2) then
	    write (10,550) '# SWH  SSB(-) Nr Var'
	 else if (mode.eq.3) then
	    write (10,550) '# Date Rbias(mm) Nr Var'
	 else
	    write (10,550) '# Hour Tbias(msec) Nr Var'
	 endif
      else if (arg(1:4).eq.'dev=') then
	 dev=arg(5:)
      else if (arg(1:4).eq.'lon=') then
	 read (arg(5:),*) lonmin,lonmax
      else if (arg(1:4).eq.'lat=') then
	 read (arg(5:),*) latmin,latmax
      else if (arg(1:5).eq.'edit=') then
	 read (arg(6:),*) edit
      else if (arg(1:7).eq.'chbias=') then
	 chbias(2)=9999
	 read (arg(8:),*,iostat=ios) chbias
	 if (chbias(2).eq.9999) chbias(2)=chbias(1)
      else if (arg(1:6).eq.'tbias=') then
         tbias(2)=1d30
         read (arg(7:),*,iostat=ios) tbias
         if (tbias(2).gt.1d20) tbias(2)=tbias(1)
* tbias in milliseconds. Positive when tags are late.
* Will be multiplied by orbit rate in mm/s to give micrometres
	 fdaux=0
      else if (arg(1:5).eq.'tmod=') then
         tmod(2)=-1
	 read (arg(6:),*,iostat=ios) tmod
	 if (tmod(2).eq.-1) tmod(2)=tmod(1)
      else if (arg(1:4).eq.'ssb=') then
         read (arg(5:),*) ssb
* Addition ssb in percentage. Positive when ssb should be higher.
         ssb=ssb*1d1
      else if (arg(1:2).eq.'-v') then
	 read (arg(3:),*,iostat=ios) vmax
	 vmax=vmax/100
	 dovar=.true.
      else if (arg(1:5).eq.'xmax=') then
	 read (arg(6:),*) xmax
	 xmax=xmax*1d4	! cm to micrometers
      else if (arg(:6).eq.'-tpalt') then
         tpalt=.true.
      else if (arg(1:2).eq.'-o') then
	 write (*,550) '(Ocean only)'
	 ocean=.true.
	 shallow=.false.
      else if (arg(1:2).eq.'-c') then
	 combi=.true.
      else if (arg(1:2).eq.'-y') then
	 years=.true.
      else if (arg(:4).eq.'fld=') then
	 read (arg(5:),*) fld
         write (*,550) '(auxiliary field statistics)'
	 fld=(fld-1)*2
	 fdaux=0
      else if (arg(1:2).eq.'-s') then
	 write (*,550) '(Shallow waters)'
	 shallow=.true.
	 ocean=.false.
      else if (arg(:2).eq.'-a' .or. arg(:2).eq.'-e') then
	 write (*,550) '(asc - des, ers - t/p)'
	 if (arg(1:3).eq.'-eA') ersasc=.true.
	 if (arg(1:3).eq.'-eD') ersdes=.true.
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
      else if (arg(1:6).eq.'title=') then
	 title=arg(7:)
      else if (arg(1:4).eq.'mrk=') then
	 read (arg(5:),*) mrk
      else if (arg(:2).eq.'-z') then
	 zero=.true.
      else
	 nf=nf+1
	 if (nf.gt.50) call fin("xttag: too many files")
         filenm(nf)=arg
      endif
      enddo
      if (nf.eq.0) goto 1300

* If variability grid is requested: load grid

      if (dovar) then
         nx=0
         ny=maxgrd
	 varnm='/user/altim'
	 call checkenv('ALTIM',varnm,l)
	 varnm(l+1:)='/data/world/ssvar95a.grd'
         call gridrd4(varnm,nx,ny,var,x0,x1,y0,y1,z0,z1)
         fx=(x1-x0)/(nx-1)
         fy=(y1-y0)/(ny-1)
      endif

* Initialise 

      call matsy1(maxt,h)
      pi=4*atan(1d0)
      rad=pi/180
      sin_i=sin(98.5d0*rad)
      sin_e=sin(23.5d0*rad)
      u_k=90d0*rad
      do j=1,2
         u1(j)=u1(j)*rad
	 u2(j)=u2(j)*rad
      enddo

      if=0

* Process all xover files

10    continue
      n1=0
      rtr0=0

* Intialise periods

      if (mode.eq.4) then
         slo=0
	 shi=day
	 bin=bin/24 ! in mode=4 bin size is in hours
      else if (slo.lt.shi) then
         smin=max(smin,slo)
         smax=min(smax,shi)
      endif
      tmin=slo/bin/day
      tmax=shi/bin/day

* First determine file type (XXO, XXS, XXE, XXF and XXB supported)

11    if=if+1
      arg=filenm(if)
      fddlt=-1
12    fd=openf(arg,'r')
      if (fd.lt.0) call fin("Error opening xover file")
      ios=readf(fd,4,spec)
      if (spec.eq.'@D_H') then
         fddlt=fd
	 call seekf(fddlt,24,0)
         l=index(arg,'.xxo')
         if (l.le.0) call fin("Error in DLT filename")
         arg=arg(:l+3)
         goto 12
      endif
      ios=readf(fd,4,nrec)
      ios=readf(fd,4,npar)
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
	 if (fdaux.lt.0) call fin('error reading AUX file')
	 ios=readf(fdaux,4,spec)
	 if (spec.ne.'@AUX') call fin ('illegal file type on aux file')
	 ios=readf(fdaux,4,i)
	 if (i.ne.nrec) call fin('aux file does not match main file')
	 ios=readf(fdaux,4,naux)
	 if (naux.gt.maux) call fin('too many parameters in aux file')
	 call seekf(fdaux,naux*2,0)
	 fldn=naux-2
      endif

* Read all xovers and store height difference in xover()

      do i=1,nrec

* Read different data types

	 ios=readf(fd,recl,xxf4)
	 if (fdaux.ge.0) ios=readf(fdaux,naux*2,aux)
	    
	 if (ftype.eq.3) then
	    txo(1)=xxf4(3)
	    txo(2)=xxf4(4)
	    dtrk=xxf2(9)-xxf2(10)
	    dssh=xxf4(8)-xxf4(9)
	    lat=xxf4(1)
	    lon=xxf4(2)
         else
	    txo(1)=xxf4(3)
	    txo(2)=xxf4(5)
	    dtrk=xxf2(13)-xxf2(14)
	    dssh=xxf4(8)-xxf4(9)
	    lat=xxf4(1)
	    lon=xxf4(2)
	    if (ftype.eq.1) dorb=xxf4(10)-xxf4(11)
         endif
         if (fddlt.ge.0) then
            ios=readf(fddlt,4,dlt)
            if (dlt(1).eq.-32768 .or. dlt(2).eq.-32768) goto 20
            dssh=dssh+(dlt(1)-dlt(2))*1d3
         endif

	 dt=abs(txo(1)-txo(2))
	 rlat=lat/1e6
	 rlon=lon/1e6

* Model T/P altitude rate if requested

	 if (tpalt) then
	    aux(fldn+1)=nint(
     |		tppar(1)*sin(tppar(2)*rlat)+tppar(3)*sin(2*tppar(2)*rlat))
	    aux(fldn+2)=-aux(fldn+1)
	 endif

* Compute timing bias according to model, if requested

	 do j=1,2
	    if (tmod(j).gt.0) then
	       k=tmod(j)
	       u=asin(sin(lat/1d6*rad)/sin_i)
	       if (j.eq.2) u=pi-u
	       mjd=txo(j)/86400d0+46066
	       u_e=-asin(sin_e*sin(epharg(3,mjd)*rad))
	       tbias(j)=b0(k)+b1(k)*(sin(u-u_e-u1(k))-sin(u_k-u_e-u1(k)))+
     |			b2(k)*(sin(2*(u-u_e-u2(k)))-sin(2*(u_k-u_e-u2(k))))
	    endif
	 enddo

* Which of the two measurements goes first? In other words, do we use
* (1)-(2) or (2)-(1)?

	 if ((sort.eq.+1 .and. txo(1).gt.txo(2)) .or.
     |	     (sort.eq.-1 .and. txo(1).lt.txo(2)) .or.
     |	     (sort.eq.+2 .and. dtrk.gt.0) .or.
     |	     (sort.eq.-2 .and. dtrk.lt.0) .or.
     |	     (dorb.gt.100 000 000)) then
	    ! Do (2)-(1) only when time/track sorting is requested or (1) is T/P
	    j1=2
	    j2=1
	    dssh=-dssh
	 else
	    ! The normal case: (1)-(2)
	    j1=1
	    j2=2
	 endif

* Zero out SSH (for test purposes only)

	 if (zero) dssh=0d0

*	 if (xxf4(1).lt.-34000000 .and. xxf4(1).gt.-35000000)
*     |		write (*,*) xxf4(1),aux(fldn+j1),aux(fldn+j2)

	 if (fld.ge.0) then
* Replace SSH difference with auxiliary field if requested
	    dssh=(aux(j1+fld)-aux(j2+fld))*1000
	 else
* Apply timing bias and/or sea state bias correction
	    dssh=dssh+ssb*(aux(j1)-aux(j2))-
     |		   (tbias(1)*aux(j1+fldn)-tbias(2)*aux(j2+fldn))
         endif

* Change bias if requested

	 if (chbias(1).ne.0) then
	    call altbias(-chbias(1),txo(j1),cor0,cor1,cor2)
	    dssh=dssh-(cor0-cor1-cor2)*1d6
	    call altbias(+chbias(1),txo(j1),cor0,cor1,cor2)
	    dssh=dssh+(cor0-cor1-cor2)*1d6
         endif
	 if (chbias(2).ne.0) then
	    call altbias(-chbias(2),txo(j2),cor0,cor1,cor2)
	    dssh=dssh+(cor0-cor1-cor2)*1d6
	    call altbias(+chbias(2),txo(j2),cor0,cor1,cor2)
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
	 else if (abs(dssh).gt.xmax) then
	 else if (ersasc .and. dorb.gt. 100000000) then
	 else if (ersdes .and. dorb.lt.-100000000) then
	 else
            if (n1.gt.mxover) goto 1310
	    if (dovar) then
	       l=nint((rlat-y0)/fy)*nx+nint((rlon-x0)/fx)+1
	       if (var(l).gt.vmax) goto 20
	       n1=n1+1
	       v(n1)=var(l)*100
	    else
	       n1=n1+1
	    endif
	    xover(n1)=dssh*1d-6
	    t(n1,1)=txo(j1)/bin/day
	    t(n1,2)=txo(j2)/bin/day
	    if (mode.eq.1) then
	       dxdt(n1,1)=+aux(j1+fldn)/1d3
	       dxdt(n1,2)=-aux(j2+fldn)/1d3
	    else if (mode.eq.2) then
	       dxdt(n1,1)=-aux(j1)/1d2
	       dxdt(n1,2)=+aux(j2)/1d2
	       t(n1,1)=aux(j1)/bin/1d3
	       t(n1,2)=aux(j2)/bin/1d3
	    else if (mode.eq.3) then
	       dxdt(n1,1)=1
	       dxdt(n1,2)=1
	    else
	       dxdt(n1,1)=+aux(j1+fldn)/1d3
	       dxdt(n1,2)=-aux(j2+fldn)/1d3
	       t(n1,1)=mod(txo(j1)/day,1.0)/bin
	       t(n1,2)=mod(txo(j2)/day,1.0)/bin
	    endif
	    tmin=min(tmin,t(n1,1),t(n1,2))
	    tmax=max(tmax,t(n1,1),t(n1,2))
	    ! Do not use second time tag when sorting or T/P
	    if (sort.ne.0 .or. abs(dorb).gt.100000000) t(n1,2)=1e30
	 endif
20	 continue
      enddo
      if (fd.gt.0) ios=closef(fd)
      if (fdaux.gt.0) ios=closef(fdaux)
      if (fddlt.gt.0) ios=closef(fddlt)
      if (combi .and. if.ne.nf) goto 11

* Initialise matrix

      itmin=int(tmin)
      itmax=int(tmax+1)
      nt=itmax-itmin
      if (nt.gt.maxt) call fin("xttag: too many days")
      do i=1,nt*(nt+1)/2
	 ata(i)=0d0
      enddo
      do i=1,nt
	 tnr(i)=0
	 atr(i)=0d0
      enddo

      do i=1,n1
	 it1=int(t(i,1)-itmin+1)
	 it2=int(t(i,2)-itmin+1)
	 if (dxdt(i,1).ne.0) tnr(it1)=tnr(it1)+1
	 if (dxdt(i,2).ne.0) tnr(it2)=tnr(it2)+1
	 j=h(it1)+it1
	 ata(j)=ata(j)+dxdt(i,1)*dxdt(i,1)
	 j=h(it2)+it2
	 ata(j)=ata(j)+dxdt(i,2)*dxdt(i,2)
	 if (it1.lt.it2) then
	    j=h(it2)+it1
	 else
	    j=h(it1)+it2
	 endif
	 ata(j)=ata(j)+dxdt(i,1)*dxdt(i,2)
	 if (it1.eq.it2) ata(j)=ata(j)+dxdt(i,1)*dxdt(i,2)
	 atr(it1)=atr(it1)+dxdt(i,1)*xover(i)
	 atr(it2)=atr(it2)+dxdt(i,2)*xover(i)
	 rtr0=rtr0+xover(i)**2
      enddo

* Add (mild) constraint

      do i=1,nt
	 j=h(i)+i
	 ata(j)=ata(j)+1d-20
      enddo

* Save AtA and AtR matrices

      do i=1,nt
         atr_save(i)=atr(i)
      enddo
      do i=1,nt*(nt+1)/2
         ata_save(i)=ata(i)
      enddo

* Solve

      call dpptrf('U',nt,ata,ios)
      if (ios.ne.0) call fin("xttag: error solving normal equation")
      call dpptrs('U',nt,1,ata,atr,nt,ios)

* Compute xover variance per parameter. This is Xt*AtA*X/n1=Xt*AtR/n1

      rtr2=rtr0
      do i=1,nt
	 rtr1(i)=atr_save(i)*atr(i)
	 rtr2=rtr2-rtr1(i)
	 if (rtr1(i).ge.0) then
	    rtr1(i)=sqrt(rtr1(i)/n1)
	 else
	    rtr1(i)=-sqrt(-rtr1(i)/n1)
         endif
      enddo

* Run through all values, do editing

      bmin=-1d40
      bmax=1d40
      do j=1,3
      sum=0
      rms=0
      cor=0
      mean=0
      do i=1,nt
         tnr(i)=abs(tnr(i))
 	 if (tnr(i).ne.0) then
	    if (mode.eq.2) then
	       s=tnr(i)*(itmin+i-0.5)*bin
	    else
	       s=tnr(i)
	    endif
	    if (atr(i).lt.bmin .or. atr(i).gt.bmax) then
	       tnr(i)=-tnr(i)
	    else
	       mean=mean+s*atr(i)
	       rms=rms+s*atr(i)**2
	       cor=cor+s*atr(i)**2*(itmin+i-0.5)
	       sum=sum+s
	    endif
 	 endif
      enddo
      mean=mean/sum
      cor=cor/sum
      rms=sqrt(rms/sum-mean**2)
      bmin=mean-3.5*rms
      bmax=mean+3.5*rms
      enddo
      write (*,31) sum,mean*1d3,rms*1d3,cor*1d6,
     |             n1,sqrt(rtr0/n1)*1d3,sqrt(rtr2/n1)*1d3

* If table is requested

      if (table.ne.' ') then
      do i=1,nt
         tnr(i)=abs(tnr(i))
 	 if (tnr(i).ne.0) then
	    if (mode.eq.2) then
	       write (10,30) (itmin+i-0.5)*bin,atr(i)*1d3,
     |		tnr(i),rtr1(i)*1d3
	    else if (mode.eq.4) then
	       write (10,30) (itmin+i-0.5)*bin*24,atr(i)*1d3,
     |		tnr(i),rtr1(i)*1d3
	    else
	       write (10,29) mdate(1,int(46066+(tmin+i-1)*bin)),atr(i)*1d3,
     |		tnr(i),rtr1(i)*1d3
	    endif
	 endif
      enddo
      write (10,31) sum,mean*1d3,rms*1d3,cor*1d6,
     |              n1,sqrt(rtr0/n1)*1d3,sqrt(rtr2/n1)*1d3
      endif
29    format (i6.6,f12.6,i10,f9.3)
30    format (f9.2,f12.6,i10,f9.3)
31    format (
     |'#           Nr        Mean    Rms@Mean        WRms'/
     |'# ',f12.1,3f12.6/
     |'# Number of xovers            = ',i9/
     |'# A priori xover rms (mm)     = ',f9.3/
     |'# A posteriori xover rms (mm) = ',f9.3)

* If plotting is requested: plot xover height differences

      if (dev.ne.'-') then

* First determine reference day ( = beginning of the year )

      yref=mdate(1,int(tmin*bin+46066))/10000
      tref=mdate(2,yref*10000+0101)-1-46066
      if (zmin.ne.zmax) then
      else if (mode.eq.1) then
         zmin=-3.
         zmax=+1.
      else if (mode.eq.2) then
         zmin=-0.5
         zmax=+1.0
      else
	 zmin=-30.
	 zmax=+30.
      endif
      if (.not.pgactive) then
         if (pgbeg(0,dev,1,1).ne.1)
     |		call fin("Error opening plot device")
         call pgscr(0,0.8,0.8,0.8)
         call pgscr(1,0.0,0.0,0.0)
         call pgscr(2,1.0,0.0,0.0)
         call pgscr(3,0.0,0.5,0.0)
	 pgactive=.true.
      else
         call pgpage
      endif
      if (mode.eq.2) then
	 call pgswin(int(tmin*bin)*1.0,int(tmax*bin)*1.0,real(zmin),real(zmax))
	 call pgbox('abcnsti',0.0,0,'abcnsti',0.0,0)
      else
         call pgsci(2)
         call pgslw(1)
         call pgsch(0.7)
         call pgswin
     |   (real((tmin*bin-2381)/3.5d0),real((tmax*bin-2381)/3.5d0),0.0,1.0)
         call pgbox('csti',1.0,0,' ',0.0,0)
         if (tmin.gt.3761) then
            call pgswin(real((tmin*bin-2381)/3.5d0-394.5d0),
     |	    real((tmax*bin-2381)/3.5d0-394.5d0),0.0,1.0)
         else
            call pgswin(real((tmin*bin-2381)/3.5d0-0.5d0),
     |		real((tmax*bin-2381)/3.5d0-0.5d0),0.0,1.0)
         endif
         call pgbox('m',0.0,0,' ',0.0,0)
         call pgsch(1.0)
         call pgsci(1)
	 if (years) then
            call pgswin(real(tmin*bin/365.25d0),real(tmax*bin/365.25d0),
     |		real(zmin),real(zmax))
            call pgbox('cst',1.0,12,' ',0.0,0)
            call pgbox('bsti',1.0,12,' ',0.0,0)
            call pgswin(real(tmin*bin/365.25d0+1984.5d0),
     |		real(tmax*bin/365.25d0+1984.5d0),0.0,1.0)
            call pgbox('n',1.0,0,' ',0.0,0)
            text='year'
         else
            call pgswin(real(tmin*bin-tref),real(tmax*bin-tref),real(zmin),real(zmax))
            call pgbox('abnti',0.0,0,'bcnsti',1.0,10)
	    if (yref.ge.50) then
               write (text,700) 1900+yref
	    else
               write (text,700) 2000+yref
	    endif
700         format ('time (days of the year ',i4,')')
         endif
         call pgswin(real(tmin),real(tmax),real(zmin),real(zmax))
	 call pgbox(' ',0.0,0,'bcnsti',1.0,10)
      endif
      arg=filenm(if)
      if (title.ne.'-') then
      else if (combi) then
         title='* combi *'
      else
         title=arg
      endif
      if (mode.eq.1 .or. mode.eq.4) then
         call pglab(text,'resolved time-tag bias (ms)',title)
      else if (mode.eq.2) then
         call pglab('SWH (m)',
     |		'resolved sea-state bias (% of SWH)',title)
      else
         call pglab(text,'resolved range bias (mm)',title)
      endif

* Plot dayly statistics

      if (mode.eq.2) then
	 call pgslw(3)
	 call pgsch(2.)
	 do i=1,nt
	    call pgpt(1,real((itmin+i-0.5)*bin),real(atr(i)*1d3),mrk)
	 enddo
      else
         call pgsci(1)
         call pgslw(lw)
         call pgmove(real(itmin),real(atr(1)*1d3))
         do i=1,nt
	    if (tnr(i-1).le.0) then
	       call pgmove(real(itmin+i-1),real(atr(i)*1d3))
	    else if (tnr(i).gt.0) then
	       call pgdraw(real(itmin+i-1),real(atr(i)*1d3))
	    endif
	    if (tnr(i).gt.0) then
	       call pgdraw(real(itmin+i  ),real(atr(i)*1d3))
	    endif
         enddo
      endif

      endif ! END plotting

      if (if.ne.nf) goto 10

  550 format (a)

      goto 9999

 1300 write (*,1301)
 1301 format ('usage: xttag/xssb [options] xoverfile ...'/,/
     |'where [options] are:'/,/
     |'  edit=mult       : use mult*sigma as edit level (def=3.5)'/
     |'  step=step       : step through file by "step"'/
     |'  xmax=xmax       : use xmax as edit level (cm)'/
     |'  t=tmin,tmax     : specify time interval (yymmdd,mjd)'/
     |'                ... or use mjd=, doy=, ymd=, sec='/
     |'  T=tmin,tmax     : specify plotting interval (yymmdd,mjd)'/
     |'                ... or use MJD=, DOY=, YMD=, SEC='/
     |'  dt=dtmin,dtmax  : only xovers with dtmin<dt<dtmax days'/
     |'  lat=lat0,lat1   : only xovers with lat0<lat<lat1'/
     |'  lon=lon0,lon1   : only xovers with lon0<lon<lon1'/
     |'  -s              : only xovers over shallow waters'/
     |'  -o              : only xovers over deep waters'/
     |'  -vVAR           : only xovers with ssvar < VAR (cm)'/
     |'  -a              : compute asc - des (default)'/
     |'  -e              : compute ERS - T/P (default)'/
     |'       -eA        : compute ERS (asc) - T/P'/
     |'       -eD        : compute ERS (des) - T/P'/
     |'  -t              : compute first - last (time sort)'/
     |'  +t              : compute last - first (time sort)'/
     |'  -T              : compute sat A - sat B (tracknr sort)'/
     |'  +T              : compute sat B - sat A (tracknr sort)'/
     |'  -c              : combine statistics of all xover files'/
     |'  -y              : label time axis with years'/
     |'  -z              : zero out SSH (for test purposes only)'/
     |'  -tpalt          : model T/P altitude rate'/
     |'  tmod=n          : use tbias model n'/
     |'  bin=dt          : bin size for plotting (days or metres)'/
     |'  z=zmin,zmax     : vertical range of window (cm or ms)'/
     |'  tbias=tb1(,tb2) : remove time tag bias (ms)'/
     |'  ssb=ssb         : apply extra ssb (%)'/
     |'  tab=table       : make table of bias estimated'/
     |'  title=title     : put title on top of plot'/
     |'  mrk=number      : use marker number (xssb only)'/
     |'  dev=device      : plot xover differences on device')
      goto 9999

 1310 write (*,1311) n1,mxover
 1311 format('xttag: too many xovers: ',i7,'>',i7)
 9999 if (pgactive) then
	 call pgask(.false.)
	 call pgend
      endif
      end
