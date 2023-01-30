      program odrdiff

      integer npnt,nhist,mrev
      parameter(npnt=70000,nhist=200,mrev=1000)
      integer nsmth,limit,nstep,nstat,narg,iargc,iarg,
     |  irep,iarc,irem,i,nrec1,nrec2,it10,it11,it20,it21,itst1,itst2,
     |  ibegin,iend,it0,it1,itstat0,itstat1,ioff1,ioff2,nrec,layout,
     |  nrev,ix,iy,ih,j,itime,nstat0,krev,krev0,lnblnk,nout,jmax,n,
     |  mjd92
      real*4 dummy,rms_sin,rms_cos,d,rmin,rmax
      integer rec(4),lon(npnt),lat(npnt)
      real*4 t(npnt),h(npnt),x(npnt),y(npnt)
      real*4 tmin,tmax,hmin,hmax,smax,xmin,xmax,ymin,ymax,dx,dy
      character*80 file1/' '/,file2/' '/,arg,dev,title,nmlfile
      character satel*8,date0*15,date1*15,text*160
      real*8 hmean,hrms,hsigma,delta,xrms,xmean,yrms,ymean,
     |  dt(npnt),dh(npnt),wk1(npnt),wk2(npnt),prob,day,rev,
     |  ofac,hifac,pi,rad,realtemp(2),amp,
     |  tnode(mrev),phase,tnode0,dlat,dlon,hgt,step
      complex*16 wk3(npnt),temp,expph0
      logical test/.false./,mid/.false./,out/.false./,xgf/.false./,
     |		getorbit/.false./,ltlend
      equivalence (wk1,wk3),(temp,realtemp)
      real*4 vx(4),vy(8),ch,hifreq
      integer getorb,fd1,fd2,scale1/1/,scale2/1/,odrinfo,pgbeg,l
      namelist /nml/ vx,vy,ch,nsmth,hifreq,nstep,ofac,layout

      nmlfile='/user/altim'
      call checkenv('ALTIM',nmlfile,l)
      nmlfile(l+1:)='/nml/odrdiff.nml'
      open (7,file=nmlfile,status='old',err=1)
      read (7,nml)
1     close (7)
      open (7,file='odrdiff.nml',status='old',err=2)
      read (7,nml)
2     close (7)

      if (layout.le.1) then
         vy(5)=vy(2)
         vy(8)=vy(3)
      endif
*
* Layout = 1:
* vx() 1  2  3  4
* vy() 1  7  8  4
*
* Layout = 2:
* vx() 1  2  3  4
* vy() 1      2 3      4
*      1  5  6   7  8  4
*
      mjd92=48622-46066
      day=86400d0
      limit=999999999
      nstat=0
      tmin=1e10
      tmax=-1e10
      hmin=1e10
      hmax=-1e10
      xmin=1e10
      xmax=-1e10
      ymin=1e10
      ymax=-1e10
      title='NONE'
      dev=' '
      pi=4*atan(1d0)
      rad=pi/180

      narg=iargc()
      if (narg.lt.2) goto 1310

      do iarg=1,narg
      call getarg(iarg,arg)
      if (arg(1:4).eq.'lim=') then
	 read (arg(5:),*) delta
	 limit=nint(delta*day)
      else if (arg.eq.'-out') then
         out=.true.
      else if (arg.eq.'-mid') then
	 mid=.true.
      else if (arg.eq.'-test') then
	 test=.true.
      else if (arg(1:5).eq.'ofac=') then
	 read (arg(6:),*) ofac
      else if (arg(1:7).eq.'hifreq=') then
	 read (arg(8:),*) hifreq
      else if (arg(1:6).eq.'title=') then
	 title=arg(7:)
      else if (arg(1:4).eq.'xgf=') then
	 open (36,file=arg(5:),status='new',form='unformatted',
     |      access='direct',recl=18)
	 xgf=.true.
      else if (arg(1:5).eq.'step=') then
	 read (arg(6:),*) nstep
      else if (arg(1:7).eq.'layout=') then
	 read (arg(8:),*) layout
      else if (arg(1:4).eq.'dev=') then
	 dev=arg(5:)
      else if (file1.eq.' ') then
	 file1=arg
      else if (file2.eq.' ') then
	 file2=arg
      else
	 dev=arg
      endif
      enddo

      i=odrinfo(fd1,file1(2:),satel,irep,iarc,irem,nrec1,
     |             it10,it11,itst1,ibegin,iend,rev)
      if (nrec1.eq.0) goto 1300      
      if (i.eq.1) scale1=10

      if (file2(1:4).eq.'dir=') then
         it20=it10
         it21=it11
         itst2=itst1
         getorbit=.true.
	 file2=file2(5:)
      else
         i=odrinfo(fd2,file2(2:),satel,irep,iarc,irem,nrec2,
     |             it20,it21,itst2,ibegin,iend,rev)
         if (nrec2.eq.0) goto 1300
         if (i.eq.1) scale2=10
      endif

      if (itst1.ne.itst2) goto 1330
*     if (mid) then
*        it0=(ibegin/itst1+1)*itst1
*        it1=(iend/itst1)*itst1
*     else
*        it0=max(it10,it20)
*        it1=min(it11,it21)
*        it1=min(it1,it10+limit,it20+limit)
*     endif
      it0=max(it10,it20)
      it1=min(it11,it21)
      it1=min(it1,it10+limit,it20+limit)
      if (mid) then
 	 itstat0=(ibegin/itst1+1)*itst1
 	 itstat1=(iend/itst1)*itst1
      else
	 itstat0=it0
	 itstat1=it1
      endif
      ioff1=(it0-it10)/itst1+2
      ioff2=(it0-it20)/itst2+2
      nrec=(it1-it0)/itst1+1
      if (nrec.gt.npnt) goto 1350
      step=max(itst1,itst2)

      nrev=0
      do i=1,nrec
         read (fd1,rec=i+ioff1) rec
	 if (ltlend()) call i4swap(4,rec)
         ix=rec(3)/scale1
         iy=rec(2)/scale1
         ih=rec(4)
         if (getorbit) then
            j=getorb(DBLE(rec(1)),dlat,dlon,hgt,file2,.true.)
	    if (j.ne.0) then
	       nrec=i-1
	       itstat1=itime
	       goto 100
	    endif
            itime=rec(1)
            dh(i)=ih/1d3-hgt
	    lon(i)=nint(dlon*1d6)
	    lat(i)=nint(dlat*1d6)
         else
            itime=rec(1)
            read (fd2,rec=i+ioff2) rec
	    if (ltlend()) call i4swap(4,rec)
            if (itime.ne.rec(1)) stop
            dh(i)=(ih-rec(4))/1d3
            lon(i)=rec(3)/scale2
            lat(i)=rec(2)/scale2
         endif
         dt(i)=itime
         t(i)=itime/day-mjd92
         if (test) dh(i)=cos(2*pi*itime/rev)
         h(i)=dh(i)
	 if (itime.ge.itstat0 .and. itime.le.itstat1) then
	    if (nstat.eq.0) nstat0=i
	    tmax=max(t(i),tmax)
	    tmin=min(t(i),tmin)
	    hmax=max(h(i),hmax)
	    hmin=min(h(i),hmin)
	    hrms=hrms+h(i)**2
	    hmean=hmean+h(i)
	    dx=(ix-lon(i))/1e6
	    if (dx.gt.180) dx=dx-360
	    if (dx.lt.-180) dx=dx+360
	    xmax=max(dx,xmax)
	    xmin=min(dx,xmin)
	    xrms=xrms+dx**2
	    xmean=xmean+dx
	    dy=(iy-lat(i))/1e6
	    ymax=max(dy,ymax)
	    ymin=min(dy,ymin)
	    yrms=yrms+dy**2
	    ymean=ymean+dy
	    nstat=nstat+1
	 endif
	 if (out) write (*,660) itime,dh(i)
         if (i.gt.1 .and. lat(i).gt.0 .and. lat(i-1).le.0 .and.
     |          dt(i)-dt(i-1).le.120d0) then
            nrev=nrev+1
            tnode(nrev)=dt(i)-
     |          lat(i)/dble(lat(i)-lat(i-1))*(dt(i)-dt(i-1))
         endif
      enddo

  100 continue
      krev=nrev-1
  105 krev0=krev
      rev=(tnode(nrev)-tnode(1))/krev
      if (rev.lt.6000) rev=6000
      if (rev.gt.7200) rev=7200
      tnode0=0
      krev=0
      do i=2,nrev
         krev=krev+nint((tnode(i)-tnode(i-1))/rev)
         tnode0=tnode0+tnode(i)-krev*rev
      enddo
      tnode0=tnode0/(nrev-1)
      if (krev.ne.krev0) goto 105
*     write (*,*) nrev,krev,rev,tnode0
      hmean=hmean/nstat
      hrms=sqrt(hrms/nstat)
      hsigma=sqrt(hrms**2-hmean**2)
      xmean=xmean/nstat
      xrms=sqrt(xrms/nstat)
      ymean=ymean/nstat
      yrms=sqrt(yrms/nstat)
      call strf1985(date0,'%y%m%d %H:%M:%S',itstat0)
      call strf1985(date1,'%y%m%d %H:%M:%S',itstat1)
      if (.not.out)
     |   write (*,610) nstat,date0,date1,xmin,xmax,xmean,xrms,
     |   ymin,ymax,ymean,yrms,hmin,hmax,hmean,hrms
      l=lnblnk(file1)

      if (dev.eq.' ') goto 300

      if (pgbeg(0,dev,1,1).ne.1)
     |stop "odrdiff: error opening plot device"
      call pgsvp(0.0,1.0,0.0,1.0)
      call pgsch(0.8)
      if (title.eq.'NONE') title=file1(1:l)//' - '//file2
      call pgmtxt('T',-1.0,0.5,0.5,title)
      call pgsch(ch)
      call pgsvp(vx(1),vx(2),vy(8),vy(4))
      dy=hmax-hmin
      hmin=hmin-dy/20
      hmax=hmax+dy/20
      write (text,620) hmean*100,hrms*100
      if (layout.eq.0) then
         tmin=90
	 tmax=-90
         do i=nstat0,nstat0+nstat-1,nstep
	    t(i)=lat(i)/1e6
	    tmax=max(tmax,t(i))
	    tmin=min(tmin,t(i))
	 enddo
         call label('latitude (deg)','\\gDh (cm)',' ',
     |'\\(0791) Radial orbit difference    \\(0794)'//
     |' 1\\gs error    \\(0792) 2\\gs error',text)

         l=index(file1,'ODR')
	 file1(l:)='radsig.out'
         l=index(file2,'ODR')
	 file2(l:)='radsig.out'
	 open (11,file=file1,status='old')
	 open (12,file=file2,status='old')
	 i=0
111      i=i+1
         read (11,*,end=112) dummy,dummy,dummy,x(i)
         read (12,*,end=112) dummy,y(i),dummy,dummy
	 x(i)=sqrt(x(i)**2+dummy**2)
	 goto 111
112      continue
	 close (11)
	 close (12)

	 call pgsls(4)
         call pgswin(tmin,tmax,hmin,hmax)
         call line(nstat0,nstat,nstep,y,x,1e20,1e20)
         call pgswin(tmin,tmax,-hmin,-hmax)
         call line(nstat0,nstat,nstep,y,x,1e20,1e20)
	 call pgsls(2)
         call pgswin(tmin,tmax,hmin/2,hmax/2)
         call line(nstat0,nstat,nstep,y,x,1e20,1e20)
         call pgswin(tmin,tmax,-hmin/2,-hmax/2)
         call line(nstat0,nstat,nstep,y,x,1e20,1e20)
	 call pgsls(1)
         call pgswin(tmin,tmax,100*hmin,100*hmax)
         call pgbox('ABCNST',0.0,0,'ABCNST',0.0,0)
         call nulldot(2,text)
         call pgswin(tmin,tmax,hmin,hmax)
	 call pgslw(3)
         call line(nstat0,nstat,nstep,t,h,1e20,1e20)
	 call pgslw(1)
      else
         call pgswin(0.0,tmax-tmin,100*hmin,100*hmax)
         call label('latitude (deg)','\\gDh (cm)',' ',
     |'\\(0791) Radial orbit difference    \\(0794) 1-cpr modulation',
     |   text)
         call pgbox('ABCNST',0.0,0,'ABCNST',0.0,0)
         call nulldot(2,text)
         call pgswin(tmin,tmax,hmin,hmax)
         call line(nstat0,nstat,nstep,t,h,1e20,1e20)
      endif

      if (out) write (*,640) hrms,hrms,1
      if (layout.eq.0) goto 9999

      hifac=hifreq*(2*itst1/rev)
      call spfper(dt(nstat0),dh(nstat0),nstat,hmean,hsigma,
     |ofac,hifac,wk1,wk2,npnt,nout,jmax,prob)
      smax=0
      do i=1,nout
         x(i)=wk1(i)*rev
	 amp=2*hsigma*sqrt(wk2(i)/nstat)
	 h(i)=amp
         if (out) write (*,630) i,wk1(i)*rev,amp
	 smax=max(smax,h(i))
      enddo
      smax=smax*1.05
      call pgsvp(vx(3),vx(4),vy(3),vy(4))
      call pgswin(x(1),x(nout),0.0,100*smax)
      call pgbox('ACST',0.0,0,'ABCNST',0.0,0)
      call pgbox('BNSTI',0.0,0,' ',0.0,0)
      call label('frequency (cycl/rev)','amplitude (cm)',' ',
     |     'Periodogram',' ')
      call pgswin(x(1),x(nout),0.0,smax)
      call pgline(nout,x,h)

  300 call spcdem(nrec,dh,wk3,step/rev,nsmth)
      n=0
      if (xgf) then
	 do i=1,nrec
	    itime=nint(dt(i))
	    if (itime.ge.itstat0 .and. itime.le.itstat1) then
	       n=n+1
	       write (36,rec=n+1) itime,lat(i),lon(i),nint(dh(i)*1d6)
	    endif
	 enddo
	 write (36,rec=1) '@XGF',n
      endif
      if (dev.eq.' ') goto 9999

      hmean=0
      hrms=0
      expph0=exp(CMPLX(0d0,tnode0/rev*2*pi))
      xmin=1e20
      xmax=-1e20
      rms_sin=0
      rms_cos=0
      do i=1,nrec
         h(i)=dh(i)
         temp=wk3(i)/expph0
         if (layout.eq.3) then
            x(i)=realtemp(1)
            y(i)=realtemp(2)
            xmin=min(xmin,x(i),y(i))
            xmax=max(xmax,x(i),y(i))
            phase=(dt(i)-tnode0)/rev*2*pi
            d=realtemp(1)*cos(phase)
            rms_cos=rms_cos+d*d
            d=realtemp(2)*sin(phase)
            rms_sin=rms_sin+d*d
         else
            phase=atan2(realtemp(2),realtemp(1))/rad
*           phase=phase-tnode0/rev*360
            phase=phase-nint(phase/360)*360
            if (phase.lt.0) phase=phase+360
            y(i)=phase
            x(i)=abs(temp)
            xmax=max(xmax,x(i))
         endif
	 itime=nint(dt(i))
	 if (itime.ge.itstat0 .and. itime.le.itstat1) then
	    rmin=min(rmin,h(i))
	    rmax=max(rmax,h(i))
	    hmean=hmean+dh(i)
            hrms=hrms+dh(i)**2
	 endif
      enddo
      d=rmax-rmin
      rmin=rmin-d/20
      rmax=rmax+d/20
      hmean=hmean/nstat
      hrms=sqrt(hrms/nstat)
      if (layout.eq.1) then
         call pgsvp(vx(1),vx(2),vy(3),vy(4))
         call pgswin(tmin,tmax,hmin,hmax)
         call pgsls(4)
         call line(nstat0,nstat,nstep,t,y,1e20,1e20)
         call pgsls(1)
      else if (layout.eq.3) then
         call pgsvp(vx(1),vx(2),vy(6),vy(7))
         d=xmax-xmin
         xmin=xmin-d/20
         xmax=xmax+d/20
         call pgswin(0.0,tmax-tmin,100*xmin,100*xmax)
         call pgbox('ABCNST',0.0,0,'ABCNST',0.0,0)
         call pgswin(tmin,tmax,xmin,xmax)
         call line(1,nrec,nstep,t,y,1e20,1e20)
         call pgsls(3)
         call line(1,nrec,nstep,t,x,1e20,1e20)
         call pgsls(1)
         rms_sin=sqrt(rms_sin/nrec)
         rms_cos=sqrt(rms_cos/nrec)
         write (text,621) rms_sin*100,rms_cos*100
	 call nulldot(2,text)
*        write (*,550) text
         call label('time (days)','amplitude (cm)',' ',
     |'\\(0791) sine part       1-cpr signal       \\(0793)'//
     |' cosine part',text)
      else
         call pgsvp(vx(1),vx(2),vy(6),vy(7))
         xmin=0
         xmax=1.05*xmax
         call pgswin(0.0,tmax-tmin,100*xmin,100*xmax)
         call pgbox('BCNST',0.0,0,'ABNST',0.0,0)
         call pgswin(tmin,tmax,xmin,xmax)
         call line(nstat0,nstat,nstep,t,x,1e20,1e20)
         call pgswin(0.0,tmax-tmin,0.,360.)
         call pgbox('AST',0.0,0,'CMST',90.0,3)
         call pgswin(tmin,tmax,0.,360.)
	 call pgsls(4)
         call line(nstat0,nstat,nstep,t,y,1e20,180.)
	 call pgsls(1)
         call label('time (days)','amplitude (cm)','phase (deg)',
     |'\\(0791) amplitude       1-cpr signal       \\(0794) phase',
     |' ')
      endif

      call pgsvp(vx(3),vx(4),vy(1),vy(2))
      call hist(nstat,h(nstat0),rmin,rmax,nhist,x,y,smax)
      call pgswin(100*rmin,100*rmax,0.0,105*smax)
      call pgbox('ACST',0.0,0,'ABCNST',0.0,0)
      call pgbox('BNSTI',0.0,0,' ',0.0,0)
      call label('\\gDh (cm)','Occurance (%)',' ',
     |     'Histogram (1-cpr removed)',' ')
      call pgswin(rmin,rmax,0.0,1.05*smax)
      call pgline(nhist,x,y)

      call pgsvp(vx(1),vx(2),vy(1),vy(5))
      call pgswin(0.0,tmax-tmin,100*rmin,100*rmax)
      call pgbox('ABCNST',0.0,0,'ABCNST',0.0,0)
      write (text,620) hmean*100,hrms*100
      call nulldot(2,text)
      call label('time (days)','\\gDh (cm)',' ',
     |     'Demodulated signal (1-cpr removed)',text)
      call pgswin(tmin,tmax,rmin,rmax)
      call line(nstat0,nstat,nstep,t,h,1e20,1e20)
      goto 9999
*
* Formats
*
  550 format (a)
  610 format ('Recs   : ',i10,t40,a,' - ',a/
     |13x,'Minimum   Maximum      Mean       RMS'/
     |'Longitude:',4f10.6/'Latitude :',4f10.6/'Altitude :',4f10.3)
  620 format ('Mean =',f5.1,' cm   RMS =',f5.1,' cm')
  621 format ('sine RMS =',f5.1,' cm   cosine RMS =',f5.1,' cm')
  630 format (i5,f10.6,2f13.6)
  640 format (2f10.6,i5)
  660 format (i9,f8.3)
  670 format ('usage: odrdiff [options] ODR1 ODR2'//
     |'Options are:'/
     |'-mid          : Use only middle section of arc.'/
     |'lim=days      : Use only "days" from begin of arc.'/
     |'title="title" : Use "title" as heading of the plot,',
     |' default is "ODR1 - ODR2".'/
     |'hifreq=freq   : Use "freq" as highest frequency (cycl/rev) in',
     |' periodogram (def=3).'/
     |'step=nstep    : Plot signal sampled at interval "nstep" ',
     |'(def=10)'/
     |'dev=dev       : Use plotdevice "dev"'/
     |'layout=ind    : Use layout type "ind"'/
     |'xgf=filenm    : Write demodulated signal to XGF file')
*
* Errors
*
 1300 write (0,550) 'odrdiff: improper ODR file.'
      goto 9999
*
 1310 write (0,670)
      goto 9999
*
 1330 write (0,550) 'odrdiff: incompatible stepsizes'
      goto 9999
*
 1350 write (0,550) 'odrdiff: too many records in overlap: ',nrec
 9999 call pgend
      end

      subroutine label(text1,textl,textr,text3,text4)
      character*(*) text1,textl,textr,text3,text4
      call pgmtxt('B',2.5,0.95,1.0,text1)
      call pgmtxt('L',2.0,0.50,0.5,textl)
      call pgmtxt('R',2.5,0.50,0.5,textr)
      call pgmtxt('T',1.0,0.50,0.5,text3)
      call pgmtxt('B',2.5,0.10,0.0,text4)
      end

      subroutine hist(n,x,x0,x1,m,xhist,h,hmax)
      integer n,m,j,i,jmin,jmax
      real x(n),x0,x1,xhist(m),h(m),hmax,cutoff,dx
      cutoff=0.05*n/m
   10 hmax=0
      dx=(x1-x0)/m
      do j=1,m
	 h(j)=0
      enddo
      do i=1,n
	 j=(x(i)-x0)/dx+1
	 if (j.ge.1 .and. j.le.m) h(j)=h(j)+1
      enddo
      jmin=m
      jmax=1
      do j=1,m
	 if (h(j).gt.cutoff) then
	    jmin=min(jmin,j-15)
	    jmax=max(jmax,j+15)
	 endif
      enddo
      jmin=max(1,jmin)
      jmax=min(m,jmax)
      if (jmin.ne.1 .or. jmax.ne.m) then
	 x1=x0+jmax*dx
	 x0=x0+(jmin-1)*dx
	 goto 10
      endif
      do j=1,m
	 h(j)=h(j)/n
         xhist(j)=x0+(j-0.5)*dx
	 hmax=max(h(j),hmax)
      enddo
      end

      subroutine line(n0,n,nstep,x,y,xgap,ygap)
      integer n0,n,nstep,i
      real x(*),y(*),xgap,ygap,x0,y0
      call pgbbuf
      x0=-1e30
      y0=-1e30
      do i=n0,n0+n-1,nstep
         if (x(i)-x0.gt.xgap .or. abs(y(i)-y0).gt.ygap) then
            call pgmove(x(i),y(i))
         else
            call pgdraw(x(i),y(i))
         endif
         x0=x(i)
         y0=y(i)
      enddo
      call pgebuf
      end

      subroutine nulldot(n,text)
      integer n,i,l
      character*(*) text
      do i=1,n
         l=index(text,' -.')
         if (l.gt.0) text(l:l+2)='-0.'
         l=index(text,' .')
         if (l.gt.0) text(l:l+1)='0.'
      enddo
      end
