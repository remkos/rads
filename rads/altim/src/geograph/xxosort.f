      program xxosort

* Sort xovers (in xxo-file) by latitude and longitude, and compute
* mean and variance at each location, to be stored in output XGF file.
*
* Syntax:
* xxosort [ options ]  xxo-file(s)  xgf-file
*
*-
* 30-Jul-1996 - Created by Remko Scharroo
* 22-Apr-1997 - Add XGF support
* 15-Aug-1997 - Add t= flag
* 30-Mar-1998 - Complete new design. Make latitude seperation based on
*    time difference, rather than latitude. Store data more efficiently.
* 12-Apr-2000 - Added DATEARG function
* 22-Dec-2000 - Extended with with trend estimation
*  9-Jul-2002 - Ensured initialisation at zero
*-----------------------------------------------------------------------
      implicit none

      character arg*80,spec*4,filenm*80,date(3)*15
      integer*4	i,ii,j,k,l,m,nmax,nmin/2/,ntot,nread,nrev,
     |		irec,nrec,nused,npar,recl,iarg,iargc,nargs,kx,ky
      integer*4 mrev,maxloc
      parameter (mrev=501,maxloc=100000)
      integer*4 nloc,loc(mrev,mrev),ny(mrev)
      integer*4 n(maxloc)
      integer*2 loci(maxloc),locj(maxloc)
      ! mpar = number of parameters to solve per grid point
      ! mata = size of AtA matrix per grid point
      ! mstat= number of statistic entries
      ! mtot = total size of work array per grid point
      integer*4 mpar,mata,mstat,mtot
      parameter (mpar=6,mata=(mpar+1)*mpar/2,mstat=3,
     |		mtot=mpar+mata+mstat)
      integer*4 tmp(5),tmin,tmax
      real*8	pi,rad,lat,latc,r,d,rev,sini,u,tani,lon,lon0,dlon/1d0/,
     |		z,z1,z2/0d0/,period,dum,tmid,
     |		year,tmean,trms,tvar,xmax/60d-2/,edit/3.5d0/,incl,
     |		zlo,zhi,vhi,tbias(2)/2*0d0/,phase,bias(2)/2*0d0/
      real*8	y(mrev),zmean(maxloc),zvar(maxloc),
     |		zmin(maxloc),zmax(maxloc),
     |		b(mtot,0:359,-90:89),v(mpar),w(mpar),a(mtot,maxloc),wgt,
     |		t0/-1d40/,t1/1d40/
      integer*4 h(mpar),if,delt0/0/,delt1/0/
      logical   tsort/.false./,flip/.false./,datearg,diffstat/.false./
      namelist  /xxosort_nml/ period,nrev,incl,lon0

* File handling

      integer*4 openf,readf,closef,ios,fd/-1/,fdaux/-1/,fddlt/-1/
      integer*4 i4(6),ssh(2),orb(2)/2*800000000/,xgf4(4),xxf4(16)
      integer*4 maux,naux/0/,fldn/0/,fld/-1/
      parameter (maux=20)
      integer*2 aux(maux)/maux*0/,dlt(2)/2*0/,xgf2(9),it(2)
      equivalence (xgf2,xgf4)
      equivalence (i4,xxf4(1)),(it,xxf4(7)),(ssh,xxf4(8)),(orb,xxf4(10))

* Initialise

      pi=4*atan(1d0)
      rad=pi/180
      year=86400d0*365.25d0

* Read namelists

      arg='/user/altim'
      call checkenv('ALTIM',arg,l)
      arg(l+1:)='/nml/xxosort.nml'
      open (7,file=arg,status='old')
      read (7,xxosort_nml)
      close (7)
      open (7,file="xxosort.nml",status='old',err=2)
      read (7,xxosort_nml)
      close (7)
2     continue

* More initialisation

      nread=0
      nmax=0
      ntot=0
      nloc=0
      tmid=0
      tmin=2147483647
      tmax=-2147483647-1

      nargs=iargc()
      if (nargs.le.1) then
         write (*,600)
         stop
      endif

      call matsy1(mpar,h)

      do l=1,maxloc
         n(l)=0
	 zmean(l)=0
	 zvar(l)=0
	 do k=1,mtot
	    a(k,l)=0
	 enddo
      enddo
      do j=1,mrev
	 ny(j)=0
	 y(j)=0
         do ii=1,mrev
	    loc(ii,j)=0
	 enddo
      enddo
      do ky=-90,89
         do kx=0,359
	    do k=1,mtot
	       b(k,kx,ky)=0
	    enddo
	 enddo
      enddo
      v(1)=1d0
      w(1)=1d0
      open(90,status='scratch',
     |  recl=20,access='direct',form='unformatted')

* Scan command line

      call getarg(nargs,filenm)

      do iarg=1,nargs-1

      call getarg(iarg,arg)
      if (arg(:4).eq.'nml=') then
	 open (7,file=arg(5:))
	 read (7,xxosort_nml)
	 close (7)
      else if (arg(:5).eq.'nmin=') then
         read (arg(6:),*) nmin
      else if (arg(:5).eq.'xmax=') then
         read (arg(6:),*) xmax
         xmax=xmax/100
      else if (arg(:5).eq.'edit=') then
         read (arg(6:),*) edit
      else if (arg(:4).eq.'-noD') then
	 delt0=284060346
	 delt1=292607918
      else if (arg(:2).eq.'-t') then
	 tsort=.true.
      else if (arg(:2).eq.'-d') then
	 diffstat=.true.
      else if (arg.eq.'tbias=0') then
	 tbias(1)=0d0
	 tbias(2)=0d0
	 fdaux=0
      else if (arg(:6).eq.'tbias=') then
	 tbias(2)=1d30
         read (arg(7:),*,iostat=ios) tbias
	 if (tbias(2).gt.1d20) tbias(2)=tbias(1)
* tbias in milliseconds. Positive when tags are late.
         fdaux=0
      else if (arg(:5).eq.'bias=') then
	 bias(2)=1d30
         read (arg(6:),*,iostat=ios) bias
	 if (bias(2).gt.1d20) bias(2)=bias(1)
      else if (arg(:4).eq.'fld=') then
	 read (arg(5:),*) fld
         write (*,550) '(auxiliary field statistics) '
	 fld=(fld-1)*2
	 fdaux=0
      else if (datearg(arg,t0,t1,dum)) then
      else

* Initialize orbit related parameters

	 rev=period*86400d0/nrev
	 dlon=180d0/nrev
	 sini=sin(incl*rad)
	 tani=tan(incl*rad)

* Process input XXO one by one
	 
	 fddlt=-1
12       fd=openf(arg,'r')
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
	    npar=0
	    recl=44
         else if (spec.eq.'xXXO') then
            recl=44+npar*4
         else
	    write (*,*) 'xxosort: unknown file type: ',spec
	    stop
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

* Read data records

         do irec=1,nrec
	    ios=readf(fd,recl,xxf4)
            if (fdaux.ge.0) ios=readf(fdaux,naux*2,aux)
	    if (fld.ge.0) then
	       ssh(1)=aux(fld+1)*1000
	       ssh(2)=aux(fld+2)*1000
	    else if (fddlt.ge.0) then
               ios=readf(fddlt,4,dlt)
	       if (dlt(1).eq.-32768 .or. dlt(2).eq.-32768) goto 100
	       ssh(1)=ssh(1)+dlt(1)*1000
	       ssh(2)=ssh(2)+dlt(2)*1000
	    endif

* Do some initial editing on epoch

	    if (i4(3).gt.delt0 .and. i4(3).lt.delt1) goto 100
	    if (i4(5).gt.delt0 .and. i4(5).lt.delt1) goto 100
	    if (i4(3).lt.t0 .or. i4(3).gt.t1) goto 100
	    if (i4(5).lt.t0 .or. i4(5).gt.t1) goto 100
	    lat=i4(1)/1d6
	    lon=i4(2)/1d6

* Flip means des-asc is needed because of the track-wise sorting
	    
	    flip=(tsort .and. it(1).gt.it(2))

* Convert to metres and remove excessively high xover differences

	    if (flip) then
	       z1=(ssh(1)-tbias(2)*aux(fldn+1))/1d6+bias(2)/1d3
	       z2=(ssh(2)-tbias(1)*aux(fldn+2))/1d6+bias(1)/1d3
	       z=z2-z1
	    else
	       z1=(ssh(1)-tbias(1)*aux(fldn+1))/1d6+bias(1)/1d3
	       z2=(ssh(2)-tbias(2)*aux(fldn+2))/1d6+bias(2)/1d3
	       z=z1-z2
	    endif
	    if (abs(z).gt.xmax) goto 100

* Track-sort is usually used for dual-satellite xovers. Thus we compute the
* latitude index j from the latitude.
 
	    if (tsort) then
	       call geocen(lat*rad,orb(1)/1d3,latc,r)
 	       u=asin(sin(latc)/sini)
	       d=asin(tan(latc)/tani)-u/501*35
*              write (*,'(2f12.6,2i12)') lat,d/rad*501/180,
*    |		nint((lon-lon0)/dlon)
	       j=-nint(d/rad*501/180)
	    else

* For single-satellite xovers it is much more efficient to look at the time
* difference to resolve the latitude index j.

	       j=nint((i4(3)-i4(5))/rev-0.5)
	       j=mod(j,nrev)
	    endif
	    if (j.lt.1) j=j+nrev
	    y(j)=y(j)+lat
	    ny(j)=ny(j)+1

* The longitude index i comes from the seperation from the first nodal longitude
* lon0. It is then mapped to ii ranging from 1 to nrev.

            i=nint((lon-lon0)/dlon)
	    if (i.lt.-nrev) i=i+nrev*2
	    if (i.ge. nrev) i=i-nrev*2
	    if (i.lt.-nrev) goto 100
	    ii=(i+nrev)/2+1
*	    write (*,*) i,ii,j

* loc is a memory mapping function. If gives the location of the xover position
* (ii,j) in one-dimentional arrays n, zvar, zmean, etc.
* loci and locj map back to i and j

	    if (loc(ii,j).eq.0) then
	       nloc=nloc+1
	       if (nloc.gt.maxloc) stop "xxosort: too many xover locations"
	       loc(ii,j)=nloc
	       loci(nloc)=i
	       locj(nloc)=j
*	    else
*	       write (*,*) nread,ii,j,loc(ii,j),lat,lon
	    endif
	    l=loc(ii,j)
	    n(l)=n(l)+1
	    zmean(l)=zmean(l)+z
	    zvar (l)=zvar (l)+z*z
	    nmax=max(nmax,n(l))
*	    if (nread.gt.maxpnt) stop "xxosort: too many xovers"
	    if (flip) l=-l

* Location, xover heights and epochs are saved
	    
	    nread=nread+1
	    tmp(1)=l
	    tmp(2)=nint(z1*1d6)
	    tmp(3)=nint(z2*1d6)
	    tmp(4)=i4(3)
	    tmp(5)=i4(5)
	    tmid=tmid+tmp(4)+tmp(5)
	    tmin=min(tmin,tmp(4),tmp(5))
	    tmax=max(tmax,tmp(4),tmp(5))
	    write(90,rec=nread)tmp
100	    continue
         enddo
         ntot=ntot+nrec

         if (fd.gt.0) ios=closef(fd)
         if (fdaux.gt.0) ios=closef(fdaux)
         if (fddlt.gt.0) ios=closef(fddlt)

      endif

      enddo

      tmid=tmid/nread/2
      call strf1985(date(1),'%y%m%d %H:%M:%S',tmin)
      call strf1985(date(2),'%y%m%d %H:%M:%S',tmax)
      call strf1985(date(3),'%y%m%d %H:%M:%S',nint(tmid))
      write (*,610) ntot,nread,nloc,date,nmin,xmax,edit

* All xovers are read
* Compute statistics per location

      nused=0
      nrec =0
      tmean=0
      trms =0
      tvar =0
      do j=1,nrev
	 if (ny(j).gt.0) then
	    y(j)=y(j)/ny(j)
	 else
	    y(j)=1d30
	 endif   
      enddo
      do l=1,nloc
	 if (n(l).ge.nmin) then

* Determine location, mean and variance

	    zmean(l)=zmean(l)/n(l)
	    zvar (l)=sqrt((zvar(l)-n(l)*zmean(l)**2)/(n(l)-1))
	    zmin (l)=zmean(l)-edit*zvar(l)
	    zmax (l)=zmean(l)+edit*zvar(l)

* Accumulate statistics

	    nrec =nrec +1
	    nused=nused+n(l)
	    tmean=tmean+zmean(l)
	    trms =trms +zmean(l)**2
	    tvar =tvar +zvar (l)**2

	 else

* Make sure points with n < nmin are removed

	    zmin (l)=1
	    zmax (l)=0

	 endif

* Reset statistics for editing pass

	 zmean(l)=0
	 zvar (l)=0
	 n    (l)=0

      enddo

      tmean=tmean/nrec
      trms =sqrt(trms/nrec-tmean**2)
      tvar =sqrt(tvar/nrec)
      write (*,620) 'Before editing',nused,nrec,
     |nused/dble(nrec),nmax,
     |tmean*100,trms*100,tvar*100

* Determine bounds on mean height and variance

      zhi  =tmean+edit*trms
      zlo  =tmean-edit*trms
      vhi  =      edit*tvar

* Start editing. Loop through all measurements.

      do irec=1,nused
*	 l =xol (irec)
*	 z1=xoz1(irec)/1d6
*	 z2=xoz2(irec)/1d6
	 read(90,rec=irec)tmp
	 l=tmp(1)
	 z1=tmp(2)/1d6
	 z2=tmp(3)/1d6
	 z=z1-z2
	 if (l.lt.0) then
	    l=-l
	    z=-z
	 endif

	 if (z.ge.zmin(l) .and. z.le.zmax(l)) then

	    n(l)=n(l)+1
	    zmean(l)=zmean(l)+z
	    zvar (l)=zvar (l)+z*z

* Build matrix for computing annual and semi-annual cycle.
* Store AtA and AtR in one vector of size mata+mpar.
*      1 to mata     : AtA
* mata+1 to mata+mpar: AtR
*         mata+mpar+1: height difference (**2) = zvar
*         mata+mpar+2: height residuals (**2)
*         mata+mpar+3: height difference (mean) = zmean
*
* Distinguish between time-sorted xovers (tsort=true) and
* ascending-descendig xovers (tsort=false).
* However, to compute (semi)annual cycles in asc-des differences
* use diffstat=true.

	    if (tsort .or. diffstat) then

	    phase=(tmp(4)+tmp(5))/year*pi
	    v(2)=cos(  phase)
	    v(3)=sin(  phase)
	    v(4)=cos(2*phase)
	    v(5)=sin(2*phase)
	    v(6)=((tmp(4)+tmp(5))/2d0-tmid)/year
	    do i=1,mpar
	       do k=1,i
		  m=h(i)+k
		  a(m,l)=a(m,l)+v(i)*v(k)
	       enddo
	    enddo
	    do i=1,mpar
	       a(i+mata,l)=a(i+mata,l)+v(i)*z
	    enddo
	    a(mtot-2,l)=a(mtot-2,l)+z*z
	    a(mtot-1,l)=a(mtot-1,l)+z*z
	    a(mtot  ,l)=a(mtot  ,l)+z

	    else

* Ascending pass

	    phase=tmp(4)/year*2*pi
	    v(2)=cos(  phase)
	    v(3)=sin(  phase)
	    v(4)=cos(2*phase)
	    v(5)=sin(2*phase)
	    v(6)=(tmp(4)-tmid)/year

* Descending pass

	    phase=tmp(5)/year*2*pi
	    w(2)=cos(  phase)
	    w(3)=sin(  phase)
	    w(4)=cos(2*phase)
	    w(5)=sin(2*phase)
	    w(6)=(tmp(5)-tmid)/year

	    do i=1,mpar
	       do k=1,i
		  m=h(i)+k
		  a(m,l)=a(m,l)+v(i)*v(k)+w(i)*w(k)
	       enddo
	    enddo
	    do i=1,mpar
	       a(i+mata,l)=a(i+mata,l)+v(i)*z1+w(i)*z2
	    enddo
	    a(mtot-2,l)=a(mtot-2,l)+z*z*2
	    a(mtot-1,l)=a(mtot-1,l)+z1*z1+z2*z2
	    a(mtot  ,l)=a(mtot  ,l)+z*2

	    endif
	 endif
      enddo
      close(90)

* Spread the AtA/AtR matrices over a regular grid

      do l=1,nloc
	 if (a(1,l).gt.0) then
	    lat=y(locj(l))
	    lon=lon0+loci(l)*dlon
*	    write (*,*) l,locj(l),loci(l),ny(locj(l)),a(1,l),lat,lon
	    call spreadn(lon-0.5d0,lat-0.5d0,mtot,a(1,l),b)
	 endif
      enddo

* Determine edit level on total weight

      if (tsort .or. diffstat) then
         edit=nmin
      else
         edit=nmin*2
      endif

* Solve annual and semi annual cycle and reduce variability.

      do kx=0,359
	 do ky=-90,89
	    wgt=b(1,kx,ky)
	    if (wgt.gt.edit) then

* Save total weight for this point (wgt), and the AtR matrix (w)

	       do k=1,mpar
		  w(k)=b(k+mata,kx,ky)
	       enddo

* Solve the system AtWA x = AtWR

	       call dpptrf('U',mpar,b(1,kx,ky),if)
	       call dpptrs('U',mpar,1,b(1,kx,ky),b(mata+1,kx,ky),mpar,if)

* mtot-2 becomes mesoscale variability = rms xover diff / sqrt(2)
* mtot-1 becomes apriori variance = sqrt (apriori rms **2 - apost mean**2)
* mtot   becomes mean xover difference
* mata   becomes aposteriori RMS ( computed from apriori )
* mata-1 becomes weight

	       b(mtot-2,kx,ky)=sqrt(b(mtot-2,kx,ky)/wgt)
	       r=b(mtot-1,kx,ky)
	       b(mtot-1,kx,ky)=sqrt(r/wgt-b(mata+1,kx,ky)**2)
	       b(mtot  ,kx,ky)=b(mtot  ,kx,ky)/wgt
	       do k=1,mpar
		  r=r-w(k)*b(k+mata,kx,ky)
	       enddo
	       b(mata,kx,ky)=sqrt(r/wgt)
	       b(mata-1,kx,ky)=wgt
	    else
	       do k=1,mtot
		  b(k,kx,ky)=1d35
	       enddo
	    endif
	 enddo
      enddo

* Write cycles as grids

      write (*,550) '** Grid statistics **'
      l=index(filenm,' ')-1
      call gridout(filenm(:l)//'_mean.grd',mtot,b(mata+1,0,-90),1d1)
      call gridout(filenm(:l)//'_an_c.grd',mtot,b(mata+2,0,-90),1d1)
      call gridout(filenm(:l)//'_an_s.grd',mtot,b(mata+3,0,-90),1d1)
      call gridout(filenm(:l)//'_sa_c.grd',mtot,b(mata+4,0,-90),1d1)
      call gridout(filenm(:l)//'_sa_s.grd',mtot,b(mata+5,0,-90),1d1)
      call gridout(filenm(:l)//'_trnd.grd',mtot,b(mata+6,0,-90),1d1)
      call gridout(filenm(:l)//'_mvar.grd',mtot,b(mtot-2,0,-90),1d1)
      call gridout(filenm(:l)//'_avar.grd',mtot,b(mtot-1,0,-90),1d1)
      call gridout(filenm(:l)//'_diff.grd',mtot,b(mtot  ,0,-90),1d1)
      call gridout(filenm(:l)//'_pvar.grd',mtot,b(mata  ,0,-90),1d1)
      call gridout(filenm(:l)//'_wght.grd',mtot,b(mata-1,0,-90),1d6)

* Open output XGF for writing

      open (20,file=filenm(:l)//'.xgf',
     |		status='new',form='unformatted',
     |		access='direct',recl=18)

* Write out XGF file
* Lines are sorted by latitude

      nmax =0
      nrec =0
      nused=0
      tmean=0
      trms =0
      tvar =0
      do j=1,nrev
         ny(j)=j
      enddo
      call dqsort(ny,y,nrev)
      do k=1,nrev
	 j=ny(k)
	 lat=y(j)
         do ii=1,nrev
	    l=loc(ii,j)
*	    write(*,*) ii,j,l
	    if (l.gt.0) then
	    lon=lon0+loci(l)*dlon

* Determine mean and variance

	    zmean(l)=zmean(l)/n(l)
	    zvar (l)=sqrt((zvar(l)-n(l)*zmean(l)**2)/(n(l)-1))

	    if (n(l).ge.nmin .and.
     |	       zmean(l).ge.zlo .and. zmean(l).le.zhi .and.
     |	       zvar (l).le.vhi) then

*		write (*,'(4i6,4f12.6)')
*     |		ii,i,j,n(l),lat,lon,
*     |		zmean(l),zvar(l)

* Write XGF record

	       nrec=nrec+1
	       nmax=max(nmax,n(l))
	       xgf2(1)=j
	       xgf2(2)=n(l)
	       xgf4(2)=nint(lat*1d6)
	       xgf4(3)=nint(lon*1d6)
	       xgf4(4)=nint(zmean(l)*1d6)
	       xgf2(9)=nint(zvar(l)*1d3)
	       write (20,rec=nrec+1) xgf2

* Accumulate statistics

   	       nused=nused+n(l)
	       tmean=tmean+zmean(l)
	       trms =trms +zmean(l)**2
	       tvar =tvar +zvar (l)**2

	    endif
	    endif
         enddo
      enddo

* Write XGF header and close file

      write (20,rec=1) '@XGF',nrec,0,0,'Mn'
      close (20)

* Write overall statistics

      tmean=tmean/nrec
      trms =sqrt(trms/nrec-tmean**2)
      tvar =sqrt(tvar/nrec)
      write (*,620) 'After editing',nused,nrec,
     |nused/dble(nrec),nmax,
     |tmean*100,trms*100,tvar*100

550   format(a)
600   format(
     |'xxosort: Sort xover file (XXO), compute mean and var at'/
     |'         each location and write XGF'//
     |'Syntax:  xxosort [ options ] xxo-file(s) prefix'//
     |'with [ options ] :'/
     |'nmin=n   : specify minimum number of xovers per location',
     |' (def:2)'/
     |'xmax=x   : specify maximum xover difference (cm) (def:60)'/
     |'tbias=tb1(,tb2) : apply time tag bias (ms) (requires AUX file)'/
     |'bias=basc(,bdes): add bias (cm) to as/descending ssh'/
     |'edit=x   : specify edit multiplier (def:3.5)'/
     |'t=t0,t1  : select time interval (t0,t1)'/
     |'         ... or use mjd=, doy=, ymd=, sec='/
     |'-t       : do a low-high tracknr difference (ERS-TOPEX)'/
     |'-d       : compute annual variations in xover differences'/
     |'-noD	: exclude ERS Phase D'//
     |'Output files will start with ''prefix''')
610   format(
     |'Number of input xovers      :',i9/
     |'Number of xovers processed  :',i9/
     |'Number of xover locations   :',i9/
     |'Time of first measurement   : ',a/
     |'Time of last measurement    : ',a/
     |'Mean epoch of measurements  : ',a/
     |'** Edit criteria **'/
     |'Min nr of xovers/location   :',i9/
     |'Max  xover difference  (cm) :',f9.3/
     |'Edit level                  :',f9.3)
620   format(
     |'** ',a,' **'/
     |'Number of xovers used       :',i9/
     |'Number of xover locations   :',i9/
     |'Avg nr of xovers/location   :',f9.3/
     |'Max nr of xovers/location   :',i9/
     |'Mean of mean xover dif (cm) :',f9.3/
     |'Var  of mean xover dif (cm) :',f9.3/
     |'Xover variance         (cm) :',f9.3)
      end

      subroutine gridout(filenm,mtot,z,zlim)
      implicit none
      integer mtot
      character*(*) filenm
      integer*4 kx,ky,n
      real*8    z(mtot,0:359,-90:89),zlim
      real*4    a(-180:179,-90:89)
      real*4    zmin,zmax,zmean,zrms,zsigma

      do ky=-90,89
         do kx=0,179
            a(kx,ky)=z(1,kx,ky)
         enddo
         do kx=-180,-1
            a(kx,ky)=z(1,kx+360,ky)
         enddo
	 do kx=-180,179
	    if (abs(a(kx,ky)).gt.zlim) a(kx,ky)=1e30
	 enddo
      enddo
      call gridwr4(filenm,360,180,a,360,-179.5,179.5,-89.5,89.5)
      call grstat4(a,360,180,360,-89.5,89.5,3.5,
     |n,zmin,zmax,zmean,zrms,zsigma)
      write (*,600) filenm,n,zmin*100,zmax*100,
     |zmean*100,zrms*100,zsigma*100
600   format(a30,i9,5f9.3)
      end

      subroutine spreadn(x,y,mtot,a,b)
      implicit none
      integer mtot
      real*8 x,y,a(mtot),b(mtot,0:359,-90:89),w,sig2
      parameter (sig2=1.5d0**2)
      integer dx,dy,kx,ky,k
      do dy=-3,3
         ky=nint(y+dy)
         do dx=-3,3
            kx=nint(x+dx)
            w=exp(-((x-kx)**2+(y-ky)**2)/sig2)
            if (kx.lt.0) kx=kx+360
            if (kx.gt.359) kx=kx-360
            do k=1,mtot
               b(k,kx,ky)=b(k,kx,ky)+a(k)*w
            enddo
         enddo
      enddo
      end
