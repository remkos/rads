      program adr2xgf
*
* Programs adr2xgf and adr2adr convert ADR files to XGF or makes changes
* to ADR files.
*-
*  3-Nov-1999 - Added DATEARG function
* 17-Jun-2000 - Introduced asc2adr
* 28-Jul-2000 - Bug fixed that caused zero-records when using asc2adr;
*               Introduced fit= flag
* 21-Sep-2000 - Apply tbias correctly when no orbit change is done
* 15-Nov-2000 - Add -topo option to place depth in orbit field
* 14-Mar-2001 - Do not print min,max,etc for files without selected data
*  3-Dec-2001 - Introduced field= flag
* 10-May-2002 - Included satids for T/P, Jason-1 and Envisat (ED)
* 20-Jan-2003 - Included additional initializations for Intel compiler (ED)
* 17-Jun-2003 - Special provisions for stdin/stdout taken over by FASTIO
* 22-Nov-2004 - Removed all ocean tide routines
*----------------------------------------------------------------------
      implicit none
      integer*4    mn,mf
      parameter    (mn=61,mf=400)
      character*80 filenm(mf),arg,spec*4,satel*8/' '/,dir,gridnm,
     |		   header*28,line*80
      integer*4    rec(4),adr4(7),latpnt(mn),lonpnt(mn),orbpnt(mn),
     |		   if,nf,getorb,mjd,mjd85,mrec,nrec,irec,irec0,nr,nt,
     |             narg,iargc,iarg,mode,i,j,l,lnblnk,ix,iy,tbias,
     |		   nx,ny,mx,istep,nrmpnt,ipnt/0/,nrnrmpnt,satid,field/1/,
     |             openf,readf,fdin,fdout,xgflen/18/,
     |		   iadrlen,oadrlen,ierr/0/
      parameter    (mx=4321*2161)
      integer*2    xgf2(9),adr2(14),bound(4),naux,maux/0/,null(14)
      real*4       grid(mx),glon0,glon1,glat0,glat1,x,y,gr3int4
      real*8	   hpnt(mn),rpnt(mn),tpnt(0:mn),lat,lon,maxvar/0/,
     |             r,rmin,rmax,rrms,rmean,ssscale/1d0/,
     |		   t,tmin,tmax,trms,tmean,hrms,ssvar,sigma0,
     |		   rmsmax,signrm,orb,orbdiff/0d0/,dx/1d0/,dy/1d0/,dum,tpass,
     |		   sec,fact,t0,t1,h0/-1d40/,h1/1d40/,
     |		   lon0/-180d0/,lon1/180d0/,lat0/-90d0/,lat1/90d0/,
     |		   hfact/1d0/,gfact/-1d0/,bias/0d0/,tmjd,
     |		   topo0/1d30/,topo1/-1d30/,fit(4)/4*0d0/,rad
      logical      ref/.false./,error,datearg,headerskipped/.false./,
     |		   orbin/.false./,orbout/.false./,orbdiffout/.false./,
     |		   useptide/.false./,noaux/.false./,
     |		   usefit/.false./,addtopo/.false./
      equivalence (rec,xgf2),(adr4,adr2)
      equivalence (spec,header(1:4)),(nrec,header(21:24))
      equivalence (bound,header(13:20)),(satel,header(5:12))
      equivalence (naux,header(25:26))
      parameter   (mjd85=46066)
      integer*4	   minint4,maxint4
      parameter   (maxint4=2147483647,minint4=-maxint4-1)

* Check output mode

      narg=iargc()
      call getarg(0,arg)
      if (index(arg,'adr2adr').gt.0) then
	 mode=1
      else if (index(arg,'adr2xgf').gt.0) then
	 mode=2
      else if (index(arg,'adrlist').gt.0) then
	 mode=3
      else if (index(arg,'asc2adr').gt.0) then
         mode=4
      else
	 mode=0
      endif

* Initialize

      tbias=0
      nrec=0
      mrec=0
      istep=1
      nrmpnt=0
      tmin=1d40
      tmax=-1d40
      t0=minint4
      t1=maxint4
      trms=0
      tmean=0
      nt=0
      nf=0
      irec0=1
      rad=atan(1d0)/45
      do i=1,7
         adr4(i)=0
      enddo

* Scan argument list

      do iarg=1,narg
         call getarg(iarg,arg)
	 if (arg(:5).eq.'step=') then
	    read (arg(6:),*) istep
	    write (0,*) 'Step = ',istep
	 else if (arg(:6).eq.'field=') then
	    read (arg(7:),*) field
	 else if (arg(:4).eq.'nrm=') then
	    read (arg(5:),*) nrmpnt,rmsmax,signrm
	    nrmpnt=min(nrmpnt/2*2+1,mn)
	    write (0,610) nrmpnt,rmsmax,signrm
	 else if (arg(:4).eq.'ref=') then
	    gridnm=arg(5:)
	    j=lnblnk(gridnm)
	    write (0,*) 'Ref grid = '//gridnm(:j)
	    nx=0
	    ny=mx
	    call gridrd4(gridnm,nx,ny,grid,
     |		glon0,glon1,glat0,glat1,rmin,rmax)
	    dx=(glon1-glon0)/(nx-1)
	    dy=(glat1-glat0)/(ny-1)
	    ref=.true.
	 else if (arg(:5).eq.'-topo') then
	    addtopo=.true.
	    gridnm='/user/altim/data/world/etopo5.grd'
	    nx=0
	    ny=mx
	    call gridrd4(gridnm,nx,ny,grid,
     |		glon0,glon1,glat0,glat1,rmin,rmax)
	    dx=(glon1-glon0)/(nx-1)
	    dy=(glat1-glat0)/(ny-1)
	 else if (arg(:6).eq.'depth=') then
	    read (arg(7:),*) topo1,topo0
	    write (0,*) 'Depth = ',topo1,topo0
	    topo0=-topo0
	    topo1=-topo1
	    gridnm='/user/altim/data/world/etopo5.grd'
	    nx=0
	    ny=mx
	    call gridrd4(gridnm,nx,ny,grid,
     |		glon0,glon1,glat0,glat1,rmin,rmax)
	    dx=(glon1-glon0)/(nx-1)
	    dy=(glat1-glat0)/(ny-1)

	 else if (arg(:5).eq.'topo=') then
	    read (arg(6:),*) topo0,topo1
	    write (0,*) 'Topo heights = ',topo0,topo1
	    gridnm='/user/altim/data/world/etopo5.grd'
	    nx=0
	    ny=mx
	    call gridrd4(gridnm,nx,ny,grid,
     |		glon0,glon1,glat0,glat1,rmin,rmax)
	    dx=(glon1-glon0)/(nx-1)
	    dy=(glat1-glat0)/(ny-1)
	 else if (arg(:2).eq.'-a') then
	    ierr=-2
	 else if (arg(:2).eq.'-x') then
	    noaux=.true.
	 else if (arg(:2).eq.'-v') then
	    nx=0
	    ny=mx
	    gridnm='/user/altim'
	    call checkenv('ALTIM',gridnm,j)
	    gridnm(j+1:)='/data/world/ssvar95a.grd'
	    call gridrd4(gridnm,nx,ny,grid,
     |		glon0,glon1,glat0,glat1,rmin,rmax)
	    dx=(glon1-glon0)/(nx-1)
	    dy=(glat1-glat0)/(ny-1)
	    read (arg(3:),*,iostat=i) maxvar
	    if (maxvar.eq.0) maxvar=1d20
	 else if (arg(:6).eq.'ssvar=') then
	    nx=0
	    ny=mx
	    i=index(arg,'*')
	    if (i.gt.0) then
	       read (arg(7:i-1),*) ssscale
	    else
	       i=6
	    endif
	    call gridrd4(arg(i+1:),nx,ny,grid,
     |		glon0,glon1,glat0,glat1,rmin,rmax)
	    dx=(glon1-glon0)/(nx-1)
	    dy=(glat1-glat0)/(ny-1)
	    maxvar=1d20
	 else if (arg(:7).eq.'sigma0=') then
	    read (arg(8:),*) sigma0
	    sigma0=(sigma0/1d2)**2
	 else if (arg(:6).eq.'gfact=') then
	    read (arg(7:),*) gfact
	 else if (arg(:6).eq.'hfact=') then
	    read (arg(7:),*) hfact
	 else if (datearg(arg,t0,t1,dum)) then
	 else if (arg(:2).eq.'h=') then
	    read (arg(3:),*) h0,h1
	    h0=h0/100
	    h1=h1/100
	 else if (arg(:4).eq.'orb=') then
	    dir=arg(5:)
	    orbin=.true.
	 else if (arg(:5).eq.'bias=') then
	    read (arg(6:),*) bias
	    bias=bias/100
	 else if (arg(:4).eq.'fit=') then
	    read (arg(5:),*) fit
	    usefit=.true.
	 else if (arg(:6).eq.'tbias=') then
	    read (arg(7:),*) t
	    tbias=nint(t*1d3)
	 else if (arg(:4).eq.'lon=') then
	    read (arg(5:),*) lon0,lon1
	 else if (arg(:4).eq.'lat=') then
	    read (arg(5:),*) lat0,lat1
	 else if (arg(:8).eq.'-orbdiff') then
	    orbdiffout=.true.
	 else if (arg(:4).eq.'-orb') then
	    orbout=.true.
	 else if (arg(:3).eq.'-pt') then
	    useptide=.true.
	 else if (arg(:4).eq.'sat=') then
	    satel=arg(5:)
	    if (satel.eq.'e1' .or. satel.eq.'ers1') satel='ERS-1'
	    if (satel.eq.'e2' .or. satel.eq.'ers2') satel='ERS-2'
	    if (satel.eq.'tx' .or. satel.eq.'tpx' .or.
     |		satel.eq.'tp') satel='TOPEX'
	    if (satel.eq.'pn' .or. satel.eq.'pos') satel='POSEIDON'
	    if (satel.eq.'j1' .or. satel.eq.'jason1') satel='JASON-1'
	    if (satel.eq.'n1' .or. satel(:4).eq.'envi') satel='ENVISAT1'
	 else
	    nf=nf+1
	    if (nf.gt.mf) call fin('too many input files')
	    filenm(nf)=arg
	 endif
      enddo

* Help !

      if (mode.eq.4 .and. satel.eq.' ') then
         write (0,626)
	 nf=0
      endif
      if (nf.lt.1) then
	 if (mode.eq.1) then
	    write (0,620)
	    write (0,623)
	 else if (mode.eq.2) then
	    write (0,621)
	    write (0,623)
	 else if (mode.eq.3) then
	    write (0,624)
	    write (0,623)
	 else if (mode.eq.4) then
	    write (0,625)
	    write (0,623)
	 else
	    write (0,622)
	    write (0,623)
	 endif
	 goto 9999
      endif

* Open output file

      if (mode.eq.1 .or. mode.eq.2 .or. mode.eq.4) then
         fdout=openf(filenm(nf),'w')
         nf=nf-1
      endif

* Process ADR or ASCII files

      do if=1,nf

	 fdin=openf(filenm(if),'r')
*         if(fdin.lt.0) call perrorf("FASTIO ")
         if(fdin.lt.0) stop
	 if (mode.ne.4) l=readf(fdin,24,header)

         if (satel.eq.'ERS-1') then
            satid=9105001
         else if (satel.eq.'ERS-2') then
            satid=9502101
         else if (satel.eq.'TOPEX') then
            satid=9205201
         else if (satel.eq.'POSEIDON') then
            satid=9205201
         else if (satel.eq.'JASON-1') then
            satid=0105501
         else if (satel.eq.'ENVISAT1') then
            satid=0200901
         endif

* Check header
	    
50	 if (mode.eq.4) then
	    l=readf(fdin,0,line)-1
	    if (l.lt.0) goto 200
	    read (line(:l),410) tpass,nrec
	    spec='aADR'
	    istep=1
	 endif
	 if (spec.eq.'@ADR') then
	    fact=1d2
	    iadrlen=24
	    oadrlen=iadrlen
	    naux=0
	    maux=naux
	 else if (spec.eq.'aADR') then
	    fact=1d3
	    iadrlen=24
	    oadrlen=iadrlen
	    naux=0
	    maux=naux
	 else if (spec.eq.'xADR') then
	    fact=1d3
	    l=readf(fdin,2,naux)
	    iadrlen=22+naux*2
	    oadrlen=iadrlen
	    maux=naux
	    l=readf(fdin,iadrlen-26,header(27:))
	    if (noaux) then
	       oadrlen=24
	       maux=0
	       spec='aADR'
	    endif
	 else
	    i=lnblnk(filenm(if))
	    write (0,680) filenm(if)(:i)
	    goto 9999
	 endif

	 mrec=mrec+nrec
	 rmin=1d40
	 rmax=-1d40
	 rrms=0
	 rmean=0
	 nr=0

* Before writing any record to the output XGF or ADR file we need
* to reserve space for the header that is written before the file is
* close. This is only done ONCE.

	 if (.not.headerskipped) then
	    if (mode.eq.1 .or. mode.eq.4) then
	       call writef(fdout,oadrlen,null)
	    else if (mode.eq.2) then
	       call writef(fdout,xgflen,null)
	    endif
	    headerskipped=.true.
	 endif

* Check for time limits if given

	 if (t0.ne.minint4 .and. t1.ne.maxint4) then
	    l=readf(fdin,iadrlen,adr2)
	    if (adr4(1).gt.t1) goto 120
	    call seekf(fdin,nrec*iadrlen,0)
	    l=readf(fdin,iadrlen,adr2)
	    if (adr4(1).lt.t0) goto 120
	    call seekf(fdin,iadrlen,0)
	 endif

* Process data

	 do irec=irec0,nrec,istep
	    if (mode.eq.4) then
	       l=readf(fdin,0,line)-1
	       read (line(:l),420) t,lat,lon,r
	       if (r.lt.-9000) goto 100
	       t=tpass+t-tbias/1d6
	       adr4(1)=int(t)
	       adr4(2)=nint((t-adr4(1))*1d6)
	       adr4(3)=nint(lat*1d6)
	       if (lon.gt.180) lon=lon-360
	       adr4(4)=nint(lon*1d6)
	       adr4(5)=800000000
	       adr2(11)=nint(r*fact)
	       adr2(12)=20
	    else
	       if (istep.gt.1) call seekf(fdin,iadrlen*irec,0)
	       l=readf(fdin,iadrlen,adr2)
	       adr4(2)=adr4(2)-tbias
	       t=adr4(1)+adr4(2)/1d6
	    endif
	    if (t.lt.t0 .or. t.gt.t1) goto 100
	    lat=adr4(3)/1d6
	    lon=adr4(4)/1d6
	    ssvar=adr2(12)/1d3
* For xADR convert SWH to sigSSH
	    if (naux.ne.0 .and. maux.eq.0) then
	       if (ssvar.gt.1.5d0) then
		  ssvar=1.6d-2+0.67d-2*ssvar
	       else
		  ssvar=2.3d-2+0.20d-2*ssvar
	       endif
	    endif
		  
	    if (lat.lt.lat0 .or. lat.gt.lat1) goto 100
	    if (lon.lt.lon0) lon=lon+360d0
	    if (lon.gt.lon1) lon=lon-360d0
	    if (lon.lt.lon0) goto 100
	    if (maxvar.gt.0) then
	       iy=nint((lat-glat0)/dy)
	       ix=nint((lon-glon0)/dx)
	       ssvar=grid(iy*nx+ix+1)*ssscale
	       if (ssvar.gt.maxvar) goto 100
	    endif
	    if (sigma0.gt.0) ssvar=sqrt(ssvar**2+sigma0)
	    if (orbin) then
	       j=getorb(t,lat,lon,orb,dir,.true.)
	       if (j.gt.0 .or. j.lt.ierr) goto 100
	       if (lon.lt.lon0) lon=lon+360d0
	       if (lon.gt.lon1) lon=lon-360d0
	       orbdiff=orb-adr4(5)/1d3
	       adr4(3)=nint(lat*1d6)
	       adr4(4)=nint(lon*1d6)
	       adr4(5)=nint(orb*1d3)
	       adr2(11)=nint(adr2(11)+orbdiff*fact)
	    else if (tbias.ne.0) then
	       orbdiff=-(tbias/1d6)*(adr2(14)/1d3)
	       adr4(5)=nint(adr4(5)+orbdiff*1d3)
	       adr2(11)=nint(adr2(11)+orbdiff*fact)
	    endif
	    if (orbout) then
	       r=hfact*adr4(5)/1d6
	    else if (orbdiffout) then
	       r=hfact*orbdiff
	    else if (addtopo) then
	       y=(lat-glat0)/dy+1
	       x=(lon-glon0)/dx+1
	       r=gr3int4(grid,x,y,nx,ny,nx)
	       adr4(5)=nint(r)
	       r=hfact*adr2(10+field)/fact
	    else if (topo1.gt.topo0) then
	       y=(lat-glat0)/dy+1
	       x=(lon-glon0)/dx+1
	       r=gr3int4(grid,x,y,nx,ny,nx)
	       if (r.lt.topo0 .or. r.gt.topo1) goto 100
	       r=hfact*adr2(10+field)/fact
	    else if (ref) then
	       y=(lat-glat0)/dy+1
	       x=(lon-glon0)/dx+1
	       r=hfact*adr2(10+field)/fact+gfact*gr3int4(grid,x,y,nx,ny,nx)
	       if (abs(r).gt.1d20) goto 100
	    else
	       r=hfact*adr2(10+field)/fact
	    endif
	    tmjd=t/86400+mjd85
	    if (usefit) then
	       r=r-(fit(1)+
     |		 cos(lat*rad)*(fit(2)*cos(lon*rad)+fit(3)*sin(lon*rad))+
     |		 fit(4)*sin(lat*rad))
	    endif
	    r=r+bias
	    if (nrmpnt.gt.0) then
	       if (t-tpnt(ipnt).gt.1.5d0) ipnt=0
	       ipnt=ipnt+1
	       tpnt(ipnt)=t
	       hpnt(ipnt)=r
	       latpnt(ipnt)=adr4(3)
	       lonpnt(ipnt)=adr4(4)
	       orbpnt(ipnt)=adr4(5)
	       if (ipnt.eq.nrmpnt) then
                  call normal(nrmpnt,hpnt,rpnt,rmsmax,signrm,hrms,error)
		  if (.not.error) then
		     i=nrmpnt/2+1
		     nrnrmpnt=nrnrmpnt+1
		     r=rpnt(i)
		     if (r.lt.h0 .or. r.gt.h1) goto 100
	             rmin=min(rmin,r)
	             rmax=max(rmax,r)
	             rmean=rmean+r
	             rrms=rrms+r*r
	             nr=nr+1
		     if (maxvar.le.0) ssvar=hrms
		     if (mode.eq.1 .or. mode.eq.4) then
			adr4(1)=int(tpnt(i))
			adr4(2)=nint((tpnt(i)-adr4(1))*1d6)
			adr4(3)=latpnt(i)
			adr4(4)=lonpnt(i)
			adr4(5)=orbpnt(i)
			adr2(11)=nint(r*fact)
			adr2(12)=nint(ssvar*1d3)
			call writef(fdout,oadrlen,adr2)
		     else if (mode.eq.2) then
               	        rec(1)=nint(tpnt(i))
               	        rec(2)=latpnt(i)
                        rec(3)=lonpnt(i)
                        rec(4)=nint(r*1d6)
               	        xgf2(9)=nint(ssvar*1d3)
			call writef(fdout,xgflen,xgf2)
		     else if (mode.eq.3) then
			write (6,665) tpnt(i),
     |  latpnt(i)/1d6,lonpnt(i)/1d6,r,ssvar
		     else
		        mjd=tpnt(i)/86400
		        sec=tpnt(i)-mjd*86400
			write (6,666) mjd+mjd85,sec,
     |  latpnt(i)/1d6,lonpnt(i)/1d6,orbpnt(i)/1d3-r,0,ssvar,satid
		     endif
		  endif
		  ipnt=0
	       endif
	    else
	       if (r.lt.h0 .or. r.gt.h1) goto 100
	       rmin=min(rmin,r)
	       rmax=max(rmax,r)
	       rmean=rmean+r
	       rrms=rrms+r*r
	       nr=nr+1
	       if (mode.eq.1 .or. mode.eq.4) then
		  adr2(11)=nint(r*fact)
		  adr2(12)=nint(ssvar*1d3)
		  call writef(fdout,oadrlen,adr2)
	       else if (mode.eq.2) then
	          rec(1)=nint(t)
	          rec(2)=adr4(3)
	          rec(3)=adr4(4)
	          rec(4)=nint(r*1d6)
		  xgf2(9)=nint(ssvar*1d3)
		  call writef(fdout,xgflen,xgf2)
	       else if (mode.eq.3) then
		  write (6,665) t,
     |  adr4(3)/1d6,adr4(4)/1d6,r,ssvar
	       else
		  mjd=t/86400
		  sec=t-mjd*86400
		  write (6,666) mjd+mjd85,sec,
     |  adr4(3)/1d6,adr4(4)/1d6,adr4(5)/1d3-r,0,ssvar,satid
	       endif
	    endif
100         continue
   	 enddo
	 tmin=min(tmin,rmin)
	 tmax=max(tmax,rmax)
	 tmean=tmean+rmean
	 trms=trms+rrms
	 nt=nt+nr
	 rmean=rmean/nr
	 rrms=sqrt(rrms/nr)
	 if (istep.gt.1) irec0=irec-nrec+istep-1
120	 continue
	 l=lnblnk(filenm(if))
	 if (l.gt.29) then
	    arg='<'//filenm(if)(l-27:l)
	 else
	    arg=filenm(if)
	 endif
	 if (nr.gt.0) then
	    write (0,600) arg,nrec,nr,rmin,rmax,rmean,rrms,
     |		sqrt(rrms**2-rmean**2)
	 else
	    write (0,600) arg,nrec,nr
	 endif
	 if (mode.eq.4) goto 50
         call closef(fdin)
200      continue
      enddo
      call seekf(fdout,0,0)
      if (mode.eq.1 .or. mode.eq.4) then
         bound(1)=nint(lon0)
         bound(2)=nint(lon1)
         bound(3)=nint(lat0)
         bound(4)=nint(lat1)
	 nrec=nt
	 naux=maux
	 call writef(fdout,oadrlen,header)
      else if (mode.eq.2) then
	 call writef(fdout,4,'@XGF')
	 call writef(fdout,4,nt)
      endif
      call closef(fdout)
      tmean=tmean/nt
      trms=sqrt(trms/nt)
      write (0,600) 'Total',mrec,nt,tmin,tmax,tmean,trms,
     |		sqrt(trms**2-tmean**2)
      goto 9999

600   format (a30,2i9,5f8.3)
610   format ('Normal points = ',i4,2f7.3)
620   format ('adr2adr converts ADR files to ADR files'//
     |'usage: adr2adr [options] ADR(s) ADR'//
     |'where [options] are:')
621   format ('adr2xgf converts ADR files to XGF files'/
     |'usage: adr2xgf [options] ADR(s) XGF'//
     |'where [options] are:')
622   format ('adr2asc converts ADR files to ASCII printout (stdout)'/
     |'usage: adr2asc [options] ADR(s)'//
     |'where [options] are:')
624   format ('adrlist converts ADR files to listing (stdout)'/
     |'usage: adrlist [options] ADR(s)'//
     |'where [options] are:')
625   format ('asc2adr converts Michiel''s ASCII files to ADR files'/
     |'usage: asc2adr [options] ADR(s)'//
     |'where [options] are:'/
     |'sat=SATNAME     : specify satellite name for ADR file (8 char',
     |' max., required!)')
626   format ('asc2adr: option sat=SATNAME is required!'/)
623   format (
     |'ref=gridnm      : use grid as reference'/
     |'gfact=FACTOR    : multiply reference by FACTOR (def:-1)'/
     |'hfact=FACTOR    : multiply input data by FACTOR (def:1)'/
     |'step=STEP       : take each out of STEP points'/
     |'nrm=NUM,RMS,SIG : create normal points outoff NUM points,'/
     |'                  editing at rms=RMS (m) or outlier=SIG*rms'/
     |'orb=orbit       : replace orbit'/
     |'field=nr        : use data field nr instead of SSH (nr=1)'/
     |'depth=dep0,dep1 : use data in depth-range (dep0,dep1) (m)'/
     |'topo=top0,top1  : use data in topo-range (top0,top1) (m)'/
     |'-a              : use entire arc in orbit replacement'/
     |'-x              : convert extended ADR to normal ADR'/
     |'ssvar=gridnm    : replace sigma field by value from grid'/
     |'-v              : remove data in high variability areas'/
     |'sigma0=sig      : add sig (cm) to ssvar (in RSS sense)'/
     |'-ptide          : correct for pole tide'/
     |'lon=lon0,lon1   : specify longitude boundaries'/
     |'lat=lat0,lat1   : specify latitude  boundaries'/
     |'t=t0,t1         : specify time interval',
     |' (real:[yy]yymmdd[hhmmss],mjd)'/
     |'              ... or use mjd=, doy=, ymd=, sec='/
     |'h=h0,h1         : specify height interval (cm)'/
     |'bias=bias       : add bias (cm) to sea surface heights'/
     |'fit=dR,dX,dY,dZ : remove effect of coordinate shift (m) from'/
     |'                  sea surface heights (output of ''earthfit'')'/
     |'tbias=dt        : subtract time tag bias (msec)'/
     |'-topo           : write topographic height in orbit field'/
     |'-orbdiff        : write out orbdiff i.s.o. sea height'/
     |'-orb            : write out orbit i.s.o. sea height')
665   format (f17.6,f11.6,f12.6,f9.3,f6.3)
666   format (i5,f10.3,2f12.6,f12.3,i7,f8.3,i9)
680   format ('adr2adr/adr2asc/adr2xgf: input file ',a,' is not ADR')
410   format (f17.6,i4)
420   format (f12.6,f11.6,f12.6,f9.6)

9999  end
