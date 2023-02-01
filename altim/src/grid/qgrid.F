**QGRID -- Make grid of ASCII or BINARY data

      program qgrid

* A Quick'n'dirty gridding program for lots of altimeter formats
* XAF, XXF, XGF, ADR and ASCII data
*
*-
* Last modifications:
* 23-Oct-1997 - 9710.0 - Quick entry into XAF files when timeinterval is
* selected by using XTF file
* 24-Oct-1997 - 9710.1 - Time selection for ADR implemented
* 13-Apr-1999 - 9904.0 - Add option num= and -w
*  3-Nov-1999 - 9911.0 - Added DATEARG function
* 10-Jan-2000 - 0001.0 - Circumvent F90 compiler for Linux
* 17-Jun-2002 - 0206.0 - Increase maxfiles
* 15-Mar-2003 - 0303.0 - Included grid= options
*-----------------------------------------------------------------------
      implicit none

      character*6 version
      parameter (version="0206.0")
      integer*4 maxfiles
      parameter (maxfiles=100)
      character*80 filenm(maxfiles),meannm,varnm,rmsnm,numnm,arg
      integer narg,iargc,iarg,nx,ny,ifill,k,jx,minnr,nf,if
      integer ios,project,l,lnblnk
      real*8  xa,xb,ya,yb,dum
      logical wrapit,datearg

      integer*4 minint4,maxint4
      parameter (maxint4=2147483647,minint4=-maxint4-1)

#ifdef f90
      real*8, allocatable :: pntmean(:),pntrms(:)
      integer, allocatable :: pntnr(:)
#else
      integer nrgridpoints
      parameter (nrgridpoints=5000000)
      real*8 pntmean(nrgridpoints),pntrms(nrgridpoints)
      integer pntnr(nrgridpoints)
#endif

      include "qgrid.inc"

      minnr=1
      meannm=' '; varnm=' '; rmsnm=' '; numnm=' '
      x0=-178d0; x1=+178d0; dx=4d0
      y0= -80d0; y1= +80d0; dy=4d0
      range0=-1e35; range1=+1e35
      t0=minint4; t1=maxint4
      trk0=0; trk1=999999999
      field=0; ifill=0
      tsort=.false.; sum=.false.; wrapit=.true.; useaux=.false.
      ocean=.false.
      ptype=1; project=1
      nf=0
      tbias(1)=0; tbias(2)=0; ssb=0
      
      narg=iargc()
      if (narg.eq.0) goto 1300
	 
      do iarg=1,narg
	 call getarg(iarg,arg)
	 if (arg(:5).eq.'proj=') then
	    read (arg(6:),*) project
	 else if (arg(:4).eq.'lat=') then
	    read (arg(5:),*) y0,y1,dy
	 else if (arg(:4).eq.'Lat=') then
	    read (arg(5:),*) y0,y1,dy
	    y0=y0+dy/2
	    y1=y1-dy/2
	 else if (arg(:4).eq.'lon=') then
	    read (arg(5:),*) x0,x1,dx
	 else if (arg(:4).eq.'Lon=') then
	    read (arg(5:),*) x0,x1,dx
	    x0=x0+dx/2
	    x1=x1-dx/2
	 else if (datearg(arg,t0,t1,dum)) then
	 else if (arg(:4).eq.'trk=') then
	    read (arg(3:),*) trk0,trk1
	 else if (arg(:4).eq.'fil=') then
	    read (arg(5:),*) ifill
	 else if (arg(:4).eq.'rng=') then
	    read (arg(5:),*) range0,range1
	    range0=range0/100
	    range1=range1/100
	 else if (arg(:4).eq.'min=') then
	    read (arg(5:),*) minnr
	 else if (arg(:5).eq.'mean=') then
	    meannm=arg(6:)
	 else if (arg(:4).eq.'var=') then
	    varnm=arg(5:)
	 else if (arg(:4).eq.'rms=') then
	    rmsnm=arg(5:)
	 else if (arg(:4).eq.'num=') then
	    numnm=arg(5:)
	 else if (arg(:5).eq.'grid=') then
	    l=lnblnk(arg)
	    meannm=arg(6:l)//'_mean.grd'
	    varnm=arg(6:l)//'_var.grd'
	    numnm=arg(6:l)//'_num.grd'
	 else if (arg(:6).eq.'tbias=') then
            tbias(2)=1d30
            read (arg(7:),*,iostat=ios) tbias
            if (tbias(2).gt.1d20) tbias(2)=tbias(1)
* tbias in milliseconds. Positive when tags are late.
* Will be multiplied by orbit rate in mm/s to give micrometres
            useaux=.true.
      else if (arg(1:4).eq.'ssb=') then
            read (arg(5:),*) ssb
* Addition ssb in percentage. Positive when ssb should be higher.
            ssb=ssb*1d1
* ssb in permille will be multiplied by SWH in mm to give micrometres
            useaux=.true.
	 else if (arg(:2).eq.'-h') then
	    goto 1300
	 else if (arg(:2).eq.'-o') then
	    ocean=.true.
         else if (arg(1:2).eq.'-a') then
            field=0
         else if (arg(1:2).eq.'-p') then
            field=1
         else if (arg(1:2).eq.'-d') then
            field=2
         else if (arg(1:2).eq.'-m') then
            sum=.true.
         else if (arg(1:2).eq.'-x') then
            sum=.false.
         else if (arg(1:2).eq.'-w') then
            wrapit=.false.
         else if (arg(1:2).eq.'-t' .or. arg(1:2).eq.'-e') then
            tsort=.true.
         else if (arg(1:2).eq.'-s') then
            field=1
	    read (arg(3:),*,iostat=ios) field
	 else
	    nf=nf+1
	    if (nf.gt.maxfiles) call fin('too many files')
	    filenm(nf)=arg
	 endif
      enddo

* Compute number of cells

      nx=nint((x1-x0)/dx+1)
      ny=nint((y1-y0)/dy+1)

* Convert the boudary y-coordinates for use with cilindrical projections.
* This has no effect for straight rectangular projections.
* Keep original values

      ya=y0; yb=y1; xa=x0; xb=x1
      call setproj(project)
      if (project.eq.41 .or. project.eq.42) then
         call convxy(0d0,y0)
	 y1=-y0; x0=y0; x1=y1
	 xa=x0; xb=x1; ya=y0; yb=y1
	 ny=2*ny-1; nx=ny
      else
         call convxy(x0,y0)
         call convxy(x1,y1)
      endif

* Recompute cell size

      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)

#ifndef f90
      if (nx*ny.gt.nrgridpoints) goto 1320
      do k=1,nx*ny
         pntmean(k)=0; pntrms(k)=0; pntnr(k)=0
      enddo
#else

* Free memory for mean, rms, and nr grids

      allocate (pntmean(nx*ny))
      allocate (pntrms(nx*ny))
      allocate (pntnr(nx*ny))

      pntmean(1:nx*ny)=0
      pntrms (1:nx*ny)=0
      pntnr  (1:nx*ny)=0
#endif

      do if=1,nf
         call gridding(nx,ny,pntmean,pntrms,pntnr,filenm(if))
      enddo

      jx=nint(360d0/dx)+1

      if (jx.le.nx .and. wrapit)
     |	call wrap(jx,nx,ny,pntmean,pntrms,pntnr)
      if (ifill.gt.0)
     |	call fillup(nx,ny,pntmean,pntrms,pntnr,ifill)

      call toolittle(nx,ny,pntmean,pntrms,pntnr,minnr)

      if (meannm.ne.' ') then
	 call gridwr8(meannm,nx,ny,pntmean,nx,xa,xb,ya,yb)
      endif
      if (rmsnm.ne.' ') then
	 call gridwr8(rmsnm,nx,ny,pntrms,nx,xa,xb,ya,yb)
      endif
      if (varnm.ne.' ') then
         call sigma(nx,ny,pntmean,pntrms,pntnr,minnr)
         call gridwr8(varnm,nx,ny,pntrms,nx,xa,xb,ya,yb)
      endif
      if (numnm.ne.' ') then
         do jx=1,nx*ny
	    if (pntnr(jx).ge.minnr) then
	       pntrms(jx)=pntnr(jx)
	    else
	       pntrms(jx)=1e35
	    endif
	 enddo
         call gridwr8(numnm,nx,ny,pntrms,nx,xa,xb,ya,yb)
      endif
      goto 9999

  550 format (a)
  600 format (
     |'qgrid (',a,
     |') : Quick grid height values in ASC/ADR/XGF/XAF/XXF/XXO'//
     |'usage: qgrid [ options ] filename(s)'//
     |'with required fields:'//
     |'filename(s)   : name of the input file(s) (- for stdin)'/
     |16x,'These files may be in ADR/XGF/XAF/XXF/XXO format or ASCII'//
     |'and options:'//
     |'-o            : use only values over deep oceans (def: off)'/
     |'-a            : use a priori value (def, XAF/XXF only)'/
     |'-d            : use a correction value (XAF/XXF only)'/
     |'-p            : use a posteriori value (XAF/XXF only)'/
     |'-x            : compute asc-des (def, XXF/XXO only)'/
     |'-m            : compute mean value (XXF/XXO only)'/
     |'-t            : compute first-last (XXF/XXO only)'/
     |'-e            : compute ERS-T/P (def, XXO only)'/
     |'-s            : use sigma value (ADR/XGF only)'/
     |'-w            : do not wrap around 360 degrees'/
     |'proj=type     : use projection'/
     |'lat=y0,y1,dy  : latitude  boundaries and cellsize (deg)',
     |' (def: -80,80,4)'/
     |'lon=x0,x1,dx  : longitude boundaries and cellsize (deg)',
     |' (def: -178,178,4)'/
     |'Lat=y0,y1,dy  : latitude  cellboundaries and cellsize (deg)',
     |' (def: -82,82,4)'/
     |'Lon=x0,x1,dx  : longitude cellboundaries and cellsize (deg)',
     |' (def: -180,180,4)'/
     |'t=t0,t1       : time selection ([yy]yymmdd[hhmmss],mjd)'/
     |'            ... or use mjd=, doy=, ymd=, sec='/
     |'fil=number    : fill gaps if enclosed by 2 or 4 valid cells'/
     |'rng=r0,r1     : use only values in range (r0,r1 [cm])'/
     |'min=n         : minimum nr of points per cell (def: 1)'/
     |'mean=gridname : create mean grid            | use 1 or'/
     |'var=gridname  : create variance grid        | more of'/
     |'rms=gridname  : create rms grid             | these'/
     |'num=gridname  : create grid with nr of pnts | fields'//
     |'ASCII files are free format but must contain 4 fields:'/
     |'time  latitude(deg)  longitude(deg)  value(m?)'/
     |'with time in YYMMDDHHMMSS.SS, YYMMDD.DD, MJD.DD, or SEC85.'/
     |'Lines starting with # are disregarded')

 1300 write (0,600) version
      goto 9999

 1320 write (0,550) 'qgrid: too many gridpoints'
 9999 continue
      end
