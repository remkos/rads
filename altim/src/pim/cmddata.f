      subroutine cmddata
      implicit none
      include "pim.inc"
      real mrk,bufsize,colnum(3)
      integer i,col(3)

* Read and plot data file (as points)

      call pgsave
10    if (pop1cmd('DATA',argum)) then
	 colnum(1)=1
	 colnum(2)=2
	 colnum(3)=0
	 ci=-999
	 ls=1
	 lw=1
	 ch=1
	 mrk=-1
	 bufsize=maxgrd/256/3
         xname=''
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 call pimopt('ls=',argum,ls,dum,dum,dum)
	 call pimopt('lw=',argum,lw,dum,dum,dum)
	 call pimopt('ch=',argum,ch,dum,dum,dum)
	 call pimopt('mrk=',argum,mrk,dum,dum,dum)
	 call pimopt('rng=',argum,rmins,rmaxs,dum,dum)
	 call pimopt('buf=',argum,bufsize,dum,dum,dum)
         call pimopt('col=',argum,colnum(1),colnum(2),colnum(3),dum)
         call pimopt('fld=',argum,colnum(1),colnum(2),colnum(3),dum)
	 call strip(argum,xname)
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgsch(ch)
	 do i=1,3
	    col(i)=nint(colnum(i))
	 enddo
	 call dataplot(xname,nint(ci),nint(mrk),
     |		min(nint(bufsize),maxgrd/256/3),work1,col)
	 goto 10
      endif
      call pgunsa
      end

**DATAPLOT -- Plot data set.
*+
      SUBROUTINE DATAPLOT (FILENM, IDX, MARKER, NWORK, WORK, COL)
      implicit none
      CHARACTER*(*) FILENM
      INTEGER IDX, MARKER, NWORK, COL(3)
      REAL WORK(NWORK,3,0:255)
*
* Routine plots a data set FILENM as markers in a plot.
* It buffers a number of NWORK points in the buffer WORK first.
* Markers depend on variable MARKER.
*
* Arguments:
*  FILENM (input): File name of ASCII data set
*  IDX    (input): Color index for plotting the markers. Use -999 to
*                  plot colors according to third specified column COL(3)
*  MARKER (input): Marker number or -999 for reading the marker number from
*                  third specified column COL(3)
*  NWORK  (input): Size of the workspace WORK.
*  WORK          : Working space.
*  COL    (input): Specifies input columns (lat, lon, hgt)
*                  If COL(3)=0 markers use IDX for colour.
*-
      include "pim.inc"

      integer i,j,k,n(0:255),mark,unit,freeunit
      integer columns
      real value(30)
      real v,vfact,cmax,lon,lat,hgt

* Open data file on new unit

      i=index(filenm,' ')-1
      write (0,600) filenm(:i)
600   format ('Plotting file ',a,' ...')
      unit=freeunit()
      open (unit,file=filenm,status='old')

* Clean counter

      do j=0,255
	 n(j)=0
      enddo

      j=idx
      mark=marker
      vfact=(nc1+1)/(rmaxs-rmins)
      cmax=nc1+0.999
      columns=max(col(1),col(2),col(3))

* Sort markers by colour

      i=0
40    read (unit,*,iostat=ios,end=51) (value(k),k=1,columns)
      i=i+1
      lat=value(col(1))
      lon=value(col(2))
      hgt=value(col(3))
*      write (*,*) i,ios,lat,lon,hgt,col(1),col(2),col(3),columns

* Convert longitude to requested range. Check longitude and latitude.

      if (lon.lt.xw0) then
         lon=lon+360
      else if (lon.gt.xw1) then
         lon=lon-360
      endif
      if (lon.ge.xw0 .and. lat.ge.yw0 .and. lat.le.yw1) then
         if (idx.lt.0) then
            v=(hgt-rmins)*vfact
            j=c_0+max(0.,min(v,cmax))
         endif
         n(j)=n(j)+1
         work(n(j),1,j)=lon
         work(n(j),2,j)=lat
         work(n(j),3,j)=mark
         if (n(j).eq.nwork) call xgfplot1(n(j),
     |		work(1,1,j),work(1,2,j),work(1,3,j),marker,j)
      endif
      goto 40
51    continue

      close (unit)
      do j=0,255
	 call xgfplot1(n(j),
     |		work(1,1,j),work(1,2,j),work(1,3,j),marker,j)
      enddo
      end
