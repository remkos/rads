**PMWDB -- draw coastlines, rivers, boundaries or land masses.
*+
      SUBROUTINE PMWDB (NAME, I1, I2)
      CHARACTER*(*) NAME
      INTEGER I1, I2
*
* Draw coastlines/rivers/boundaries/landmasses/lakes from World Data
* Bank (WDB) in the present PGPLOT window. The map boundaries and projection
* type must have been specified using the proper PMPLOT routines (PMDEF and
* PMWINDOW). The map drawn by this routine may extend to any longitude
* (x-coordinate) and in latitude (y-coordinate) from -90 to +90 degrees.
*
* ** In case a line-segment data set (.bth, .cil, .bdy, .pby, .mnt, .bth, .plt)
*    is used:
*
* The line segments stored in the WDB all have rank-numbers. The
* most important line-segments have rank 1, then comes rank 2, etc.
* You may specify the ranks to be drawn by selecting proper values for
* I1 (minimum rnak number) and I2 (maximum rank number). Use I2<I1 if all
* ranks must be plotted.
* All relevant line segments are drawn with line properties set by
* PGSCI, PGSCR, PGSLS, PGSLW.
*
* Arguments:
*  NAME (input) : Name of the WDB dataset to be drawn, without the
*                 extensions '.DAT' or '.TAB'. If the dataset cannot
*                 be found in the current directory, the default directory
*                 specified by the environment variable WDB_DIR is searched.
*  I1   (input) : Minimum rank.
*  I2   (input) : Maximum rank.
*
* ** In case a landmasses data set (.lnd) is used:
*
* Landmasses and lakes are drawn as polygons with color indices I1 and I2,
* respectively. The color index selection before the call to PMWDB is
* restored after execution.
*
* Arguments:
*  NAME (input) : Name of the WDB dataset to be drawn, without the
*                 extensions '.DAT' or '.TAB'. If the dataset cannot
*                 be found in the current directory, the default directory
*                 specified by the environment variable WDB_DIR is searched.
*  I1   (input) : Color index for landmasses.
*  I2   (input) : Color index for lakes.
*--
*  9-Jan-1991 - created [Remko Scharroo].
* 14-Jan-1991 - use GRVCT0 to draw dot and implement conic projections.
* 17-Jan-1991 - import outside working space.
* 20-Feb-1991 - support negative ranks.
* 10-Apr-1991 - compressed data sets to use on dutlru2/3 with increased speed.
* 15-Apr-1991 - remove use of outside working space.
* 31-Oct-1991 - shift window if there is space west of -180 or east of 180.
* 29-Nov-1991 - include a default WDB directory.
* 19-Mar-1992 - Standardize PMPLOT include landmasses plotting.
* 23-Apr-1992 - search default directory if not in current one.
* 16-Jun-1993 - Faster version using fastio routines
* 21-Jul-1993 - PGBBUF and PGEBUF removed from PGLND1
* 21-Oct-1993 - Including new format.
* 11-Jul-1996 - New PMCPOLY.
* 17-Jun-2003 - Removed sysdep.h
*-----------------------------------------------------------------------
* fastio version
*
      include 'pmplot.inc'
      integer fd_tab,fd_dat
      character*132 filenm,wdbdir
      character*4 spec
      integer l,m,lnblnk,ires,ind,n,GRIFIL
      parameter (n=18000)
      real delta,x(n),y(n),res
      integer*1 work(n*4)
      logical exist,new,swap,ltlend

      common /cpmwdb/ fd_tab,fd_dat,swap

      if (pmopen.lt.4) return
*
* Substitute the default directory name if it is not given.
*
      filenm=name
      l=lnblnk(filenm)
      filenm(l+1:)='.TAB'
      inquire (file=filenm,exist=exist)
      if (.not.exist) then
         call getenv('WDB_DIR',wdbdir)
         m=lnblnk(wdbdir)
         filenm=wdbdir
	 filenm(m+1:)='/'
	 filenm(m+2:)=name
         l=l+m+1
         filenm(l+1:)='.TAB'
      endif
      swap=ltlend()
*
* Open table file:
* NAME.TAB = 36-bytes binary direct access table containing
*            positions of line segments and ranks
*
* - If this is an old wdb data set open the data file
* NAME.DAT = 2-bytes binary direct access file containing lat/lon
* - If this is a new wdb data set open the data file
* NAME.DAT = 4-bytes binary direct access file containing lat/lon
* - For land data sets, open
* NAME.DAT = 4-bytes binary direct access file containing lat/lon
*
      fd_tab=GRIFIL(filenm)
      call GRRFIL(fd_tab,4,spec)
      call GRRFIL(fd_tab,4,ires)
      if (swap) call i4swap(1,ires)
      res=ires
      filenm(l+1:)='.DAT'
      if (spec.eq.'@WDB' .or. spec.eq.'@LND' .or. spec.eq.'@WDN') then
	 fd_dat=GRIFIL(filenm)
      else
	 call grwarn('pmwdb: not a WDB data set')
	 call GRCFIL(fd_tab)
	 return
      endif
      new=(spec.eq.'@WDN')

      delta=0.
      call pgqci(ind)
   10 if (lonmin.lt.delta+180.) then
	 call GRSEEK(fd_tab,36,0)
	 if (spec.eq.'@LND') then
            call pmlnd1(res,i1,i2,delta,n,x,y,work)
	 else
            call pmwdb1(res,i1,i2,delta,n,x,y,work,work,new)
	 endif
         delta=delta-360.
         goto 10
      endif
      delta=360.
   20 if (lonmax.gt.delta-180.) then
	 call GRSEEK(fd_tab,36,0)
	 if (spec.eq.'@LND') then
            call pmlnd1(res,i1,i2,delta,n,x,y,work)
	 else
            call pmwdb1(res,i1,i2,delta,n,x,y,work,work,new)
	 endif
         delta=delta+360.
         goto 20
      endif
      call pgsci(ind)
      call GRCFIL(fd_tab)
      call GRCFIL(fd_dat)
      end
      
      subroutine pmwdb1(res,ranklo,rankhi,delta,n,x,y,wk1,wk2,new)
*
      include 'pmplot.inc'
      integer fd_tab,fd_dat
      logical new,swap
      character*1 wk1(2,n)
      integer*2   wk2(2,n)
      integer rank,rank1,rank0,tab(9),ranklo,rankhi,n,pmconv,ios
      real x(n),y(n),res,delta
      integer ixhi,ixlo,iyhi,iylo,npnt,nrec,lon0,lon1,lat0,lat1
      integer ix0,iy0,j
      equivalence (tab(1),rank),(tab(2),npnt),(tab(3),nrec),
     .(tab(4),ixlo),(tab(5),ixhi),(tab(6),iylo),(tab(7),iyhi),
     .(tab(8),ix0),(tab(9),iy0)

      common /cpmwdb/ fd_tab,fd_dat,swap

      if (rankhi.lt.ranklo) then
         rank0=-2147483647-1
         rank1=+2147483647
      else
         rank0=ranklo
         rank1=rankhi
      endif
*
* determine window boundaries in resolution units. note that the window
* may be shifted by delta degrees in longitude.
*
      lon0=nint((lonmin-delta)*res-.5)
      lon1=nint((lonmax-delta)*res+.5)
      lat0=nint(latmin*res-.5)
      lat1=nint(latmax*res+.5)
*
* read table entries and check if line segment is within the window
*
   30 call GRRFIL(fd_tab,36,tab)
      if (swap) call i4swap(9,tab)
      if (nrec.le.0) return
      if (ixhi.lt.lon0 .or. ixlo.gt.lon1 .or.
     .    iyhi.lt.lat0 .or. iylo.gt.lat1 .or.
     .    rank.lt.rank0 .or. rank.gt.rank1) goto 30
*
* read and draw line segments
*
      x(1)=ix0/res+delta
      y(1)=iy0/res
      if (npnt.eq.1) then
         ios=pmconv(1,x,y)
         call grvct0(3,.false.,1,x,y)
         goto 30
      else if (new) then
         call wdbdatrd(4*nrec,4*npnt-4,wk2(1,2))
         if (swap) call i2swap(2*npnt-2,wk2(1,2))
         do j=2,npnt
	    ix0=ix0+wk2(1,j)
	    iy0=iy0+wk2(2,j)
	    x(j)=ix0/res+delta
            y(j)=iy0/res
         enddo
      else
         call wdbdatrd(2*nrec,2*npnt-2,wk1(1,2))
         do j=2,npnt
            ix0=ix0+ichar(wk1(1,j))-127
            iy0=iy0+ichar(wk1(2,j))-127
            x(j)=ix0/res+delta
            y(j)=iy0/res
         enddo
      endif
      ios=pmconv(npnt,x,y)
      call pgbbuf
      call grmova(x(1),y(1))
      do j=2,npnt
         call grlina(x(j),y(j))
      enddo
      call pgebuf
      goto 30
      end

      subroutine pmlnd1(res,land,lake,delta,n,x,y,work)
*
      include 'pmplot.inc'
      integer fd_tab,fd_dat
      logical swap
      integer pmcpoly,pmconv
      integer tab(9),rank,land,lake,n
      real x(n),y(n),res,delta
      integer ixhi,ixlo,iyhi,iylo,npnt,nrec,lon0,lon1,lat0,lat1
      integer ix0,iy0,j,lonu0,lonu1,latu0,latu1
      integer*2 work(2,n)
      equivalence (tab(1),rank),(tab(2),npnt),(tab(3),nrec),
     .(tab(4),ixlo),(tab(5),ixhi),(tab(6),iylo),(tab(7),iyhi),
     .(tab(8),ix0),(tab(9),iy0)

      common /cpmwdb/ fd_tab,fd_dat,swap
*
* determine window boundaries in resolution units. note that the window
* may be shifted by delta degrees in longitude.
*
      lon0=nint((lonmin-delta)*res-.5)
      lon1=nint((lonmax-delta)*res+.5)
      lat0=nint(latmin*res-.5)
      lat1=nint(latmax*res+.5)
      lonu0=nint((lonumin-delta)*res-.5)
      lonu1=nint((lonumax-delta)*res+.5)
      latu0=nint(latumin*res-.5)
      latu1=nint(latumax*res+.5)
*
* read table entries and check if line segment is within the window
*
   30 call GRRFIL(fd_tab,36,tab)
      if (swap) call i4swap(9,tab)
      if (nrec.le.0) return
      if (ixhi.lt.lon0 .or. ixlo.gt.lon1 .or.
     .    iyhi.lt.lat0 .or. iylo.gt.lat1 .or. npnt.lt.3) goto 30
      call wdbdatrd(4*nrec,4*npnt,work)
      if (swap) call i2swap(2*npnt,work)
      do j=1,npnt
	 x(j)=(ix0+work(1,j))/res+delta
         y(j)=(iy0+work(2,j))/res
      enddo
      if (rank.eq.2) then
         call pgsci(lake)
      else
         call pgsci(land)
      endif
      if (.not.azmtal .and.
     .   ixlo.ge.lonu0 .and. ixhi.le.lonu1 .and.
     .   iylo.ge.latu0 .and. iyhi.le.latu1) then
	 if (pmconv(npnt,x,y).gt.0) call pgpoly(npnt,x,y)
      else
	 if (pmcpoly(npnt,x,y).gt.0) call pgpoly(npnt,x,y)
      endif
      goto 30
      end

      subroutine wdbdatrd(skip,bytes,buffer)
*
* read data record of 2*2 bytes.
*
      integer skip,bytes
      integer*1 buffer(*)
      integer fd_tab,fd_dat
      logical swap
      common /cpmwdb/ fd_tab,fd_dat,swap
      call GRSEEK(fd_dat,skip,0)
      call GRRFIL(fd_dat,bytes,buffer)
      end
