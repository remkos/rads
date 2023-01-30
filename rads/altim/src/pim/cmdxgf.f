      subroutine cmdxgf
      include "pim.inc"
      real mrk,t0,t1,bufsize,fieldnum
      real*8 sec85
      integer it0,it1,itrk0,itrk1,field
      logical pimopt,l,bad

* Call XGFPLOT if requested

      call pgsave
10    if (pop1cmd('XGF',argum)) then
	 field=0
         fieldnum=1.0
	 ci=-999
	 ls=1
	 lw=1
	 ch=1
	 mrk=-1
	 bufsize=maxgrd/256/3
         xname='mean.xgf'
	 itrk0=minint4
	 itrk1=maxint4
	 it0=minint4
	 it1=maxint4
	 bad=pimopt('-bad',argum,dum,dum,dum,dum)
	 l=pimopt('ci=',argum,ci,dum,dum,dum)
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('ch=',argum,ch,dum,dum,dum)
	 l=pimopt('mrk=',argum,mrk,dum,dum,dum)
	 l=pimopt('rng=',argum,rmins,rmaxs,dum,dum)
	 l=pimopt('buf=',argum,bufsize,dum,dum,dum)
         l=pimopt('fld=',argum,fieldnum,dum,dum,dum)
	 t1=0
	 if (pimopt('trk=',argum,t0,t1,dum,dum)) then
	    if (t1.eq.0) t1=t0
	    itrk0=nint(t0)
	    itrk1=nint(t1)
	 endif
	 t1=0
	 if (pimopt('t=',argum,t0,t1,dum,dum)) then
	    it0=nint(sec85(0,dble(t0)))
	    if (t1.eq.0) then
	       it1=it0+86400
	    else
	       it1=nint(sec85(0,dble(t1)))
	    endif
	 endif
	 if (pimopt('-sig',argum,dum,dum,dum,dum)) then
	    field=2	! plot sigma field
	 else if (pimopt('-mean',argum,dum,dum,dum,dum)) then
	    field=0	! plot mean field
	 else if (pimopt('-trend',argum,dum,dum,dum,dum)) then
	    field=3	! plot mean field
	 else
	    field=int(fieldnum)	! plot field=fieldnum or data field by default
	 endif
	 xgftrnd=pimopt('-show',argum,dum,dum,dum,dum)
	 call strip(argum,xname)
	 read (argum,*,iostat=ios) ci,mrk,lw,ch
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgsch(ch)
	 call xgfplot(xname,nint(ci),nint(mrk),bad,
     |		itrk0,itrk1,it0,it1,
     |		min(nint(bufsize),maxgrd/256/3),work1,field)
	 goto 10
      endif
      call pgunsa
      end

**XGFPLOT -- Plot XGF data set.
*+
      SUBROUTINE XGFPLOT (FILENM, IDX, MARKER, BAD,
     |      TRK0, TRK1, T0, T1, NWORK, WORK, FIELD)
      CHARACTER*(*) FILENM
      LOGICAL BAD
      INTEGER IDX, MARKER, TRK0, TRK1, T0, T1, NWORK, FIELD
      REAL WORK(NWORK,3,0:255)
*
* Routine plots an XGF or XXO or ADR data set FILENM as markers in a plot.
* It buffers
* a number of NWORK points in the buffer WORK first. Markers depend on
* variable MARKER.
*
* Arguments:
*  FILENM  (input): XGF filenm
*  IDX     (input): Color index for plotting the markers. Use -999 to
*                   plot colors according to height field.
*  MARKER  (input): Marker number or -999 for reading the marker number from
*                   the XGF file (sigma field).
*  BAD     (input): Plot bad points in stead of good
*  TRK0, TRK1 (in): Limits on tracknumber.
*  T0, T1  (input): Limits on time (Sec85).
*  NWORK   (input): Size of the workspace WORK.
*  WORK           : Working space.
*  FIELD   (input): Specifies input field.
*-
      include "pim.inc"

      integer i,j,k,n(0:255),lat,lon,nrec,itime,mark,unit,hgt,l0,l1
      integer naux,nrec2,itime0,freeunit
      integer ssh(2),itime2(2),itime4(4),u(2)
      integer*2 sig,trk,trk2(2),naux2,hgt2
      integer*2 columns
      real rhgt(30)
      real v,vfact,cmax,a,b,r,rlon,rlat,rdum
      character spec*4,spec2*2
      character line*80
      equivalence (trk2,itime)

* Open XGF file on new unit

      i=index(filenm,' ')-1
      write (0,600) filenm(:i)
550   format (a)
600   format ('Plotting file ',a,' ...')
      unit=freeunit()
      open (unit,file=filenm,status='old',form='unformatted',
     |   access='direct',recl=26)
      read (unit,rec=1) spec,nrec,naux,dum,dum,nrec2,naux2
      l0=1
      l1=nrec
      if (spec.eq.'@XXO') then
	 close (unit)
         open (unit,file=filenm,status='old',form='unformatted',
     |   access='direct',recl=44)
      else if (spec.eq.'@XAB') then
	 close (unit)
         open (unit,file=filenm,status='old',form='unformatted',
     |   access='direct',recl=28)
         call xtflimits(filenm,t0,t1,trk0,trk1,l0,l1)
      else if (spec.eq.'@XXB') then
	 close (unit)
         open (unit,file=filenm,status='old',form='unformatted',
     |   access='direct',recl=48)
      else if (spec.eq.'aADR') then
	 close (unit)
         open (unit,file=filenm,status='old',form='unformatted',
     |   access='direct',recl=24)
	 nrec=nrec2
	 l1=nrec
      else if (spec.eq.'xADR') then
	 close (unit)
         open (unit,file=filenm,status='old',form='unformatted',
     |   access='direct',recl=22+naux2*2)
	 nrec=nrec2
	 l1=nrec
      else if (spec.eq.'@XGF') then
	 close (unit)
         open (unit,file=filenm,status='old',form='unformatted',
     |   access='direct',recl=18)
         read (unit,rec=1) spec,nrec,dum,dum,spec2
      else if (spec.eq.'# RA') then
         close(unit)
         open (unit,file=filenm,status='old',access='sequential',
     |         form='formatted')
         spec = 'RASC'
         columns = 0
         do while (.true.)
             read (unit,550,err=13,end=13) line
             if(line(1:1).ne.'#') then
               backspace unit
               goto 14
	     else if (line(:10).eq.'# RADS_ASC') then
	       columns = 0
             else if (index(line,"Col").gt.0) then
               columns = columns + 1
             endif
         enddo 
      endif
      goto 14

  13  continue 
      write(0,550) '... wrong file type. Exit'
      return

  14  continue 
      if (spec.eq.'RASC'.and.(columns.lt.1.or.field.gt.columns-3)) then
        write(0,550) '... not enough columns. '
        goto 13
      endif

* Clean counter

      do j=0,255
	 n(j)=0
      enddo

      j=idx
      mark=marker
      vfact=(nc1+1)/(rmaxs-rmins)
      cmax=nc1+0.999

*      write(0,*) j, nc1, rmaxs, rmins, vfact, cmax

* Sort markers by colour

      i=l0-1
40    i=i+1
	 if (spec.eq.'@XXO' .or. spec.eq.'xXXO') then
	    read (unit,rec=i+1,err=50) lat,lon,itime4,trk2,ssh
	    if (ci.eq.-888 .and. trk2(1).lt.trk2(2)) then
	       hgt=ssh(1)-ssh(2)
	    else
	       hgt=ssh(2)-ssh(1)
	    endif
	    if (itime4(1).lt.t0 .or. itime4(1).gt.t1) goto 50
	    if (itime4(3).lt.t0 .or. itime4(3).gt.t1) goto 50
	    if (trk2(1).lt.trk0 .or. trk2(1).gt.trk1) goto 50
	    if (trk2(2).lt.trk0 .or. trk2(2).gt.trk1) goto 50
	 else if (spec.eq.'@XAB') then
	    read (unit,rec=i+1,err=50) itime,lat,lon,ssh,u(1),sig,trk
	    if (itime.lt.t0 .or. itime.gt.t1) goto 50
	    if (trk.lt.trk0 .or. trk.gt.trk1) goto 50
	    if ((sig.lt.0).neqv.bad) goto 50
	    hgt=ssh(2)
	 else if (spec.eq.'@XXB') then
	    read (unit,rec=i+1,err=50) lat,lon,itime2,trk2,ssh,ssh,u,sig
	    if (itime2(1).lt.t0 .or. itime2(1).gt.t1) goto 50
	    if (itime2(2).lt.t0 .or. itime2(2).gt.t1) goto 50
	    if (trk2(1).lt.trk0 .or. trk2(1).gt.trk1) goto 50
	    if (trk2(2).lt.trk0 .or. trk2(2).gt.trk1) goto 50
	    if ((sig.lt.0).neqv.bad) goto 50
	    hgt=(ssh(1)+ssh(2))/2
	 else if (spec.eq.'aADR' .or. spec.eq.'xADR') then
	    read (unit,rec=i+1,err=50) itime,dum,lat,lon,dum,hgt2,sig
	    if (itime.lt.t0 .or. itime.gt.t1 .or. sig.lt.0) goto 50
	    if (field.eq.2) then
	       hgt=sig*1000
	    else
	       hgt=hgt2*1000
	    endif
         else if (spec.eq.'RASC') then
            read (unit,*,err=50,end=51) rdum,rlat,rlon,
     |                                     (rhgt(k),k=1,columns-3)
            hgt=nint(rhgt(field)*1e6)
            lat=nint(rlat*1e6)
            lon=nint(rlon*1e6)
	 else if (spec.eq.'@XGF' .and. spec2.ne.'Tm') then
	    read (unit,rec=i+1,err=50) itime,lat,lon,hgt,sig
	    if (itime.lt.t0 .or. itime.gt.t1 .or. sig.lt.0) goto 50
	    if (field.eq.2) hgt=sig*1000
	    if (marker.eq.-999) mark=sig
	 else if (spec.eq.'@XGF' .and. spec2.eq.'Tm') then
	    read (unit,rec=i+1,err=50) itime,lat,lon,hgt,sig
	    if (field.eq.0) then
	       if (trk.lt.trk0 .or. trk.gt.trk1 .or. sig.lt.0) goto 50
	    else if (field.eq.1) then
	       if (itime.lt.t0 .or. itime.gt.t1 .or. sig.ge.0) goto 50
	    else if (field.eq.2) then
	       if (trk2(1).lt.trk0 .or. trk2(1).gt.trk1 .or.
     |			sig.lt.0) goto 50
	       hgt=sig*1000
	    else
	       if (trk2(1).lt.trk0 .or. trk2(1).gt.trk1 .or.
     |			sig.lt.0) goto 50
	       nstack=0
10	       read (unit,rec=i+nstack+2,err=50) itime,lat,lon,hgt,sig
	       if (nstack.eq.0) itime0=itime
	       if (sig.lt.0) then
		  nstack=nstack+1
	          tx(nstack)=(itime-itime0)/86400/36525d0
		  ty(nstack)=hgt
		  goto 10
	       endif
	       call regres(nstack,tx,ty,a,b,r)
	       hgt=nint(b)
	       i=i+nstack
	    endif
	    if (marker.eq.-999) mark=sig
	 endif

* Convert longitude to requested range. Check longitude and latitude.

	 rlon=lon/1e6
	 rlat=lat/1e6
	 if (rlon.lt.xw0) then
	    rlon=rlon+360
	 else if (rlon.gt.xw1) then
	    rlon=rlon-360
	 endif
	 if (rlon.ge.xw0 .and. rlat.ge.yw0 .and. rlat.le.yw1) then
            if (idx.lt.0) then
	       v=(hgt/1e6-rmins)*vfact
	       j=c_0+max(0.,min(v,cmax))
	    endif
	    n(j)=n(j)+1
	    work(n(j),1,j)=rlon
	    work(n(j),2,j)=rlat
	    work(n(j),3,j)=mark
	    if (n(j).eq.nwork) call xgfplot1(n(j),
     |		work(1,1,j),work(1,2,j),work(1,3,j),marker,j)
         endif
50    if (i.lt.l1) goto 40
51    continue

      close (unit)

      do j=0,255
	 call xgfplot1(n(j),
     |		work(1,1,j),work(1,2,j),work(1,3,j),marker,j)
      enddo

      end

      subroutine xgfplot1(n,x,y,z,marker,ic)
      integer n,marker,i,ic
      real x(*),y(*),z(*)

      if (n.eq.0) return
      call pgsci(ic)
      call pmconv(n,x,y)
      if (marker.eq.-999) then
         do i=1,n
	    call pgpt(1,x(i),y(i),nint(z(i)))
	 enddo
      else
	 call pgpt(n,x,y,marker)
      endif

      n=0
      end

**REGRES -- Regression analysis
*+
      SUBROUTINE REGRES (N, X, Y, A, B, R)
      INTEGER*4 N
      REAL  X(N), Y(N), A, B, R
*
* This subroutine computes the best fitting straight line trough a number
* of points with coordinates (X,Y). Upon return, A and B will be the
* coefficients of the line
*
*     Y = A + B * X
*
* that is the best fit through the data points. R is the correlation
* between the fit and the data.
*
* Arguments:
*  N  (input) : Number of data points
*  X  (input) : Array of X coordinates of the data points
*  Y  (input) : Array of Y coordinates of the data points
*  A (output) : Coefficient of linear fit (Y = A + B * X)
*  B (output) : Coefficient of linear fit (Y = A + B * X)
*  R (output) : Correlation
*-
*  2-Jul-1993 - New manual
*-----------------------------------------------------------------------
      integer*4 i
      real*8 sumx,sumy,sumxx,sumxy,sumyy,uxx,uxy,uyy
      sumx=0d0
      sumy=0d0
      sumxy=0d0
      sumxx=0d0
      sumyy=0d0
      do i=1,n
         sumx =sumx +x(i)
         sumy =sumy +y(i)
         sumxy=sumxy+x(i)*y(i)
         sumxx=sumxx+x(i)*x(i)
         sumyy=sumyy+y(i)*y(i)
      enddo
      uxx=n*sumxx-sumx*sumx
      uxy=n*sumxy-sumx*sumy
      uyy=n*sumyy-sumy*sumy
      b=uxy/uxx
      a=(sumy-b*sumx)/n
      r=uxy/sqrt(uxx*uyy)
      end
