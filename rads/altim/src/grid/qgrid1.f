      subroutine gridding(nx,ny,mean,rms,nr,filenm)
      implicit none
      character*80 filenm,name
      integer nx,ny,nr(nx,ny),l0,l1
      real*8 mean(nx,ny),rms(nx,ny),t
      integer kx,ky,irec,istat
      real*8 value,dlon,dlat,range/0d0/
      logical sland
      integer openf,ios,seekf,readf,fdin,fdaux,fddlt

      character spec*4,line*80
      integer*4 nrec,type,l
      integer*2 xgf2(9),xaf2(14),xxf2(24),adr2(14),aux(6),dlt(2)
      integer*4 xgf4(4),xaf4(7),xxf4(12),adr4(7)
      equivalence (xgf2,xgf4),(xgf2,xaf2),(xgf2,xaf4),
     |(xgf2,xxf2),(xgf2,xxf4),(xgf2,adr2),(xgf2,adr4)

      data istat /0/

      include "qgrid.inc"

      fddlt=-1
      fdaux=-1

      name=filenm
10    fdin=openf(filenm,'r')
      ios=readf(fdin,4,spec)
      ios=readf(fdin,4,nrec)
      if (spec.eq.'@D_H') then
	 fddlt=fdin
	 ios=readf(fddlt,16,xxf4)
	 l=index(filenm,'.xxo')
	 if (l.le.0) stop "Error in DLT filename"
	 filenm=filenm(:l+3)
	 goto 10
      else if (spec.eq.'@XGF') then
         type=1
         ios=seekf(fdin,18,0)
      else if (spec.eq.'@XAB') then
         type=2
	 l0=1
	 l1=nrec
	 call xtflimits(filenm,nint(t0),nint(t1),trk0,trk1,l0,l1)
	 ios=seekf(fdin,28*l0,0)
	 nrec=l1-l0+1
      else if (spec.eq.'@XXB') then
         type=3
         ios=seekf(fdin,48,0)
      else if (spec.eq.'@XXO') then
	 type=4
	 ios=seekf(fdin,44,0)
      else if (spec.eq.'aADR') then
	 type=5
	 ios=seekf(fdin,20,0)
	 ios=readf(fdin,4,nrec)
      else if (spec.eq.'xADR') then
	 type=6
	 ios=seekf(fdin,20,0)
	 ios=readf(fdin,4,nrec)
	 ios=seekf(fdin,28,0)
      else
         type=7
	 nrec=999999999
      endif

* If needed, open auxiliary file

      if (useaux) then
         l=index(filenm,' ')-1
         fdaux=openf(filenm(:l)//'.aux','r')
         if (fdaux.lt.0) stop "qgrid: no AUX file available"
         ios=readf(fdaux,12,xxf4)
      endif

      l=index(name,' ')-1
      if(nrec.ne.999999999) call statbar(istat,nrec,name(max(1,l-24):l))

      do irec=1,nrec
	 if (nrec.ne.999999999) call statbar(1,irec,' ')
         if (type.eq.1) then
            ios=readf(fdin,18,xgf2)
	    if (xgf2(9).lt.0) goto 100
	    if (xgf4(1).lt.t0 .or. xgf4(1).gt.t1) goto 100
            dlat=xgf4(2)/1d6
            dlon=xgf4(3)/1d6
	    if (field.eq.0) then
	       value=xgf4(4)/1d6
	    else
	       value=xgf2(9)/1d3
	    endif
	    range=value
         else if (type.eq.2) then
            ios=readf(fdin,28,xaf2)
	    if (xaf4(1).lt.t0 .or. xaf4(1).gt.t1) goto 100
            dlat=xaf4(2)/1d6
            dlon=xaf4(3)/1d6
	    if (field.eq.0) then
	       value=xaf4(4)/1d6
            else if (field.eq.1) then
	       if (xaf2(13).lt.0) goto 100
               value=xaf4(5)/1d6
            else
	       if (xaf2(13).lt.0) goto 100
               value=(xaf4(4)-xaf4(5))/1d6
            endif
	    range=value
         else if (type.eq.3) then
            ios=readf(fdin,48,xxf2)
	    if (xxf4(3).lt.t0 .or. xxf4(3).gt.t1) goto 100
	    if (xxf4(4).lt.t0 .or. xxf4(4).gt.t1) goto 100
            if (field.eq.0) then
               if (sum) then
                  value=(xxf4(6)+xxf4(7))/2d6
               else
                  value=(xxf4(6)-xxf4(7))/1d6
               endif
            else if (field.eq.1) then
               if (xxf2(23).lt.0 .or. xxf2(24).lt.0) goto 100
               if (sum) then
                  value=(xxf4(8)+xxf4(9))/2d6
               else
                  value=(xxf4(8)-xxf4(9))/1d6
               endif
	    else
               if (xxf2(23).lt.0 .or. xxf2(24).lt.0) goto 100
               if (sum) then
		  value=(xxf4(6)-xxf4(8)+xxf4(7)-xxf4(9))/2d6
               else
		  value=(xxf4(6)-xxf4(8)-xxf4(7)+xxf4(9))/1d6
               endif
            endif
	    if (tsort .and. .not.sum .and. xxf2(9).gt.xxf2(10))
     |			    value=-value
            dlat=xxf4(1)/1d6
            dlon=xxf4(2)/1d6
	    range=(xxf4(8)-xxf4(9))/1d6
	 else if (type.eq.4) then
	    ios=readf(fdin,44,xxf2)
	    if (xxf4(3).lt.t0 .or. xxf4(3).gt.t1) goto 100
	    if (xxf4(5).lt.t0 .or. xxf4(5).gt.t1) goto 100
	    if (useaux) then
	       ios=readf(fdaux,12,aux)
	       xxf4(8)=xxf4(8)+nint(ssb*aux(1)-tbias(1)*aux(5))
	       xxf4(9)=xxf4(9)+nint(ssb*aux(2)-tbias(2)*aux(6))
	    endif
	    if (fddlt.ge.0) then
	       ios=readf(fddlt,4,dlt)
	       if (dlt(1).eq.-32768 .or. dlt(2).eq.-32768) goto 100
	       xxf4(8)=xxf4(8)+dlt(1)*1000
	       xxf4(9)=xxf4(9)+dlt(2)*1000
	    endif
	    if (sum) then
	       value=(xxf4(8)+xxf4(9))/2d6
	    else
	       value=(xxf4(8)-xxf4(9))/1d6
	       if ((tsort .and. xxf2(13).gt.xxf2(14)) .or.
     |			xxf4(10)-xxf4(11).gt.100000000) value=-value
	    endif
            dlat=xxf4(1)/1d6
            dlon=xxf4(2)/1d6
	    range=(xxf4(8)-xxf4(9))/1d6
	 else if (type.eq.5) then
	    ios=readf(fdin,24,adr2)
	    if (adr4(1).lt.t0 .or. adr4(2).gt.t1) goto 100
	    dlat=adr4(3)/1d6
	    dlon=adr4(4)/1d6
	    value=adr2(11+field)/1d3
	    range=value
	 else if (type.eq.6) then
	    ios=readf(fdin,28,adr2)
	    if (adr4(1).lt.t0 .or. adr4(2).gt.t1) goto 100
	    dlat=adr4(3)/1d6
	    dlon=adr4(4)/1d6
	    value=adr2(11+field)/1d3
	    range=value
	 else if (type.eq.7) then
	    l=readf(fdin,0,line)-1
	    if (l.lt.0) goto 200
	    if (line(1:1).eq.'#' .or. line(:l).eq.' ') goto 100
	    value=0
	    read (line(:l),*,iostat=ios) t,dlat,dlon,value
	    if (t.lt.t0 .or. t.gt.t1) goto 100
	    range=value
         endif

	 if (ocean.and.sland(dlat,dlon)) then
	 else if (range0.le.range .and. range.le.range1) then
	    if (dlon.lt.x0) dlon=dlon+360
	    if (dlon.gt.x1) dlon=dlon-360
	    call convxy(dlon,dlat)
	    kx=nint((dlon-x0)/dx+1)
	    ky=nint((dlat-y0)/dy+1)
	    if (kx.ge.1 .and. kx.le.nx .and. ky.ge.1 .and. ky.le.ny) then
	       mean(kx,ky)=mean(kx,ky)+value
	       rms(kx,ky)=rms(kx,ky)+value**2
	       nr(kx,ky)=nr(kx,ky)+1
	    endif
	 endif
100      continue
      enddo
200   if (fdin .gt.0) call closef(fdin)
      if (fdaux.gt.0) call closef(fdaux)
      if (fddlt.gt.0) call closef(fddlt)
      return

      end

      subroutine wrap(jx,nx,ny,mean,rms,nr)
      implicit none
      integer nx,ny,nr(nx,ny),kx,ky,jx,k
      real*8 mean(nx,ny),rms(nx,ny)
      do kx=jx,nx
      do ky=1,ny
	 k=kx-jx+1
	 mean(k,ky)=mean(k,ky)+mean(kx,ky)
	 mean(kx,ky)=mean(k,ky)
	 rms(k,ky)=rms(k,ky)+rms(kx,ky)
	 rms(kx,ky)=rms(k,ky)
	 nr(k,ky)=nr(k,ky)+nr(kx,ky)
	 nr(kx,ky)=nr(k,ky)
      enddo
      enddo
      end

      subroutine fillup(nx,ny,mean,rms,nr,ifill)
      implicit none
      integer nx,ny,nr(nx,ny),ifill,kx,ky,n,m
      real*8 mean(nx,ny),rms(nx,ny),a,b
      do kx=1,nx
	 do ky=1,ny
	    if (nr(kx,ky).eq.0) then
	       n=0
	       m=0
	       a=0
	       b=0
	       if (nr(kx-1,ky  ).gt.0) then
	          a=a+mean(kx-1,ky  )
		  b=b+rms (kx-1,ky  )
		  m=m+nr  (kx-1,ky  )
		  n=n+1
	       endif
	       if (nr(kx+1,ky  ).gt.0) then
	          a=a+mean(kx+1,ky  )
		  b=b+rms (kx+1,ky  )
		  m=m+nr  (kx+1,ky  )
		  n=n+1
	       endif
	       if (nr(kx  ,ky-1).gt.0) then
	          a=a+mean(kx  ,ky-1)
		  b=b+rms (kx  ,ky-1)
		  m=m+nr  (kx  ,ky-1)
		  n=n+1
	       endif
	       if (nr(kx  ,ky+1).gt.0) then
	          a=a+mean(kx  ,ky+1)
		  b=b+rms (kx  ,ky+1)
		  m=m+nr  (kx  ,ky+1)
		  n=n+1
	       endif
	       if (n.ge.ifill) then
	          mean(kx,ky)=a/m
	          rms(kx,ky)=b/m
 	          nr(kx,ky)=-1
	       endif
	    endif
	 enddo
      enddo
      end

      subroutine toolittle (nx,ny,mean,rms,nr,minnr)
      integer nx,ny,nr(*),minnr,k
      real*8 mean(*),rms(*)
      do k=1,nx*ny
	 nr(k)=abs(nr(k))
	 if (nr(k).ge.minnr) then
	    rms(k)=sqrt(rms(k)/nr(k))
	    mean(k)=mean(k)/nr(k)
	 else
	    mean(k)=1d35
	    rms(k)=1d35
	 endif
      enddo
      end

      subroutine sigma(nx,ny,mean,rms,nr,minnr)
      integer nx,ny,nr(*),minnr,k
      real*8 mean(*),rms(*)
      do k=1,nx*ny
	 if (nr(k).ge.minnr) then
	    rms(k)=sqrt(rms(k)**2-mean(k)**2)
	 else
	    rms(k)=1d35
	 endif
      enddo
      end

      subroutine convxy(x,y)

* Adjusted version of PMCONV

      include 'qgrid.inc'
      real*8 x,y,loncen,a,f
      parameter (loncen=0)
*
      if (ptype.le.1) then
         return
      else if (ptype.eq.2) then
         y=factk*sin(y*rad)
      else if (ptype.eq.3 .or. ptype.eq.4) then
         y=factk*log(tan(qpi+y*factl))
      else if (ptype.eq.41) then
         f=(90-y)*factk
         a=(x-loncen)*rad
         x= f*sin(a)
         y=-f*cos(a)
      else if (ptype.eq.42) then
         f=(y+90)*factk
         a=(x-loncen)*rad
         x=f*sin(a)
         y=f*cos(a)
      endif
      return
      end

      subroutine setproj(project)
      include "qgrid.inc"
      integer project
      real*8 cor1,ppara1

* Adjusted version of PMSWIN

      ptype=project
      ppara1=(y0+y1)/2
      cor1=cos(ppara1*rad)
      qpi=atan(1d0)
      rad=qpi/45
      if (ptype.eq.2) then
         factk=1/cor1**2/rad
      else if (ptype.eq.3) then
	 factl=rad/2
         factk=1/rad
      else if (ptype.eq.4) then
         factl=0.4d0*rad
	 factk=1.25d0/rad
      else if (ptype.eq.41 .or. ptype.eq.42) then
         factk=cor1/(abs(ppara1)*rad)
      else
         ptype=1
      endif
      end
