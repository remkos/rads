      subroutine cmdgpstr
      include "pim.inc"
      logical pimopt,l

* Call GPSTRACK if requested

      call pgsave
10    if (pop1cmd('GPSTR',argum)) then
	 ci=-999
	 ls=1
	 lw=1
         xname='file.plt'
	 l=pimopt('ci=',argum,ci,dum,dum,dum)
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 call strip(argum,xname)
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call gpstrack(xname,nint(ci))
	 goto 10
      endif
      call pgunsa
      end

**GPSTRACK -- Plot Track Point file (OziExplorer)
*+
      SUBROUTINE GPSTRACK (FILENM, IDX)
      implicit none
      CHARACTER*(*) FILENM
      INTEGER IDX
*
* Routine plots the track from a GPS track file. For the time being it only
* supports the OziExplorer track file format. (.plt)
* Tracks can take colour according to the velocity (IDX=-999) or a specified
* colour (IDX>=0).
*
* Arguments:
*  FILENM  (input): Name of the track point file
*  IDX     (input): Color index for plotting the track. Use -999 to
*                   plot colors according to track velocity
*-
      include "pim.inc"

      integer i,j,unit,freeunit,comma(10),ndata,pen
      real v,vfact,cmax
      real*8 rad,dlon,dlon0,dlat,dlat0,time,time0,re,sfdist,height,vel,
     |		dist
      character*80 line
      parameter (re=6375d0)

* Initialize

      rad=atan(1d0)/45
      vfact=(nc1+1)/(rmaxs-rmins)
      cmax=nc1+0.999

* Open GPS track file on new unit

      i=index(filenm,' ')-1
      write (0,600) filenm(:i)
550   format (a)
600   format ('Plotting GPS track ',a,' ...')
610   format ('... ',a)
      unit=freeunit()
      open (unit,file=filenm,status='old')
      !
      ! Check file format
      ! And read header
      !
      read (unit,550) line
      if (line(1:40).eq.'OziExplorer Track Point File Version 2.0') then
         read (unit,*)
         read (unit,*)
         read (unit,*)
	 read (unit,610) line
	 j=0
	 do i=1,80
	    if (line(i:i).eq.',') then
	       j=j+1
	       comma(j)=i
	    endif
	 enddo
	 if (j.ge.4) then
	    write (0,550) line(comma(3)+1:comma(4)-1)
	 endif
	 read (unit,*) ndata
      else
	 write (0,550) '... wrong file type. Exit'
	 return
      endif

* Set colour index if specified

      if (idx.ge.0) call pgsci(idx)

* Plot track

      do i=1,ndata
         read (unit,*,end=9999) dlat,dlon,pen,height,time
	 x=dlon
	 y=dlat
	 call pmconv(1,x,y)
	 if (pen.eq.1) then	! pick up pen
	    call pgmove(x,y)
	 else
	    if (idx.lt.0) then
	       dist=re*sfdist(dlat0*rad,dlon0*rad,dlat*rad,dlon*rad)
	       vel=dist/(time-time0)/24
	       ! write (*,*) i,dlat0,dlat,dlon0,dlon,dist,vel
	       v=(vel-rmins)*vfact
	       j=c_0+max(0.,min(v,cmax))
	       call pgsci(j)
	       time0=time
	       dlon0=dlon
	       dlat0=dlat
	    endif
	    call pgdraw(x,y)
	 endif
      enddo
9999  end
