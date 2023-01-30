      subroutine cmdwaypo
      include "pim.inc"
      real style,mrk,ci1,ci2
      logical pimopt,l

* Call WAYPOINT if requested

      call pgsave
10    if (pop1cmd('WAYPO',argum)) then
	 ci1=1
	 ci2=0
	 ls=1
	 lw=1
	 ch=def_ch
	 style=1
	 mrk=3
	 xname='file.wpt'
	 l=pimopt('ci=',argum,ci1,ci2,dum,dum)
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('ch=',argum,ch,dum,dum,dum)
	 l=pimopt('mrk=',argum,mrk,dum,dum,dum)
	 l=pimopt('style=',argum,style,dum,dum,dum)
	 call strip(argum,xname)
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgsch(ch)
	 call waypoint(xname,nint(style),nint(ci1),nint(ci2),nint(mrk))
	 goto 10
      endif
      call pgunsa
      end

**WAYPOINT -- Plot Waypoint file (OziExplorer)
*+
      SUBROUTINE WAYPOINT (FILENM, STYLE, CI1, CI2, MARKER)
      CHARACTER*(*) FILENM
      INTEGER STYLE, CI1, CI2, MARKER
*
* Routine plots markers and names for waypoints from a GPS waypoint file.
* For the time being it only supports the OziExplorer waypoint file format. (.wpt)
* Waypoints can be plotted with the standard 16 GPS markers, or with a fixed
* marker type.
*
* Input arguments:
*  FILENM  : Name of the waypoint file
*  STYLE   : Specify type of name to be plotted:
*            0 = No name (only marker)
*            1 = Abbreviated (6-character) name
*            2 = Full name
*  CI1     : Colour index for the text
*  CI2     : Colour index for the background box (if CI2 < 0 no box is drawn)
*  MARKER  : Marker type or -999 for standard GPS markers
*-
      include "pim.inc"

      integer i,j,unit,freeunit,mark
      real*8 dlon,dlat,time
      real*4 angle/0.0/,just/0.5/,xbox(4),ybox(4),dy
      character*160 line,text

* Open GPS track file on new unit

      i=index(filenm,' ')-1
      write (0,600) filenm(:i)
550   format (a)
600   format ('Plotting GPS waypoints ',a,' ...')
      unit=freeunit()
      open (unit,file=filenm,status='old')
      !
      ! Check file format
      ! And read header
      !
      read (unit,550) line
      if (line(1:37).eq.'OziExplorer Waypoint File Version 1.0') then
         read (unit,*)
         read (unit,*)
         read (unit,*)
      else
	 write (0,550) '... wrong file type. Exit'
	 return
      endif

* Set foreground colour index

      call pgsci(ci1)

* Plot track

100   read (unit,550,end=9999) line
      read (line(21:58),*) dlat,dlon,time,mark
      x=dlon
      y=dlat
      call pmconv(1,x,y)
      if (marker.eq.-999) then
         mark=mark+0
      else
         mark=marker
      endif
      if (ci2.ge.99999) then
         call pgqtxt(x,y,angle,just,'AAAAA',xbox,ybox)
	 dy=ybox(3)-ybox(1)
	 call pgsci(ci2)
	 call pgsfs(1)
	 call pgrect(x-0.5*dy,x+0.5*dy,
     |		     y-0.5*dy,y+0.5*dy)
	 call pgsci(ci1)
	 call pgsfs(0)
	 call pgrect(x-0.5*dy,x+0.5*dy,
     |		     y-0.5*dy,y+0.5*dy)
      endif
      call pgpt(1,x,y,mark)
      if (style.eq.1) then
         text=line(6:19)
      else if (style.eq.2) then
         text=line(88:127)
      else
         text=' '
      endif
      if (text.ne.' ') then
         call pgqtxt(x,y,angle,just,text,xbox,ybox)
	 dy=ybox(3)-ybox(1)
	 if (ci2.ge.0) then
	    call pgsci(ci2)
	    call pgsfs(1)
	    call pgrect(xbox(1)-0.1*dy,xbox(3)+0.1*dy,
     |		        ybox(1)+0.9*dy,ybox(3)+1.1*dy)
	    call pgsci(ci1)
	    call pgsfs(0)
	    call pgrect(xbox(1)-0.1*dy,xbox(3)+0.1*dy,
     |		        ybox(1)+0.9*dy,ybox(3)+1.1*dy)
         endif
         call pgptxt(x,y+dy,angle,just,text)
      endif
      goto 100
9999  end
