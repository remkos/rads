      program ft7merge

* This program will merge a Geodyn input file (giis.in.c) with a
* parameter solution file (fort.7). Solutions for EPOCH, ELEMS, DRAG, SOLRAD,
* and ACCEL9 will be picked from the last solution fort.7 and will replace
* the values in the giis.in.c file.
*
*-
*        1998 - Created by Remko Scharroo
*  4-Apr-2000 - Last change: include -a and -g flags
*-----------------------------------------------------------------------

      character*80 giis,ft7,line,sline(1000,6),command(6)*6
      data command /'EPOCH','ELEMS','DRAG','SOLRAD','ACCEL9','STAPOS'/
      integer i,l,n,type,lnblnk,unit,glob,adjusted,mline(6),nline(6),
     |		iargc
      logical error,global_update,arc_update
      real*8  sigma,cd,cd0,cdmin,cdmax,value
      data mline/6*0/,nline/6*0/,glob/0/,adjusted/0/
      data giis/' '/,ft7/' '/,global_update/.true./arc_update/.true./

* Read arguments

      do i=1,iargc()
         call getarg(i,line)
	 if (line(:2).eq.'-h') then
	    goto 1300
	 else if (line(:2).eq.'-a') then
	    global_update=.false.
	 else if (line(:2).eq.'-g') then
	    arc_update=.false.
	 else if (giis.eq.' ') then
	    giis=line
	 else
	    ft7=line
	 endif
      enddo

* Open files

      if (giis.eq.' ') giis='giis.in.c'
      if (ft7.eq.' ') ft7='fort.7'

      open (11,file=giis,status='old',err=1300)
      open (12,file=ft7,status='old',err=1300)

* Scan fort.7 for ELEM, DRAG, SOLRAD, ACCEL9 cards and keep them in memory.

100   read (12,550,end=200) line
      type=0
      if (line(1:5).eq.'/ARC/') then
         glob=0
	 do i=1,5
	    nline(i)=0
	 enddo
	 read (12,550)
	 read (12,550)
	 read (12,550)
      else if (line(1:6).eq.'/GLOB/') then
         nline(6)=0
	 glob=1
	 read (12,550)
	 read (12,550)
	 read (12,550)
      else if (line(1:5).eq.'EPOCH') then
         type=1
      else if (line(1:5).eq.'ELEMS') then
         type=2
      else if (line(1:4).eq.'DRAG') then
         type=3
      else if (line(1:6).eq.'SOLRAD') then
         type=4
      else if (line(1:6).eq.'ACCEL9') then
         type=5
         read (line(60:72),'(d13.1)') sigma
	 if (sigma.le.0d0) goto 100 ! Skip the fully constrained ACCEL9 cards.
      else if (line(1:6).eq.'TIMVEL') then
      else if (line(1:6).eq.'STAVEL') then
      else if (glob.ge.1) then
	 type=6
      endif

      if (type.gt.0) then
         nline(type)=nline(type)+1
         sline(nline(type),type)=line
      endif

      goto 100
200   continue
      close (12)

* Compute the average of the DRAG cards; only the ones that have been
* changed, so not the ones that are equal to the first DRAG card.
* Also, we remove the maximum and minimum number to compute the average.

      n=0
      cd0=0
      cdmin=1d40
      cdmax=-1d40
      do i=1,nline(3)
         read (sline(i,3),"(24x,d20.8,15x,d13.1)") value,sigma
	 if (sigma.eq.0d0) then
	    cd0=value
	 else if (abs(value-cd0).gt.1d-8) then
	    n=n+1
	    cd=cd+value
	    cdmin=min(cdmin,value)
	    cdmax=max(cdmax,value)
	 endif
      enddo
      if (n.gt.4) then
         cd=(cd-cdmin-cdmax)/(n-2)
      else
         cd=cd/n
      endif

* Go through the giis.in.c file and scan for lines that could be
* replaced by lines from fort.7 stored in memory (sline). The parameter
* 'type' indicates the type of parameter. The default is '0', meaning
* 'no parameter substitution for this line'. A stting of '.false.' for
* either 'arc_update' or 'global_update' may prevent updating of arc or
* global parameters.

      unit=11
210   read (unit,550,end=300) line
      type=0
      if (line(1:5).eq.'EPOCH') then
         if (arc_update) type=1
      else if (line(1:5).eq.'ELEMS') then
         if (arc_update) type=2
      else if (line(1:4).eq.'DRAG') then
         if (arc_update) type=3
      else if (line(1:6).eq.'SOLRAD') then
         if (arc_update) type=4
      else if (line(1:6).eq.'ACCEL9') then
         if (arc_update) type=5
         read (line(60:72),'(d13.1)') sigma
	 if (sigma.le.0d0) type=0 ! Do not replace fully constrained ACCEL9 cards.
      else if (line.eq.'#include <giis.common>') then
         mline(6)=nline(6)
      else if (line.eq.'#include <stapos.in>') then
         mline(6)=nline(6)
      else if (line.eq.'#include "accel9.ers1"') then
	 unit=unit+1
         open (unit,file='accel9.ers1',status='old')
	 type=-1
      else if (line.eq.'#include "accel9.ers2"') then
	 unit=unit+1
         open (unit,file='accel9.ers2',status='old')
	 type=-1
      else if (line.eq.'#include "accel9.in"') then
	 unit=unit+1
         open (unit,file='accel9.in',status='old')
	 type=-1
      else if (line(1:8).eq.'ADJUSTED') then
         adjusted=1
      else if (line(1:6).eq.'ENDSTA') then
         adjusted=0
      else if (adjusted.gt.0 .and. line(1:6).ne.'STAVEL'
     |		.and. line(1:6).ne.'TIMVEL') then
         do i=1,nline(6)
	    if (line(17:20).eq.sline(i,6)(17:20)) then
	       mline(6)=i
	       if (global_update) type=6
	    endif
	 enddo
      endif

* If 'type' is not '0' then a line has been found that ought to be
* replaced by a line from memory (stored in 'sline').

      if (type.ge.1 .and. type.le.4) then
         mline(type)=mline(type)+1
         line=sline(mline(type),type)
      else if (type.eq.5) then
         mline(type)=mline(type)+1
         line(11:)=sline(mline(type),type)(11:)
      else if (type.eq.6) then
         line=sline(mline(type),type)
      endif
      if (type.ge.0) then
         l=lnblnk(line)
         write (*,550) line(:l)
      endif
      goto 210

* Close the current open unit. If the unit number is greater than 11, we
* will have decended into an #include file. If this is the case,
* continue with the rest of the file from which the #include was called.

300   continue
      close (unit)
      if (unit.gt.11) then
         unit=unit-1
	 goto 210
      endif

* Check whether all replacements have been done

      do type=1,6
         if (type.ge.5) then
	 else if (arc_update .and. nline(type).ne.mline(type)) then
	    error=.true.
	    write (0,600) command(type),mline(type),nline(type)
	 endif
      enddo
550   format(a)
600   format('ft7merge: incompatible number of ',a,' cards in ',
     |'giis and ft7 file: ',2i3)

* End of the program execution

      goto 9999

* Errors

1300  write (*,1301)
1301  format ('FT7MERGE -- Update giis.in with fort.7 file'//
     |'usage:  ft7merge [options] [giis.in.c fort.7]'//
     |'where [options] are:'/
     |' -h        : print this help message'/
     |' -a        : perform update of ARC parameters ONLY'/
     |' -g        : perform update of GLOBAL parameters ONLY'/
     |' giis.in.c : name of the GEODYN II setup file'/
     |' fort.7    : name of the GEODYN II parameter solution file'//
     |'If filenames are ommitted the defaults (giis.in.c and fort.7)',
     |' will be used.'/'This program replaces only STAPOS, EPOCH,',
     |' ELEMS, DRAG, SOLRAD, and ACCEL9 cards.'/
     |'Output (a new giis.in.c file) comes to standard output.')
9999  end
