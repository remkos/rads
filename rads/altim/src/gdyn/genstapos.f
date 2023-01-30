      program genstapos

* This program generates Geodyn station position cards for a GIIS setups
* based on input provided about:
* - which stations
* - epoch
* - station XYZ coordinates
* - station aliases
* - eccentricities
* - adjusted coordinates
* - Nuvel1A plate model
*
*-
*  9-Sep-1998 - Created by Remko Scharroo
*  3-Nov-1999 - Added DATEARG function
*  1-Feb-2000 - Allow comment lines in coordinate set
* 15-May-2000 - No fatal error on missing station. Just warning.
* 14-Mar-2001 - Ignore empty lines in input. If t= not given, take
*               epoch from the first input line
* 13-Apr-2001 - Use new subroutine LOADCOORD
* 10-May-2001 - Repeat station id in columns 67-70
*  8-Aug-2001 - Include time limits in adjustments file
* 13-Aug-2002 - Read alias file (Eelco Doornbos)
* 28-Feb-2004 - Print warnings on standard error not standard out
*-----------------------------------------------------------------------
      implicit	none
      integer*4	msta
      parameter (msta=1000)
      integer*2 pnt(9999)
      logical	datearg,lees,l
      real*8	s_pos(6,msta),s_sig(6,msta)
      integer*4	s_id(msta),s_plate(msta),noise,s_wave(msta)
      character s_code(msta)*6
      real*8	ecc(3),lat,lon,height,r,epoch,epoch_ref,dum,sec85,year,
     |		adjust,date
      parameter (year=365.25d0*86400d0)
      integer*4	i,j,k,p,nsta,iargc,statinfo,mjdref,mdate
      character*80 text,xyz,adj,aliaslst
      integer*4 maxalias, nalias, findalias
      parameter (maxalias=50)
      integer*4 alias(maxalias,4)

* Initialise

      do i=1,9999
         pnt(i)=0
      enddo
      xyz=' '
      adj=' '
      aliaslst=' '
      epoch=-1d20

* Scan arguments

600   format('genstapos - Generate GIIS STAPOS cards'//
     |'syntax: genstapos [t=epoch] xyz adj [alias=aliaslst] '//
     |'< list > cards where'/
     |'  t=epoch  : epoch of the setup (MJD, [YY]YYMMDD[HHMMSS])'/
     |'             (where t= is optional)'/
     |'         ... or use mjd=, doy=, ymd=, sec='/
     |'             (if omitted: epoch must be first line of input)'/
     |'  xyz      : file containing XYZ coordinates and velocities'/
     |'  adj      : list of station adjustments'/
     |'  aliaslst : list of station aliases'/
     |'  list     : list of station numbers (read from stdin)'/
     |'  cards    : STAPOS cards (written to stdout)')

      do i=1,iargc()
         call getarg(i,text)
         if (text(:2).eq.'-h') then
            write (0,600)
	    goto 9999
         else if (datearg(text,epoch,dum,dum)) then
         else if (xyz.eq.' ') then
            xyz=text
         else if (adj.eq.' ') then
            adj=text
         else
           aliaslst = text
           call readalias(aliaslst,alias,nalias)
         endif
      enddo

* If xyz or adj not given, write syntax

      if (adj.eq.' ') then
         write (0,600)
	 goto 9999
      endif

* If time not given, read it from the first line of standard input

      if (epoch.lt.0) then
         l=lees(5,text)
         read (text,*) epoch
	 epoch=sec85(0,epoch)
      endif
      epoch=epoch/86400d0+46066d0

* Read the XYZ coordinates and velocities from the requested file

      call loadcoord (xyz,mjdref,nsta,s_id,s_pos,s_sig)
      epoch_ref=mdate(1,mjdref)

* Set pointers (only set those that are zero)
* Reduce velocities from m/year to m/s.
* Default adjustment to zero.

      do i=1,nsta
         if (pnt(s_id(i)).eq.0) pnt(s_id(i))=i
	 do j=4,6
	    s_pos(j,i)=s_pos(j,i)/year
	 enddo
	 s_sig(1,i)=0
	 s_plate(i)=9876
      enddo

* Read the adjustment information and store it in s_sig(1,pnt(k))

      open (11,file=adj,status='old')
60    if (lees(11,text)) goto 70
      date=0
      read (text,*,err=65) k,adjust,date
      date=sec85(0,date)/86400d0+46066d0
65    if (pnt(k).ne.0 .and. date.le.epoch) s_sig(1,pnt(k))=adjust
      goto 60
70    continue

* Read the list of stations.
* Check if station is in coordinate file and system.data.
* Compute the optical centre coordinates by applying the eccentricities
* Flip pointer from positive to negative for successful stations.

10    if (lees(5,text)) goto 20
      read (text,*) k
      p=k
      if(aliaslst(1:1).ne.' ') 
     | k=findalias(k,nint(epoch),alias,nalias)
      i=pnt(k)
      if (i.eq.0) then
         write (0,1302) k
      else
	 j=statinfo(epoch,p,text,s_code(i),dum,
     |		noise,s_wave(i),s_plate(i),ecc)
	 if (j.ne.0) then
	    write (0,1301) k
	 else
	    call xyzgeo(s_pos(1,i),r,lat,lon,height)
	    s_pos(1,i)=s_pos(1,i)
     |		-ecc(1)*sin(lat)*cos(lon)-ecc(2)*sin(lon)
     |		+ecc(3)*cos(lat)*cos(lon)
	    s_pos(2,i)=s_pos(2,i)
     |		-ecc(1)*sin(lat)*sin(lon)+ecc(2)*cos(lon)
     |		+ecc(3)*cos(lat)*sin(lon)
	    s_pos(3,i)=s_pos(3,i)
     |		+ecc(1)*cos(lat)+ecc(3)*sin(lat)
         endif
         pnt(k)=-i
      endif
      goto 10
20    continue
1301  format ('WARNING: Station ',i4.4,' not found in system.data. ',
     |'Station ignored. Please edit system.data')
1302  format ('WARNING: Station ',i4.4,' not found in station ',
     |'coordinate file. Station ignored.')

* Write the STAPOS cards. First the fixed sites, then the adjusted.
* 10-May-2001: Repeat station id in columns 67-70

      write (*,120) 'STAPOS 1'
      write (*,120) 'ELCUTOFF            10.0'
      write (*,120) 'FIXED'
      do k=1,9999
	 i=-pnt(k)
         if (i.gt.0 .and. s_sig(1,i).eq.0) then
           p=findalias(-k,nint(epoch),alias,nalias)
           write (*,100) s_wave(i)/1d3,s_code(i),s_plate(i),
     |		p,(s_pos(j,i),j=1,3),p,p,(s_pos(j,i),j=4,6),
     |		p,epoch_ref
         endif
      enddo
100   format ('INSTRMNT       0',t36,f5.3/
     |a6,'  0',i3,i8,3f15.6,1x,i4/'STAVEL',6x,i8,3d15.6/
     |'TIMVEL',6x,i8,f15.6)
110   format ('ADJUSTED2',t21,3d15.6)
120   format (a)
      do k=1,9999
	 i=-pnt(k)
         if (i.gt.0 .and. s_sig(1,i).ne.0) then
            p=findalias(-k,nint(epoch),alias,nalias)
	    write (*,110) (s_sig(1,i),j=1,3)
	    write (*,100) s_wave(i)/1d3,s_code(i),s_plate(i),
     |		p,(s_pos(j,i),j=1,3),p,p,(s_pos(j,i),j=4,6),
     |		p,epoch_ref
         endif
      enddo
      write (*,120) 'ENDSTA'
         
9999  end

********************************************************************************

      function lees(unit,text)
      integer unit
      character*(*) text
      logical lees
      lees=.true.
10    read (unit,'(a)',end=9999) text
      if (text(:1).eq.'#' .or. text.eq.' ') goto 10
      lees=.false.
9999  end

********************************************************************************

      subroutine readalias(filename,alias,nalias)

      character*(*) filename
      integer*4 maxalias 
      parameter (maxalias = 50)
      integer*4 alias(maxalias,4)
      integer*4 i, nalias, mdate
      logical lees, l
      character*80 line

      nalias=0
      open (11,file=filename,status='old')

   10 l=lees(11,line)
      if(.not.l) then
        nalias = nalias+1
        read(line(6:),*) (alias(nalias,i), i=1,4)
        alias(nalias,3) = mdate(2,alias(nalias,3))
        alias(nalias,4) = mdate(2,alias(nalias,4))
      else
        goto 9999
      endif
      goto 10

 9999 close(11)

      end

********************************************************************************

      function findalias(id,epoch,list,n)

* Use the true station id to find the alias id 
* Use the negative value of the station alias id to get the true station id
* If no alias is found, findalias returns id

      integer*4 maxalias, n, a, b
      parameter (maxalias = 50)
      integer*4 findalias,id,list(maxalias,4), epoch, i

      if(id.lt.0) then
        a=2
        b=1
        id = -id
      else
        a=1
        b=2
      end if

      findalias = id
 
      do i=1,n
        if(id.eq.list(i,a).and.epoch.gt.list(i,3).and.
     |                         epoch.lt.list(i,4)) then
          findalias = list(i,b)
        endif
      enddo
 
      end
