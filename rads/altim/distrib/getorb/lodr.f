**LODR -- Program to list contents of ODR file
*
* This program prints the contents of an orbital data records file (ODR)
* in ASCII to standard output. For a full description of the syntax
* type 'lodr' at the command prompt.
*-
* $Log: lodr.f,v $
* Revision 1.5  2015/07/10 14:30:24  rads
* Cut lines before or at 72nd character.
*
* Revision 1.4  2010/11/15 18:10:05  rads
* - Expanded to 4-digit ODR file numbering (arc_????.odr)
*
* Revision 1.3  2006/01/10 14:19:15  rads
* - Removed FASTIO; reverted to FORTRAN I/O
*
* Revision 1.2  2004/11/22 14:32:19  remko
* - Allow arc numbers higher than 999
*
* 10-Jan-2006 - Removed fastio
*  3-Nov-1999 - Added DATEARG function
* 20-Oct-1996 - Supports new ODR format.
* 16-Nov-1992 - Created by Remko Scharroo, DUT/DEOS
*-----------------------------------------------------------------------
      program	lodr
      integer	i,iargc,unit,ios,n0,n1
      integer	rep,begin,end,time0,time1,tstep,ver,arc,nrec
      integer	rec(4),odrinfo
      character filenm*80,arg*80,satel*8,spec*4
      character date*15,date0*15,date1*15
      real*8	t0,t1,mu,rev,dum
      logical   summary,datearg,swap,ltlend
      data	filenm/' '/,n0/0/,n1/0/,t0/-1d40/,t1/1d40/,
     |		summary/.false./

* Scan arguments

      do i=1,iargc()
         call getarg(i,arg)
         if (arg(1:2).eq.'n=') then
            read (arg(3:),*,err=1300) n0,n1
         else if (datearg(arg,t0,t1,dum)) then
	 else if (arg(1:2).eq.'-s') then
	    summary=.true.
         else
            filenm=arg
         endif
      enddo
      if (filenm.eq.' ') goto 1300

* Open ODR and get header information

      if (odrinfo(unit,filenm,satel,rep,arc,ver,nrec,
     |time0,time1,tstep,begin,end,rev).eq.0) then
	 spec='@ODR'
	 mu=1d-6
      else
	 spec='xODR'
	 mu=1d-7
      endif

* Convert dates

      call chrdat(time0,date0)
      call chrdat(time1,date1)
      call chrdat(begin,date)

* Write summary information

      if (.not.summary) then
      else if (index(filenm,'ODR.').gt.0) then
	 write (*,610) arc/100,mod(arc,100),date0,date1,satel,rep/1d3,
     |	 	ver,date
	 goto 9999
      else
	 write (*,611) arc,date0,date1,satel,rep/1d3,ver,date
	 goto 9999
      endif

* Write output header

      write (*,620) spec,satel,rep/1d3,arc,nrec,ver,date
      if (n0.le.0 .or. n0.gt.nrec) n0=1
      if (n1.le.0 .or. n1.gt.nrec) n1=nrec
      swap=ltlend()

* Write data records

      do i=n0,n1
	 read (unit,rec=i+2,iostat=ios) rec
         if (ios.ne.0) call fin("lodr: premature end of file")
         if (swap) call i4swap(4,rec)
         if (rec(1).ge.t0 .and. rec(1).le.t1) then
            call chrdat(rec(1),date)
            write (*,630) date,rec(2)*mu,rec(3)*mu,rec(4)/1d3
         endif
      enddo
      goto 9999

* Formats

610   format (z1.1,i2.2,2x,a12,' - ',a12,2x,a8,2x,f7.3,i4,2x,a15)
611   format (i4.4,1x,a12,' - ',a12,2x,a8,2x,f7.3,i4,2x,a15)
620   format (a4,1x,a8,2x,'Rep:',f7.3,2x,'Arc:',i6,2x,'Rec:',i5,
     |2x,'Ver:',i5,2x,'Beg:',a15)
630   format (a15,2f13.7,f12.3)

* Errors

1300  write (*,1301)
1301  format ('lodr: list contents of Orbital Data Records (ODR).'//
     |'usage: lodr [options] ODR-filename'//
     |'where [options] are'/
     |' -s      : give only a summary'/
     |' n=n0,n1 : list only record numbers n0 through n1'/
     |' t=t0,t1 : list only record of period t0 through t1, where'/
     |'           t0 and t1 may be in YYMMDD, MJD, or YYDDD'/
     |'       ... or use mjd=, doy=, ymd=, sec=')
9999  end
