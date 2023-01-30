**MAKETEX -- Strip manual from C or Fortran source
*+
      program maketex
*
* This program strips the manual out of a C or Fortran source. The output
* is in LaTeX (when called as manstrip) or HTML (when called as htmlstrip).
* Put the following characters in the second column of your C or Fortran
* code:
* * Indicates the (sub)program title
* + Is the line just above the (sub)program call sequence. All text
*   thereafter is included in the output until:
* - Indicates the end of the verbose.
* = Stops scanning altogether
* & Restarts verbose
*
*-
* 15-Mar-2008 - Made compatible with Fortran90
*  9-Aug-2000 - Initialise count to 0
* 16-Feb-1998 - Remko Scharroo
*-----------------------------------------------------------------------
      character*80 text,filenm
      character*1 backslash
      logical verb,copy,tex
      integer count/0/,iargc,i,m,iarg,lnblnk,l

      backslash=char(92)

      call getarg(0,filenm)
      tex=index(filenm,'maketex').gt.0

      do iarg=1,iargc()
	 call getarg(iarg,filenm)
	 open (10,file=filenm,status='old')
	 rewind (10)
	 verb=.false.
	 copy=.false.
   10    read (10,550,end=100) text
	 l=lnblnk(text)
	 m=3
	 if (text(2:2).eq.'*') then
	    count=count+1
            if (.not.verb) then
	    else if (tex) then
	       write (6,620) backslash
	    else
	       write (6,820)
	    endif
	    verb=.false.
	    if (tex) then
	       i=index(text,'$')
	       if (i.ne.0) text(i:i)='.'
	       i=index(text,'_')
	       if (i.ne.0) text(i:i)='.'
	       i=index(text,'_')
	       if (i.ne.0) text(i:i)='.'
	       write (6,600) backslash,text(m:l)
	    else
	       write (6,800) count,text(m:l)
	    endif
	 else if (text(2:2).eq.'&') then
	    count=count+1
            if (.not.verb) then
	    else if (tex) then
	       write (6,620) backslash
	    else
	       write (6,820)
	    endif
	    verb=.false.
	    if (tex) then
	       i=index(text,'$')
	       if (i.ne.0) text(i:i)='.'
	       i=index(text,'_')
	       if (i.ne.0) text(i:i)='.'
	       i=index(text,'_')
	       if (i.ne.0) text(i:i)='.'
	       write (6,601) backslash,text(m:l)
	    else
	       write (6,801) count,text(m:l)
	    endif
	 else if (text(2:2).eq.'+') then
	    if (tex) then
	       write (6,610) backslash
	    else
	       write (6,810)
	    endif
	    verb=.true.
	    copy=.true.
	 else if (text(2:2).eq.'-') then
	    copy=.false.
	 else if (text(2:2).eq.'=') then
	    goto 100
	 else if (copy) then
	    if (text(2:2).ne.' ') m=1
	    if (l.lt.m) then
	       write (6,550)
	    else
	       write (6,550) text(m:l)
	    endif
	 endif
	 goto 10
  100    if (.not.verb) then
	 else if (tex) then
	    write (6,620) backslash
	 else
	    write (6,820)
	 endif
      enddo
  550 format (a)
  600 format (a1,'SECTION{',a,'}')
  601 format (a1,'SUBSECTION{',a,'}')
  610 format (a1,'begin{verbatim}')
  620 format (a1,'end{verbatim}'/)
  800 format ('<a name=',i4.4,'><h2>',a,'</h2>')
  801 format ('<a name=',i4.4,'><h3>',a,'</h3></a>')
  810 format ('<pre>')
  820 format ('</pre><hr>')
      end
