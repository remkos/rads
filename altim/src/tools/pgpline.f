      program pgpline

* Compresses LP commands in PGPLOT PostScript output.
* - Removes "0 0 LP"
* - Joins "1 -1 LP 1 -1 LP" to "2 -2 LP", etc.
*-
* 27-Aug-2001 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      character*256 line
      integer l,lnblnk

10    read (*,'(a)',end=900) line
      l=lnblnk(line)
      if (line(:1).eq.'%') then
         call pr(line,l)
      else if (index(line,'LP').le.0) then
         call pr(line,l)
      else
         call scanline(line,l)
      endif
      goto 10

900   end

      subroutine pr(line,l)
      character*(*) line
      integer l
      write (*,'(a)') line(:l)
      end

      subroutine scanline(line,l)
      character*(*) line
      integer i,last,l,n,x,x0,y,y0,sp
      real*8  f,g
      x0=0
      y0=0
      i=0
      n=0
      last=0
      sp=0
10    i=i+1
      if (i.gt.l) then
         goto 100
      else if (line(i:i+1).eq.'LP') then
         if (i-last.ge.5 .and. (sp.ge.3 .or.(sp.eq.2.and.last.eq.0)))
     |   then
            read (line(last+1:i-1),*) x,y
	    f=atan2(dble(y),dble(x))
	    g=atan2(dble(y0),dble(x0))
	    if (x.eq.0 .and. y.eq.0) then
	    else if (f.eq.g) then
	       x0=x0+x
	       y0=y0+y
	    else
	       call prlp(x0,y0,n)
	       x0=x
	       y0=y
	    endif
         else
            call prlp(x0,y0,n)
	    write (*,'(a,$)') line(last+1:i+1)
	    n=1
	 endif
	 i=i+1
	 last=i
      else if (line(i:i).eq.' ') then
         sp=sp+1
      else if (line(i:i).eq.'-' .or. line(i:i).eq.'+') then
      else if (line(i:i).ge.'0' .and. line (i:i).le.'9') then
      else
         call prlp(x0,y0,n)
	 write (*,'(a,$)') line(last+1:i)
         last=i
	 n=1
      endif
      goto 10
100   continue
      call prlp(x0,y0,n)
      write (*,*)
      end

      subroutine prlp(x,y,n)
      integer x,y,n,lnblnk
      character*25 line
      if (x.eq.0 .and. y.eq.0) return
      n=n+1
      write (line,*) x,y,' LP'
      if (n.eq.1) then
         write (*,'(a,$)') line(2:lnblnk(line))
      else
         write (*,'(a,$)') line(:lnblnk(line))
      endif
      x=0
      y=0
      end
