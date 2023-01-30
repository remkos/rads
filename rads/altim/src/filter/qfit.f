**QFIT -- Fit a linear and seasonal trend
*+
      program qfit

* This program fits a linear trend and seasonal cycle through a set of
* data points.
*-
* $Log: qfit.f,v $
* Revision 1.4  2005/12/12 15:26:29  rads
* - Added option offset=
*
* Revision 1.3  2004/12/09 17:31:08  remko
* - Added f7.2 format
*
* Revision 1.2  2004/09/24 14:30:29  remko
* - Cosmetic changes only
*
* Revision 1.1  2004/09/21 21:05:17  remko
* - Filter routines added to altim/src suite
*
*-----------------------------------------------------------------------
      implicit none
      integer mcol,mrec
      parameter (mcol=200,mrec=10000)
      integer	i,j,k,m,xcol/1/,ycol/2/,n/0/,unit/0/,iargc
      real*8	scale/1d0/,offset/0d0/,twopi,ata(10)/10*0d0/,a(4),
     |		b(5)/5*0d0/,array(mcol),x0/-1d30/,x1/1d30/,y0/-1d30/,y1/1d30/,
     |		t(mrec),y(mrec),year/31557600d0/
      character*1024 line
      logical	csh/.false./

      do i=1,iargc()
         call getarg(i,line)
	 if (line(:5).eq.'xcol=') then
	    read (line(6:),*) xcol
	 else if (line(:5).eq.'ycol=') then
	    read (line(6:),*) ycol
	 else if (line(:2).eq.'x=') then
	    read (line(3:),*,iostat=k) x0,x1
	 else if (line(:2).eq.'y=') then
	    read (line(3:),*,iostat=k) y0,y1
	 else if (line(:6).eq.'scale=') then
	    read (line(7:),*) scale
	 else if (line(:7).eq.'offset=') then
	    read (line(8:),*) offset
	 else if (line(:2).eq.'-s') then
	    csh=.false.
	 else if (line(:2).eq.'-c') then
	    csh=.true.
	 else if (line.eq.'-') then
	    unit=6
	 else
	    unit=10
	    open (unit,file=line)
         endif
      enddo

      m=max(xcol,ycol)
      if (m.gt.mcol) call fin("too many columns")
      twopi=8*atan(1d0)

10    read (*,'(a)',end=100) line
      if (line(:1).eq.'#') goto 10
      read (line,*) (array(i),i=1,m)
      n=n+1
      if (n.gt.mrec) call fin("too many records")
      t(n)=array(xcol)/year-15d0
      y(n)=array(ycol)*scale+offset
      if (array(xcol).lt.x0 .or. array(xcol).gt.x1) goto 10
      if (array(ycol).lt.y0 .or. array(ycol).gt.y1) goto 10
      a(1)=1d0
      a(2)=t(n)
      a(3)=cos(t(n)*twopi)
      a(4)=sin(t(n)*twopi)

      k=0
      do i=1,4
         do j=1,i
	    k=k+1
	    ata(k)=ata(k)+a(i)*a(j)
	 enddo
	 b(i)=b(i)+a(i)*y(n)
      enddo

      goto 10

100   continue
      call dppsv('U',4,1,ata,b,4,i)
      b(5)=sqrt(b(3)**2+b(4)**2)

      do i=1,5
	 if (i.gt.1) write (*,'("; ",$)')
	 if (csh) write (*,'("set ",$)')
	 write (*,'(a1,"=",$)') char(96+i)
         if (b(i).ge.999.995d0 .or. b(i).le.-99.995d0) then
	    write (*,'(f7.2,$)') b(i)
         else if (b(i).ge.99.995d0 .or. b(i).le.-9.995d0) then
	    write (*,'(f6.2,$)') b(i)
         else if (b(i).ge.9.995d0 .or. b(i).le.-0.005d0) then
	    write (*,'(f5.2,$)') b(i)
         else
	    write (*,'(f4.2,$)') b(i)
	 endif
      enddo
      write (*,*)
      if (unit.gt.0) then
         do i=1,n
	    write (unit,*) (t(i)+15d0)*year,y(i),b(1)+b(2)*t(i),
     |		b(3)*cos(t(i)*twopi)+b(4)*sin(t(i)*twopi)
         enddo
	 close (unit)
      endif
      end
