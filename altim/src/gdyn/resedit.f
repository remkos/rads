* resedit: Subroutine belonging to the program lres.
*
* This routine uses the PGPLOT interface to display and edit
* one pass of SLR (or other range residual) data.
*
* Some plotting defaults are set through the namelist file
* 'lres.nml'. Note that additional interfacing is done through
* the common block cedit

      subroutine resedit(nobs,resid,sigma,deriv,elev,line,name)
      integer nobs
      real*8 resid(*),sigma(*),deriv(*),elev(*)
      character line*(*),name*32

      integer i,ci,l,lnblnk,j,n,pgcurs
      real x,y,xcurs,ycurs,fact
      character text*80,ch*4

      logical trend

      real*8 rmean,rrms,erms,fit,a,b

      include "lres.inc"

      save xcurs,ycurs,fact
      data trend/.false./,fact/1.0/

      j=(nobs+1)/2

300   call resfit(nobs,deriv,resid,sigma,n,rmean,rrms,erms,a,b,fit)

400   call pgpage
      call pgswin(xmin,xmax,0.0,90.0)
      call pgbox(' ',0.0,0,'cmst',30.0,3)
      call pgslw(5)
      do i=1,nobs
	 x=deriv(i)/1000
	 y=elev(i)
	 call pgpt(1,x,y,1)
      enddo
      call pgslw(1)
      call pgswin(xmin,xmax,ymin*fact,ymax*fact)
      call pgbox('abcnst',0.0,0,'abnst',0.0,0)
      l=lnblnk(line)
      call pglab('range rate (km/s)','range residual (cm)',' ')
      call pgmtxt('T',1.5,0.5,0.5,line(:40)//'  '//name)
      call pgmtxt('R',2.5,0.5,0.5,'elevation (deg)')
      call pgslw(20)
      do i=1,nobs
         x=deriv(i)/1000
         if (trend) then
            y=(resid(i)-a-b*deriv(i))*100
         else
            y=resid(i)*100
         endif
         ci=3
	 if (sigma(i).eq.0) ci=5
	 if (sigma(i).lt.0) ci=2
         if (y.gt.ymax*fact) then
            y=0.999*ymax*fact
            if (ci.eq.3) ci=4
         else if (y.lt.ymin*fact) then
            y=0.999*ymin*fact
            if (ci.eq.3) ci=4
         endif
	 call pgsci(ci)
         call pgpt(1,x,y,-1)
      enddo
      call pgsci(1)
      call pgslw(1)
      write (text,402) 'Obs = ',n,nobs-n
      call pgmtxt('T',-3.0,0.1,0.0,text)
      write (text,403) 'RMS = ',rrms*100,erms*100,' cm'
      call pgmtxt('T',-4.5,0.1,0.0,text)
      write (text,401) 'Range bias = ',a*100,' cm'
      call pgmtxt('T',-6.0,0.1,0.0,text)
      write (text,401) 'Time bias  = ',b*1d6,' \\gms'
      call pgmtxt('T',-7.5,0.1,0.0,text)
      write (text,401) 'RMS of fit = ',fit*100,' cm'
      call pgmtxt('T',-9.0,0.1,0.0,text)
      write (text,402) 'Pass',ii
      call pgsci(2)
      call pgmtxt('T',-1.5,0.1,0.0,text)
      call pgsci(5)
401	 format (a,f6.1,a)
402	 format (a,2i6,a)
403	 format (a,2f6.1,a)
      if (.not.trend) then
         y=a*1e2+b*1e5*xmin
         call pgmove(xmin,y)
	 y=a*1e2+b*1e5*xmax
         call pgdraw(xmax,y)
      else
	 call pgmtxt('T',-10.5,0.1,0.0,'DETRENDED')
      endif
      if (line(l:l).eq.'!') then
	 call pgsci(2)
	 call pgmtxt('T',-12.0,0.1,0.0,'WARNING !')
      endif
      call pgsci(1)

420   continue
      xcurs=deriv(j)/1000
      if (trend) then
         ycurs=(resid(j)-a-b*deriv(j))*100
      else
         ycurs=resid(j)*100
      endif
      ycurs=min(ycurs,ymax*fact)
      ycurs=max(ycurs,ymin*fact)
      curs=pgcurs(xcurs,ycurs,ch)
      i=ichar(ch(1:1))

      if (curs.eq.0) then		! No cursor device
         di=+1
      else if (i.eq.111. .or. i.eq.57) then	! o / Page Up  (numpad)
         di=-10
      else if (i.eq.44 .or. i.eq.51) then	! , / Page Down
         di=+10
      else if (i.eq.105 .or. i.eq.56) then	! i / Up
	 di=-1
      else if (i.eq.109 .or. i.eq.50 .or. i.eq.88) then	! m / Down / Right Mouse
	 di=+1
      else if (i.eq.117 .or. i.eq.55) then	! u / Home
         di=-1000000
      else if (i.eq.110 .or. i.eq.49) then	! n / End
         di=1000000
      else if (i.eq.32 .or. i.eq.68) then	! Space / Middle Mouse
	 call find(nobs,xcurs,ycurs,
     |		xmin,xmax,ymin,ymax,deriv,resid,trend,a,b,j)
	 sigma(j)=-sigma(j)
	 goto 300
      else if (i.eq.27) then		! Esc
	 di=0
      else if (i.eq.113) then		! q
	 di=0
      else if (i.eq.120) then		! x
	 di=0
	 returncode=999
      else if (i.eq.43) then		! +
	 fact=fact*2
	 goto 400
      else if (i.eq.45) then		! -
	 fact=fact/2
	 goto 400
      else if (i.eq.65) then		! Left mouse
	 call find(nobs,xcurs,ycurs,
     |		xmin,xmax,ymin,ymax,deriv,resid,trend,a,b,j)
	 goto 420
      else if (i.eq.106 .or. i.eq.52) then		! j / Left
	 if (j.gt.1) j=j-1
	 goto 420
      else if (i.eq.108 .or. i.eq.54) then		! l / Right
	 if (j.lt.nobs) j=j+1
	 goto 420
      else if (i.eq.47) then		! Slash
	 trend=.not.trend
	 goto 400
      else if (i.eq.98 .or. i.eq.48) then		! b / Insert (Numpad)
	 do i=1,nobs
	    sigma(i)=+abs(sigma(i))
	 enddo
	 goto 300
      else if (i.eq.46) then		! . / Del (Numpad)
	 do i=1,nobs
	    sigma(i)=-abs(sigma(i))
	 enddo
	 goto 300
      else if (i.eq.119) then		! w
         select=.true.
	 if (line(l:l).eq.'!') goto 420
      else if (i.eq.97) then		! a
         select=.false.
	 goto 420
      else
         write (6,551) char(7)
         goto 420
      endif

551   format(a,$)
      end

      subroutine find(n,x0,y0,
     |		xmin,xmax,ymin,ymax,x,y,trend,a,b,j)
      integer n,i,j
      real x0,y0,x1,y1,xmin,xmax,ymin,ymax,z,zmin
      real*8 x(*),y(*),a,b
      logical trend

      zmin=1e30
      do i=1,n
	 x1=x(i)/1000
	 if (trend) then
	    y1=(y(i)-a-b*x(i))*100
	 else
	    y1=y(i)*100
	 endif
	 if (y1.lt.ymin) y1=ymin
	 if (y1.gt.ymax) y1=ymax
	 z=((x1-x0)/(xmax-xmin))**2+((y1-y0)/(ymax-ymin))**2
	 if (z.lt.zmin) then
	    j=i
	    zmin=z
	 endif
      enddo
      end
