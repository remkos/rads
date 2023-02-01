      program elemplot

      parameter (n=999)

      real imin,imax
      real start(n),end(n),a(n),e(n),incl(n),og0(n),ogdot(n),
     .  w0(n),wdot(n),wm0(n),wmdot(n),arc(n)
      real xv0/0.03/,xv1/0.98/,yv0/0.1/,yv1/0.9/
      
      character*80 dev,line*256

      namelist /nml/ amin,amax,emin,emax,imin,imax,ogmin,ogmax,
     .wmin,wmax,dev,ch

      open (7,file='elemplot.nml')
      read (7,nml)
      close (7)

      i=0
10    continue
	 read (*,'(a)',end=20) line
	 if (line(:1).eq.'#') goto 10
	 i=i+1
         read (line,*) arc(i),start(i),end(i),step,a(i),e(i),
     |		incl(i),og0(i),ogdot(i),w0(i),wdot(i),wm0(i),wmdot(i)
	 start(i)=start(i)-48622
	 end(i)=end(i)-48622
	 a(i)=a(i)/1000
	 e(i)=e(i)*1000
	 d=start(i)-nint(start(i))
	 og0(i)=og0(i)-d*360
	 if (og0(i).lt.0) og0(i)=og0(i)+360
 	 og0(i)=(og0(i)/360*24-22)*60
	 ogdot(i)=(ogdot(i)-360)/360*24*60
      goto 10

20    nrec=i
      call pgbegin(0,dev,1,5)

      call pgsch(ch)
      call pgpage
      call pgvport(xv0,xv1,yv0,yv1)
      call pgwindow(arc(1),arc(n)+1,0.0,1.0)
      call pgbox('CMST',10.0,10,' ',0.0,0)
      call pgwindow(start(1),end(nrec),amin,amax)
      call pgbox('BNST',0.0,0,'BCNST',0.0,0)
      call plot(nrec,start,end,a)
      call label ('a (km)')

      call pgpage
      call pgvport(xv0,xv1,yv0,yv1)
      call pgwindow(arc(1),arc(nrec)+1,0.0,1.0)
      call pgbox('CMST',10.0,10,' ',0.0,0)
      call pgwindow(start(1),end(nrec),emin,emax)
      call pgbox('BNST',0.0,0,'BCNST',0.0,0)
      call plot(nrec,start,end,e)
      call label ('e (10\\u-3\\d)')

      call pgpage
      call pgvport(xv0,xv1,yv0,yv1)
      call pgwindow(arc(1),arc(nrec)+1,0.0,1.0)
      call pgbox('CMST',10.0,10,' ',0.0,0)
      call pgwindow(start(1),end(nrec),imin,imax)
      call pgbox('BNST',0.0,0,'BCNST',0.0,0)
      call plot(nrec,start,end,incl)
      call label ('i (deg)')

      call pgpage
      call pgvport(xv0,xv1,yv0,yv1)
      call pgwindow(arc(1),arc(nrec)+1,0.0,1.0)
      call pgbox('CMST',10.0,10,' ',0.0,0)
      call pgwindow(start(1),end(nrec),wmin,wmax)
      call pgbox('BNST',0.0,0,'BCNST',0.0,0)
      call plot(nrec,start,end,w0)
      call label ('\\gw (deg)')

      call pgpage
      call pgvport(xv0,xv1,yv0,yv1)
      call pgwindow(arc(1),arc(nrec)+1,0.0,1.0)
      call pgbox('CMST',10.0,10,' ',0.0,0)
      call pgwindow(start(1),end(nrec),ogmin,ogmax)
      call pgbox('BNST',0.0,0,'BCNST',0.0,0)
      call plot(nrec,start,end,og0)
      call label ('LST 22:00 + (min)')

      call pgend
      end

      subroutine plot(n,a,b,c)
      real a(n),b(n),c(n)
      do i=1,n
	 call pgmove(a(i),c(i))
	 call pgdraw(b(i),c(i))
      enddo
      end

      subroutine label(a)
      character*(*) a
      call pgmtext ('B',2.4,1.0,1.0,'days since 2.0 Jan 1992')
      call pgmtext ('T',2.0,0.0,0.0,'arcnr')
      call pgmtext ('L',2.2,0.5,0.5,a)
      end
