      parameter (npnt=100000)
      real x(npnt),y(npnt)
      integer rec(4)
      character*80 filenm

      xmin=1e40
      xmax=-1e40
      ymin=1e40
      ymax=-1e40
      day=86400

      call getarg(1,filenm)
      open (10,file=filenm,status='old',form='unformatted',
     .access='direct',recl=18)
      read (10,rec=1) spec,n
      do i=1,n
	 read (10,rec=i+1) rec
	 x(i)=rec(1)/day
         y(i)=rec(4)/1e4
	 xmin=min(x(i),xmin)
	 xmax=max(x(i),xmax)
	 ymin=min(y(i),ymin)
	 ymax=max(y(i),ymax)
      enddo
      dx=(xmax-xmin)/3

      call getarg(2,filenm)
      call pgbegin(0,filenm,1,3)
      call pgpage
      call pgwindow(xmin,xmin+dx,ymin,ymax)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
      call pgline(n,x,y)
      call pgpage
      call pgwindow(xmin+dx,xmin+2*dx,ymin,ymax)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
      call pgline(n,x,y)
      call pgpage
      call pgwindow(xmin+2*dx,xmax,ymin,ymax)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
      call pgline(n,x,y)
      call pgend
      end
