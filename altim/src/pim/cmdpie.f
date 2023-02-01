      subroutine cmdpie
      include "pim.inc"   
      real a,a0,a1,da,rad
      integer i,n

* Plot a pie chunk

      rad=atan(1.)/45
      call pgsave
  320 if (pop1cmd('PIE',argum)) then
	 ci=c_fg
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 read (argum,*,iostat=ios) xv0,xv1,yv0,yv1,a0,a1
	 call pgsvp(xv0,xv1,yv0,yv1)
	 call pgwnad(-1.,1.,-1.,1.)
	 call pgsci(nint(ci))
	 da=3
	 n=int((a1-a0)/da+1)
	 da=(a1-a0)/(n-1)
	 do i=1,n
	    a=a0+da*(i-1)
	    work1(i)=sin(a*rad)
	    work2(i)=cos(a*rad)
	 enddo
	 work1(n+1)=0
	 work2(n+1)=0
	 call pgpoly(n+1,work1,work2)
	 write (0,'(a)') 'Drawing pie chunk ...'
	 goto 320
      endif
      call pgunsa
      end
