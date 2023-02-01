      subroutine cmdlogo(cmdname,atend)
      character cmdname*(*)
      include "pim.inc"
      real xp0,xp1,yp0,yp1
      integer ciconv(0:255),i
      logical l,pimopt,atend
     
* Plot logo. Get PGM file.

      call pgsave
      if (popcmd(cmdname,argum)) then
	 write (0,'(a)') 'Plotting logo ...'
	 if (postscript) then
	 write (0,'(a)') '... logo for postscript not implemented'
         else if (pimopt('atend',argum,x,dum,dum,dum).eqv.atend) then
         x=0.01
         y=0.01
	 name='logo_col.pgm'
	 l=pimopt('x=',argum,x,dum,dum,dum)
	 l=pimopt('y=',argum,y,dum,dum,dum)
	 l=pimopt('pos=',argum,x,y,dum,dum)
	 call strip(argum,name)
	 read (argum,*,iostat=ios) x,y
	 do i=0,255
	    ciconv(i)=c_bg
	 enddo
	 do i=0,15
	    ciconv(i*16)=i
	 enddo
	 call readpgm(name,work1,work2,ix,iy,ciconv)
	 call pgsvp (0.,1.,0.,1.)
	 call pgqvp(3,xp0,xp1,yp0,yp1)
	 x=xp0+(xp1-xp0)*x
	 y=yp0+(yp1-yp0)*y
	 call grpxpx(work2,ix,iy,1,ix,1,iy,x,y)
	 call pgsvp (xv0,xv1,yv0,yv1)
	 endif
      endif
      call pgunsa
      end
