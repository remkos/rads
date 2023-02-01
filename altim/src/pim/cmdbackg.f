      subroutine cmdbackg
      include "pim.inc"
      real pixx0,pixx1,pixy0,pixy1
      real*4    r,g,b,hue,lgt,sat
      real halftone/0.85/
      integer i,j,ciconv(0:255)

* Handle background

      if (popcmd('BACKG',argum)) then
	 write (0,'(a)') 'Setting up background ...'
         if (postscript) then
	 write (0,'(a)') '... background for postscript not implemented'
	 else
	 name='logo_bw.pgm'
	 call pgsvp (0.,1.,0.,1.)
	 call pgqvp (3,pixx0,pixx1,pixy0,pixy1)
	 call pimopt('tone=',argum,halftone,dum,dum,dum)
	 call strip(argum,name)
	 read (argum,*,iostat=ios) halftone
	 i=c_bg
	 call grxhls(rgb(1,i),rgb(2,i),rgb(3,i),hue,lgt,sat)
	 if (lgt.gt.0.5) then
	    lgt=lgt*halftone
	 else if (lgt/halftone.lt.1-halftone) then
	    lgt=1-halftone
	 else
	    lgt=lgt/halftone
	 endif
	 call grxrgb(hue,lgt,sat,r,g,b)
	 call pgscr(c_half,r,g,b)
	 do i=0,255
	    ciconv(i)=c_bg
	 enddo
	 do i=1,15
	    ciconv(i*16)=c_half
	 enddo
	 call readpgm(name,work1,work2,ix,iy,ciconv)
	 do j=1,(pixy1-pixy0+1)/iy
	    do i=1,(pixx1-pixx0+1)/ix
	       call grpxpx(work2,ix,iy,1,ix,1,iy,(i-1)*real(ix),
     .                  (j-1)*real(iy))
	    enddo
	 enddo
	 call pgsvp (xv0,xv1,yv0,yv1)
	 endif
      endif
      end
