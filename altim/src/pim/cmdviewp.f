      subroutine cmdviewp
      include "pim.inc"

* Set viewport

      call pgsvp(0.,1.,0.,1.)
      call pgqvp(3,xv0,xv1,yv0,yv1)
      if (debug) write (*,*) 'fullsize:',xv0,xv1,yv0,yv1
      if (popcmd('VIEWP',argum)) then
	 call grtoup(argum,argum)
	 if (argum.eq.'SMALL') then
	    call pgsvp(0.10,0.90,0.10,0.90)
	 else if (argum.eq.'BIG') then
	    call pgsvp(0.05,0.95,0.05,0.95)
	 else if (argum.eq.'FULLSCREEN') then
	    call pgsvp(0.00,1.00,0.00,1.00)
	 else
	    read (argum,*,iostat=ios) xv0,xv1,yv0,yv1
	    call pgsvp(xv0,xv1,yv0,yv1)
	 endif
      else
	 call pgsvp(0.05,0.95,0.05,0.95)
      endif
      call pgqvp(0,xv0,xv1,yv0,yv1)
      end
