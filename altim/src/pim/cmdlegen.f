      subroutine cmdlegen1
      include "pim.inc"
      real xblc,xtrc,yblc,ytrc
      logical horleg,pimopt,l
      common /legend/ horleg

* Adjust viewport to fit in legend

      if (popcmd('LEGEN',argum)) then
         if (project/10.eq.1) then
            call pmswin((xw0+xw1)/2,(yw0+yw1)/2,0.0,0.0)
         else
            call pmswin(xw0,xw1,yw0,yw1)
         endif
	 call pgqvp(0,xblc,xtrc,yblc,ytrc)
	 horleg=(xblc.le.1-ytrc)
	 if (pimopt('-h',argum,dum,dum,dum,dum)) horleg=.true.
	 if (pimopt('-v',argum,dum,dum,dum,dum)) horleg=.false.
	 if (horleg) then
	    yv0=yv0+0.10
	 else
	    xv1=xv1-0.08
	 endif
	 call pgsvp(xv0,xv1,yv0,yv1)
      endif
      if (project/10.eq.1) then
         call pmswin((xw0+xw1)/2,(yw0+yw1)/2,0.0,0.0)
      else
         call pmswin(xw0,xw1,yw0,yw1)
      endif
      call pgqvp(0,xv0,xv1,yv0,yv1)
      end

      subroutine cmdlegen2
      include "pim.inc"
      real rint,sub,xp0,xp1,yp0,yp1,xblc,xtrc,yblc,ytrc,ver,
     |leglength,legwidth
      logical horleg,pimopt,l
      common /legend/ horleg

* Plot legend if requested

*     if (sname.eq.' ') return
      call pgsave
      if (popcmd('LEGEN',argum)) then
	 rmin=rmins
	 rmax=rmaxs	 
	 rint=0
	 sub=0
	 ver=0
	 ch=def_ch
	 lw=1
	 leglength=0.50
	 legwidth=0.02
	 name='meters'
	 if (pimopt('-h',argum,dum,dum,dum,dum)) horleg=.true.
	 if (pimopt('-v',argum,dum,dum,dum,dum)) horleg=.false.
	 l=pimopt('ver=',argum,ver,dum,dum,dum)
	 if (horleg) then
	    x=0.5
	    y=yv0-0.08
	    if (nint(ver).eq.1) y=yv0-0.04
	 else
	    x=xv1+0.10
	    y=0.5
	 endif
	 l=pimopt('rng=',argum,rmin,rmax,dum,dum)
	 l=pimopt('min=',argum,rmin,dum,dum,dum)
	 l=pimopt('max=',argum,rmax,dum,dum,dum)
	 l=pimopt('tix=',argum,rint,sub,dum,dum)
	 l=pimopt('ch=',argum,ch,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('pos=',argum,x,y,dum,dum)
	 l=pimopt('size=',argum,leglength,legwidth,dum,dum)
	 call strip(argum,name)
	 read (argum,*,iostat=ios) rmin,rmax,rint,sub
	 call pgsch(ch)
	 call pgslw(nint(lw))
	 call pgqvp(0,xp0,xp1,yp0,yp1)
	 call pgqwin(xblc,xtrc,yblc,ytrc)
	 if (horleg) then
	    write (0,'(a)') 'Plotting legend, horizontal bar ...'
	    call pgclh(rmin,rmax,x,y,name,work1,rint,nint(sub),nint(ver),
     |		leglength,legwidth)
	 else
	    write (0,'(a)') 'Plotting legend, vertical bar ...'
	    call pgclv(rmin,rmax,x,y,name,work1,rint,nint(sub),nint(ver),
     |		leglength,legwidth)
	 endif
	 call pgsvp(xp0,xp1,yp0,yp1)
	 call pgswin(xblc,xtrc,yblc,ytrc)
      endif
      call pgunsa
      end
