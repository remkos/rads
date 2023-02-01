      subroutine cmdcircl
      include "pim.inc"
      real s,c,ci1,ci2
      logical l,xy,pimopt

* Plot symbols at requested position

      call pgsave

360   if (pop1cmd('CIRCL',argum)) then
	 ci1=c_fg
	 ci2=-1
	 ls=1
	 lw=1
	 fs=2
	 ch=def_ch
	 s=1000
	 c=-1
	 xy=.true.
	 if (pimopt('-xy',argum,dum,dum,dum,dum)) then
	    xy=.true.
	 endif
	 if (pimopt('-yx',argum,dum,dum,dum,dum)) then
	    xy=.false.
	 endif
	 l=pimopt('ci=',argum,ci,dum,dum,dum)
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('fs=',argum,fs,dum,dum,dum)
	 l=pimopt('r=',argum,s,c,dum,dum)
	 call pgsch(ch)
	 call pgsci(nint(ci))
	 call pgsfs(nint(fs))
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 write (0,'(a)') 'Plotting range circles ...'
	 if (argum(1:5).eq.'file=') then
	    open (20,file=argum(6:),status='old')
400	    read (20,*,end=410) x,y
	    call circle(x,y,s,c,xy)
	    goto 400
410	    close (20)
	 else
	    read (argum,*,iostat=ios) x,y
	    call circle(x,y,s,c,xy)
	 endif
	 goto 360
      endif
      call pgunsa
      end

      subroutine circle(lon,lat,s,c,xy)
      logical xy
      real s,c,r,rad,re,lat,lon
      include "pim.inc"

      rad=atan(1.)/45
      re=6371.0

* If cutoff-elevation c > 0, then convert satellite altitude into range

      if (c.ge.0) then
         r = re * ( acos ( re/(re+s) * cos(c*rad) ) - c*rad )
      else
         r = s
      endif

      if (xy) then
         call rngcir(lon,lat,r,360,1.0,360.0,work1,work2)
      else
         call rngcir(lat,lon,r,360,1.0,360.0,work1,work2)
      endif

      call pmconv(360,work1,work2)
      call pgpoly(360,work1,work2)
      end
