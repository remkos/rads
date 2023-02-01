      subroutine cmdveloc
      include "pim.inc"
      real dx,dy,xg0,xg1,yg0,yg1,vecfact,vecch
      common /cvec/ vecfact,vecch

* Call DYNTOPO for plotting velocity vectors from dynamic topography

      call pgsave
340   if (pop1cmd('VELOC',argum)) then
         ci=1
	 ls=1
	 lw=1
         fact=1e-2
         ch=0.2
         rmin=0
         rmax=10
         dx=1
         dy=1
	 ix=0
	 iy=maxgrd
	 xg0=0
	 xg1=0
	 yg0=0
	 yg1=0
         call pimopt('ci=',argum,ci,dum,dum,dum)
         call pimopt('ls=',argum,lw,dum,dum,dum)
         call pimopt('lw=',argum,lw,dum,dum,dum)
         call pimopt('ch=',argum,ch,dum,dum,dum)
         call pimopt('fact=',argum,fact,dum,dum,dum)
         call pimopt('rng=',argum,rmin,rmax,dum,dum)
         call pimopt('step=',argum,dx,dy,dum,dum)
	 vecch=ch
	 vecfact=fact
         call strip(argum,name)
	 write (0,'(a)') 'Reading X-velocity grid ...'
	 if (popnextcmd('GRIDA',sargum))
     |		read (argum,*,iostat=ios) xg0,xg1,yg0,yg1
	 call gridread(work1,work3,work4,ix,iy,name,
     |		xg0,xg1,yg0,yg1)
         call strip(argum,name)
	 write (0,'(a)') 'Reading Y-velocity grid ...'
	 call gridread(work2,work3,work4,ix,iy,name,
     |		xg0,xg1,yg0,yg1)
         call pgsci(nint(ci))
         call pgsls(nint(ls))
         call pgslw(nint(lw))
         call pgsch(ch)
	 write (0,'(a)') 'Plotting velocity vectors ...'
	 call velocity(work1,work2,ix,iy,
     |		nint(dx),nint(dy),xg0,xg1,yg0,yg1,rmin,rmax,fact)
         goto 340
      endif
      call pgunsa
      end
