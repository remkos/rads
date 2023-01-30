      subroutine cmddynto
      include "pim.inc"
      real dx,dy,xg0,xg1,yg0,yg1,vecfact,vecch
      common /cvec/ vecfact,vecch

* Call DYNTOPO for plotting velocity vectors from dynamic topography

      call pgsave
340   if (pop1cmd('DYNTO',argum)) then
	 name=sname
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
         call pimopt('ls=',argum,ls,dum,dum,dum)
         call pimopt('lw=',argum,lw,dum,dum,dum)
         call pimopt('ch=',argum,ch,dum,dum,dum)
         call pimopt('fact=',argum,fact,dum,dum,dum)
         call pimopt('rng=',argum,rmin,rmax,dum,dum)
         call pimopt('step=',argum,dx,dy,dum,dum)
	 vecch=ch
	 vecfact=fact
         call strip(argum,name)
	 if (name.eq.' ' .or. name.eq.'=') name=sname
	 write (0,'(a)') 'Reading dynamic topography ...'
	 if (popnextcmd('GRIDA',argum))
     |		read (argum,*,iostat=ios) xg0,xg1,yg0,yg1
	 call gridread(work1,work2,work4,ix,iy,name,
     |		xg0,xg1,yg0,yg1)
	 write (0,'(a)') 'Computing current velocity ...'
	 call dyntopo(work1,work2,work3,ix,iy,ix,xg0,xg1,yg0,yg1)
         call pgsci(nint(ci))
         call pgsls(nint(ls))
         call pgslw(nint(lw))
         call pgsch(ch)
	 write (0,'(a)') 'Plotting velocity vectors ...'
	 call velocity(work2,work3,ix,iy,
     |		nint(dx),nint(dy),xg0,xg1,yg0,yg1,rmin,rmax,fact)
         goto 340
      endif
      call pgunsa
      end
