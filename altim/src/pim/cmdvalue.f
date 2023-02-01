      subroutine cmdvalue
      include "pim.inc"
      logical l,pimopt
      real just,angle

* Draw surface grid mesh

      if (sname.eq.' ') return
      call pgsave
      if (popcmd('VALUE',argum)) then
	 write (0,550) 'Plotting mesh ...'
	 ch=def_ch
         lw=1
	 ls=1
         ci=c_cont
	 angle=0
	 just=0.5
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('ch=',argum,ch,dum,dum,dum)
	 l=pimopt('ci=',argum,ci,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('just=',argum,just,dum,dum,dum)
	 l=pimopt('angle=',argum,angle,dum,dum,dum)
	 call pgsch(ch)
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgsci(nint(ci))
	 call values(xs0,xs1,ys0,ys1,nsx,nsy,work1,angle,just)
      endif
  550 format (a)
      call pgunsa
      end

      subroutine values(x0,x1,y0,y1,nx,ny,grid,angle,just)
      real x0,x1,y0,y1,dx,dy,x,y,z,grid(nx,ny),angle,just
      integer nx,ny,kx,ky,mm,pp,nc
      character*80 text

      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)
      do ky=1,ny
	 do kx=1,nx
	    z=grid(kx,ky)
	    if (z.lt.1e20) then
	       x=x0+(kx-1)*dx
	       y=y0+(ky-1)*dy
	       pp=int(log10(abs(z)))-3
	       mm=nint(z/10.**pp)
	       call pmconv(1,x,y)
	       call pgnumb(mm,pp,0,text,nc)
	       call pgptxt(x,y,angle,just,text(:nc))
	    endif
	 enddo
      enddo
      end
