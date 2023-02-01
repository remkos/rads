      subroutine cmdmesh
      include "pim.inc"
      logical l,pimopt
      real mark

* Draw surface grid mesh

      if (sname.eq.' ') return
      call pgsave
      if (popcmd('MESH',argum)) then
	 write (0,550) 'Plotting mesh ...'
	 ch=def_ch
         lw=1
	 ls=1
         ci=c_cont
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('ch=',argum,ch,dum,dum,dum)
	 l=pimopt('ci=',argum,ci,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('mark=',argum,mark,dum,dum,dum)
	 call pgsch(ch)
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgsci(nint(ci))
	 if (pimopt('-points',argum,dum,dum,dum,dum)) then
	    call mesh1(xs0,xs1,ys0,ys1,nsx,nsy,nint(mark))
	 else if (pimopt('-crosses',argum,dum,dum,dum,dum)) then
	    call mesh1(xs0,xs1,ys0,ys1,nsx,nsy,1)
	 else if (pimopt('-x',argum,dum,dum,dum,dum)) then
	    call mesh1(xs0,xs1,ys0,ys1,nsx,nsy,1)
	 else if (pimopt('-grid',argum,dum,dum,dum,dum)) then
	    call mesh2(xs0,xs1,ys0,ys1,nsx,nsy,cells,project)
	 else
	    call mesh1(xs0,xs1,ys0,ys1,nsx,nsy,-1)
	 endif
      endif
  550 format (a)
      call pgunsa
      end

      subroutine mesh1(x0,x1,y0,y1,nx,ny,mark)
      real x0,x1,y0,y1,dx,dy,x,y
      integer nx,ny,mark,kx,ky

      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)
      do ky=1,ny
	 do kx=1,nx
	    x=x0+(kx-1)*dx
	    y=y0+(ky-1)*dy
	    call pmconv(1,x,y)
	    call pgpoint(1,x,y,mark)
	 enddo
      enddo
      end

      subroutine mesh2(x0,x1,y0,y1,nx,ny,cell,project)
      real x0,x1,y0,y1,dx,dy,x,y,x2,y2,x3,y3
      integer nx,ny,kx,ky,mx,my,project
      logical cell

      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)
      if (cell) then
	 x2=x0-dx/2
	 y2=y0-dy/2
	 x3=x1+dx/2
	 y3=y1+dy/2
	 mx=nx+1
	 my=ny+1
      else
	 x2=x0
	 y2=y0
	 x3=x1
	 y3=y1
	 mx=nx
	 my=ny
      endif

      if (project.lt.10) then
      do kx=1,mx
	 x=x2+(kx-1)*dx
	 y=y2
	 call pmconv(1,x,y)
	 call pgmove(x,y)
	 x=x2+(kx-1)*dx
	 y=y3
	 call pmconv(1,x,y)
	 call pgdraw(x,y)
      enddo
      do ky=1,my
	 x=x2
	 y=y2+(ky-1)*dy
	 call pmconv(1,x,y)
	 call pgmove(x,y)
	 x=x3
	 y=y2+(ky-1)*dy
	 call pmconv(1,x,y)
	 call pgdraw(x,y)
      enddo
      else
      do kx=1,mx
	 do ky=1,my-1
	    x=x2+(kx-1)*dx
	    y=y2+(ky-1)*dy
	    call pmconv(1,x,y)
	    call pgmove(x,y)
	    x=x2+(kx-1)*dx
	    y=y2+ky*dy
	    call pmconv(1,x,y)
	    call pgdraw(x,y)
	 enddo
      enddo
      do ky=1,my
	 do kx=1,mx-1
	    x=x2+(kx-1)*dx
	    y=y2+(ky-1)*dy
	    call pmconv(1,x,y)
	    call pgmove(x,y)
	    x=x2+kx*dx
	    y=y2+(ky-1)*dy
	    call pmconv(1,x,y)
	    call pgdraw(x,y)
	 enddo
      enddo
      endif
      end
