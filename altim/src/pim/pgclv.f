**PGCLV -- vertical color legend bar
*
      subroutine pgclv(leg0,leg1,xv,yv,units,ileg,rint,nsubs,ver,
     |leglength,legwidth)

      include 'pim.inc'
      integer ic2,nsubs,ver
      real    leg0,leg1,xv,yv,rint
      real    leglength,legwidth
      integer ileg(0:maxcol,0:maxcol),i,k
      character*(*) units
      real xp0,xp1,yp0,yp1,vecfact,vecch
      common /cvec/ vecfact,vecch

      leglength=.50
      legwidth=.02

      if (nc2.eq.0) then
	 ic2=0
         do k=0,nc1
	    ileg(0,k)=c_0+k
	 enddo
      else
	 ic2=min(31,nc2)
	 do k=0,nc1
	    do i=0,ic2
	       ileg(ic2-i,k)=sh_mat(c_0+k,nint(i/real(ic2)*nc2))
	    enddo
	 enddo
      endif
      call pgsvp(xv-legwidth,xv,
     .		   yv-leglength/2,yv+leglength/2)
      call pgswin(0.,real(ic2+1),0.,real(nc1+1))
      call pgpixl(ileg(0,0),maxcol+1,maxcol+1,1,ic2+1,1,nc1+1,
     &		    0.,real(ic2+1),0.,real(nc1+1))
      call pgswin(0.0,1.0,leg0,leg1)
      if (ver.eq.1) then
         call pgbox('BC',0.0,0,'C',rint,nsubs)
         call pgbox(' ',0.0,0,'BNSTIV',rint,nsubs)
      else
         call pgbox('BC',0.0,0,'BCNSTIV',rint,nsubs)
      endif

* Now draw unit text at top of legenda

      if (ver.eq.1) then
         call pgmtext('T',1.2,1.0,1.0,units)
      else
         call pgmtext('T',1.2,1.0,1.0,units)
      endif

* If requested, plot vector below the legend in separate box

      if (vecfact.le.0) return
      call pgqvp(2,xp0,xp1,yp0,yp1)
      if (vecfact.le.1e-2) then
         call pgsvpx(2,xp1-vecfact*1e3,xp1,yp0-2*(xp1-xp0),yp0)
         call pgswin(0.,1.,-1.,3.)
         call pgmtxt('B',1.2,0.5,0.5,'1 m/s')
      else
         call pgsvpx(2,xp1-vecfact*1e2,xp1,yp0-2*(xp1-xp0),yp0)
         call pgswin(0.,1.,-1.,3.)
         call pgmtxt('B',1.2,0.5,0.5,'10 cm/s')
      endif
      call pgsch(vecch)
      call pgsah(2,60.0,1.0)
      call pgarro(0.,0.,1.,0.)

      end
