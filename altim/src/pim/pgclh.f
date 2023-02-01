**PGCLH -- horizontal color legend bar
*
      subroutine pgclh(leg0,leg1,xv,yv,units,ileg,rint,nsubs,ver,
     |leglength,legwidth)

      include 'pim.inc'
      integer ic2,nsubs,ver
      real    leg0,leg1,xv,yv,rint
      real    leglength,legwidth
      integer ileg(0:maxcol,0:maxcol),i,k
      character*(*) units
      real xp0,xp1,yp0,yp1,vecfact,vecch
      common /cvec/ vecfact,vecch

      if (nc2.eq.0) then
	 ic2=0
         do k=0,nc1
	    ileg(k,0)=c_0+k
	 enddo
      else
	 ic2=min(31,nc2)
	 do k=0,nc1
	    do i=0,ic2
	       ileg(k,i)=sh_mat(c_0+k,nint(i/real(ic2)*nc2))
	    enddo
	 enddo
      endif
      call pgsvp(xv-leglength/2,xv+leglength/2,
     &             yv-legwidth,yv)
      call pgswin(0.,real(nc1+1),0.,real(ic2+1))
      call pgpixl(ileg(0,0),maxcol+1,maxcol+1,1,nc1+1,1,ic2+1,
     &		    0.,real(nc1+1),0.,real(ic2+1))
      call pgswin(leg0,leg1,0.0,1.0)
      if (ver.eq.1) then
         call pgbox('C',rint,nsubs,'BC',0.0,0)
         call pgbox('BNSIT',rint,nsubs,' ',0.0,0)
      else
         call pgbox('BCMSIT',rint,nsubs,'BC',0.0,0)
      endif

* Now draw unit text at top of legend

      if (ver.eq.1) then
         call pgmtext('RV',1.2,0.5,0.0,units)
      else
         call pgmtext('B',1.2,0.5,0.5,units)
      endif

* If requested, plot vector next to the legend in separate box

      if (vecfact.le.0) return
      call pgqvp(2,xp0,xp1,yp0,yp1)
      if (vecfact.le.1e-2) then
         call pgsvpx(2,xp1,xp1+2*vecfact*1e3,yp0,yp1)
         call pgswin(-1.,1.,-1.,1.)
         call pgmtxt('B',1.2,0.75,0.5,'1 m/s')
      else
         call pgsvpx(2,xp1,xp1+2*vecfact*1e2,yp0,yp1)
         call pgswin(-1.,1.,-1.,1.)
         call pgmtxt('B',1.2,0.5,0.0,'10 cm/s')
      endif
      call pgsch(vecch)
      call pgsah(2,60.0,1.0)
      call pgarro(0.,0.,1.,0.)

      end
