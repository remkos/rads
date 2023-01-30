**PMSWIN -- set geographical projection window and adjust viewport
*+
      SUBROUTINE PMSWIN (A, B, C, D)
      REAL A, B, C, D
*
* Change the window in real-world coordinate space that is to be mapped on
* to the viewport, and simultaneously adjust the viewport so that the
* projection type and scale selected in PMDEF are satisfied. If the
* scale was set to 0.0 in PMDEF, the new viewport is the largest one
* that can fit within the previously set viewport while retaining the
* required aspect ratio, and the scale will be computed.
* The window's coordinates will be map coordinates (x,y). However,
* the input is in real-world (longitude,latitude) coordinates (degrees).
*
* Arguments:
* For all projections others than Azimuthal:
*  A (input) : The most western longitude of the area to be mapped (degrees).
*  B (input) : The most eastern longitude of the area to be mapped (degrees).
*  C (input) : The most southern latitude of the area to be mapped (degrees).
*  D (input) : The most northern latitude of the area to be mapped (degrees).
*
* For Azimuthal projections:
*  A (input) : The longitude of the central point in the map (degrees).
*  B (input) : The latitude of the central point in the map (degrees).
*  C (input) : (for Perspective projection) The height of the observer
*              above the earth's surface (kilometers).
*              (for Azimuthal Equi-distant projection) The maximum range
*              to be plotted (kilometers). If C=0.0 the whole world will be
*              plotted.
*  D         : (not used).
*--
*  9-Jan-1991 - created [Remko Scharroo].
* 14-Jan-1991 - include conic projections.
* 16-Jan-1991 - support azimuthal projections.
*  4-Jul-1991 - bug fixed.
* 30-Oct-1991 - support more projections.
* 28-Mar-1992 - Standardize PMPLOT.
*  2-Apr-1993 - Include tilted rectangular projection. LATSP -> PPARA
* 30-Jun-1994 - Include polar projections.
* 19-Jul-1996 - Renamed to PMSWIN from PMWINDOW.
* 14-Aug-1997 - Add ERS projection.
*-----------------------------------------------------------------------
      INCLUDE 'pmplot.inc'
      character alt*10    
      REAL X(8),Y(8),KMPERX
      real xv0,xv1,xvs,xvm,xw0,xw1
      real yv0,yv1,yvs,yvm,yw0,yw1
      real cor1,cor2,temp,scale,rc
      integer nc
*
      if (PMOPEN.LT.2) then
	 CALL GRWARN('PMSWIN: Use PMDEF first')
	 RETURN
      endif
      xcurve=.false.
      ycurve=.false.
      azmtal=.false.
*
* Determine plot scale
*
      kmperx=rmean*rad
      rview=rmean+C
      cor2=0
      if (ptype.eq.22) then
         if (ppara1.lt.-90..or.ppara1.gt.+90.) ppara1=(3*C+D)/4
         if (ppara2.lt.-90..or.ppara2.gt.+90.) ppara2=(C+3*D)/4
         cor1=cos(ppara1*rad)
         cor2=cos(ppara2*rad)
      else if (ptype.eq.33) then
         ppara1=45*rad
         cor1=cos(ppara1*rad)
      else if (ptype.eq.34) then
         if (ppara1.lt.0..or.ppara1.gt.1.) ppara1=1.
         if (ppara2.lt.-90..or.ppara2.gt.+90.) ppara2=0.
	 cor1=cos((C+D)*rad/2)
      else
         if (ppara1.lt.-90..or.ppara1.gt.+90.) ppara1=(C+D)/2
         cor1=cos(ppara1*rad)
      endif
*
* Set projection parameters
*
      if (PTYPE.le.0) then
         PTYPEC='No projection'
         kmperx=0
      else if (PTYPE.eq.1) then
         PTYPEC='Equi-rectangular projection'
         kmperx=kmperx*cor1  
      else if (PTYPE.eq.2) then
         PTYPEC='Peters projection'
         factk=1/cor1**2/rad
         kmperx=kmperx*cor1
      else if (PTYPE.eq.3) then
         PTYPEC='Mercator projection'
         factl=rad/2
         factk=1/rad
         kmperx=kmperx*cor1
      else if (PTYPE.eq.4) then
         PTYPEC='Miller projection'
         factl=.4*rad
         factk=1.25/rad
         kmperx=kmperx*cor1
      else if (PTYPE.eq.5) then
         PTYPEC='Gall''s stereographic projection'
         kmperx=kmperx*cor1
         factk=(1/cor1+1)/rad
      else if (PTYPE.eq.6) then
         PTYPEC='ERS projection'
         kmperx=kmperx*cor1
*	 factk=2/cor1/rad
	 factk=2/cor1*90
	 factl=tan(82*rad)/(pi/2)
      else if (PTYPE.eq.11) then
         PTYPEC='Orthographic projection'
         azmtal=.true.
      else if (PTYPE.eq.12) then
         rview=rmean+C
         call pgnumb(nint(rview-rmean),0,1,alt,nc)
         PTYPEC='Perspective projection ('//alt(1:nc)//' km alt.)'
         azmtal=.true.
      else if (PTYPE.eq.13) then
         PTYPEC='Azimuthal equal-area projection'
         azmtal=.true.
      else if (PTYPE.eq.14) then
         PTYPEC='Azimuthal equi-distant projection'
         rview=C
         if (C.le.0.0) rview=pi*rmean
         azmtal=.true.
      else if (PTYPE.eq.21) then
         PTYPEC='Ptolemy projection'
         if (ppara1.eq.0)
     .call grwarn('PMDEF: Standard parallel must not be the equator.')
         factl=sin(ppara1*rad)*rad
         factk=1/tan(ppara1*rad)/rad+ppara1
         xcurve=.true.
      else if (PTYPE.eq.22) then
         PTYPEC='Kavraiskiy IV projection'
         if (ppara1.eq.ppara2)
     .call grwarn('PMDEF: Standard parallel must not be the identical.')
         temp=min(ppara1,ppara2)
         ppara2=max(ppara1,ppara2)
         ppara1=temp
         factl=(cor1-cor2)/(ppara2-ppara1)
         factk=(cor1/factl+ppara1+cor2/factl+ppara2)/2
         xcurve=.true.
      else if (PTYPE.eq.31) then
         PTYPEC='Sinusoidal projection'
         ycurve=.true.
      else if (PTYPE.eq.32) then
         PTYPEC='Mollweide projection'
         ycurve=.true.
      else if (PTYPE.eq.33) then
         PTYPEC='Bartholomew''s ''The Times'' projection'
         factl=rad/2.5
         factk=(1/cor1+1)/rad*cos(factl*ppara1)
         ycurve=.true.
      else if (PTYPE.eq.34) then
	 PTYPEC='Tilted rectangular projection'
	 factk=ppara1
	 factl=sin(ppara2*rad)
         kmperx=kmperx*cor1
      else if (PTYPE.eq.41 .or. PTYPE.eq.42) then
	 PTYPEC='Polar projection'
	 factk=cor1/(abs(ppara1)*rad)
	 kmperx=kmperx/factk
         xcurve=.true.
      else
         ptype=0
         PTYPEC='Illegal projection type.'
         call grwarn('PMSWIN: projection type not recognized.')
      endif
      if (azmtal) then
         kmperx=rmean
         xcurve=.true.
         ycurve=.true.
      endif
*
* Set loncen to middle longitude
*
      loncen=(A+B)/2
      latcen=(C+D)/2

      lonumin=A
      lonumax=B
      latumin=C
      latumax=D
*
* For polar projections and less than 90 degrees longitude coverage:
* set loncen to left hand side of the area
*
      if (ptype.eq.41 .and. b-a.le.90) loncen=A
      if (ptype.eq.42 .and. b-a.le.90) loncen=A
*
* Set window and adjust viewport
*
      if (ptype.eq.-1) then
         call pgswin(A,B,C,D)
      else if (ptype.eq.0) then
         call pgwnad(A,B,C,D)
      else if (ptype.eq.1) then
         call pgwnad(A*cor1,B*cor1,C,D)
         call pgswin(A,B,C,D)
      else if (ptype.eq.11) then
         call pgwnad(-1.,1.,-1.,1.)
      else if (ptype.eq.12) then
         rc=sqrt((rview-rmean)/(rview+rmean))
         call pgwnad(-rc,rc,-rc,rc)
      else if (ptype.eq.13) then
         rc=sqrt(2.)
         call pgwnad(-rc,rc,-rc,rc)
      else if (ptype.eq.14) then
         rc=rview/rmean
         call pgwnad(-rc,rc,-rc,rc)
      else
*
* Set edges in (lon,lat)
*  4  5  6 --- most Northern edge (D)
*  7     8 --- at Equator         (if appropriate)
*  1  2  3 --- most Soutern edge  (C)
*  |  |  `---- Eastern edge       (B)
*  |  `------- Mid-longitude
*  `---------- Western edge       (A)
*
      x(1)=A
      y(1)=C
      x(2)=loncen
      y(2)=C
      x(3)=B
      y(3)=C
      x(4)=A
      y(4)=D
      x(5)=loncen
      y(5)=D
      x(6)=B
      y(6)=D
      X(7)=A
      Y(7)=0.0
      X(8)=B
      Y(8)=0.0
*
* If equator not in latitude range, set 7 to 1 and 8 to 3
*
      if (C*D.GE.0.0) then
         Y(7)=Y(1)
         Y(8)=Y(3)
      endif
*
* For polar projections and more than 180 degrees longitude coverage:
* set 7 and 8 to +/- 90 degrees wrt center longitude
*
      if (ptype.eq.41 .and. b-a.ge.180) then
	 x(7)=loncen-90.
	 y(7)=c
	 x(8)=loncen+90.
	 y(8)=c
      endif
      if (ptype.eq.42 .and. b-a.ge.180) then
	 x(7)=loncen-90.
	 y(7)=d
	 x(8)=loncen+90.
	 y(8)=d
      endif
*
* Convert corners (1-8) to (X,Y)
*
      call pmconv(8,x,y)
*
* Set window in (x,y).
*
      xw0=min(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8))
      xw1=max(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8))
      yw0=min(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8))
      yw1=max(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8))
      call pgwnad(xw0,xw1,yw0,yw1)
      endif
*
* Adjust SIZE of viewport (if requested)
*
      call pgqwin(xw0,xw1,yw0,yw1)
      call pgqvp(1,xv0,xv1,yv0,yv1)
      scale=kmperx*(xw1-xw0)/(xv1-xv0)/25.4e-6
      if (pscale.gt.0) then
         xvm=(xv0+xv1)/2
         xvs=(xv1-xv0)/2*scale/pscale
         yvm=(yv0+yv1)/2
         yvs=(yv1-yv0)/2*scale/pscale
         call pgvsiz(xvm-xvs,xvm+xvs,yvm-yvs,yvm+yvs)
      else
         pscale=scale
      endif
*
* Determine actual map boundaries
*
      if (ptype.eq.1) then
         lonmin=a
         lonmax=b
         latmin=c
         latmax=d
      else if (azmtal) then
         loncen=a
         latcen=b
         lonmin=a-180
         lonmax=a+180
         latmin=-90
         latmax=+90
      else if (ptype.eq.32 .or. ptype.eq.41 .or. ptype.eq.42) then
         lonmin=loncen-180
         lonmax=loncen+180
         latmin=-90
         latmax=+90
      else
         x(1)=xw0
         y(1)=yw0
	 x(2)=(xw0+xw1)/2
         y(2)=yw0
	 x(3)=xw1
         y(3)=yw0
         x(4)=xw0
         y(4)=yw1
	 x(5)=x(2)
         y(5)=yw1
         x(6)=xw1
         y(6)=yw1
         x(7)=xw0
	 y(7)=(yw0+yw1)/2
         x(8)=xw1
         y(8)=y(7)
         call pmcinv(8,x,y)
         lonmin=min(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8))
         lonmax=max(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8))
         latmin=min(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8))
         latmax=max(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8))
      endif
      sinb0=sin(latcen*rad)
      cosb0=cos(latcen*rad)
      pmopen=4

      return
      end
