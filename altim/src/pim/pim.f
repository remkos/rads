**PIM -- Perfect Image Maker
*
* This is PIM -- DUT/SSR&T's Perfect Image Maker.
*
* This graphics package requires input in the form of 'input cards'
* consisting of a command with arguments. Read more about it
* in pim.tex
*
* 1994:
* 21-Jul - v2.0 - Resurrection of PIM
* 29-Jul - v2.1 - Better scanning of command lines. Less memory needed
*        (only three major buffers)
*  9-Aug - v2.2 - Quantisation of colours. More comments. Lots of changes in
*        subroutines. No logo for postscript. VIEWPORT BIG default.
* 12-Aug - v2.3 - Better shading (first slopes, than interpolation).
* 16-Sep - v2.4 - Proper (black) colour for TITLE in PostScript mode.
* 27-Sep - v2.5 - New XGFPLOT possibilities. Quit when no grid.
* 11-Oct - v2.6 - XGFPLOT revised. Bugs removed. Change of input.
*                 RANGE card included.
* 28-Oct - v2.7 - Bug in GRIDREAD removed that prevented adding of grids
* 16-Nov - v2.8 - PROFILE card included. New GR1INT.
* 23-Nov - 9411.2 - Bug removed in interpolation mode. Default interpolation for
*		  Postscript devices limited to cells 36 per inch.
*  2-Dec - 9412.1 - Multiple CONTOURS. Free choise of grid dimensions.
*  5-Dec - 9412.2 - Bug removed that prevented contours in PS mode.
*  9-Dec - 9412.3 - Use external gr1int4 routine from gridlib
* 13-Dec - 9412.4 - Include TEXT and SYMBOLS command
* 20-Dec - 9412.5 - New colourmap format. INCLUDE command. Command file
*                 name on command line.
* 1995:
*  3-Jan - 9501.1 - Include velocity plotting with VELOCITY command.
*  4-Jan - 9501.2 - Include velocity plotting with DYNTOPO command. VELOCITY
*                 removed. Improved setup. PGSAVE around each command. Separate
*                 command routines. SSRT logo. Include -h/-v flag in
*                 LEGEND. Properly handle GRIDAREA.
*  5-Jan - 9501.3 - Include velocity vector in legend. Unit above vertical
*                 legend moved up and left.
*  6-Jan - 9501.4 - ANNOTATE included. Replaces NOTATE.
* 13-Jan - 9501.5 - New vector sizes in legend. Masking vectors around equator
* 13-Jan - 9501.6 - New gridrd routine that makes it possible to select
*                 an area from a very large grid. Also, it automatically
*                 wraps the grid around the 'date-line' if needed.
* 17-Jan - 9501.7 - Include selection on time and tracknumber in XGF command
* 24-Jan - 9501.8 - Improvement of the subgrid selection
* 27-Jan - 9501.9 - Adding grids of unequal resolution
* 30-Jan - 9501.10 - Added sig= option to surface and illum grids. Sets
*                  plot interval depending on statistics
* 31-Jan - 9501.11 - Polished out some bugs in contouring. And removed
*                  the bugs in shading and rounding which were introduced
*                  yesterday :-)  Dumping profile for steps > 0.5 deg
*  1-Feb - 9502.01 - Removed bug in sigma level for ILLUM
* 15-Feb - 9502.02 - Removed bug in determining mean/rms in Gridstat.
*		   Disallow interpolation if PROJECT=1 and INTER is off.
*		   Include velocity field plotting (VELOCITY).
* 17-Feb - 9502.03 - Properly handle non-interpolated grids (cells !)
* 21-Feb - 9502.04 - Removed bugs with INTER=off and PROJECT>1.
* 21-Feb - 9502.05 - Introduce PROJECT=0 and MESH.
* 22-Feb - 9502.06 - Change position of title. Change determination of subarea
*                  in sgridrd4.f
* 23-Feb - 9502.07 - Provide selection on tracknr and time in PROFILE.
* 28-Feb - 9502.08 - Include VALUE command.
*  8-Mar - 9503.01 - Add angle= and just= options in TEXT and VALUE.
* 10-Mar - 9503.02 - Include possibility to have only LEGEND.
* 13-Mar - 9503.03 - Include groundtrack plotting with GROUND
* 14-Mar - 9503.04 - Include all levels of the WDB file. Optional shading
*                  range in ILLUM. New illumination algorithm.
* 13-Apr - 9504.01 - Geosat orbits included in TRACKS.
* 21-Apr - 9504.02 - TRACKS properly cut off at longitude edges.
* 26-Apr - 9504.03 - Extra workspace (WORK5) for ILLUM
*  3-May - 9505.01 - WORK5 used again in ILLUM to overcome clash with CONTO
*  4-May - 9505.02 - Gridstatistics added in CMDCONTO, also when grid is
*                  already read in CMDSURFA. This will create contours
*                  throughout the entire range again.
* 19-May - 9505.03 - New logos.
* 16-Jun - 9506.01 - XXO files included in XGF command. Coast and Box moved
*                  till after XGF and TRACK
* 27-Jun - 9506.02 - cutoff removed from VELOCITY
* 28-Jun - 9506.03 - include projection 11.
*  3-Jul - 9507.01 - removed bug with true scale parallel.
* 18-Sep - 9509.01 - SYMBOL enhanced with -yx and -xy options.
* 28-Sep - 9509.02 - XGF enhanced with plotting of xXXO and aADR
* 15-Nov - 9511.01 - asc/des selection in PROFILE
*
* 1996:
* 25-Jan - 9601.01 - Aligned with PGPLOT 5.0.3. Optional labels in contouring.
*  8-Apr - 9604.01 - Area check in CMDSYMBO included.
*  5-May - 9605.01 - New gridstat (weighted).
*  7-Nov - 9611.00 - gridstat replaced by grstat4.
*
* 1997:
* 22-Apr - 9704.00 - CMDTRACK changed to have tracks up to the edge of the
*                    plotarea.
* 15-Jul - 9707.00 - buf= flag added to CMDXGF. ci= flag added to CMDCONTO.
*  5-Aug - 9708.00 - Added CMDCIRCL.
* 23-Oct - 9710.00 - Changes CMDXGF to do area checking and longitude conversion.
*
* 1998:
* 15-Feb - 9802.0 - Multiplatform version - Use $ALTIM/pim
* 11-Mar - 9803.0 - Bugs in selections removed in CMDXGF
* 23-Mar - 9803.1 - Provide window info when -debug is used
* 17-Apr - 9804.0 - Swap order of text and symbol (symbol first)
* 26-Nov - 9811.0 - Extend card lines to 1024 characters. Slopes in /day.
*
* 1999:
* 23-Feb - 9902.0 - Add posibility of HLS values in colour maps
* 13-Apr - 9904.0 - Add LINE command
* 18-Apr - 9904.1 - Add COG and BAD commands
*
* 2000:
* 24-Sep - 0009.0 - Added GPS track and waypoint plotting
*
* 2001:
* 29-Sep - 0109.0 - Fixed bugs related to gridread
* 30-Sep - 0109.1 - Introduced reading of complex grids
*
* 2003:
*  6-Jun - 0306.0 - Add option to plot logo at end
*-----------------------------------------------------------------------
      program pim

      character*6	version
      parameter (version='0506.0')

      include 'pim.inc'

      logical log,pimopt
      real*4  scale

      integer iarg,iargc

      real xblc,xtrc,yblc,ytrc
      real pixx0,pixx1,pixy0,pixy1

      real angle,dx,dy

* Start the job

      write (0,'("This is PIM, version ",a/)') version
      write (0,550) 'Copying input cards to scratch file ...'
      open (7,status='scratch')

* Copy command file(s) to scratch

      call prepro1('default.pim')
      log=.false.
      do iarg=1,iargc()
	 call getarg(iarg,argum)
	 if (argum.eq.'-debug') then
	    debug=.true.
	 else if (argum.eq.'-quiet') then
	 else if (argum(1:2).eq.'-h') then
	    write (0,3050) version
	    goto 9010
	 else
	    call prepro1(argum)
	    log=.true.
	 endif
      enddo
      if (.not.log) call prepro1('-')

* Set defaults

      sname=' '
      cname=' '
      lname=' '
      dumped=.true.
      nsx=0
      nsy=maxgrd
      nlx=0
      nly=maxgrd
      ncx=0
      ncy=maxgrd
      shade0=0.20
      shade1=0.95

* Open device, viewport, and get projection

      call cmddevic
      call cmdviewp
      call cmdproje
      call cmdrange
      call cmdbad
      call cmddithe

* Get the plotarea. This will come from the PLOTAREA card, GRIDAREA card,
* or one of the grids. Get also gridinfo.

      if (popcmd('PLOTA',argum)) then
	 read (argum,*,iostat=ios) xw0,xw1,yw0,yw1
         call setarea(2,nsx,nsy,xs0,xs1,ys0,ys1,.false.)
         call setarea(2,nlx,nly,xl0,xl1,yl0,yl1,.false.)
         call setarea(2,ncx,ncy,xc0,xc1,yc0,yc1,.false.)
      endif
      if (popcmd('PLOTC',argum)) then
	 dx=360
	 dy=180
	 read (argum,*,iostat=ios) xw0,yw0,dx,dy
	 xw0=xw0-dx/2
	 yw0=yw0-dy/2
	 xw1=xw0+dx
	 yw1=yw0+dy
         call setarea(2,nsx,nsy,xs0,xs1,ys0,ys1,.false.)
         call setarea(2,nlx,nly,xl0,xl1,yl0,yl1,.false.)
         call setarea(2,ncx,ncy,xc0,xc1,yc0,yc1,.false.)
      endif
      if (pop1cmd('SURFA',sargum)) then
	 call strip1(sargum,sname)
	 cells=pimopt('-c',sargum,dum,dum,dum,dum)
	 if (popnextcmd('GRIDA',argum))
     |		read (argum,*,iostat=ios) xs0,xs1,ys0,ys1
	 call gridinfo(sname,nsx,nsy,xs0,xs1,ys0,ys1,rmins,rmaxs)
	 call setarea(1,nsx,nsy,xs0,xs1,ys0,ys1,cells)
	 rewind (7)
      endif
      if (pop1cmd('ILLUM',largum)) then
	 call strip1(largum,lname)
	 celll=pimopt('-c',largum,dum,dum,dum,dum)
	 log=pimopt('shade=',largum,shade0,shade1,dum,dum)
	 if (lname.eq.' '.or.lname.eq.'=') lname=sname
	 if (popnextcmd('GRIDA',argum))
     |		read (argum,*,iostat=ios) xl0,xl1,yl0,yl1
	 call gridinfo(lname,nlx,nly,xl0,xl1,yl0,yl1,rminl,rmaxl)
	 call setarea(1,nlx,nly,xl0,xl1,yl0,yl1,celll)
	 rewind (7)
      endif
105   if (pop1cmd('CONTO',cargum)) then
 	 call strip1(cargum,cname)
	 cellc=pimopt('-c',cargum,dum,dum,dum,dum)
 	 if (cname.eq.' '.or.cname.eq.'=') cname=sname
	 if (popnextcmd('GRIDA',argum))
     |		read (argum,*,iostat=ios) xc0,xc1,yc0,yc1
 	 call gridinfo(cname,ncx,ncy,xc0,xc1,yc0,yc1,rminc,rmaxc)
	 call setarea(1,ncx,ncy,xc0,xc1,yc0,yc1,cellc)
	 goto 105
      endif
      if (pop1cmd('AMPHA',cargum)) then
 	 call strip1(cargum,cname)
	 if (popnextcmd('GRIDA',argum))
     |		read (argum,*,iostat=ios) xc0,xc1,yc0,yc1
 	 call gridinfo(cname,ncx,ncy,xc0,xc1,yc0,yc1,rminc,rmaxc)
	 call setarea(1,ncx,ncy,xc0,xc1,yc0,yc1,cellc)
 	 call strip1(cargum,sname)
	 if (popnextcmd('GRIDA',argum))
     |		read (argum,*,iostat=ios) xs0,xs1,ys0,ys1
 	 call gridinfo(sname,nsx,nsy,xs0,xs1,ys0,ys1,rmins,rmaxs)
	 call setarea(1,nsx,nsy,xs0,xs1,ys0,ys1,cells)
      endif
      if (xw0.eq.xw1 .or. yw0.eq.yw1) stop "pim: no PLOTAREA"

* Set all colours

      write (0,550) 'Setting colours ...'
      call pimcol(lname.ne.' ')

* Print background and logo

      call cmdbackg
      call cmdlogo('LOGO',.false.)
      call cmdlogo('DEOS',.false.)

      call cmdtitle
      call cmdlegen1

* Now set the gridarea boundaries for those that are not already set

      call setarea(2,nsx,nsy,xs0,xs1,ys0,ys1,cells)
      call setarea(2,nlx,nly,xl0,xl1,yl0,yl1,celll)
      call setarea(2,ncx,ncy,xc0,xc1,yc0,yc1,cellc)

* Set interpolation resolution. If the INTER card is provided, use the
* number of cells given on the card. If no cells are specified, determine
* the device resolution though pgqvp (3,..). For PostScript devices this is
* 1000 cells per inch, which is beyond the actual (visible) resolution.
* Hence, reduce this to 50 cells per inch. If interpolation is performed
* to the highest resolution, fast contouring can be used.
*
* Use INTER -1 to prevent any interpolation.

      call pgqvp(3,pixx0,pixx1,pixy0,pixy1)
      call pgqwin(xblc,xtrc,yblc,ytrc)
      if (debug) then
         write (*,*) 'window:',xblc,xtrc,yblc,ytrc
	 write (*,*) 'pixels:',pixx0,pixx1,pixy0,pixy1
      endif
      if (project.ne.0) call pmqdef(argum,scale)
      if (popcmd('INTER',argum)) then
	 log=pimopt('off',argum,dum,dum,dum,dum)
	 read (argum,*,iostat=ios) px,py
	 if (log) then
	    inter=-1
	    if (nsx.eq.0) nsy=0
	    if (nlx.eq.0) nly=0
	    if (ncx.eq.0) ncy=0
	    px=max(nsx,nlx,ncx)
	    py=max(nsy,nly,ncy)
	 else if (px*py.eq.0) then
	    px=nint(pixx1)-nint(pixx0)+1
	    py=nint(pixy1)-nint(pixy0)+1
	    px=nint(pixx1-pixx0)+1
	    py=nint(pixy1-pixy0)+1
*	    write (*,*) 'pixx0,pixx1,px',pixx0,pixx1,px
*	    write (*,*) 'pixy0,pixy1,py',pixy0,pixy1,py
	    if (postscript) then
	       px=px/20
	       py=py/20
	       inter=1
	    else
	       inter=2
	    endif
	 else
	    inter=1
	 endif
      else
	 if (nsx*nlx.ne.0 .and. project.ne.1)
     |		stop "pim: need interpolation for projection type not 1"
	 if (nsx.eq.0) nsy=0
	 if (nlx.eq.0) nly=0
	 if (ncx.eq.0) ncy=0
	 px=max(nsx,nlx,ncx)
	 py=max(nsy,nly,ncy)
	 if ((nsx.ne.0 .and. (nsx.ne.px .or. nsy.ne.py)) .or.
     |       (nlx.ne.0 .and. (nlx.ne.px .or. nly.ne.py)) .or.
     |       (ncx.ne.0 .and. (ncx.ne.px .or. ncy.ne.py)))
     |		stop "pim: cannot use incompatible grids without INTER"
	 inter=0
	 if (project.gt.1) inter=1
      endif
      if (debug) write (*,*) 'inter:',inter,px,py

* Process surface, illumination, and contour grids

      call cmdland(0)
      call cmdcoast(0)
      call cmdsurfa
      call cmdillum
      call cmdconto

* If pixel map not dumped, dump it now

      if (.not.dumped) call dumppix

* Print info

      if (popcmd('INFO',argum)) then
	 write (0,3001) xw0,xw1,yw0,yw1
	 if (sname.ne.' ') write (0,3010) sname,rmins,rmaxs
	 if (cname.ne.' ') write (0,3020) cname,rminc,rmaxc,rintc
	 if (lname.ne.' ') write (0,3030) lname,rminl,rmaxl,angle
      endif

      call cmdmesh
      call cmdvalue
      call cmdannot(rintc)
      call cmdbathy
      call cmddynto
      call cmdveloc
      call cmdland(1)
      call cmddata
      call cmdxgf
      call cmdprofi
      call cmdcircl
      call cmdcoast(1)
      call cmdtrack
      call cmdgpstr
      call cmdwaypo
      call cmdrect
      call cmdbox
      call cmdcog
      call cmdline
      call cmdsymbo
      call cmdtext
      call cmdlegen2
      call cmdpie

* Update screen

      call pgupdt
      call showtrnd
      call cmdlogo('LOGO',.true.)
      call cmdlogo('DEOS',.true.)
      call cmdcurso(rintc)

* Closing plot

      write (0,550) 'Closing plot ...'
      call pgend

* Formats

550   format(a)
3001  format (' Window: Longitude: from : ',f6.1,10x,'to : ',f6.1/
     |        '         Latitude : from : ',f6.1,10x,'to : ',f6.1)
3010  format (' Surface file name      : ',a40/
     |        ' Background level       : ',f6.2,
     |        '  foreground level      : ',f6.2)
3020  format (' Contour file name      : ',a40/
     |        ' Minimum level          : ',f8.2,
     |        '  maximum level         : ',f8.2/
     |        ' 1st Contour interval   : ',f8.2/
     |        ' 2nd Contour interval   : ',f8.2)
3030  format (' Illumination file name : ',a40/
     |        ' Minimum level          : ',f8.2,
     |        '  maximum level         : ',f8.2/
     |        ' Illumination angle     : ',f8.2)
3040  format (/25x,'All is well that ends well')
3050  format ('PIM ',a,' -- The Perfect Image Manager'//
     |'usage: pim commandfile(s)'/'   or: pim < commandfile'//
     |'where the commandfile(s) may contain the following commands:'/
     |'DEVICe     : specify plot device name'/
     |'SURFAce    : specify surface file and range'/
     |'ILLUMinate : specify illumination'/
     |'CONTOur    : specify contour file and levels'/
     |'ANNOTate   : annotate contours'/
     |'GRIDArea   : specify grid boundaries'/
     |'PLOTArea   : specify plotting boundaries'/
     |'VIEWPort   : specify viewport boundaries'/
     |'PROJEction : specify projection'/
     |'COLOUrmap  : load colourmap'/
     |'INTERpolate: interpolate input grids before plotting'/
     |'DITHER     : specify dithering'/
     |'C_????     : specify single colours'/
     |'INCLUde    : include other commandfile'/
     |'TEXT       : plot text'/
     |'SYMBOls    : plot symbols'/
     |'LEGENd     : plot a legend'/
     |'TITLE      : plot title'/
     |'LAND       : plot landmask'/
     |'COASTlines : plot coastlines'/
     |'BATHYmetry : plot bathymetric lines'/
     |'DATA       : plot ASCII dataset'/
     |'XGF        : plot XGF dataset'/
     |'PROFIle    : plot profile from XGF file'/
     |'DYNTOpo    : plot velocity vectors according to dynamic',
     |' topography'/
     |'BOX        : specify PMBOX parameters'/
     |'BACKGround : plot DEOS logo on background'/
     |'LOGO       : plot DEOS logo')
*
      write (0,3040)
9010  close (7)
      end
