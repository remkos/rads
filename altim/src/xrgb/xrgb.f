**XRGB -- Colourmap managing program
*+
      program xrgb
*
* Program to view and change colourmaps for PIM
*
* usage: xrgb [colormap(s)]
*
* 9412.1 - First operational version
* 9902.0 - Better initialisation - Printing
*-
      real version
      parameter (version=9902.0)

      integer idx,i,kx,ky
      integer iarg,narg,iargc
      real x0,x1,y0,y1
      character*80 text
      include 'xrgb.inc'

* Initialize position of buttons 1-12

      dbx=0.15
      dby=0.06
      do i=0,10,2
	 bx0(i+1)=(1-dbx)*i/10
	 bx1(i+1)=bx0(i+1)+dbx
	 bx0(i+2)=bx0(i+1)
	 bx1(i+2)=bx0(i+2)+dbx
	 by1(i+1)=0.23
	 by0(i+1)=by1(i+1)-dby
	 by1(i+2)=0.14
	 by0(i+2)=by1(i+2)-dby
      enddo

* Initialize text of buttons

      ctxt(0)='BG'
      ctxt(1)='FG'
      ctxt(2)='LAN'
      ctxt(3)='BAD'
      ctxt(4)='COA'
      ctxt(5)='CON'
      do i=16,maxind
         write (ctxt(i),610) i-15
	 used(i)=.false.
      enddo
      
* Initialize text for RGB and HLS

      textrgb(1)='R'
      textrgb(2)='G'
      textrgb(3)='B'
      texthls(1)='H'
      texthls(2)='L'
      texthls(3)='S'
      write (text,600) version

* Initialize colour fields

      xblksize=1./16
      yblksize=.75/9
      do kx=1,16
         x0=(kx-1)*xblksize
         x1=x0+xblksize
         do ky=9,1,-1
            y0=0.25+(ky-1)*yblksize
            y1=y0+yblksize
            i=(9-ky)*16+(kx-1)
            xl(i)=x0
            xr(i)=x1
            yb(i)=y0
            yt(i)=y1
	 enddo
      enddo

* Initialize arrays

      do i=0,255
         cmap(1,i)=0
	 cmap(2,i)=0
	 cmap(3,i)=0
      enddo

  550 format (a)
  551 format (a,$)
  600 format ('XRGB ',f6.1)
  610 format (i3.3)

* Open window

      call pgbeg(0,'/xw',1,1)
      call pgsvp(0.,1.,0.,1.)
      call pgswin(0.,1.,0.,1.)
      call pgask(.false.)
      call boxes

* Scan colourmap names from the command line

      narg=iargc()
      iarg=1

* Main menu

      call getarg(iarg,cmfile)
      call cmload
30    call menu(0,' ')
      call menu(3,'Reload')
      call menu(4,'Load ...')
      call menu(5,'Save')
      call menu(6,'Save as ...')
      call menu(7,'Print')
      call menu(8,text)
      call menu(9,'Change')
      call menu(10,'Copy')
      call menu(11,'Linear')
      call menu(12,'Exit')
35    call menu(13,'File: '//cmfile)
      if (iarg.gt.1) call menu(1,'Previous')
      if (iarg.lt.narg) call menu(2,'Next')
40    call pgcurse(x,y,ch)
      ich=ichar(ch)
      call whatindx(idx)
      if (ich.ne.65) then
	 goto 40
      else if (idx.ge.0) then
	 call cmchange(idx)
	 goto 30
      else if (idx.eq.-1 .and. iarg.gt.1) then
	 iarg=iarg-1
	 call getarg(iarg,cmfile)
	 call cmload
         goto 35
      else if (idx.eq.-2 .and. iarg.lt.narg) then
	 iarg=iarg+1
	 call getarg(iarg,cmfile)
	 call cmload
         goto 35
      else if (idx.eq.-3) then
	 call cmload
         goto 35
      else if (idx.eq.-4) then
	 write (6,551) 'Load colourmap: enter name -> '
	 read (5,550) cmfile
	 call cmload
         goto 35
      else if (idx.eq.-5) then
         call cmsave
         call menu(13,'Saved: '//cmfile)
	 goto 40
      else if (idx.eq.-6) then
	 write (6,551) 'Save colourmap: enter name -> '
	 read (5,550) cmfile
         call cmsave
         call menu(13,'Saved: '//cmfile)
	 goto 40
      else if (idx.eq.-7) then
         call pgopen('xrgb.ps/cps')
         call pgsvp(0.,1.,0.,1.)
         call pgswin(0.,1.,0.,1.)
         call pgask(.false.)
         call cmload
         call boxes
	 call pgend
	 call pgclos
	 goto 9999
      else if (idx.eq.-8) then
         call menu(13,'This is '//text)
	 goto 40
      else if (idx.eq.-9) then
         call cmchange(idx)
	 goto 30
      else if (idx.eq.-10) then
         call cmcopy
	 goto 30
      else if (idx.eq.-11) then
         call cmlinear
	 goto 30
      else if (changed) then
	 call menu(0,' ')
	 call menu(13,'Changes are not saved. Exit anyway ?')
	 call menu(11,'Do not exit')
	 call menu(12,'Exit now')
	 call pgcurse(x,y,ch)
	 call whatindx(idx)
	 if (idx.ne.-12) goto 30
      endif
      call menu(-1,' ')
      call menu(13,'Exiting ...')
      call pgend
9999  end
