      subroutine cmload
      include 'xrgb.inc'

      logical exist
      character*80 filenm,cmapdir
      character*5 command
      integer l,ic,ios,ir,ig,ib,idx
      real r,g,b

      if (changed) then
	 call menu(0,' ')
	 call menu(13,'Changed are not saved. Prepared to lose changes ?')
	 call menu(11,'Quit loading')
	 call menu(12,'Lose changes')
	 call pgcurse(x,y,ch)
	 call whatindx(idx)
	 if (idx.ne.-12) return
      endif
      do ic=0,maxind
	 cmap(1,ic)=-1
      enddo

      filenm=cmfile
      inquire (file=filenm,exist=exist)
      if (.not.exist) then
	 call getenv('CMAP_DIR',cmapdir)
	 l=index(cmapdir,' ')-1
	 filenm=cmapdir(:l)//'/'//cmfile
      endif
      open (92,file=filenm,status='old')
      read (92,'(a)') filenm
      if (filenm(1:4).ne.'#COL') stop 'illegal colourmap'

20	 read (92,'(a)',end=30) filenm
	 l=index(filenm,' ')-1
	 command=filenm(1:l)
	 call grtoup(command,command)
	 read (filenm(l+1:),*,iostat=ios) ir,ig,ib
	 if (command(1:1).eq.'#') then
	 else if (command.eq.'C_BG') then
	    ic=0
	 else if (command.eq.'C_FG') then
	    ic=1
	 else if (command.eq.'C_LAN') then
	    ic=2
	 else if (command.eq.'C_BAD') then
	    ic=3
	 else if (command.eq.'C_COA') then
	    ic=4
	 else if (command.eq.'C_CON') then
	    ic=5
	 else if (command.eq.'C_RNG') then
	    ic=ic+1
	    if (ic.lt.16) ic=16
	 endif
	 if (index(filenm,'HLS').gt.0) then
	    call grxrgb(ir*1.,ig/100.,ib/100.,r,g,b)
	    cmap(1,ic)=nint(r*255)
	    cmap(2,ic)=nint(g*255)
	    cmap(3,ic)=nint(b*255)
	 else
	    cmap(1,ic)=ir
	    cmap(2,ic)=ig
	    cmap(3,ic)=ib
	 endif
	 goto 20

30    do ic=0,maxind
	 call colupdt(ic)
      enddo
      changed=.false.
      close (92)
      end
