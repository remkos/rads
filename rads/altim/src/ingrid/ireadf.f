      subroutine ireadf(iform,quit)
*
* Read file(s) in one of three formats
* Files are specified in unit# 11
*
      implicit none
      include "ingrid.inc"
      logical sea(MA),error,quit
      character*80 filenm
      equivalence (a,sea)
      integer nxp,nyp,nxc,nyc,ipnt,number,kb,i,iform,gridrd4
      real x0,x1,y0,y1,z0,z1
*
  550 format (a)
  551 format (a,$)
  570 format (a35,' -> ',$)
  575 format (a22,' -> Try again -> ',$)
  580 format ('      Target buffer ->',i3)
 1000 format ((6e12.5))
*
      quit=.false.
      ipnt=0
      number=0
*
* If function file type selected, ask for size first
*
      if (iform.eq.3) then
         write (0,550) 'Function file selected'
    5    write (0,551) 'Enter number of cells in x- and y-direction -> '
         read (5,*,end=190,err=5) nxc,nyc
         if (nxc.le.0 .or. nyc.le.0) then
            write (0,550) 'ingrid: at least one dimension too small'
            goto 20
         endif
      endif
*
   20 read (11,550,err=200,end=200) filenm
      write (0,570) filenm
      nxp=0
      nyp=MA-idnext+1
      if (iform.eq.1) then
*
* Read BINARY (UNFORMATTED) grid
*
         error=(gridrd4(filenm,nxp,nyp,a(idnext),
     1		x0,x1,y0,y1,z0,z1).ne.0)
      else if (iform.eq.3) then
*
* Function file type selected
*
         nxp=nxc+1
         nyp=nyc+1
         open (10,file=filenm,status='old',form='formatted',err=140)
         rewind (10)
         read (10,1000,end=100,err=100) (a(idnext+i-1),i=1,nxp*nyp)
         close (10)
      else
*
* Mask file type selected
*
*        call maskrd(sea(idnext),nxp,nyp,filenm)
*        do 50 i=idnext,idnext+nxp*nyp-1
*           if (sea(i)) then
*              a(i)=1
*           else
*              a(i)=0
*           endif
*  50    continue
      write (0,550) 'ingrid: maskrd function no longer available'
      endif
      if (nxp*nyp.eq.0) goto 20
*
* If no more target buffer numbers are memorized, ask for more
*
   60 ipnt=ipnt+1
      if (ipnt.gt.number) then
         call multin('Enter target buffer',MBUF,number,idata,error,quit)
         if (quit) goto 190
         if (error) goto 20
         ipnt=1
         kb=idata(1)
      else
         kb=idata(ipnt)
         write (0,580) kb
         error=.false.
      endif
      if (quit) goto 190
      if (error) goto 20
*
* Check for unwritten changes
* If illegal buffer number was selected, read next file
* If modified file was kept or buffer was blocked, ask for other buffer number
*
      call noclob(kb,error,quit)
      if (quit) goto 190
      if (kb.eq.0) goto 20
      if (error) goto 130
      if (blkd(kb)) goto 120
*
      nx(kb)=nxp
      ny(kb)=nyp
      xmin(kb)=x0
      xmax(kb)=x1
      ymin(kb)=y0
      ymax(kb)=y1
      id(kb)=idnext
      idnext=idnext+nxp*nyp
      proc(kb)=.false.
      used(kb)=.true.
      fname(kb)=filenm
      blkd(kb)=.true.
      goto 20
*
  100 close (10)
      write (0,550) 'ingrid: error during read; no data stored'
      goto 20
*
  120 write (0,550) 'ingrid: can not use this buffer number'//
     .' more than once during this action'
  130 write (0,575) filenm
      goto 60
*
*  140 call perror('ingrid')
  140 write(0,*) 'ingrid: error opening ' // filenm
      goto 20
*
* Error return after EOD on user input
*
  190 quit=.true.
  200 end
