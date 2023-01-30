      subroutine iwritef(filenm,iform,kb,nxp,nyp,quit)
*
* Write file in one of three formats
*
      implicit none
      include "ingrid.inc"
      logical sea(MA),error,quit
      character*80 filenm
      integer nxp,nyp,iform,i,nspace,kb,gridwr4
      equivalence (a,sea)
*
  550 format (a)
 1000 format ((6e12.5))
*
      quit=.false.
      error=.true.
      nspace=nxp*nyp
      if (iform.eq.1) then
*
* BINARY (UNFORMATTED) GRID format
*
         error=(gridwr4(filenm,nxp,nyp,a(id(kb)),nxp,
     1		xmin(kb),xmax(kb),ymin(kb),ymax(kb)).ne.0)
      else if (iform.eq.3) then
*
* function file
*
         open (10,file=filenm,form='formatted',status='new',err=130)
         write (10,1000,err=100) (a(id(kb)+i-1),i=1,nspace)
         close (10)
         error=.false.
      else
*
* mask file
*
         write (0,550) 'ingrid: maskwr function no longer available'
      endif
  100 if (error) then
         write (0,550) 'ingrid: write unsuccessful'
      else
         fname(kb)=filenm
         proc(kb)=.false.
         do i=1,MBUF
            if (i.ne.kb .and. fname(i).eq.fname(kb)) fname(i)=' '
         enddo
      endif
      return
*
*  130 call perror('ingrid')
  130 write(0,*) 'ingrid: error opening ' // filenm
      return
*
      quit=.true.
      end
