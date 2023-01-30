      implicit none
      integer*4 adr4a(7),adr4b(7)
      integer*2 adr2a(14),adr2b(14),extra(2),bound(4)
      equivalence (adr4a,adr2a),(adr4b,adr2b)
      character*80 file1,file2,spec*4,satel*8
      integer nrec1,nrec2,mrec1,mrec2,mrec3,mrec4,mrec5,i

      call getarg(1,file1)
      call getarg(2,file2)

      open (11,file=file1,form='unformatted',status='old',
     | access='direct',recl=28)
      read (11,rec=1) spec,satel,bound,nrec1,extra
      open (12,file=file2,form='unformatted',status='old',
     | access='direct',recl=28)
      read (12,rec=1) spec,satel,bound,nrec2,extra

      open (13,file='OnlyIn1.adr',form='unformatted',status='new',
     | access='direct',recl=28)
      open (14,file='OnlyIn2.adr',form='unformatted',status='new',
     | access='direct',recl=28)
      open (15,file='Diff1-2.adr',form='unformatted',status='new',
     | access='direct',recl=28)
      
      mrec1=1
      mrec2=1
      mrec3=0
      mrec4=0
      mrec5=0

10    continue
      if (mrec1.gt.nrec1-1000 .or. mrec2.gt.nrec2-1000)
     |write (*,*) mrec1,nrec1,mrec2,nrec2
      read (11,rec=mrec1+1) adr4a
      read (12,rec=mrec2+1) adr4b

      if (adr4a(1).eq.adr4b(1) .and. adr4a(2).eq.adr4b(2)) then
         adr4a(5)=adr4a(5)-adr4b(5)
	 do i=11,14
            adr2a(i)=adr2a(i)-adr2b(i)
	 enddo
	 mrec5=mrec5+1
         write (15,rec=mrec5+1) adr4a
         mrec1=mrec1+1
         mrec2=mrec2+1
      else if (adr4a(1)+adr4a(2)/1d6.lt.adr4b(1)+adr4b(2)/1d6) then
	 mrec3=mrec3+1
         write (13,rec=mrec3+1) adr4a
         mrec1=mrec1+1
      else
	 mrec4=mrec4+1
         write (14,rec=mrec4+1) adr4b
         mrec2=mrec2+1
      endif

      if (mrec1.le.nrec1 .and. mrec2.le.nrec2) goto 10

      write (13,rec=1) spec,satel,bound,mrec3,extra
      write (14,rec=1) spec,satel,bound,mrec4,extra
      write (15,rec=1) spec,satel,bound,mrec5,extra

      close (11)
      close (12)
      close (13)
      close (14)
      close (15)

      end
