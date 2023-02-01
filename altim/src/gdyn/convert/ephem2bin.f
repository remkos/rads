c
c read ASCII EPHEM table and write to binary format
c
      program epem2bin
      implicit none
      integer*4 mbuf
      parameter (mbuf=2000)
      character*80 input,output
      integer*4 ihead(5),ibuf(mbuf)
      real*8 dbuf(mbuf)
      integer*4 keypos(8)/999999999,1010,1020,1030,1040,1041,1050,1070/
      integer*4 ivers/3/,nt1041/0/,nt1050/0/,ndim/784/
      integer*4 i,iunit/10/,ounit/20/,iverss/0/,ncount

      call getarg(1,input)
      call getarg(2,output)
      if (output.eq.' ') then
	 output=input
	 input='-'
      endif
      if (input.eq.'-') then
         iunit=5
      else
         open(unit=iunit,status='old',form='formatted',file=input)
      endif
      open(unit=ounit,status='unknown',form='unformatted',file=output)

c header
  100 read(iunit,'(5i15)') ihead
      write(ounit) ihead
      do i=1,8
         if (ihead(4).ne.keypos(i)) cycle
         goto (200,200,200,400,200,600,700,800),i
      enddo

  200 read(iunit,'(i15)',end=999) ncount
      read(iunit,'(6i15)',end=999) (ibuf(i),i=1,ncount)
      write(ounit) ncount,(ibuf(i),i=1,ncount)
      if (ncount.eq.1) goto 100
      goto 200

  400 read(iunit,'(i15)',end=999) ncount
      read(iunit,'(6d24.17)',end=999) (dbuf(i),i=1,ncount)
      write(ounit) ncount,(dbuf(i),i=1,ncount)
      if (ncount.eq.1) goto 100
      ivers=nint(dbuf(6))
      if (ivers.ne.3.and.ivers.ne.4.and.ivers.ne.6) ivers=5
      iverss=ivers-2
      goto 400

  600 read(iunit,'(i15)',end=999) ncount
      read(iunit,'(6d24.17)',end=999) (dbuf(i),i=1,ncount)
      write(ounit) ncount,(dbuf(i),i=1,ncount)
      nt1041=nt1041+1
      if (ncount.eq.1) goto 100
      if (iverss.gt.2) goto 600
      if (nt1041.ne.1) goto 610
      goto 600
  610 if (nt1041.ne.2) goto 600
      goto 600

  700 read(iunit,'(i15)',end=999) ncount
      read(iunit,'(6i15)',end=999) (ibuf(i),i=1,ncount)
      write(ounit) ncount,(ibuf(i),i=1,ncount)
      if (ncount.eq.1) goto 100
      nt1050=nt1050+1
      if (nt1050.ne.2) goto 700
      ndim=ibuf(34)+ibuf(35)*ibuf(36)*2-1-2
      if (ndim.gt.mbuf) stop 'out of buffer space'
      goto 700

  800 read(iunit,'(i15)',end=999) ncount
      read(iunit,'(6d24.17)',end=999) (dbuf(i),i=1,ncount)
      write(ounit) ncount,(dbuf(i),i=1,ncount)
      goto 800

  999 end
