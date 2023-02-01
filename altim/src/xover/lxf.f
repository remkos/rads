      program lxf

      character arg*80,filenm*80/' '/,spec*4
      integer*4 trk/0/
      integer*4 i,iargc,nrec,npar
      logical swap/.false./

      do i=1,iargc()
         call getarg(i,arg)
	 if (arg(1:4).eq.'trk=') then
	    read (arg(5:),*) trk
	 else if (arg(:2).eq.'-S') then
	    swap=.true.
	 else if (arg(:2).eq.'-h') then
	    goto 1300
	 else
	    filenm=arg
	 endif
      enddo
	 
      if (filenm.eq.' ') goto 1300

      open (10,file=filenm,recl=12,access='direct',status='old',
     .form='unformatted')
      read (10,rec=1) spec,nrec,npar
      if (swap) then
      	 call i4swap(1,nrec)
	 call i4swap(1,npar)
      endif
      close (10)
      if (spec.eq.'@XGF') then
	 call xgf(filenm,nrec,swap)
      else if (spec.eq.'@XTF') then
	 call xtf(filenm,nrec,3,trk,swap)
      else if (spec.eq.'@XTE') then
	 call xtf(filenm,nrec,npar,trk,swap)
      else if (spec.eq.'@XTB') then
	 call xtf(filenm,nrec,npar,trk,swap)
      else if (spec.eq.'@XAF') then
	 call xaf(filenm,nrec,3,trk,swap)
      else if (spec.eq.'@XAE') then
	 call xaf(filenm,nrec,npar,trk,swap)
      else if (spec.eq.'@XAB') then
	 call xaf(filenm,nrec,1,trk,swap)
      else if (spec.eq.'@XXF') then
	 call xxf(filenm,nrec,3,trk,swap)
      else if (spec.eq.'@XXS') then
	 call xxs(filenm,nrec,trk,swap)
      else if (spec.eq.'@XXO') then
	 call xxo(filenm,nrec,   0,trk,swap)
      else if (spec.eq.'xXXO') then
	 call xxo(filenm,nrec,npar,trk,swap)
      else if (spec.eq.'@XXE') then
	 call xxf(filenm,nrec,npar,trk,swap)
      else if (spec.eq.'@XXB') then
	 call xxb(filenm,nrec,trk,swap)
      else if (spec.eq.'@STK') then
	 call stk(filenm,nrec,swap)
      else if (spec.eq.'@AUX') then
	 call aux(filenm,nrec,npar,swap)
      else
	 write (6,*) 'lxf: unknown file type: ',spec
      endif
      goto 9999

1300  write (0,1301)
1301  format
     .('lxf: list any xover file format: XAF, XTF, XXF, XGF, ...'//
     .'usage: lxf [options] filename'//
     .'where [options] are:'/
     .'  trk=nnn : select only tracknumber nnn')
9999  end

      subroutine xgf(filenm,nrec,swap)
      character filenm*(*),spec*4,spec2*2
      integer*2 i2,i4_1(2)
      integer*4 i4(4),i4_2(3)
      integer*4 irec,nrec,it0,it1
      logical norm,swap
      equivalence (i4_1,i4(1)),(i4_2,i4(2))
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=18)
      read (10,rec=1) spec,nrec,it0,it1,spec2
      if (swap) then
      	 call i4swap(1,nrec)
      	 call i4swap(1,it0)
      	 call i4swap(1,it1)
      endif
      norm=(spec2.eq.'Tm'.or.spec2.eq.'Mn')
      do irec=2,nrec+1
	 read (10,rec=irec) i4,i2
	 if (swap) then
	    call i4swap(4,i4)
	    call i4swap(1,i2)
	 endif
	 if (i2.gt.0.and.norm) then
	    write (6,*) i4_1,i4_2,i2
	 else
	    write (6,*) i4,i2
	 endif
      enddo
      end

      subroutine xtf(filenm,nrec,npar,itrk,swap)
      character*(*) filenm
      integer*2 i2(4),flags
      integer*4 i4(6),par(10)
      integer*4 irec,nrec,npar,length,i,itrk
      logical swap
      length=34+npar*8
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=length)
      do irec=2,nrec+1
	 read (10,rec=irec) i2,i4,(par(i),i=1,npar*2),flags
	 if (swap) then
	    call i2swap(4,i2)
	    call i4swap(6,i4)
	    call i4swap(npar*2,par)
	 endif
	 if (itrk.eq.0 .or. itrk.eq.i2(1))
     1		write (6,*) i2,i4,(par(i),i=1,npar*2),flags
      enddo
      end

      subroutine xxf(filenm,nrec,npar,itrk,swap)
      character*(*) filenm
      integer*2 is(2),it(2)
      integer*4 i4(4),ssh(2),fun(10)
      integer*4 irec,nrec,npar,length,itrk,i
      logical swap
      length=24+npar*8
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=length)
      do irec=2,nrec+1
	 read (10,rec=irec) i4,it,ssh,(fun(i),i=3,npar*2),is
	 if (swap) then
	    call i4swap(4,i4)
	    call i2swap(2,it)
	    call i4swap(2,ssh)
	    call i4swap(npar*2,fun)
	    call i2swap(2,is)
	 endif
	 if (itrk.eq.0 .or. itrk.eq.it(1) .or. itrk.eq.it(2))
     1		 write (6,*) i4,it,ssh,(fun(i),i=3,npar*2),is
      enddo
      end

      subroutine xxs(filenm,nrec,itrk,swap)
      character*(*) filenm
      integer*2 it(2)
      integer*4 i4(6),ssh(2)
      integer*4 irec,nrec,length,itrk
      logical swap
      length=36
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=length)
      do irec=2,nrec+1
	 read (10,rec=irec) i4,it,ssh
	 if (swap) then
	    call i4swap(6,i4)
	    call i2swap(2,it)
	    call i4swap(2,ssh)
	 endif
	 if (itrk.eq.0 .or. itrk.eq.it(1) .or. itrk.eq.it(2))
     1	 write (6,*) i4,it,ssh
      enddo
      end

      subroutine xxb(filenm,nrec,itrk,swap)
      character*(*) filenm
      integer*2 it(2),sigma(2)
      integer*4 i4(4),ssh(4),omega(2)
      integer*4 irec,nrec,length,itrk
      logical swap
      length=48
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=length)
      do irec=2,nrec+1
         read (10,rec=irec) i4,it,ssh,omega,sigma
	 if (swap) then
	    call i4swap(4,i4)
	    call i2swap(2,it)
	    call i4swap(4,ssh)
	    call i4swap(2,omega)
	    call i2swap(2,sigma)
	 endif
	 if (itrk.eq.0 .or. itrk.eq.it(1) .or. itrk.eq.it(2))
     1      write (6,*) i4,it,ssh,omega,sigma
      enddo
      end

      subroutine xxo(filenm,nrec,npar,itrk,swap)
      character*(*) filenm
      integer*2 it(2),aux(8)
      integer*4 i4(6),ssh(2),orb(2)
      integer*4 irec,nrec,npar,length,itrk,i
      logical swap
      length=44+npar*4
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=length)
      do irec=2,nrec+1
	 read (10,rec=irec) i4,it,ssh,orb,(aux(i),i=1,npar*2)
	 if (swap) then
	    call i4swap(6,i4)
	    call i2swap(2,it)
	    call i4swap(2,ssh)
	    call i4swap(2,orb)
	    call i2swap(npar*2,aux)
	 endif
	 if (itrk.eq.0 .or. itrk.eq.it(1) .or. itrk.eq.it(2))
     |		 write (*,'(2i12,2(i12,1x,i6.6),2i7,4i12,10i7)')
     |           i4,it,ssh,orb,(aux(i),i=1,npar*2)
      enddo
      end

      subroutine xaf(filenm,nrec,npar,itrk,swap)
      character*(*) filenm
      integer*4 i4(6)
      integer*2 i2(2)
      integer*4 irec,nrec,npar,length,itrk
      logical swap
      length=28+(npar-1)*4
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=length)
      do irec=2,nrec+1
	 read (10,rec=irec) i4,i2
	 if (swap) then
	    call i4swap(6,i4)
	    call i2swap(2,i2)
	 endif
	 if (itrk.eq.0 .or. itrk.eq.i2(2))
     1		 write (*,'(6i12,2i6)') i4,i2
      enddo
      end

      subroutine stk(filenm,nrec,swap)
      character filenm*(*),spec*4
      integer*2 i2(10),i2_1(4)
      integer*4 i2_5(3)
      equivalence (i2_1,i2(1)),(i2_5,i2(5))
      integer*4 irec,nrec,i
      logical swap
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=20)
      read (10,rec=1) spec,nrec
      irec=1
10    irec=irec+1
      read (10,rec=irec) i2
      if (swap) call i2swap(10,i2)
      write (6,*) i2_1,i2_5
      do i=1,int(i2(4)),10
	 irec=irec+1
         read (10,rec=irec) i2
	 if (swap) call i2swap(10,i2)
         write (*,'(10i7)') i2
      enddo
      if (irec.lt.nrec+1) goto 10
      end

      subroutine aux(filenm,nrec,npar,swap)
      character*(*) filenm
      integer*2 i2(10)
      integer*4 nrec,npar,length,irec0,irec,i
      logical swap
      length=2*npar
      irec0=11/length+1
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=length)
      do irec=irec0+1,irec0+nrec
	 read (10,rec=irec) (i2(i),i=1,npar)
	 if (swap) call i2swap(npar,i2)
     	 write (*,'(20i7)') (i2(i),i=1,npar)
      enddo
      end
