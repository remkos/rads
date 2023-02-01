      program rewrnor

      integer nparam,npar,i,j,n,m
      parameter (npar=5088,nparam=5088)
      integer ics(nparam),iord(nparam),ideg(nparam)
      integer par(npar)
      integer h(nparam)
      real*8 jnorm(npar)
      real*8 dum
      real*4 ata(nparam*(nparam+1)/2),atb(nparam)
      real*4 solve(nparam),sigma(nparam),scale(nparam)
      character*80 filenm

      call getarg(1,filenm)
      open (20,file=filenm,status='old')

      h(1)=0
      do i=2,npar
         h(i)=h(i-1)+(i-1)
      enddo

      read (20,*) dum,dum,dum
      m=nint(dum)
      if (m.gt.npar) stop 'm>npar'
      if (m.lt.n) stop 'not possible'
      read (20,*) dum
      write (*,*) dum
      read (20,'(7(i5,i3,2i4,1x))')
     |		(par(i),ics(i),ideg(i),iord(i),i=1,m)
      read (20,*) dum
      write (*,*) dum
      read (20,220) (solve(i),i=1,m)
      read (20,*) dum
      write (*,*) dum
      read (20,220) (sigma(i),i=1,m)

      call getarg(2,filenm)
      open (30,file=filenm,status='new',form='unformatted')
      n=3717
      write (30) n
      do i=1,n
	 scale(i)=1
         write (30) ics(i),ideg(i),iord(i),solve(i),sigma(i),scale(i)
      enddo

      do i=1,n
         read (20,*) dum
         write (*,*) dum
         read (20,220) (jnorm(j),j=i,m),dum
	    atb(i)=-dum ! inverse sign of coefficient corrections !
	    do j=i,n
	       ata(h(j)+i)=jnorm(j)
	    enddo
	    write (*,*) dum,ata(h(i)+i)
      enddo
      close (20)

      do i=1,n
         write (30) (ata(h(i)+j),j=1,i)
      enddo
      write (30) (atb(i),i=1,n)
      close (30)
220   format (5d24.14)
      end
