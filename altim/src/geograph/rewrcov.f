      program rewrcov

      integer nparam,i,j,n,m
      parameter (nparam=5035)
      integer ics(nparam),iord(nparam),ideg(nparam)
      integer jcs(nparam),jord(nparam),jdeg(nparam)
      real*4 solve(nparam),sigma(nparam),scale(nparam),norm(nparam)
      real*4 jsolve(nparam),jsigma(nparam)
      real*4 dum
      character*80 filenm

      call getarg(1,filenm)
      open (10,file=filenm,status='old',form='unformatted')

      read (10) n
      do i=1,n
	 read (10) ics(i),ideg(i),iord(i),solve(i),sigma(i),scale(i)
	 write (*,*) ics(i),ideg(i),iord(i)
      enddo

      call getarg(2,filenm)
      open (20,file=filenm,status='old')
      call getarg(3,filenm)
      open (30,file=filenm,status='new',form='unformatted')

      read (20,*) dum,dum,dum
      m=nint(dum)
      if (m.lt.n) stop 'not possible'
      read (20,*) dum
      write (*,*) dum
      read (20,'(7(5x,i3,2i4,1x))') (jcs(i),jdeg(i),jord(i),i=1,m)
      read (20,*) dum
      write (*,*) dum
      read (20,*) (jsolve(i),i=1,m)
      read (20,*) dum
      write (*,*) dum
      read (20,*) (jsigma(i),i=1,m)

      do i=1,n
         do j=1,m
	    if (ics(i).eq.jcs(j) .and. ideg(i).eq.jdeg(j)
     |		.and. iord(i).eq.jord(j)) then
     		if (solve(i)/scale(i).ne.jsolve(j))
     |	write (*,*) ics(i),ideg(i),iord(i),solve(i),jsolve(j)
	       sigma(i)=sqrt(jsigma(j))/scale(i)
	    endif
         enddo
      enddo

      write (30) n
      do i=1,n
         write (30) ics(i),ideg(i),iord(i),solve(i),sigma(i),scale(i)
      enddo
      do i=1,n
         read (10) (norm(j),j=1,i)
         write (30) (norm(j),j=1,i)
      enddo
      end
