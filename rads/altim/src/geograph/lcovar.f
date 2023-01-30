      program lcovar

* List covariance matrix
*-
* 11-Apr-1998 - Remko Scharroo
*-
      integer nparam,i,j,n
      parameter (nparam=5035)
      integer ics(nparam),iord(nparam),ideg(nparam)
      real*4 solve(nparam),sigma(nparam),scale(nparam),norm(nparam)
      character*80 filenm
      logical normal/.true./

      call getarg(1,filenm)
      open (10,file=filenm,status='old',form='unformatted')
      call getarg(2,filenm)
      if (filenm(1:2).eq.'-n') normal=.false.

      read (10) n
      write (*,*) n
      do i=1,n
	 read (10) ics(i),ideg(i),iord(i),solve(i),sigma(i),scale(i)
	 write (*,'(4i5,3d16.8)') i,ics(i),ideg(i),iord(i),
     |		solve(i)*scale(i),sigma(i)*scale(i),scale(i)
      enddo

      do i=1,n
         if (normal) write (*,*) i
         read (10) (norm(j),j=1,i)
	 if (normal) write (*,200) (norm(j),j=1,i)
      enddo

      read (10,end=90) (norm(i),i=1,n)
      write (*,*) n
      write (*,200) (norm(i),i=1,n)
90    close (10)
200   format ((5d16.8))
      end
