      program readcov

      integer nparam,i,j,n,lmax
      parameter (nparam=5035)
      integer ics(nparam),iord(nparam),ideg(nparam)
      real*4 solve(nparam),sigma(nparam),scale(nparam),norm(nparam)
      real*4 corr,cormin
      character*80 filenm

      call getarg(1,filenm)
      open (10,file=filenm,status='old',form='unformatted')
      read (10) n
      call getarg(2,filenm)
      read (filenm,*) cormin
      call getarg(3,filenm)
      read (filenm,*) lmax

      do i=1,n
	 read (10) ics(i),ideg(i),iord(i),solve(i),sigma(i),scale(i)
	 if (ideg(i).le.lmax) write (*,*) ics(i),ideg(i),iord(i),
     |		solve(i)*scale(i),sigma(i)*scale(i)
      enddo
      write (*,600)
      do i=1,n
	 read (10) (norm(j),j=1,i)
	 do j=1,i
            corr=norm(j)/scale(i)/scale(j)
	    if (abs(corr).ge.cormin.and.
     |		ideg(i).le.lmax.and.ideg(j).le.lmax) then
	       write (*,630)
     |         ideg(i),iord(i),ideg(j),iord(j),ics(i)-1,ics(j)-1,corr
	    endif
	 enddo
      enddo
600   format('COVARIANCES GRAVITY MODEL - JGM-2'/
     |'EARTH 000      70 70    0.39860044150000D+150.637813630D+07',
     |'0.2982564D+03')
630   format (4i3,2i2,d17.9)
      end
