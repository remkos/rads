c
c write Geodyn ASCII Residual File to Binary
c
      implicit none
      integer maxobs,i,j,nc,nsta,ifloc,nobs,nm,nloc,mtype,k1,kk1
      parameter(maxobs=1000)
      character*80 input,output
      real*8 rec20(20),rec8(8),recr(9*maxobs)
      real*8 resid(8*maxobs)
      character char8(10)*8,chrsta*8

      call getarg(1,input)
      call getarg(2,output)

      open ( 10 , file = input , form = 'formatted' )
      open ( 20 , file = output , form = 'unformatted' )

c read header
      read(10,'(6d22.15)') rec20
      write(20) rec20
      nc=nint(rec20(1))
      nsta=nint(rec20(10))

c reac GEODYN II unit 5 input
      do 10 i=1,nc
      read(10,'(10a8)') char8
      write(20) char8
   10 continue

c read tracking station coordinate records
      do 20 i=1,nsta
      read(10,'(a8,8d22.15)') chrsta,rec8
      write(20) chrsta,rec8
   20 continue

c read arc header
      read(10,'(6d22.15)') rec20
      write(20) rec20
      ifloc=nint(rec20(7))
      nobs=nint(rec20(8))
      if (nobs.ne.0) stop
     c 'observation data present, modification required'

c lengths record
   30 continue
      read(10,'(6d22.15)',end=999) rec20
      write(20) rec20

      nm=nint(rec20(8))
      nloc=nint(rec20(6))
      mtype=nint(rec20(4))

c location data
      if (mtype.eq.100) then
         read(10,'(d22.15)',end=999) rec20(1)
         nm=nint(rec20(1))
         if (nm.gt.1000) then
            backspace(10)
            goto 30
         endif
         write(20) rec20(1)
      endif

      if (ifloc.gt.0 .and. mtype.ne.99 .and. mtype.ne.100
     c    .and. mtype.gt.6) then
         read(10,'(6d22.15)',end=999) rec20
         write(20) rec20
         do 70 k1=1,nloc
            do 80 i=1,20
            kk1=nint(rec20(i))
            if (kk1.ne.0) then
               if (mod(kk1,2).eq.0) then
               read(10,'(6d22.15)',end=999) (recr(j),j=1,7*nm)
               write(20) (recr(j),j=1,7*nm)
               else
               read(10,'(6d22.15)',end=999) (recr(j),j=1,9*nm)
               write(20) (recr(j),j=1,9*nm)
               endif
            endif
   80       continue
   70    continue
      endif

c residual data records
      read(10,'(6d22.15)',end=999) (resid(i),i=1,8*nm)
      write(20) (resid(i),i=1,8*nm)

      goto 30

  999 end

