c
c read binary FLUX table and write to ASCII format
c
      program flux2asc
      implicit none
      character*80 input,output
      character*8 title
      real*8 r8(24)
      integer*4 i4(46),i2(96)
      integer*4 na1utc,na1ut1,npole,nsolar,nmagap,nmagkp,na1utcf,
     |          nkpday,na1utcdat,na1ut1dat,na1ut1tot,nlast,na1utctot,
     |		ndata,i,iunit/10/,ounit/20/

      call getarg(1,input)
      call getarg(2,output)
      open(unit=iunit,status='old',form='unformatted',file=input)
      if (output.eq.' ' .or. output.eq.'-') then
	 ounit=6
      else
         open(unit=ounit,status='unknown',form='formatted',file=output)
      endif

  100 format(a8,3i12/6i12/6i12/6i12/6i12/6i12/6i12/6i12/i12)
  200 format(2d24.17/2d24.17/2d24.17/2d24.17/2d24.17/2d24.17
     |  /2d24.17/2d24.17/2d24.17/2d24.17/2d24.17/2d24.17)
  300 format(12i6/12i6/12i6/12i6/12i6/12i6/12i6/12i6)

c header record for binary file
      read(iunit) title,i4
      write(ounit,100) title,i4

c number of A1-UTC values per record
      na1utc=i4(11)
c number of A1-UTC values per record
      na1ut1=i4(12)
c number of pole values per record
      npole=i4(13)
c number of solar flux values per record
      nsolar=i4(14)
c number of magnetic Ap values per record
      nmagap=i4(15)
c number of magnetic Kp values per record
      nmagkp=i4(16)
c number of A1-UTC values in file
      na1utcf=i4(17)
c number of Kp values per day
      nkpday=i4(17)

c header record A1-UTC data
      read(iunit,end=999) title,i4
      write(ounit,100) title,i4

c number of A1-UTC data records
      na1utcdat=i4(1)
c number of non-zero values in last record
      nlast=i4(2)
c total number of A1-UTC values in file
      na1utctot=i4(3)

c read and write A1-UTC data records
      do i=1,na1utcdat
         read(iunit,end=999) r8
         write(ounit,200) r8
      enddo

c header record A1-UT1 data
      read(iunit,end=999) title,i4
      write(ounit,100) title,i4

c number of A1-UT1 data records
      na1ut1dat=i4(1)
c number of non-zero values in last record
      nlast=i4(2)
c total number of A1-UT1 values in file
      na1ut1tot=i4(3)

c read and write A1-UT1 data records
      do i=1,na1ut1dat
         read(iunit,end=999) r8
         write(ounit,200) r8
      enddo

c read additional data

   30 continue

c header record
      read(iunit,end=999) title,i4
      write(ounit,100) title,i4

c number of data records
      ndata=i4(1)

c read and write data records
      do i=1,ndata
         read(iunit,end=999) i2
         write(ounit,300) i2
      enddo

      goto 30

  999 end
