*---------
* Load polar motion and time parameters from GEODYN Flux file
*---------
      subroutine geodyn_load(type,data,data2,sec,n)
      implicit none
      integer*4 na1utcf,imjdrfn,mdate
      real*8 pi, rad, rmul
      real*8 r8(30)
      real*8 a1utc(1),a1ut1(1)
      real*8 taiutc(1),taiut1(1)
      real*8 xp(1),yp(1)
      integer*4 i4(60)
      integer*4 i2(288)
      integer*4 na1utc, na1ut1,npole,nsolar,nmagap,nmagkp
      integer*4 nkpday,na1utcdat,nlast,na1utctot,imjdref
      integer*4 i,j,k,l,imjd,iymd,na1ut1dat,na1ut1tot,ndata
      integer*4 npoletot,npold

      real*8 sec85, ymdhms
      real*8 polestart, sec0
 
      character*8 c8
      integer*4 iymd, n
      real*8 data(*), data2(*), sec(*), tflux
      character*8 type

      integer month, year
      integer cal(12), ndays
      data cal /31,28,31,30,31,30,31,31,30,31,30,31/

      character*80 date
*
* Constants
*
      n = 0
      pi=4*datan(1d0)
      rad=pi/180
c BE CAREFULL WITH RMUL FACTOR: SOMETIMES IN MILLIARCSEC
c AND SOMETIMES IN 0.1 MILLIARCSEC !!!
      rmul=rad/3600d4

      open(unit=10,status='old',form='unformatted',
     | file='/uwa/pieter/Geodyn/Support/FLUX.FINAL')

c header record for binary file
      read(10) r8(1),(i4(i),i=1,46)
 
      polestart = sec85(4,i4(3)*1D4)

c Start date for A1-UT1 data
      imjdrfn=mdate(2,i4(3)/100)

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

   40 continue

      read(10,end=999) c8,(i4(i),i=1,46)
      ndata = i4(1)
*      print *, c8, ndata

c POLE data
      if (c8(1:4).eq.'POLE'.and.type(1:4).eq.'POLE') then
        sec(1) = polestart
        do i=1,ndata
          read(10,end=999) (i2(j),j=1,96)
          do j=1,48
            if(i2(j*2-1).ne.0.or.i2(j*2).ne.0) then
              n=n+1
              if(n.gt.1) sec(n)=sec(n-1) + 5 * 86400.0D0
              data(n)=dfloat(i2(j*2-1))*rmul
              data2(n)=dfloat(i2(j*2))*rmul
            endif
          enddo
        enddo

c FLUXS data
      elseif (c8(1:5).eq.'FLUXS'.and.type(1:5).eq.'FLUXS') then
        do i=1,ndata
          read(10,end=999) (i2(j),j=1,96)
          do k=1,3
            do l=1,31
              iymd = i2((k-1)*32+32)*100 + l
              tflux = (i2((k-1)*32+l)) / 10.0D0
              if (tflux.gt.0.0D0) then
                n = n + 1
                sec(n) = sec85(5,dble(iymd))
                data(n) = tflux
              endif
            enddo
          enddo
        enddo 

c FLUXAP data
      elseif (c8(1:6).eq.'FLUXAP'.and.type(1:6).eq.'FLUXAP') then
        do i=1,ndata
          read(10,end=999) (i2(j),j=1,96)
          do k=1,3
            do l=1,31
              iymd = i2((k-1)*32+32)*100 + l
              tflux = dble(i2((k-1)*32+l)) 
              if (tflux.gt.0.0D0) then
                n = n + 1
                sec(n) = sec85(5,dble(iymd))
                data(n) = tflux
*                print *, n, sec(n), data(n)
              endif
            enddo
          enddo
        enddo 

c FLUXKP8 data
      elseif (c8(1:7).eq.'FLUXKP8'.and.type(1:7).eq.'FLUXKP8') then
        do i=1,ndata
          k=0
          read(10,end=999) (i2(j),j=1,96)
          read(10,end=999) (i2(j),j=97,192)
          read(10,end=999) (i2(j),j=193,288)
          iymd = i2(288)*100+1
          sec0 = sec85(5,dble(iymd))
          month = int(ymdhms(sec0)/1d8)
          year = int(month/100)
          month = month - 100*year

          ndays = cal(month)
          if(month.eq.2.and.mod(year,4).eq.0) ndays=ndays+1
          
*          print *, iymd, year, month, ndays

          do j=1, ndays * 8
            n=n+1
            sec(n) = sec0 + k * 10800.0D0 
            data(n) = i2(j)/100.0D0
*            print '(i5,1x,f15.2,1x,f15.2,1x,f10.4)', 
*     |              n, sec(n), ymdhms(sec(n)),data(n)
            k=k+1
          enddo
        enddo

c Other data entries
      else
        do i=1,ndata
          read(10,end=999) (r8(j),j=1,24)
        enddo
      end if

      goto 40

c header record A1-UTC data
      read(10,end=999) r8(1),(i4(i),i=1,46)

c number of A1-UTC data records
      na1utcdat=i4(1)
c number of non-zero values in last record
      nlast=i4(2)
c total number of A1-UTC values in file
      na1utctot=i4(3)

c reference date is 660101
      imjdref=mdate(2,660101)
c read and store A1-UTC data records
      na1utc=0
      do 10 i=1,na1utcdat
      read(10,end=999) (r8(j),j=1,24)
      do 10 j=1,8
      imjd=imjdref+nint(r8(j*3))
      iymd=mdate(1,imjd)
c for leap seconds after 2000
      if (iymd.lt.500000) iymd=iymd+1000000
      na1utc=na1utc+1
      a1utc(na1utc)=r8((j-1)*3+1)
      taiutc(na1utc)=iymd
   10 continue

c header record A1-UT1 data
      read(10,end=999) r8(1),(i4(i),i=1,46)

c number of A1-UT1 data records
      na1ut1dat=i4(1)
c number of non-zero values in last record
      nlast=i4(2)
c total number of A1-UT1 values in file
      na1ut1tot=i4(3)

c read and store A1-UT1 data records
c increments of 5 days
      na1ut1=0
      imjd=imjdrfn
      do 20 i=1,na1ut1dat
      read(10,end=999) (r8(j),j=1,24)
      do 20 j=1,24
      iymd=mdate(1,imjd)
      na1ut1=na1ut1+1
      a1ut1(na1ut1)=r8(j)
      if (iymd.lt.500000) iymd=iymd+1000000
      taiut1(na1ut1)=iymd
      imjd=imjd+5
   20 continue

c read new header record

  999 continue

      close(10)

      return
 
      end
