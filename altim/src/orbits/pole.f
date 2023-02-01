      subroutine pole(mjd,x,y,n)

* Variables for GEODYN flux file format

      character*8 c8
      integer*4 i4(48), ii4(96)
      integer*4 x(*),y(*)
      integer*4 i,j,k,n
      real*8    mjd(*)
      real*8    interval
      integer*4 na1utc,na1ut1,npole
      integer*4 start, end, startmjd, endmjd
      integer*4 mdate

* Read the polar motion values from the GOEDYN flux file

      open(unit=10,status='old',form='unformatted',
     |     file='/uwa/pieter/Geodyn/Support/FLUX.FINAL')

*  Header record
      read(10) c8, (i4(i),i=1,46)
      start    = nint(i4(3)/1D2)
      end      = nint(i4(4)/1D2)
      interval = i4(27)/24.0D0
      startmjd = mdate(2,start)
      endmjd   = mdate(2,end)

*  A1-UTC header + data records
      read(10) c8, (i4(i),i=1,46)
      na1utc=i4(1)
      do j=1,na1utc
        read(10) (i4(i),i=1,48)
      enddo
*  A1-UT1 header + data records
      read(10) c8, (i4(i),i=1,46)
      na1ut1=i4(1)
      do j=1,na1ut1
        read(10) (i4(i),i=1,48)
      enddo
*  pole header + data records
      read(10) c8, (i4(i),i=1,46)
      npole=i4(1)
      do j=1,npole
        read(10) (ii4(i),i=1,96)
        do k=1,48
          n=n+1
          x(n)=ii4(2*k-1)
          y(n)=ii4(2*k)
          mjd(n)=startmjd+(n-1)*interval
        enddo
      enddo
      close(10)

      end
