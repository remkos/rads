**XGFSTCK -- Search location in XGF file and build stack
*+
      SUBROUTINE XGFSTCK (FILENM, RLON, RLAT, M, TX, TY, C)
      CHARACTER*(*) FILENM
      INTEGER M
      REAL RLAT, RLON, TX(*), TY(*), C(10)
*-
      integer i,lat,lon,nrec,itime,unit,freeunit,hgt
      integer irec/0/,dum,itime0/0/
      integer*2 sig
      character spec*4,spec2*2
      real*8 sumx,sumy,sumxx,sumxy,sumyy,uxx,uxy,uyy,a,b,r,f,x,y,
     |		miny,maxy,d,dmin

* Open XGF file on new unit

      i=index(filenm,' ')-1
      write (0,600) filenm(:i)
550   format (a)
600   format ('Scanning file ',a,' ...')
      unit=freeunit()
      open (unit,file=filenm,status='old',form='unformatted',
     |   access='direct',recl=18)
      read (unit,rec=1) spec,nrec,dum,dum,spec2
      if (spec.ne.'@XGF' .or. spec2.ne.'Tm') then
	 write (0,550) '... wrong file type. Exit'
	 close (unit)
	 return
      endif

* Scan for closest point

      dmin=1e30
      do i=1,nrec
	 read (unit,rec=i+1,err=50) itime,lat,lon,hgt,sig
	 if (sig.lt.0) goto 50
	 d=(lat/1e6-rlat)**2+(lon/1e6-rlon)**2
	 if (d.lt.dmin) then
	    dmin=d
	    irec=i
	 endif
50	 continue
      enddo

* Read the stack

      sumx=0d0
      sumy=0d0
      sumxy=0d0
      sumxx=0d0
      sumyy=0d0
      miny=1d35
      maxy=-1d35

      m=0
10    read (unit,rec=irec+m+2,err=20) itime,lat,lon,hgt,sig
      if (m.eq.0) itime0=itime
      if (sig.lt.0) then
         m=m+1
         tx(m)=(itime-itime0)/86400
	 ty(m)=hgt/1e6
	 x=tx(m)/36525d0
	 y=ty(m)
         sumx =sumx +x
         sumy =sumy +y
         sumxy=sumxy+x*y
         sumxx=sumxx+x*x
         sumyy=sumyy+y*y
	 miny =min(miny,y)
	 maxy =max(maxy,y)
	 goto 10
      endif
20       uxx=m*sumxx-sumx*sumx
         uxy=m*sumxy-sumx*sumy
         uyy=m*sumyy-sumy*sumy
         b=uxy/uxx
	 a=(sumy-b*sumx)/m
         r=uxy/sqrt(uxx*uyy)
	 f=sumyy+b*b*sumxx+m*a*a-2*b*sumxy-2*a*sumy+2*a*b*sumx
	 c(1)=sumy/m
	 c(2)=sqrt(sumyy/m)
	 c(3)=sqrt(sumyy/m-(sumy/m)**2)
	 c(4)=a
	 c(5)=b
	 c(6)=r
	 c(7)=sqrt(f/m)
	 c(8)=sumy
	 c(9)=miny
	 c(10)=maxy
      close (unit)

      end
