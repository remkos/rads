      program odrhead
      implicit none
      character filenm*80,satel*8,spec*4
      real*8 x,sec85
      integer rep,arc,rec,ver,yymmdd,hh,mm,ss,beg,i,j
      logical ltlend

      call getarg(1,filenm)
      open (20,file=filenm,form='unformatted',status='unknown',
     |	access='direct',recl=16)
      read (5,600) spec,satel,x,i,j,rec,ver,yymmdd,hh,mm,ss
600   format(a4,1x,a8,6x,f7.3,6x,z1,i2,2(6x,i5),6x,i6,3(1x,i2))
      arc=i*100+j
      rep=nint(x*1000)
      beg=nint(sec85(4,yymmdd*1d6+hh*1d4+mm*1d2+ss))
      if (ltlend()) then
         call i4swap(1,beg)
	 call i4swap(1,rep)
	 call i4swap(1,arc)
	 call i4swap(1,rec)
	 call i4swap(1,ver)
      endif
      write (20,rec=1) spec,satel,beg
      write (20,rec=2) rep,arc,rec,ver
      close (20)
      end
