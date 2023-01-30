      program xtfcmp

      integer nmax
      parameter (nmax=50000)
      character arg*80
      integer*2 i2(4)
      integer*4 i4(6),par1(5),par2(5),sat/4/
      integer*4 ifile,irec,j,j1,j2,i,n,npar,iargc,dt,m,mstep
      real*8 u,u0,u1,u2,d1,d2,pi,rad,rev,rms
      integer t1(nmax),t2(nmax),pnt(nmax),index(nmax)

      pi=4*atan(1d0)
      rad=pi/180

      do i=1,iargc()
	 call getarg(i,arg)
	 if (arg(1:4).eq.'sat=') then
	    read (arg(5:),*) sat
	 else if (arg(1:2).eq.'u=') then
	    read (arg(3:),*) u0,u1
	 else
	    call scanxtf(ifile,arg,n,t1,t2,pnt,sat,npar)
	    ifile=ifile+1
	 endif
      enddo

      do i=1,n
	 index(i)=i
      enddo

      call qsort(index,t1,n)

      do j=1,n-1
	 j1=index(j)
	 j2=index(j+1)
	 dt=t1(j2)-t2(j1)
	 if (dt.le.3000) then
	    ifile=pnt(j1)/40000
	    irec=mod(pnt(j1),40000)
	    read (10+ifile,rec=irec+1) i2,i4,par1
            if (sat.eq.3) then
               rev=6000.0d0
            else if (sat.eq.4) then
               rev=6035.0d0
            else if (sat.eq.5 .or. sat.eq.6) then
               rev=6745.759d0
            else
               write (6,*) i2(2)
               call fin('unknown sat')
            endif
            u1=i4(2)/1d6
            u2=u1+(i4(6)-i4(5))/rev*360

	    ifile=pnt(j2)/40000
	    irec=mod(pnt(j2),40000)
	    read (10+ifile,rec=irec+1) i2,i4,par2
            u1=i4(2)/1d6
	    mstep=1
	    if (dt.lt.0) then
	       mstep=-1
	    else if (u1.lt.u2) then
	       u1=u1+360
	    endif

	    m=0
	    rms=0
	    do i=nint(u2),nint(u1),mstep
               u=i*rad
	       d1=par1(1)+par1(2)*sin(u)+par1(3)*cos(u)+
     .			  par1(4)*sin(2*u)+par1(5)*cos(2*u)
	       d2=par2(1)+par2(2)*sin(u)+par2(3)*cos(u)+
     .			  par2(4)*sin(2*u)+par2(5)*cos(2*u)
	       m=m+1
	       rms=rms+(d1-d2)**2
	    enddo
	    rms=sqrt(rms/m)/1d6
	
	    write (6,*) t1(j1),t2(j1),t1(j2),t2(j2),dt
	    write (6,*) par1
	    write (6,*) par2
	    write (6,'(i10,f10.3)') mstep*m,rms
	 endif
      enddo

      end

      subroutine scanxtf(ifile,filenm,n,t1,t2,pnt,sat,npar)
      integer*4 ifile,unit,n,pnt(*),t1(*),t2(*)
      integer*4 sat,i4(6),par(10),length
      integer*4 npar,nrec,i,j
      integer*2 i2(4),flags
      character spec*4,filenm*80

      unit=10+ifile
      open (unit,file=filenm,form='unformatted',status='old',
     .access='direct',recl=12)
      read (unit,rec=1) spec,nrec,npar
      close (unit)
      if (spec.eq.'@XTF') then
	 npar=3
      else if (spec.eq.'@XTE' .or. spec.eq.'@XTB') then
      else
	 call fin('Input file is not XTF')
      endif

      length=34+npar*8
      open (unit,file=filenm,form='unformatted',status='old',
     .access='direct',recl=length)

      do j=1,nrec
         read (unit,rec=j+1) i2,i4,(par(i),i=1,npar*2),flags
	 if (i2(2).eq.sat.and.flags.ge.256) then
	    n=n+1
	    pnt(n)=ifile*40000+j
	    t1(n)=i4(5)
	    t2(n)=i4(6)
	 endif
      enddo
      end
