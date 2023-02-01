      program genaccel9

* Program generates ACCEL9 cards for GIIS setup.
*
*-
* 10-Sep-1998 - Created by Remko Scharroo
*  3-Nov-1999 - Added DATEARG function
*  5-Oct-2001 - Added 1cpr and const options, 
*                 maneuver file optional - Eelco Doornbos
* 21-Oct-2002 - ED: fixed num. precision bug in "if (dtn.eq.0)" line
*-----------------------------------------------------------------------
      real*8 t0,t1,dt,dint,s0,s1
      integer*4 id,i,iargc,nman,mman
      parameter (mman=300)
      real*8 m_acc(3,mman),m_t0(mman+1),m_t1(0:mman)
      integer*4 yy,mm,dd,hh,nn,mdate
      real*8 ss,dtn,dtm,acc(3),dum
      character*80 arg,man
      character*3 cpr,const
      logical datearg

600   format ('genaccel9 - Generate ACCEL9 cards for GIIS setup'//
     |'syntax: genaccel9 [options]'//
     |'Required options are:'/
     |'  t=t0,t1    : Time interval (MJD or [YY]YYMMDD[HHMMSS])'/
     |'           ... or use mjd=, doy=, ymd=, sec='/
     |'  id=id      : Satellite ID'/
     |'Optional options are:'/
     |'  man=file   : Specify filename of ESOC manoeuvre file'/
     |'  dt=dt      : Time step (days) (def: 0=arclength)'/
     |'  int=dt     : Integration step size (sec; def:60)'/
     |'  sig=s0,s1  : Specify sigmas for empiricals and maneuvres',
     |' (def: 1d-7, 1d-1)'/
     |'  1cpr=str   : Specify string to request 1-CPR accelerations'/
     |'               in a=along c=cross and r=radial directions',
     |' (def: ac)'/
     |'  const=str  : Specify string to request constant accelerations'/
     |'               in a=along c=cross and r=radial directions',
     |' (def: none)')
     
* Initialise

      t0=0
      t1=0
      dt=0
      dint=60
      s0=1d-7
      s1=1d-1
      man=' '
      id=0
      nman = 0
      cpr='ac '
      const='   '

* Scan arguments
     
      do i=1,iargc()
         call getarg(i,arg)
	 if (datearg(arg,t0,t1,dum)) then
	    t0=t0/86400d0+46066d0	! Conversion to MJD
	    t1=t1/86400d0+46066d0
	 else if (arg(:3).eq.'id=') then
	    read (arg(4:),*) id
	 else if (arg(:4).eq.'man=') then
	    man=arg(5:)
	 else if (arg(:3).eq.'dt=') then
	    read (arg(4:),*) dt
	 else if (arg(:4).eq.'int=') then
	    read (arg(5:),*) dint
	 else if (arg(:4).eq.'sig=') then
	    read (arg(5:),*) s0,s1
         else if (arg(:5).eq.'1cpr=') then
            cpr = arg(6:8)
         else if (arg(:6).eq.'const=') then
            const = arg(7:9)
	 endif
      enddo

      if (t0.eq.0 .or. t1.eq.0 .or. id.eq.0) then
         write (0,600)
	 goto 9999
      endif

      dint=dint/86400d0
      	 
* Read maneuver file. Round start and end time to nearest step size. 

      if(.not.man.eq.' ') then
      open (10,file=man,status='old')
10    read (10,15,end=50) yy,mm,dd,hh,nn,ss,acc
      nman=nman+1
      if (nman.gt.mman) call fin('genaccel9: too many maneuvers')
      m_t0(nman)=mdate(2,yy*10000+mm*100+dd)+
     |hh/24d0+nn/1440d0+ss/86400d0
15    format (i4,4i3,f11.7,3d15.8)
      read (10,15,end=50) yy,mm,dd,hh,nn,ss,acc
      m_t1(nman)=mdate(2,yy*10000+mm*100+dd)+
     |hh/24d0+nn/1440d0+ss/86400d0
      dtm=m_t1(nman)-m_t0(nman)
      m_t0(nman)=t0+nint((m_t0(nman)-t0)/dint)*dint
      m_t1(nman)=t0+nint((m_t1(nman)-t0)/dint)*dint    
      dtn=m_t1(nman)-m_t0(nman)
      if (dtn.lt.1e-8) then
         m_t0(nman)=m_t0(nman)-dint
	 dtn=dint
      endif
      m_acc(1,nman)=acc(2)*1d3*dtm/dtn
      m_acc(2,nman)=acc(3)*1d3*dtm/dtn
      m_acc(3,nman)=acc(1)*1d3*dtm/dtn
      if (m_t0(nman).le.t0+dint .or. m_t1(nman).ge.t1-dint) nman=nman-1
      goto 10
50    close (10)
      else
      nman = 0
      end if

* Set times for first and last interval

      m_t1(0)=t0
      m_t0(nman+1)=t1
      
* Cycle through the manoeuvres, print out the cards

      call empacc(id,m_t1(0),m_t0(1),dt,dint,s0,cpr,const)
      do i=1,nman
         call manacc(id,m_t1(i),m_acc(1,i),s1)
	 call empacc(id,m_t1(i),m_t0(i+1),dt,dint,s0,cpr,const)
      enddo
      
9999  end

      subroutine empacc(id,t0,t1,dt,dint,sig,cpr,const)
      real*8 t0,t1,dt,sig,t,sigma,dt1,dint
      integer*4 id,n,i
      character*12 date
      character*3 cpr, const
      sigma=sig
      if (t1-t0.lt.1d0/24) sigma=0	! set sigma to 0 for short intervals
      if (dt.eq.0) then
         n=1
      else
         n=nint((t1-t0)/dt)
	 if (n.lt.1) n=1
      endif
      dt1=(t1-t0)/n
100   format ('ACCEL9    99',5x,i7,20x,a12)
110   format ('ACCEL9  ',i2,7x,i7,d20.13,15x,d13.5)      
      do i=1,n
         t=t0+nint((i*dt1)/dint)*dint
	 call strf1985(date,'%y%m%d%H%M%S',nint((t-46066)*86400))
         write (*,100) id,date
         if(index(cpr,'a').gt.0)   write (*,110) 11,id,0d0,sigma
	 if(index(cpr,'a').gt.0)   write (*,110) 12,id,0d0,sigma
         if(index(const,'a').gt.0) write (*,110) 13,id,0d0,sigma
         if(index(cpr,'c').gt.0)   write (*,110) 21,id,0d0,sigma
	 if(index(cpr,'c').gt.0)   write (*,110) 22,id,0d0,sigma
         if(index(const,'c').gt.0) write (*,110) 23,id,0d0,sigma
         if(index(cpr,'r').gt.0)   write (*,110) 31,id,0d0,sigma
	 if(index(cpr,'r').gt.0)   write (*,110) 32,id,0d0,sigma
         if(index(const,'r').gt.0) write (*,110) 33,id,0d0,sigma
      enddo
      end 

      subroutine manacc(id,t1,acc,sig)
      real*8 t1,acc(3),sig
      integer*4 id
      character*12 date
100   format ('ACCEL9    99',5x,i7,20x,a12)
110   format ('ACCEL9  ',i2,7x,i7,d20.13,15x,d13.5)      
      call strf1985(date,'%y%m%d%H%M%S',nint((t1-46066)*86400))
      write (*,100) id,date
      write (*,110) 13,id,acc(1),sig
      write (*,110) 23,id,acc(2),sig
      write (*,110) 33,id,acc(3),sig
      end         
