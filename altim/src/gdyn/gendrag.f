      program gendrag

* Program generates DRAG cards for GIIS setup. 
* The apriori values are determined from DRAG cards in a fort.7 file
*-
* 10-Sep-1998 - Created by Remko Scharroo
*  3-Nov-1999 - Added DATEARG function
*-----------------------------------------------------------------------
      real*8 t0,t1,dt,val,sig,t,dum
      integer*4 i,id,iargc
      character arg*80
      logical datearg
      
600   format ('gendrag - Generate DRAG cards for GIIS setup'//
     |'syntax: gendrag [options]'//
     |'Required options are:'/
     |'  t=t0,t1    : Time interval (MJD or [YY]YYMMDD[HHMMSS])'/
     |'           ... or use mjd=, doy=, ymd=, sec='/
     |'  dt=dt      : Time step (days)'//
     |'Optional options are:'/
     |'  id=id      : Satellite ID'/
     |'  cd=val,sig : Specify CD values and sigma'/
     |'  ft7=name   : Get id, average and sigma from fort.7 file')
      
* Initialise

      val=0
      sig=0
      t0=0
      t1=0
      dt=0
      id=0

* Read arguments

      do i=1,iargc()
         call getarg(i,arg)
	 if (datearg(arg,t0,t1,dum)) then
	    t0=t0/86400d0+46066d0	! Conversion to MJD
	    t1=t1/86400d0+46066d0
	 else if (arg(:3).eq.'dt=') then
	    read (arg(4:),*) dt
	 else if (arg(:4).eq.'ft7=') then
	    call rdft7(arg(5:),id,val,sig)
	 else if (arg(:3).eq.'cd=') then
	    read (arg(4:),*) val,sig
	 else if (arg(:3).eq.'id=') then
	    read (arg(4:),*) id
	 endif
      enddo

      if (val.eq.0 .or. sig.eq.0 .or. t0.eq.0
     |  .or. t1.eq.0 .or. dt.eq.0 .or. id.eq.0) then
         write (0,600)
         goto 9999
      endif
      
      write (*,100) id,val
      do i=1,nint((t1-t0)/dt)
         t=t0+i*dt
	 call strf1985(arg,'%y%m%d%H%M%S',nint((t-46066)*86400))
	 write (*,100)id,val,arg,sig
100      format('DRAG',13x,i7,f18.10,2x,a12,'.00',f13.6)
      enddo	 
	 
9999  end
      
      subroutine rdft7(file,id,val,sig)
* Read fort.7 file and return average value and sigma for drag cards
      character*(*) file
      integer*4 id,n
      real*8 val,sig,vsum,vmin,vmax,t
      character*80 text
      
      open (10,file=file,status='old')
      n=0
      vsum=0d0
      vmin=1d30
      vmax=-1d30
10    read (10,'(a)',end=20) text
      if (text(1:4).eq.'DRAG') then
	 read (text(60:72),'(d13.1)') t
         if (t.gt.0) then
            read (text,'(17x,i7,d20.8,d15.3,d13.1)') id,val,t,sig
            vsum=vsum+val
            vmin=min(vmin,val)
            vmax=max(vmax,val)
	    n=n+1
	 endif
      else if (text(1:5).eq.'/ARC/') then
         n=0
	 vsum=0
	 vmin=1d30
	 vmax=-1d30
      endif
      goto 10
* Return average value, excluding the highest and lowest value
20    close (10)
      if (n.gt.4) then
         val=(vsum-vmin-vmax)/(n-2)
      else
         val=vsum/n
      endif

      end  
