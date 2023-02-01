c
c     this program computes the orbit differences in the radial,
c     along-track, and cross-track direction between two POE's
c
c     usage: orbdif POE1 POE2 DIFF
c
      implicit real*8(a-h,o-z)
      real*8 x1(6),x2(6),gtrack(6)
      real*8 pos(3),rmspos(3),rmnpos(3)
      character*100 input1,input2,output
      
      call getarg(1,input1)
      call getarg(2,input2)
      call getarg(3,output)

      write(6,*) 'PLEASE NOTE:'
      write(6,*) 'orbdif should not take longer than 5 minutes'
      write(6,*) 'for orbital arcs shorter than 1 month!!'
      write(6,*)
      
      if (input1.eq.' ') stop
     c 'usage: orbdif POE1 POE2 DIFF'
      
      open(unit=10,status='old',form='formatted',file=input1)
      open(unit=20,status='old',form='formatted',file=input2)
      if (output.ne.' ')
     c open(unit=30,status='unknown',form='formatted',file=output)

      nobs=0
      do 10 i=1,3
      rmspos(i)=0.d0
      rmnpos(i)=0.d0
   10 continue
      
c
c formats
  120 format (2f10.4,f15.4,4d15.8)
  130 format (3f10.4)

      iunit1=10
      iunit2=20

   20 continue
      
      call rdrec(iunit1,x1,gtrack,time1,iend) 
      call rdrec(iunit2,x2,gtrack,time2,iend) 
      if (iend.eq.1) goto 40
      
  210 if (time2-time1.gt.0.5) call rdrec(iunit1,x1,gtrack,time1,iend) 
      if (iend.eq.1) goto 40
      if (time2-time1.gt.0.5) goto 210
  220 if (time1-time2.gt.0.5) call rdrec(iunit2,x2,gtrack,time2,iend) 
      if (iend.eq.1) goto 40
      if (time1-time2.gt.0.5) goto 220
     
c
c     compute position differences
      call conver(x1(1),x1(2),x1(3),x1(4),x1(5),x1(6),
     c  x1(1)-x2(1),x1(2)-x2(2),x1(3)-x2(3),
     c  pos(1),pos(2),pos(3))

c
c     update statistics
      nobs=nobs+1
      do 30 i=1,3
      rmspos(i)=rmspos(i)+pos(i)**2
      rmnpos(i)=rmnpos(i)+pos(i)
   30 continue

c
c     compute geographical location
      call xyzged(gtrack(1),gtrack(2),gtrack(3),rlat,rlon,height)
     
c
c     write position differences
      if (output.ne.' ')
     c write(30,120) rlat,rlon,height,(pos(i)*1.d2,i=1,3),time1
    
      goto 20
      
   40 continue
      
c
c     write statistics
      write(6,*) 'number of diff.: ',nobs
      write(6,*) 'rms of position (rad,alo,cros): '
      write(6,130) (dsqrt(rmspos(i)/nobs)*1.d2,i=1,3)
      write(6,*) 'mean of position (rad,alo,cros): '
      write(6,130) (rmnpos(i)/nobs*1.d2,i=1,3)
      
      end
     
      subroutine conver(x,y,z,xdot,ydot,zdot,ax,ay,az,arad,aalo,acro)   
**********************************************************************  
*                                                                       
*     this subroutine converts the vector ax,ay,az (in the              
*     rectangular) frame to the vector arad,alo,acro (in the            
*     satellite local vertical, local horizontal frame)        
*                                                                       
**********************************************************************  
      implicit real*8 (a-h,o-z)                                         
*                                                                       
      range=sqrt(x**2+y**2+z**2)                                        
      vel=sqrt(xdot**2+ydot**2+zdot**2)                                 
      rinp=x*xdot+y*ydot+z*zdot                                         
      velx=(xdot-rinp*x/range**2)                                       
      vely=(ydot-rinp*y/range**2)                                       
      velz=(zdot-rinp*z/range**2)                                       
      velno=sqrt(velx**2+vely**2+velz**2)                               
      velcro=sqrt((zdot*y-ydot*z)**2+(ydot*x-xdot*y)**2+                
     c            (xdot*z-zdot*x)**2)                                   
*                                                                       
      arad=(x*ax+y*ay+z*az)/range                                       
      aalo=(velx*ax+vely*ay+velz*az)/velno                              
      acro=((ydot*z-zdot*y)*ax+(zdot*x-xdot*z)*ay+(xdot*y-ydot*x)*az)   
     c      /velcro                                                     
*                                                                       
      return
      end                                                               
      
      subroutine xyzged(x,y,z,rlat,rlon,height)
c***********************************************************************
c
c  routine to convert from body fixed x,y,z (m) to 
c  geodetic latitude (deg), longitude (deg), and height (m)
c
c  ae=6378136.3
c  flat=1/298.257
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      logical switch
      
      switch=.false.
      
c     constants relating to the reference ellipsoid
      
      ae=6378136.3d0
      flat=1d0/298.257
      fl2=1d0/(1d0-flat)**2
      ecc=dsqrt(1d0-(1d0-flat)**2)
     
      pi=4d0*datan(1d0)
      degrad=pi/180d0
      
c     compute location on reference ellipsoid (xyref,zref),
c     and height by iterative procedure
c     eps = convergence criterion
      eps=1e-6
      diff=10.0

c     compute geodetic longitude and latitude
c     rlatc : geocentric latitude, first estimate !!
      xy=dsqrt(x**2+y**2)
      rlon=datan2(y,x)/degrad
      rlatc=datan2(z,xy)/degrad
      rlat=datan(fl2*tan(rlatc*degrad))/degrad
      help=dsqrt(1d0-ecc**2*(dsin(rlat*degrad))**2)
      xyref=ae*dcos(rlat*degrad)/help
      zref=ae*(1d0-ecc**2)*dsin(rlat*degrad)/help
      height=dsqrt((xy-xyref)**2+(z-zref)**2)
      radius=dsqrt(xy**2+z**2)
      radellips=dsqrt(xyref**2+zref**2)
      if (radellips.gt.radius) switch=.true.
      if (switch) height=-height

c     iterations
   10 continue
      heightold=height
      if (diff.gt.eps) then
      xyiter=dsqrt(x**2+y**2)-height*dcos(rlat*degrad)
      ziter=z-height*dsin(rlat*degrad)
      rlatc=datan2(ziter,xyiter)/degrad
      rlat=datan(fl2*tan(rlatc*degrad))/degrad
      help=dsqrt(1d0-ecc**2*(dsin(rlat*degrad))**2)
      xyref=ae*dcos(rlat*degrad)/help
      zref=ae*(1d0-ecc**2)*dsin(rlat*degrad)/help
      height=dsqrt((xy-xyref)**2+(z-zref)**2)
      if (switch) height=-height
      diff=dabs(height-heightold)
      goto 10
      endif
      
      return
      end
      
      subroutine rdrec(iunit,x,gtrack,time,iend) 
      implicit real*8(a-h,o-z)
      real*8 x(6),gtrack(6)
      
      iend=0
      
c 
c     formats input files
  100 format (2x,i6,i2.2,i2.2,10x,5d22.16)
  110 format (6d22.16)
      
      read (iunit,100,end=10) idate,ihour,imin,rsec,
     c  grhran,xpole,ypole,et
      read (iunit,110) (x(j),j=1,6)
      read (iunit,110) gtrack
      read (iunit,*)

c
c     compute time in seconds
      imjd=mdate(2,idate)
      time=(imjd-48500)*24.0*3600.d0+ihour*3600.0+imin*60.0+rsec
      
      return
   10 iend=1
      
      return
      end
