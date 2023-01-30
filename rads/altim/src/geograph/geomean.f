      program geomean
      implicit none

      include "geograph.inc"
      real*8 phi,lambda,mean,var,x0,x1,y0,y1,dx,dy,period
      integer cycles,l
      character*80 filenm
      namelist /geograph_nml/ incl,period,cycles,deadband,test,
     | x0,x1,y0,y1,dx,dy


      write (*,*) 'Reading gravity models'
      call getarg (1,filenm)
      if (filenm.eq.'sigma') then
         call getarg(2,filenm)
         call gravrd(1.0,filenm)
      else
         call gravrd(1.0,filenm)
         call getarg(2,filenm)
         call gravrd(-1.0,filenm)
      endif

      filenm='/user/altim'
      call checkenv('ALTIM',filenm,l)
      filenm(l+1:)='/nml/geograph.nml'
      open (7,file=filenm,status='old')
      read (7,geograph_nml)
      close (7)
      open (7,file="geograph.nml",status="old",err=2)
      read (7,geograph_nml)
      close (7)
2     continue

      incl=incl*rad
      n0=(2*pi)/(period*86400d0/cycles)
      wmdot=n0
      ogdot=-nint(period)*n0/cycles
      a0=(gm/n0**2)**(1d0/3d0)

      call getarg (3,filenm)
      test=(filenm.eq.'-t')
      
      write (*,*) 'Initialize Flmp and Dlmp'
      call d_lmp
   10 write (*,551) 'phi,lambda -> '
      read (*,*,end=9999) phi,lambda
      write (*,600) phi,lambda
      call geograph(phi*rad,lambda*rad,mean,var)
      write (*,600) mean,var
      goto 10
  551 format (a,$)
  600 format (2f12.6)
 9999 end
