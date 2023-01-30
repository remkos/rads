      program geogrtst

* Program to create grid of gravity induced radial orbit errors

      include "geograph.inc"
      integer k,l,m,ip,ltot,mtot
      integer kx,ky,nx,ny
      real*8 dr_c,dr_s,h0,x0,x1,y0,y1,dx,dy,x,y,revs,period,fact,c,s

* Set area for computation of gravity induced radial orbit error

      x0=-180d0 ; x1=180d0 ; dx=1d0
      y0= -80d0 ; y1= 80d0 ; dy=1d0

* Set satellite parameters

      incl=98.56d0  ! inclination [degrees]
      period=35.0d0 ! length of repeat period [days]
      revs=501      ! nr of revolutions per repeat [-]
      deadband=2d-2 ! cutoff frequency [cpr]
      test=.false.

* Read GEODYN formatted GCOEF file

      open (31,file='jgm3',status='old')
600   format (14x,2i3,4x,d20.8,d15.8,d13.1)
      rewind (31)
      read (31,*)
      read (31,600) ltot,mtot,gm,ae,flat

* Read geopotential coefficients and store them in array cs(ip).
* ics(ip)  indicates whether element ip is C (1) or S (2).
* ideg(ip) determines the degree L of the element ip.
* iord(ip) determines the order M of the element ip.
* Note that the order of C and S coefficients is irrelevant.

      ip=0
      lmax=0
10    continue
         read (31,600,end=20) l,m,c,s
         fact=1d0
         ip=ip+1
         ics(ip)=1
         ideg(ip)=l
         iord(ip)=m
         cs(ip)=fact*c
         if (m.ne.0) then
            ip=ip+1
            ics(ip)=2
            ideg(ip)=l
            iord(ip)=m
            cs(ip)=fact*s
         endif
         lmax=max(l,lmax)
         if (ip.gt.npar) stop 'too many parameters (ip>npar)'
      goto 10
20    ipmax=ip

* Convert satellite orbit information

      incl=incl*rad
      n0=(2*pi)/(period*86400d0/revs)
      wmdot=n0
      ogdot=-nint(period)*n0/revs
      a0=(gm/n0**2)**(1d0/3d0)
      h0=a0-ae  

* Initialize F_lmp and D_lmp

      call d_lmp

* Compute the gravity induced radial orbit error on a regular
* geocentric grid

      nx=nint((x1-x0)/dx+1)
      ny=nint((y1-y0)/dy+1)

      do ky=1,ny
         y=y0+(ky-1)*dy
         do kx=1,nx
            k=k+1
            x=x0+(kx-1)*dx
            call geograph(y*rad,x*rad,dr_c,dr_s)
            write (*,610) x,y,dr_c,dr_s
         enddo
      enddo
610   format (2f6.1,2f12.4)
      end
