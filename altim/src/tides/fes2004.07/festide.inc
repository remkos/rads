      integer*4	nw,nw_prime,nw_extra,nw_long,
     |		nx,ny,px,py,ptr_ocean,ptr_load,nb_ocean,nb_load
      logical*4	haveload
      integer*2	undef,tmp_s(0:1)
      parameter	(nw=33)
      real*8	xmin,ymin,xmax,ymax,dx,dy,fx,fy,z0(nw,2),dz(nw,2)
      real*8	f(nw),v0_u(nw),freq(nw)
      real*8	pi,t_nodal,rad
      complex*16 radj
      character	wave(nw)*8
      parameter (undef=-32768)
      parameter (pi=3.141592653589793d0)
      parameter (rad=0.01745329252d0, radj=(0d0,0.01745329252d-2))

      common /cfestide/ f,v0_u,freq,t_nodal,xmin,xmax,ymin,ymax,
     |	dx,dy,fx,fy,z0,dz,wave,nw_prime,nw_extra,nw_long,
     |	nx,ny,px,py,ptr_ocean,ptr_load,nb_ocean,nb_load,haveload,
     |	tmp_s
