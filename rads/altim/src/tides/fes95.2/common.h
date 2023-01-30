c-----------------------------------------------------------------------
c extended_max : nombre max d'ondes prises en compte
c delta_max    : intervalle en heure de recalcul des corrections nodales
c-----------------------------------------------------------------------

      integer extended_max
      real*8 pi,rad,delta_max,t_nodal
      character*6 wave
      character*20 model
      real*8 freq   
      real*8 n,p,s,p1,pp,nu,xi,tt,nuprim,nusec,hp,r,iang,x1ra
	parameter (extended_max=27)
	parameter (pi = 3.141592653589793d+00,rad=pi/180d0)
	parameter (delta_max=24.D+00)

                             
      common/wave/      nbwave,f(extended_max),
     &                  v0_u(extended_max),num(extended_max),   
     &                  sreal(extended_max),
     &                  simag(extended_max),wave(extended_max)
      common/wave2/     freq(extended_max)

      common/grid/      ni,nj,xmin,ymin,dx,dy

      common/model/     wrp(0:721,361,8),
     &                  wip(0:721,361,8),spec

      common/model_name/ model

      common/nodal/     n,p,s,p1,pp,nu,xi,tt,nuprim,nusec,hp,r,iang,
     &                  x1ra,t_nodal
