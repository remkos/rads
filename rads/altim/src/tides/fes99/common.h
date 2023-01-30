c#######################################################################
c extended_max : nombre max d ondes prises en compte
c delta_max    : intervalle en heure de recalcul des corrections nodales
c iascbin      : fichier ascii iascbin=1
c              : fichier binaire iascbin=0
c#######################################################################

c---- integer*4 --------------------------------------------------------
      integer*4      extended_max
      integer*4      iascbin
      integer*4      nbwave,num,ni,nj

c---- real*4 -----------------------------------------------------------
      real*4         xmin,ymin,dx,dy
      real*4         wrp, wip, wrpload, wipload, spec
      real*4         f,v0_u,sreal,simag,srealload,simagload

c---- real*8 -----------------------------------------------------------
      real*8         pi,delta_max,t_nodal
      real*8         freq   
      real*8         n,p,s,p1,pp,nu,xi,tt,nuprim,nusec,hp,r,iang,x1ra

c---- character --------------------------------------------------------
      character*20   wave
      character*20   model
      character*1024 wave_path

c---- parameter --------------------------------------------------------
      parameter (extended_max=28)
      parameter (pi = 3.141592653589793d+00)
      parameter (delta_max=24.D+00)

c#### common ###########################################################
                             
      common/wave/      nbwave,f(extended_max),
     &                  v0_u(extended_max),num(extended_max),   
     &                  sreal(extended_max),
     &                  simag(extended_max),
     &                  srealload(extended_max),
     &                  simagload(extended_max),
     &                  wave(extended_max),
     &                  iascbin
     
      common/wave2/     freq(extended_max)

      common/grid/      ni,nj,xmin,ymin,dx,dy

      common/model/     wrp(1441,721,8), wip(1441,721,8),
     &                  wrpload(1441,721,8),wipload(1441,721,8),  spec

      common/model_name/ model, wave_path

      common/nodal/     n,p,s,p1,pp,nu,xi,tt,nuprim,nusec,hp,r,iang,
     &                  x1ra,t_nodal
