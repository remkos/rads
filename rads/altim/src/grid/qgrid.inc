      real*8 x0,x1,dx,y0,y1,dy,t0,t1,tbias(2),ssb
      real*8 range0,range1,qpi,rad,factk,factl
      logical ocean,sum,tsort,useaux
      integer field,ptype,trk0,trk1
      common /qgridc/ x0,x1,dx,y0,y1,dy,t0,t1,trk0,trk1,tbias,ssb,
     .range0,range1,qpi,rad,factk,factl,
     .ocean,sum,tsort,useaux,field,ptype
