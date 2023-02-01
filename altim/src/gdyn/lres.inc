* lres.inc: Include file shared by lres.f and resedit.f

      integer ii,di,returncode,curs
      logical select
      real rgb(3,6),xmin,xmax,ymin,ymax

      common /clres/ ii,di,select,returncode,curs,
     |		rgb,xmin,xmax,ymin,ymax
