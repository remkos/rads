      integer*4	mhead
      parameter (mhead=104)
      integer*4 ntype,nbuf,nx,ny,head(mhead/4)
      real*8	z0,dz,xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,znan
      equivalence (head(1),ntype),(head(2),nbuf),(head(3),nx),
     |		(head(4),ny),(head(5),z0),(head(7),dz),
     |		(head(9),xmin),(head(11),xmax),(head(13),dx),
     |		(head(15),ymin),(head(17),ymax),(head(19),dy),
     |		(head(21),zmin),(head(23),zmax),(head(25),znan)
      integer*1	tmp_b(0:7)
      integer*2 tmp_s(0:3)
      integer*4 tmp_i(0:1)
      real*4	tmp_f(0:1)
      real*8	tmp_d(0:0)
      equivalence(tmp_b,tmp_s),(tmp_b,tmp_i),(tmp_b,tmp_f),(tmp_b,tmp_d)
      common	/cgridbuff/ tmp_b,head
