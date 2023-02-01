* file: ingrid.inc

      integer MA,MBUF,MPOL,MVRS,MOPT
      parameter (MVRS=9905,MOPT=11)
      parameter (MA=20000000,MBUF=200,MPOL=7)
      real*4 a(MA),xmin(MBUF),xmax(MBUF),ymin(MBUF),ymax(MBUF)
      logical used(MBUF),proc(MBUF),blkd(MBUF)
      integer nx(MBUF),ny(MBUF),id(MBUF),idata(MBUF),idnext,jused
      character*80 fname(MBUF)
      common /buffer/ a,used,proc,blkd,nx,ny,fname,id,idata,idnext,
     1		xmin,xmax,ymin,ymax,jused
